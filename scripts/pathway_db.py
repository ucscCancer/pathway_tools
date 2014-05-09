#!/usr/bin/env python

import argparse
import sys
import os
import re
from copy import copy
import yaml
import csv
import shutil
import subprocess
import networkx
from pathway_tools import convert as network_convert
from glob import glob
from StringIO import StringIO
from urllib2 import urlopen

import zipfile

def log(message):
    sys.stderr.write(message + "\n")

def error(message):
    sys.stderr.write( "ERROR:" + message + "\n")



reKeyValue = re.compile(r'^(\S+)\s*=\s*(.*)$')

LOCAL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CENTRAL_REPO = "git://github.com/ucscCancer/superpathway_db.git"
LOCAL_REPO = os.path.join(LOCAL_DIR, "superpathway_db")
DATA_REPO = os.path.join(LOCAL_DIR, "data")


class Dogma:
    def __init__(self):
        self.data = None


    def read(self, handle):
        self.data = yaml.load(handle.read())

    def has_edgetype(self, src_type, edge_type, dst_type):
        if edge_type in self.data['interactions']:
            if src_type in self.data['interactions'][edge_type]['src']:
                if dst_type in self.data['interactions'][edge_type]['dst']:
                    return True
        return False

    def has_nodetype(self, node_type):
        return node_type in self.data['dogma_template']

class GraphChecker:

    def load_hugomap(self, path):
        handle = open(path)
        reader = csv.reader(handle, delimiter="\t")
        header = None
        self.hugo_map = {}
        for row in reader:
            if header is None:
                header = {}
                for i, n in enumerate(row):
                    header[n] = i
            else:
                o = {}
                for col in header:
                    o[col] = row[header[col]]
                self.hugo_map[ row[header["Approved Symbol"]] ] = o
        handle.close()


    def suggest_nodes(self, gr):
        if self.hugo_map is not None:
            for node in gr.node:
                if 'type' in gr.node[node] and gr.node[node]['type'] == 'protein':
                    if node not in self.hugo_map:
                        log("Checking for %s" % (node))
                        found = False
                        for alt in self.hugo_map:
                            for src_namespace in self.hugo_map[alt]:
                                if str(self.hugo_map[alt][src_namespace]) == str(node):
                                    yield [node, src_namespace, alt]
                                    found = True
                        if not found:
                            for alt in self.hugo_map:
                                if len(alt) > 2:
                                    if node.startswith(alt):
                                        yield [node, 'name prefix', alt]


class RepoChecker(GraphChecker):
    def __init__(self, repo):
        self.repo = repo
        self.hugo_map = None
        self.dogma = None

    def load_dogma(self, path):
        handle = open(path)
        self.dogma = Dogma()
        self.dogma.read(handle)
        handle.close()


    def check_project(self, project):

        handle = open(os.path.join(self.repo, project, "INFO"))
        conf = yaml.load(handle.read())
        handle.close()

        errors = []
        if 'PID' not in conf:
            errors.append(Exception("Missing required PID field"))

        if 'DESC' not in conf:
            errors.append(Exception("Missing required DESC description field"))

        if 'SOURCE' not in conf:
            errors.append(Exception("Missing required SOURCE field"))

        handle = open(os.path.join(self.repo, project, "graph"))
        gr = network_convert.read_spf(handle)
        handle.close()
        if self.hugo_map is not None:
            for node in gr.node:
                if gr.node[node]['type'] == 'protein':
                    if node not in self.hugo_map:
                        errors.append(Exception("Gene Symbol not found:\t%s" % (node)))

        if self.dogma is not None:
            for node in gr.node:
                ntype = gr.node[node]['type']
                if self.dogma.get_dogma(ntype) == None:
                    errors.append(Exception("No dogma for node %s type %s" % (node, ntype)))

        return errors




def main_new(args):
    parser = argparse.ArgumentParser(prog="pathway_db new")
    parser.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    parser.add_argument("-m", "--manual-id", help="Manual ID", default=None)
    args = parser.parse_args(args)

    pid = None
    if pid is None:
        pid = 0
        for p in get_project_list(args.base_dir):
            pid = max(int(p.replace("PID", "")), pid)
        pid += 1
    else:
        pid = int(args.manual_id.replace("PID", ""))

    pid_name = "PID%s" % (pid)
    dstdir = os.path.join(args.base_dir, pid_name)
    os.mkdir(dstdir)

    config = {
        'PID' : pid,
        'DESC' : ''
    }

    store_config(args.base_dir, pid_name, config)
    log(pid_name)




def add_repo(basedir, pathway_name, conf):
    dstdir = os.path.join(basedir, pathway_name)
    if not os.path.exists(dstdir):
        raise Exception("Pathway PID%s doesn't exists" % (pathway_name))
    #print "Adding", pathway_name
    subprocess.check_call("cd %s; git add %s; git commit -m 'Adding Pathway %s'" % (basedir, pathway_name, pathway_name), shell=True)


def main_add(args):
    parser = argparse.ArgumentParser(prog="pathway_db add")
    parser.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    parser.add_argument("input", help="Input Directory")
    args = parser.parse_args(args)

    if args.base_dir is None or not os.path.exists(args.base_dir):
        sys.stderr.write("Define Pathway REPO\n")
        sys.exit(1)

    try:
        pid_name = None
        dstdir = os.path.join(args.base_dir, args.input)
        if os.path.exists(dstdir):
            pid_name = args.input
        else:
            if os.path.exists(args.input):
                handle = open(os.path.join(args.input, "INFO"))
                conf = yaml.load(handle.read())
                handle.close()
                pid_name = "PID%s" % (conf['PID'])
                dstdir = os.path.join(args.base_dir, pid_name)
                os.mkdir(dstdir)
                for f in ['INFO', 'graph']:
                    shutil.copy(os.path.join(args.input, f), dstdir)
            else:
                raise Exception("Directory not found")
        repocheck = RepoChecker(args.base_dir)
        conf = repocheck.check_project(pid_name)
        add_repo(args.base_dir, pid_name, conf)
    except Exception, e:
        sys.stderr.write("Pathway Check Error: %s : %s\n" % (args.input, str(e)))
        sys.exit(1)

def main_sync(args):
    parser = argparse.ArgumentParser(prog="pathway_db sync")
    parser.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    args = parser.parse_args(args)

    if args.base_dir is None:
        sys.stderr.write("Define Pathway REPO\n")
        sys.exit(1)

    if not os.path.exists(args.base_dir):
        subprocess.check_call("git clone %s %s" %
                (CENTRAL_REPO, args.base_dir),
            shell=True)

    subprocess.check_call("cd %s ; git pull origin" %
                (args.base_dir),
            shell=True)


def scan_dir(base_dir, suffix=".xgmml"):
    if os.path.isfile(base_dir):
        return [base_dir]
    out = []
    paths = glob(os.path.join(base_dir, "*"))
    for path in paths:
        if os.path.isdir(path):
            out += scan_dir(path)
        else:
            if path.endswith(suffix):
                out.append(path)
    return out

def pathway_opener(pathway_list):

    for path_type, path in pathway_list:

        if path_type == 'xgmml-dir':
            for path_path in scan_dir(path):
                yield XGMMLOpener(path_path)

        if path_type == 'xgmml-file':
            yield XGMMLOpener(path)

        if path_type == 'spf-file':
            yield SPFOpener(path)

        if path_type == 'biopax-dir':
            for path_path in scan_dir(path, ".owl"):
                yield BioPaxOpener(path_path)

        if path_type == "library":
            for path_path in glob(os.path.join(path, "*", "graph.xgmml")):
                yield XGMMLOpener(path_path)

        if path_type == "library-biopax":
            for path_path in glob(os.path.join(path, "*", "biopax.owl")):
                yield XGMMLOpener(path_path)




class BioPaxOpener:

    def __init__(self, path):
        self.path = path
        self.name = path

    def read(self):
        raise Exception("Can't read directly from BioPAX")


class XGMMLOpener:

    def __init__(self, path):
        self.path = path
        self.name = path

    def read(self):
        handle = open(self.path)
        cur_gr = network_convert.read_xgmml(handle)
        handle.close()
        return cur_gr


class SPFOpener:

    def __init__(self, path):
        self.path = path
        self.name = path

    def read(self):
        handle = open(self.path)
        cur_gr = network_convert.read_spf(handle, strict=False)
        handle.close()
        return cur_gr


class GraphBuilder:
    def __init__(self, args):
        self.args = args
        self.merge_map = {}
        if args.merge_file is not None:
            handle = open(args.merge_file)
            for line in handle:
                tmp = line.rstrip("\r\n").split("\t")
                for i in tmp:
                    self.merge_map[i] = tmp[0]
            handle.close()

        self.exclude = {}
        if args.exclude is not None:
            handle = open(args.exclude)
            for line in handle:
                self.exclude[line.rstrip()] = True
            handle.close()

        self.dogma = None
        if args.dogma is not None:
            self.dogma = Dogma()
            handle = open(args.dogma)
            self.dogma.read(handle)
            handle.close()

    def fix_graph(self, graph):
        gr = networkx.MultiDiGraph()

        node_remap = {}

        for node in graph.node:
            skip = False

            if 'type' not in graph.node[node]:
                log("Untyped node: %s %s" % (node, graph.node[node].get('url', '')))
                if self.args.all:
                    skip = True
            else:
                if self.args.exclude_type is not None and graph.node[node]['type'] in self.args.exclude_type:
                    skip = True
                    log("Remove node: %s of type %s" % (node, graph.node[node]['type']))

            if 'label' not in graph.node[node] or graph.node[node]['label'] == 'None':
                log("Unlabeled node: %s %s" % (node, graph.node[node].get('url', '')))
                if self.args.all:
                    skip = True

            if self.dogma is not None:
                if 'type' not in graph.node[node]:
                    log("Undefined node type: %s" % (node))
                    skip = True
                elif not self.dogma.has_nodetype(graph.node[node]['type']):
                    log("Unknown node type: %s : %s" % (node, graph.node[node]['type']))
                    skip = True


            if not skip:
                data = copy(graph.node[node])
                if self.args.rename_hugo or self.args.all:
                    if 'db_xref' in data:
                        for key in data['db_xref']:
                            if key.startswith("HGNC Symbol:"):
                                log("Changing %s to %s" % (data['label'], key))
                                data['label'] = key.split(":")[1]
                                new_node = data['label']
                                node_remap[node] = new_node
                                node = new_node

                if node not in gr.node:
                    if self.args.rename_type or self.args.all:
                        if data.get('type', '') == 'complex':
                            data['label'] += " (complex)"
                        if data.get('type', '') == 'family':
                            data['label'] += " (family)"

                    if self.args.rename_space or self.args.all:
                        data['label'] = data['label'].replace(" ", "_")

                    if self.args.rename_prime or self.args.all:
                        data['label'] = data['label'].replace("5'", "5prime")
                        data['label'] = data['label'].replace("3'", "3prime")

                    if self.args.rename_ascii or self.args.all:
                        data['label'] = re.sub("&delta;", "delta", data['label'], flags=re.IGNORECASE)
                        data['label'] = re.sub("&alpha;", "alpha", data['label'], flags=re.IGNORECASE)
                        data['label'] = re.sub("&beta;", "beta", data['label'], flags=re.IGNORECASE)
                        data['label'] = re.sub("&gamma;", "gamma", data['label'], flags=re.IGNORECASE)

                    if self.args.rename_char or self.args.all:
                        data['label'] = re.sub( r'[\'\"\:\?\\\*= \,]', "_", data['label'])


                    if 'label' in data and data['label'] in self.merge_map:
                        data['label'] = self.merge_map[data['label']]

                    if 'label' not in data or data['label'] not in self.exclude:
                        gr.add_node(node, attr_dict=data)
                else:
                    if 'type' in gr.node[node] and 'type' in data and gr.node[node]['type'] != data['type']:
                        error("Mismatch Node Type: %s :%s --> %s" % (node, gr.node[node]['type'], data['type'] ))

                        if self.args.rename_nonprotein or self.args.all:
                            #because 'protein' is a default node type, if we see something not protein, then change the node to match
                            if gr.node[node]['type'] == 'protein':
                                gr.node[node]['type'] = data['type']

        for src, dst, data in graph.edges(data=True):
            src_node_type = graph.node[src].get('type', None)
            dst_node_type = graph.node[dst].get('type', None)
            if src in node_remap:
                src = node_remap[src]
            if dst in node_remap:
                dst = node_remap[dst]
            interaction = data['interaction']

            if src in gr.node and dst in gr.node:
                add_edge = True
                if self.dogma is not None:
                    if not self.dogma.has_edgetype(src_node_type, interaction, dst_node_type):
                        error("BAD_EDGETYPE: %s(%s) %s %s(%s)" % (src, src_node_type, interaction, dst, dst_node_type))
                        add_edge = False

                has_edge = False
                if dst in gr.edge[src]:
                    for i in gr.edge[src][dst]:
                        if gr.edge[src][dst][i]['interaction'] == interaction:
                            has_edge = True

                if not has_edge:
                    if add_edge:
                        if not (self.args.remove_self or self.args.all) or src != dst:
                            gr.add_edge(src, dst, attr_dict=data )
                        else:
                            log("Removing self loop: %s" % (src))
        return gr

def main_build(args):
    gr = networkx.MultiDiGraph()

    paths = pathway_opener( list( (args.pathways[i],args.pathways[i+1]) for i in range(0, len(args.pathways),2) ) )

    builer = GraphBuilder(args)

    type_count = {}
    interaction_count = {}
    duplicate_edges = 0
    for cur_path in paths:
        log("Scanning: %s" % (cur_path.name))
        cur_gr = cur_path.read()
        if args.organism is None or cur_gr.graph.get('organism', args.organism) == args.organism:
            fix_gr = builer.fix_graph(cur_gr)

            for node in fix_gr.node:
                if node in gr.node:
                    if 'type' in gr.node[node] and 'type' in fix_gr.node[node] and gr.node[node]['type'] != fix_gr.node[node]['type']:
                        error("%s failure: Mismatch Node Type: %s :%s --> %s" % (cur_path.name, node, gr.node[node]['type'], fix_gr.node[node]['type'] ))
                        if args.rename_nonprotein or args.all:
                            #because 'protein' is a default node type, if we see something not protein, then change the node to match
                            if gr.node[node]['type'] == 'protein':
                                gr.node[node]['type'] = fix_gr.node[node]['type']
                else:
                    log("Merging: %s" % node)
                    gr.add_node(node, attr_dict=fix_gr.node[node])

            for src, dst, data in fix_gr.edges(data=True):
                interaction = data['interaction']
                src_node_type = fix_gr.node[src].get('type', None)
                dst_node_type = fix_gr.node[dst].get('type', None)

                if src in gr.node and dst in gr.node:
                    has_edge = False
                    if dst in gr.edge[src]:
                        for i in gr.edge[src][dst]:
                            if gr.edge[src][dst][i]['interaction'] == interaction:
                                has_edge = True

                    if not has_edge:
                        if not (args.remove_self or args.all) or src != dst:
                            gr.add_edge(src, dst, attr_dict=data )
                        else:
                            log("Removing self loop: %s" % (src))
                    else:
                        duplicate_edges += 1

    connect_list = networkx.connected_components(networkx.Graph(gr))
    rm_list = []
    for group in connect_list:
        if len(group) < args.min_subgraph:
            rm_list.extend(group)
    for r in rm_list:
        gr.remove_node(r)


    log("+---------------------------+")
    log("|Node Count: %15d|" % (len(gr.nodes())))
    log("|Edge Count: %15d|" % (len(gr.edges())))
    log("|Duplicate Edges: %10d|" % (duplicate_edges))
    log("|Connected Components: %5d|" % (networkx.number_connected_components(networkx.Graph(gr))))
    log("+---------------------------+")
    if args.output:
        handle = open(args.output, "w")
    else:
        handle = sys.stdout

    for n_type in type_count:
        log("Node Type %s: %d" % (n_type, type_count[n_type]))
    if args.spf:
        network_convert.write_spf(gr, handle)
    elif args.sif:
        network_convert.write_sif(gr, handle)
    else:
        network_convert.write_xgmml(gr, handle)
    handle.close()


def get_modified(base_dir):
    output = subprocess.check_output("cd %s; git status" % (base_dir),
        shell=True)

    pid_list = {}
    outs = StringIO(output)
    untracked = False
    for line in outs:
        if line.startswith("# Untracked files"):
            untracked = True
        if not untracked:
            tmp = line.rstrip().split("\t")
            if len(tmp) == 2 and tmp[0] == "#" and tmp[1].startswith('modified:'):
                fname = re.sub(r'^modified:\s*', '', tmp[1])
                pid_list[os.path.dirname(fname)] = 'modified'
        else:
            tmp = line.rstrip().split("\t")
            if len(tmp) == 2 and tmp[0] == '#' and tmp[1].startswith("PID"):
                pid_list[os.path.dirname(tmp[1])] = 'untracked'
    return pid_list

def get_project_list(base_dir):
    pid_list = []
    for p in glob(os.path.join(base_dir, "*")):
        name = os.path.basename(p)
        if name.startswith("PID"):
            pid_list.append(name)
    return pid_list

def load_config(base_dir, project):
    handle = open(os.path.join(base_dir, project, "INFO"))
    conf = yaml.load(handle.read())
    handle.close()
    return conf

def load_stored_config(base_dir, project):
    output = subprocess.check_output("cd %s; git show HEAD:%s/INFO" % (base_dir, project),
        shell=True)
    outs = StringIO(output)
    conf = yaml.load(outs.read())
    return conf

def store_config(base_dir, project, config):
    text = yaml.dump(config, default_flow_style=False)
    ohandle = open(os.path.join(base_dir, project, "INFO"), "w")
    ohandle.write(text)
    ohandle.close()

def main_status(args):
    parser = argparse.ArgumentParser(prog="pathway_db status")
    parser.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    args = parser.parse_args(args)

    pid_list = get_modified(args.base_dir)

    for pid in pid_list:
        print "%s:\t%s" % (pid_list[pid], pid)

def main_hugosync(args):
    src_url = "http://www.genenames.org/cgi-bin/hgnc_downloads?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_prev_name&col=gd_aliases&col=gd_name_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_enz_ids&col=gd_pub_eg_id&col=gd_pubmed_ids&col=gd_pub_refseq_ids&col=gd_ccds_ids&col=md_eg_id&col=md_mim_id&col=md_refseq_id&col=md_prot_id&col=md_ensembl_id&col=md_ucsc_id&col=md_mgd_id&col=md_rgd_id&status=Approved&status=Entry+Withdrawn&status_opt=2&where=&order_by=gd_hgnc_id&format=text&limit=&hgnc_dbtag=on&submit=submit"
    log("Downloading HUGO")
    ihandle = urlopen(src_url)
    if not os.path.exists(DATA_REPO):
        os.mkdir(DATA_REPO)
    ohandle = open(os.path.join(DATA_REPO, "hugo.tab"), "w")
    for line in ihandle:
        ohandle.write(line)
    ihandle.close()
    ohandle.close()


def main_check(args):
    parser = argparse.ArgumentParser(prog="pathway_db check")
    parser.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    parser.add_argument("project", help="Project List", nargs="*")

    args = parser.parse_args(args)

    projects = []
    if len(args.project) == 0:
        projects = get_project_list(args.base_dir)
    else:
        projects = args.project

    checker = RepoChecker(args.base_dir)
    hugo_path = os.path.join(DATA_REPO, "hugo.tab")
    if os.path.exists(hugo_path):
        checker.load_hugomap(hugo_path)
    else:
        log("Skipping HUGO check")

    dogma_path = os.path.join(DATA_REPO, "standard_dogma.yaml")
    if os.path.exists(dogma_path):
        checker.load_dogma(dogma_path)
    else:
        log("Skipping Dogma check")

    for project in projects:
        try:
            errors = checker.check_project(project)
            if len(errors) == 0:
                log("OK: %s" % (project))
            else:
                for e in errors:
                    log("CheckError:\t%s\t%s" % (project, str(e)))
                log("BAD: %s" % (project))
        except Exception, e:
            sys.stderr.write("Pathway Check Error: %s : %s\n" % (project, str(e)))

def main_suggest(args):

    checker = GraphChecker()
    hugo_path = os.path.join(DATA_REPO, "hugo.tab")
    if os.path.exists(hugo_path):
        checker.load_hugomap(hugo_path)
    else:
        log("Can't make suggestions until hugo synced")
        sys.exit(1)


    paths = pathway_opener( list( (args.pathways[i],args.pathways[i+1]) for i in range(0, len(args.pathways),2) ) )
    for cur_path in paths:
        gr = cur_path.read()

        #try:
        suggestions = checker.suggest_nodes(gr)
        for sug in suggestions:
            log("Instead of %s (from:%s) try %s" % (sug[0], sug[1], sug[2]))
        #except Exception, e:
        #    sys.stderr.write("Pathway Check Error: %s\n" % (str(e)))


def runner(x):
    x.run()

def main_biopax_format(args):

    rename = {}
    if args.rename is not None:
        handle = open(args.rename)
        for line in handle:
            tmp = line.rstrip().split("\t")
            rename[tmp[0]] = tmp[1]
        handle.close()
    paths = pathway_opener( list( (args.pathways[i],args.pathways[i+1]) for i in range(0, len(args.pathways),2) ) )
    tasks = []
    for p in paths:
        path = p.path
        name = os.path.basename(path)
        if name in rename:
            name = rename[name]
        outdir = os.path.join(args.outdir, name)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        biopax_path = os.path.join(outdir, "biopax.owl")
        log("Copying: %s" % (biopax_path))
        shutil.copy(path, biopax_path)


def main_biopax_convert(args):
    from read_biopax import ConvertTask
    from multiprocessing import Pool
    paths = pathway_opener( [ ['library-biopax', args.library] ] )
    tasks = []
    for p in paths:
        outdir = os.path.dirname(p.path)
        tasks.append( ConvertTask(
            p.path,
            singlepath=os.path.join(outdir, "graph.xgmml")
        ) )

    p = Pool(args.cpus)
    p.map(runner, tasks)

def main_library_gmt(args):

    filter_map = None
    if args.filter is not None:
        filter_map = {}
        handle = open(args.filter)
        for line in handle:
            filter_map[line.rstrip()] = True
        handle.close()


    for path_dir in glob(os.path.join(args.library, "*")):
        if os.path.isdir(path_dir):
            path_name = os.path.basename(path_dir)
            spf_path = os.path.join(path_dir, "graph.spf")
            handle = open(spf_path)
            gr = network_convert.read_spf(handle, strict=False)
            handle.close()
            if args.organism is None or gr.graph.get('organism', args.organism) == args.organism:
                if filter_map is None:
                    out = gr.node.keys()
                else:
                    out = []
                    for a in gr.node.keys():
                        if a in filter_map:
                            out.append(a)
                if not args.skip_empty or len(out):
                    print "%s\t%s" % (path_name, "\t".join(out))


def main_library_table(args):

    for path_dir in glob(os.path.join(args.library, "*")):
        if os.path.isdir(path_dir):
            path_name = os.path.basename(path_dir)
            xgmml_path = os.path.join(path_dir, "graph.xgmml")
            handle = open(xgmml_path)
            gr = network_convert.read_xgmml(handle)
            handle.close()

            organism = gr.graph.get('organism', '')
            name = gr.graph.get('displayName', None)
            if name is None:
                name = gr.graph.get('name', None)

            print "%s\t%s\t%s" % (path_name, organism, name)

def main_library_compile(args):
    builder = GraphBuilder(args)

    for path_dir in glob(os.path.join(args.library, "*")):
        if os.path.isdir(path_dir):
            path_name = os.path.basename(path_dir)
            xgmml_path = os.path.join(path_dir, "graph.xgmml")
            handle = open(xgmml_path)
            gr = network_convert.read_xgmml(handle)
            handle.close()
            fix_gr = builder.fix_graph(gr)
            spf_path = os.path.join(path_dir, "graph.spf")
            log("Writing: %s" % (spf_path))
            handle = open(spf_path, "w")
            network_convert.write_spf(fix_gr, handle)
            handle.close()

def main_library_copy(args):
    builder = GraphBuilder(args)
    out_zip = zipfile.ZipFile(args.out, "w")
    out_base = re.sub( r'.zip$', '', os.path.basename(args.out))
    for path_dir in glob(os.path.join(args.library, "*")):
        if os.path.isdir(path_dir):
            path_name = os.path.basename(path_dir)
            xgmml_path = os.path.join(path_dir, "graph.xgmml")
            handle = open(xgmml_path)
            gr = network_convert.read_xgmml(handle)
            handle.close()
            fix_gr = builder.fix_graph(gr)
            log("Writing: %s" % (path_name))
            text = StringIO()
            network_convert.write_spf(fix_gr, text)
            if len(text.getvalue()) > 0:
                out_zip.writestr( os.path.join(out_base, path_name + ".spf"), text.getvalue())
    out_zip.close()


"""
def main_commit(args):
    parser = argparse.ArgumentParser(prog="pathway_db commit")
    parser.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    parser.add_argument("-f", "--force", help="Force Commit", action="store_true", default=False)
    parser.add_argument("-m", "--message", help="Message", default=None)
    parser.add_argument("project", help="Project")
    args = parser.parse_args(args)

    checker = RepoChecker(args.base_dir)
    hugo_path = os.path.join(DATA_REPO, "hugo.tab")
    if os.path.exists(hugo_path):
        checker.load_hugomap(hugo_path)
    else:
        log("Skipping HUGO check")

    errors = checker.check_project(args.project)
    if len(errors):
        for e in errors:
            log(str(e))
        if not args.force:
            sys.exit(1)
        else:
            log("Forcing Commit")

    config = load_config(args.base_dir, args.project)
    stored_config = load_stored_config(args.base_dir, args.project)
    version = stored_config.get('VERSION', 0)

    config['VERSION'] = version + 1

    store_config(args.base_dir, args.project, config)
    message = args.message
    if message is None:
        message = "Updating pathway %s" % (args.project)

    subprocess.check_call("cd %s ; git commit -m '%s' %s" %
                (args.base_dir, message, args.project),
            shell=True)

"""

def add_build_args(parser):
    parser.add_argument("--merge-file", default=None)
    parser.add_argument("--exclude", default=None)
    parser.add_argument("--exclude-type", action="append")
    parser.add_argument("--dogma", default=None)
    parser.add_argument("-r", "--rename-hugo", help="Rename nodes to HUGO codes if possible", action="store_true", default=False)
    parser.add_argument("--rename-type", action="store_true", default=False)
    parser.add_argument("--rename-space", action="store_true", default=False)
    parser.add_argument("--rename-char", action="store_true", default=False)
    parser.add_argument("--rename-prime", action="store_true", default=False)
    parser.add_argument("--remove-self", action="store_true", default=False)
    parser.add_argument("--rename-ascii", action="store_true", default=False)
    parser.add_argument("--rename-nonprotein", action="store_true", default=False)
    parser.add_argument('-a', '--all', help="All processing options on", action="store_true", default=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    subparsers = parser.add_subparsers(title="subcommand")

    parser_build = subparsers.add_parser('build')
    parser_build.add_argument('-p', '--spf', help="Compile SimplePathwayFormat File", action="store_true", default=False)
    parser_build.add_argument('-s', '--sif', help="Compile SIF File", action="store_true", default=False)
    parser_build.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    parser_build.add_argument("-o", "--output", default=None)
    parser_build.add_argument("--organism", default=None)
    parser_build.add_argument("--min-subgraph", type=int, default=0)
    add_build_args(parser_build)
    parser_build.add_argument("pathways", nargs="*")

    parser_build.set_defaults(func=main_build)

    parser_hugosync = subparsers.add_parser('hugosync')
    parser_hugosync.set_defaults(func=main_hugosync)

    parser_suggest = subparsers.add_parser('suggest')
    parser_suggest.add_argument("pathways", nargs="*")
    parser_suggest.set_defaults(func=main_suggest)

    parser_biopax_format = subparsers.add_parser('biopax-format')
    parser_biopax_format.set_defaults(func=main_biopax_format)
    parser_biopax_format.add_argument("--rename", default=None)
    parser_biopax_format.add_argument("--cpus", type=int, default=2)
    parser_biopax_format.add_argument("outdir")
    parser_biopax_format.add_argument("pathways", nargs="*")

    parser_biopax_convert = subparsers.add_parser('biopax-convert')
    parser_biopax_convert.set_defaults(func=main_biopax_convert)
    parser_biopax_convert.add_argument("--cpus", type=int, default=2)
    parser_biopax_convert.add_argument("library")

    parser_library_gmt = subparsers.add_parser('library-gmt')
    parser_library_gmt.set_defaults(func=main_library_gmt)
    parser_library_gmt.add_argument("--organism", default=None)
    parser_library_gmt.add_argument("--filter", default=None)
    parser_library_gmt.add_argument("--skip-empty", action="store_true", default=False)
    parser_library_gmt.add_argument("library")

    parser_library_gmt = subparsers.add_parser('library-table')
    parser_library_gmt.set_defaults(func=main_library_table)
    parser_library_gmt.add_argument("library")


    parser_library_compile = subparsers.add_parser('library-compile')
    parser_library_compile.set_defaults(func=main_library_compile)
    add_build_args(parser_library_compile)
    parser_library_compile.add_argument("library")

    parser_library_copy = subparsers.add_parser('library-copy')
    parser_library_copy.set_defaults(func=main_library_copy)
    add_build_args(parser_library_copy)
    parser_library_copy.add_argument("library")
    parser_library_copy.add_argument("out")

    args = parser.parse_args()

    args.func(args)
