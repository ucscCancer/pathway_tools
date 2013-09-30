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

def log(message):
    sys.stderr.write(message + "\n")


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

    def get_dogma(self, name):
        return self.data['dogma'].get(name, None)

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
        gr = network_convert.read_paradigm_graph(handle)
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
    print pid_name




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


def main_compile(args):

    gr = networkx.MultiDiGraph()
    base_dir = args.base_dir

    if len(args.pathways):
        paths = []
        for p in args.pathways:
            paths.append(os.path.join(base_dir, p))
    else:
        paths = glob(os.path.join(base_dir, "[A-Za-z]*"))
        if args.append:
            paths +=args.append

    type_count = {}
    interaction_count = {}
    duplicate_edges = 0
    for path in paths:
        info_path = os.path.join(path, "INFO")
        if os.path.exists(info_path):
            handle = open(info_path)
            info = yaml.load(handle.read())
            handle.close()
            handle = open(os.path.join(path, "graph"))
            for line in handle:
                tmp = line.rstrip().split("\t")
                node_type = tmp[0]
                node_name = tmp[1]
                if len(tmp) == 2:
                    if node_name not in gr.node:
                        gr.add_node( tmp[1], type=node_type )
                        type_count[node_type] = type_count.get(node_type, 0) + 1
                    else:
                        if gr.node[node_name]['type'] != node_type:
                            raise Exception("%s failure: Mismatch Node Type: %s :%s --> %s" % (path, node_name, gr.node[node_name]['type'], node_type ))
                    if 'pathway' not in gr.node[node_name]:
                        gr.node[node_name]['pathway'] = []
                        gr.node[node_name]['pid'] = []
                    gr.node[node_name]['pathway'].append(info['DESC'])
                    #gr.node[node_name]['pid'].append( "PID%s" % (info['PID']))
                elif len(tmp) == 3:
                    if tmp[0] not in gr.node:
                        raise Exception("Missing Node Declaration: %s" % (tmp[0]))
                    if tmp[1] not in gr.node:
                        raise Exception("Missing Node Declaration: %s" % (tmp[1]))
                    has_edge = False
                    if tmp[1] in gr.edge[tmp[0]]:
                        for i in gr.edge[tmp[0]][tmp[1]]:
                            if gr.edge[tmp[0]][tmp[1]][i]['interaction'] == tmp[2]:
                                has_edge = True
                    if not has_edge:
                        gr.add_edge(tmp[0], tmp[1], attr_dict={ 'interaction' : tmp[2] })
                    else:
                        duplicate_edges += 1
            handle.close()
    log("Node Count: %d" % (len(gr.nodes())))
    log("Edge Count: %d" % (len(gr.edges())))
    log("Duplicate Edges: %s" % (duplicate_edges))
    log("Connected Components: %d" % (networkx.number_connected_components(networkx.Graph(gr))))
    for n_type in type_count:
        log("Node Type %s: %d" % (n_type, type_count[n_type]))
    if args.paradigm:
        network_convert.write_paradigm_graph(gr, sys.stdout)
    elif args.sif:
        network_convert.write_sif(gr, sys.stdout)        
    else:        
        network_convert.write_xgmml(gr, sys.stdout)


def scan_dir(base_dir):
    if os.path.isfile(base_dir):
        return [base_dir]
    out = []
    paths = glob(os.path.join(base_dir, "*"))
    for path in paths:
        if os.path.isdir(path):
            out += scan_dir(path)
        else:
            if path.endswith(".xgmml"):
                out.append(path)
    return out

def main_build(args):
    gr = networkx.MultiDiGraph()

    paths = []
    for base_dir in args.pathways:
        paths += scan_dir(base_dir)

    type_count = {}
    interaction_count = {}
    duplicate_edges = 0
    for path in paths:
        log("Scanning: %s" % (path))
        handle = open(path)
        cur_gr = network_convert.read_xgmml(handle)
        handle.close()

        for node in cur_gr.node:
            if 'type' not in cur_gr.node[node]:
                log("Untyped node: %s" % (node))
            if node not in gr.node:
                data = copy(cur_gr.node[node])
                if args.rename:
                    if 'db_xref' in data:
                        for key in data['db_xref']:
                            if key.startswith("HGNC Symbol:"):
                                log("Changing %s to %s" % (data['label'], key))
                                data['label'] = key.split(":")[1]
                gr.add_node(node, attr_dict=data)
            else:
                if 'type' in gr.node[node] and gr.node[node]['type'] != cur_gr.node[node]['type']:
                    log("%s failure: Mismatch Node Type: %s :%s --> %s" % (path, node, gr.node[node]['type'], cur_gr.node[node]['type'] ))

        for src, dst, data in cur_gr.edges(data=True):

            interaction = data['interaction']
            has_edge = False
            if dst in gr.edge[src]:
                for i in gr.edge[src][dst]:
                    if gr.edge[src][dst][i]['interaction'] == interaction:
                        has_edge = True

            if not has_edge:
                gr.add_edge(src, dst, attr_dict=data )
            else:
                duplicate_edges += 1

            """

            handle = open(os.path.join(path, "graph"))
            for line in handle:
                tmp = line.rstrip().split("\t")
                node_type = tmp[0]
                node_name = tmp[1]
                if len(tmp) == 2:
                    if node_name not in gr.node:
                        gr.add_node( tmp[1], type=node_type )
                        type_count[node_type] = type_count.get(node_type, 0) + 1
                    else:
                        if gr.node[node_name]['type'] != node_type:
                            raise Exception("%s failure: Mismatch Node Type: %s :%s --> %s" % (path, node_name, gr.node[node_name]['type'], node_type ))
                    if 'pathway' not in gr.node[node_name]:
                        gr.node[node_name]['pathway'] = []
                        gr.node[node_name]['pid'] = []
                    gr.node[node_name]['pathway'].append(info['DESC'])
                    gr.node[node_name]['pid'].append( "PID%s" % (info['PID']))
                elif len(tmp) == 3:
                    if tmp[0] not in gr.node:
                        raise Exception("Missing Node Declaration: %s" % (tmp[0]))
                    if tmp[1] not in gr.node:
                        raise Exception("Missing Node Declaration: %s" % (tmp[1]))
                    has_edge = False
                    if tmp[1] in gr.edge[tmp[0]]:
                        for i in gr.edge[tmp[0]][tmp[1]]:
                            if gr.edge[tmp[0]][tmp[1]][i]['interaction'] == tmp[2]:
                                has_edge = True
                    if not has_edge:
                        gr.add_edge(tmp[0], tmp[1], attr_dict={ 'interaction' : tmp[2] })
                    else:
                        duplicate_edges += 1
            handle.close()
            """

    log("Node Count: %d" % (len(gr.nodes())))
    log("Edge Count: %d" % (len(gr.edges())))
    log("Duplicate Edges: %s" % (duplicate_edges))
    log("Connected Components: %d" % (networkx.number_connected_components(networkx.Graph(gr))))
    if args.output:
        handle = open(args.output, "w")
    else:
        handle = sys.stdout

    for n_type in type_count:
        log("Node Type %s: %d" % (n_type, type_count[n_type]))
    if args.paradigm:
        network_convert.write_paradigm_graph(gr, handle)
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
                print "OK: %s" % (project)
            else:
                for e in errors:
                    print "CheckError:\t%s\t%s" % (project, str(e))
                print "BAD: %s" % (project)
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

    handle = open(args.graph)
    gr = network_convert.read_xgmml(handle)
    handle.close()

    #try:
    suggestions = checker.suggest_nodes(gr)
    for sug in suggestions:
        print "Instead of %s (from:%s) try %s" % (sug[0], sug[1], sug[2])
    #except Exception, e:
    #    sys.stderr.write("Pathway Check Error: %s\n" % (str(e)))



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

mode_map = {
    'add' : {
        'method' : main_add
    },
    'sync' : {
        'method' : main_sync
    },
    'compile' : {
        'method' : main_compile
    },
    'status' : {
        'method' : main_status
    },
    'add' : {
        'method' : main_add
    },
    'hugo_sync' : {
        'method' : main_hugosync
    },
    'check' : {
        'method' : main_check
    },
    'suggest' : {
        'method' : main_suggest
    },
    'commit' : {
        'method' : main_commit
    },
    'new' : {
        'method' : main_new
    },
    


}

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    subparsers = parser.add_subparsers(title="subcommand")

    parser_compile = subparsers.add_parser('compile')
    parser_compile.add_argument('-a', '--append', help="Default Edge Type", action="append")
    parser_compile.add_argument('-p', '--paradigm', help="Compile Paradigm File", action="store_true", default=False)   
    parser_compile.add_argument('-s', '--sif', help="Compile SIF File", action="store_true", default=False)   
    parser_compile.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    parser_compile.add_argument("pathways", nargs="*")
    parser_compile.set_defaults(func=main_compile)

    parser_build = subparsers.add_parser('build')
    parser_build.add_argument('-p', '--paradigm', help="Compile Paradigm File", action="store_true", default=False)   
    parser_build.add_argument('-s', '--sif', help="Compile SIF File", action="store_true", default=False)   
    parser_build.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    parser_build.add_argument("-r", "--rename", help="Rename nodes to HUGO codes if possible", action="store_true", default=False)
    parser_build.add_argument("-o", "--output", default=None)
    
    parser_build.add_argument("pathways", nargs="*")
    parser_build.set_defaults(func=main_build)

    parser_hugosync = subparsers.add_parser('hugosync')
    parser_hugosync.set_defaults(func=main_hugosync)


    parser_suggest = subparsers.add_parser('suggest')
    parser_suggest.add_argument("graph")
    parser_suggest.set_defaults(func=main_suggest)

    args = parser.parse_args()

    args.func(args)

