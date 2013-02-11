#!/usr/bin/env python

import argparse
import sys
import os
import re
import yaml
import csv
import shutil
import subprocess
import networkx
import network_convert
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


class RepoChecker:
    def __init__(self, repo):
        self.repo = repo
        self.hugo_map = None

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
                        errors.append(Exception("Gene Symbol not found: %s" % (node)))

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
    print "Adding", pathway_name
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
    parser = argparse.ArgumentParser(prog="pathway_db compile")
    parser.add_argument('-a', '--append', help="Default Edge Type", action="append")
    parser.add_argument('-p', '--paradigm', help="Compile Paradigm File", action="store_true", default=False)   
    parser.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    
    args = parser.parse_args(args)

    gr = networkx.MultiDiGraph()
    base_dir = args.base_dir

    paths = glob(os.path.join(base_dir, "[A-Z]*"))
    if args.append:
        paths +=args.append

    type_count = {}
    interaction_count = {}
    duplicate_edges = 0
    for path in paths:
        info_path = os.path.join(path, "INFO")
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
                        raise Exception("Mismatch Node Type: %s :%s --> %s" % (node_name, gr.node[node_name]['type'], node_type ))
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
    log("Node Count: %d" % (len(gr.nodes())))
    log("Edge Count: %d" % (len(gr.edges())))
    log("Duplicate Edges: %s" % (duplicate_edges))
    for n_type in type_count:
        log("Node Type %s: %d" % (n_type, type_count[n_type]))
    if args.paradigm:
        network_convert.write_paradigm_graph(gr, sys.stdout)
    else:        
        network_convert.write_xgmml(gr, sys.stdout)

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
    for project in projects:
        try:
            errors = checker.check_project(project)
            if len(errors) == 0:
                print "OK: %s" % (project)
            else:
                for e in errors:
                    print "CheckError: %s : %s" % (project, str(e))
                print "BAD: %s" % (project)
        except Exception, e:
            sys.stderr.write("Pathway Check Error: %s : %s\n" % (project, str(e)))


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
    'commit' : {
        'method' : main_commit
    },
    'new' : {
        'method' : main_new
    },
    


}

if __name__ == "__main__":
    mode = None
    if len(sys.argv) > 1:
        mode = sys.argv[1]

    if mode in mode_map:
        mode_map[mode]['method'](sys.argv[2:])
    else:
        print "Modes:"
        for m in mode_map:
            print "\t%s" % (m)

