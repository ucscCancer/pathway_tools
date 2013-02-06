#!/usr/bin/env python

import argparse
import sys
import os
import re
import yaml
import shutil
import subprocess
import networkx
import network_convert
from glob import glob

def log(message):
    sys.stderr.write(message + "\n")


reKeyValue = re.compile(r'^(\S+)\s*=\s*(.*)$')


CENTRAL_REPO = "hgwdev.soe.ucsc.edu:/hive/groups/cancerGB/paradigm/superpathway_db.git"
LOCAL_REPO = "superpathway_db"

def check_repo(basedir):
    handle = open(os.path.join(basedir, "INFO"))
    conf = yaml.load(handle.read())
    handle.close()
    
    if 'PID' not in conf:
        raise Exception("Missing required PID field")

    if 'DESC' not in conf:
        raise Exception("Missing required DESC description field")
        
    if 'SOURCE' not in conf:
        raise Exception("Missing required SOURCE field")
    
    return conf

def add_repo(input, conf, basedir):
    pathway_name = "PID%s" % (conf['PID'])
    dstdir = os.path.join(basedir, pathway_name)
    if os.path.exists(dstdir):
        raise Exception("Pathway PID%s already exists" % (conf['PID']))
    print "Adding", dstdir
    os.mkdir(dstdir)
    for f in ['INFO', 'graph']:
        shutil.copy(os.path.join(input, f), dstdir)
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
        conf = check_repo(args.input)
        add_repo(args.input, conf, args.base_dir)
    except Exception, e:
        sys.stderr.write("Pathway Check Error: %s : %s\n" % (args.input, str(e)))
        sys.exit(1)

def main_sync(args):
    parser = argparse.ArgumentParser(prog="pathway_db add")
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
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--append', help="Default Edge Type", action="append")
    parser.add_argument('-p', '--paradigm', help="Compile Paradigm File", action="store_true", default=False)   
    parser.add_argument("-b", "--base-dir", help="BaseDir", default=LOCAL_REPO)
    
    args = parser.parse_args(args)

    gr = networkx.DiGraph()
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
                gr.node[node_name]['pathway'].append(info['DESC'])
            elif len(tmp) == 3:
                if tmp[0] not in gr.node:
                    raise Exception("Missing Node Declaration: %s" % (tmp[0]))
                if tmp[1] not in gr.node:
                    raise Exception("Missing Node Declaration: %s" % (tmp[1]))
                if tmp[1] not in gr.edge[tmp[0]]:
                    gr.add_edge(tmp[0], tmp[1], interaction=tmp[2])
                else:
                    #log("Duplicate Edge: %s %s (%s)" % (tmp[0], tmp[1], tmp[2]))
                    if gr.edge[tmp[0]][tmp[1]]['interaction'] != tmp[2]:
                        log("Mismatch Edge Declaration: %s %s : %s --> %s" % (tmp[0], tmp[1], gr.edge[tmp[0]][tmp[1]]['interaction'], tmp[2] ))
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


if __name__ == "__main__":
    mode = sys.argv[1]

    if mode == 'add':
        main_add(sys.argv[2:])

    if mode == 'sync':
        main_sync(sys.argv[2:])

    if mode == 'compile':
        main_compile(sys.argv[2:])


