#!/usr/bin/env python


import sys
import os
import yaml
import network_convert
import networkx
from glob import glob
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--append', help="Default Edge Type", action="append")
    parser.add_argument('-p', '--paradigm', help="Compile Paradigm File", action="store_true", default=False)
    
    parser.add_argument('basedir', help="Base Dir")
    
    args = parser.parse_args()

    gr = networkx.DiGraph()
    base_dir = args.basedir

    paths = glob(os.path.join(base_dir, "[A-Z]*"))
    if args.append:
        paths +=args.append

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
                else:
                    if gr.node[node_name]['type'] != node_type:
                        raise Exception("Mismatch Node Type: %s :%s %s" % (node_name, gr.node[node_name]['type'], node_type ))
                if 'pathway' not in gr.node[node_name]:
                    gr.node[node_name]['pathway'] = []
                gr.node[node_name]['pathway'].append(info['DESC'])
            elif len(tmp) == 3:
                gr.add_edge(tmp[0], tmp[1], interaction=tmp[2])
        handle.close()
    if args.paradigm:
        network_convert.write_paradigm_graph(gr, sys.stdout)
    else:        
        network_convert.write_xgmml(gr, sys.stdout)
