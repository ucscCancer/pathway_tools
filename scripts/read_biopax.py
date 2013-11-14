#!/usr/bin/env python 

import sys
import os
import re
import argparse
import networkx
from pathway_tools.biopax import BioPaxFile, BioPaxSparql
from pathway_tools import convert

re_clean = re.compile(r'[/\\ \n;]')
re_namesplit = re.compile(r'[/ \#\&]')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', "--output", default=None)
    parser.add_argument('-p', '--pathways', action="append")
    parser.add_argument("-s", "--sparql", action="store_true", default=False)
    parser.add_argument("-l", "--list", action="store_true", default=False)
    parser.add_argument("-d", "--out-dir", default=None)
    
    parser.add_argument("--paradigm", action="store_true", default=False)
    
    
    parser.add_argument("input")
    args = parser.parse_args()

    if args.sparql:
        b = BioPaxSparql(args.input)
    else:
        b = BioPaxFile()
        b.load(args.input)

    if args.list:
        for a,b in b.pathways().iteritems():
            print a, b.encode('ascii', errors='ignore')
    else:
        
        paths = []
        for key,value in b.pathways().iteritems():
            add = False
            if args.pathways is None or value in args.pathways:
                add = True
            elif key in args.pathways:
                add = True
            else:
                a_name = re_namesplit.split(key)[-1]
                if a_name in args.pathways:
                    add = True
            if add:
                paths.append(key)

        for subnet in b.toNet(paths):
            
            gr = networkx.MultiDiGraph()
            gr.graph['name'] = subnet.meta['name']
            gr.graph['url'] = subnet.meta['url']
            subnet.to_graph(gr)

            name = re_namesplit.split(gr.graph['url'])[-1]
            if args.out_dir:
                handle = open(os.path.join(args.out_dir, name + ".xgmml"), "w")
            else:
                if args.output:
                    if "%s" not in args.output:
                        raise Exception("Output path string not formatted correctly (use %s))")
                    handle = open(args.output % (name), "w")
                else:
                    handle = sys.stdout

            
            if args.paradigm:
                convert.write_paradigm_graph(gr, handle)
            else:
                convert.write_xgmml(gr, handle)
            
            if handle != sys.stdout:
                handle.close()