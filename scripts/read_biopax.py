#!/usr/bin/env python 

import sys
import os
import argparse
from pathway_tools.biopax import BioPaxFile, BioPaxSparql
from pathway_tools import convert

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
        for a in b.pathways():
            print a
    else:
        paths = args.pathways        
        for gr_set in b.toNet(paths):
            for gr in gr_set:
                if args.out_dir:
                    handle = open(os.path.join(args.out_dir, gr.graph['name'] + ".xgmml"), "w")
                else:
                    if args.output:
                        handle = open(args.output, "w")
                    else:
                        handle = sys.stdout

                
                if args.paradigm:
                    convert.write_paradigm_graph(gr, handle)
                else:
                    convert.write_xgmml(gr, handle)
                
                handle.close()