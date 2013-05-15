#!/usr/bin/env python 

import sys
import argparse
from pathway_tools.biopax import BioPaxFile, BioPaxSparql
from pathway_tools import convert

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', "--output", default=None)
    parser.add_argument('-p', '--pathways', action="append")
    parser.add_argument("-s", "--sparql", action="store_true")
    parser.add_argument("input")
    args = parser.parse_args()

    if args.sparql:
        b = BioPaxSparql(args.input)
    else:
        b = BioPaxFile()
        b.load(args.input)

    paths = args.pathways
    
    for gr_set in b.toNet(paths):
        for gr in gr_set:
            if args.output:
                handle = open(args.output, "w")
                convert.write_xgmml(gr, handle)
            else:
                convert.write_xgmml(gr, sys.stdout)
        