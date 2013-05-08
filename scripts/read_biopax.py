#!/usr/bin/env python 

import sys
import argparse
from pathway_tools.biopax import BioPax
from pathway_tools import convert

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', "--output", default=None)
    parser.add_argument("input")
    args = parser.parse_args()

    b = BioPax()
    b.load(args.input)
    
    for gr in b.toNet():
        #convert.write_paradigm_graph(gr, sys.stdout)
        if args.output:
            handle = open(args.output, "w")
            convert.write_xgmml(gr, handle)
        else:
            convert.write_xgmml(gr, sys.stdout)
