#!/usr/bin/env python

from argparse import ArgumentParser
import pathway_tools.convert
from networkx.algorithms import shortest_path

if __name__ == "__main__":
    parser = ArgumentParser(description="Find shortest directed path between two genes")    
    parser.add_argument("--xgmml", default=None)
    parser.add_argument("gene1")
    parser.add_argument("gene2")

    args = parser.parse_args()

    handle = open(args.xgmml)
    net = pathway_tools.convert.read_xgmml(handle)

    out = []
    path = shortest_path(net, args.gene1, args.gene2)
    for i in range(len(path)-1):
        emap = net[path[i]][path[i+1]]
        ea = []
        for e in emap:
            ea.append(emap[e]['interaction'])
        estr = ",".join(ea)
        out.append( "%s\t%s\t" % (path[i], estr) )
    print "%s%s" % ("".join(out), path[-1])