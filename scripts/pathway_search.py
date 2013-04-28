#!/usr/bin/env python

from pathway_tools import convert
import sys
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--query-paradigm', default=None)
    parser.add_argument('--query-xgmml', default=None)
    parser.add_argument('--base-paradigm', default=None)
    parser.add_argument('--base-xgmml', default=None)
    parser.add_argument("--min", "-m", type=int, default=1)

    args = parser.parse_args()

    if args.query_paradigm:
        handle = open(args.query_paradigm) 
        query_net = convert.read_paradigm_graph(handle)
    
    if args.query_xgmml:
        handle = open(args.query_xgmml) 
        query_net = convert.read_xgmml(handle)
    

    if args.base_paradigm:
        handle = open(args.base_paradigm) 
        base_net = convert.read_paradigm_graph(handle)
        
    if args.base_xgmml:
        handle = open(args.base_xgmml) 
        base_net = convert.read_xgmml(handle)
    
    pathways = {}
    for n in base_net.node:
        for pid in base_net.node[n]['pid']:
            pathways[pid] = True
    
    hits = {}
    for n in query_net.node:
        if n in base_net.node:
            for pid in base_net.node[n]['pid']:
                if pid not in hits:
                    hits[pid] = []
                hits[pid].append(n)

    k = hits.keys()
    k.sort( key=lambda x:len(hits[x]), reverse=True )
    for i in k[:10]:
        print "%s = %s" % (i, ",".join(hits[i]))
