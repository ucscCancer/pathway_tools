#!/usr/bin/env python

from pathway_tools import convert
import sys
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--query-paradigm', default=None)
    parser.add_argument('--base-paradigm', default=None)
    parser.add_argument('--base-xgmml', default=None)
    parser.add_argument("--min", "-m", type=int, default=1)

    args = parser.parse_args()

    handle = open(args.query_paradigm) 
    query_net = convert.read_paradigm_graph(handle)
    
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
    for pid in pathways:
        found = 0
        for n in query_net.node:
            if n in base_net.node:
                if pid in base_net.node[n]['pid']:
                    found += 1
        if found >= args.min:
            hits[pid] = found
    
    k = hits.keys()
    k.sort( key=lambda x:hits[x] )
    for i in k[:10]:
        print i, hits[i]
