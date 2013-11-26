#!/usr/bin/env python

from pathway_tools import convert
import sys
import argparse


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--query-spf', help="Input query in simple pathway format", default=None)
    parser.add_argument('--query-xgmml',  help="Input query in XGMML graph format", default=None)
    parser.add_argument('--target-sfp', help="Search Target in Simple Pathway Format", default=None)
    parser.add_argument('--target-xgmml', help="Search Target in XGMML graph format", default=None)
    parser.add_argument('--query-pathattr', default="pid")
    parser.add_argument('--target-pathattr', default="pid")
    parser.add_argument("--query-nodetype", default=None)
    parser.add_argument("--target-nodetype", default=None)
    parser.add_argument("--min", "-m", type=int, default=1)
    parser.add_argument("--top", type=int, default=None)

    args = parser.parse_args()

    if args.query_spf:
        handle = open(args.query_spf) 
        query_net = convert.read_spf(handle)
    
    if args.query_xgmml:
        handle = open(args.query_xgmml) 
        query_net = convert.read_xgmml(handle)
    

    if args.target_spf:
        handle = open(args.target_spf) 
        target_net = convert.read_spf(handle)
        
    if args.target_xgmml:
        handle = open(args.target_xgmml) 
        target_net = convert.read_xgmml(handle)
    
    target_pathways = {}
    for n in target_net.node:
        if args.target_nodetype is None or target_net.node[n]['type'] == args.target_nodetype:
            for pid in target_net.node[n][args.target_pathattr]:
                target_pathways[pid] = target_pathways.get(pid, 0) + 1

    query_pathways = {}
    for n in query_net.node:
        if args.query_nodetype is None or query_net.node[n]['type'] == args.query_nodetype:
            for pid in query_net.node[n][args.query_pathattr]:
                query_pathways[pid] = query_pathways.get(pid, 0) + 1

    
    hits = {}
    for n in query_net.node:
        if n in target_net.node:
            for q_pid in query_net.node[n][args.query_pathattr]:
                if q_pid not in hits:
                    hits[q_pid] = {}
                for t_pid in target_net.node[n][args.target_pathattr]:
                    if t_pid not in hits[q_pid]:
                        hits[q_pid][t_pid] = []
                    hits[q_pid][t_pid].append(n)

    for q in hits:
        hit_list = hits[q].keys()
        hit_list.sort(key=lambda x:len(hits[q][x]), reverse=True)
        if args.top is not None:
            hit_list = hit_list[:args.top]
        for t in hit_list:
            if len(hits[q][t]) >= args.min:
                print "%s = %s : %d : ( %f %f ) : (%s)" % (
                        q, t, len(hits[q][t]),
                        float(len(hits[q][t])) / float(query_pathways[q]), 
                        float(len(hits[q][t])) / float(target_pathways[t]), 
                        ",".join(hits[q][t])
                    )
