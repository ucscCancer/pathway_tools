#!/usr/bin/env python

import argparse
import rdflib
import sys
import re
import networkx
import pathway_tools.convert

def ns_convert(ns_map, url):
    for ns in ns_map:
        if url.startswith(ns):
            return re.sub("^" + ns, ns_map[ns] + ":", url)
    return url

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", nargs=2, dest="namespace", action="append")
    parser.add_argument("input")
    args = parser.parse_args()

    if args.input == "-":
        handle = sys.stdin
    else:
        handle = open(args.input)
    g = rdflib.Graph()
    result = g.parse(handle)
    
    gr = networkx.MultiDiGraph()

    ns_map = {}
    if args.namespace:
        for ns, addr in args.namespace:
            ns_map[addr] = ns

    for subj, pred, obj in g:
        subj_str = ns_convert(ns_map, str(subj))
        pred_str = ns_convert(ns_map, str(pred))
        if subj_str not in gr.node:
            gr.add_node(subj_str)
        if isinstance(obj, rdflib.Literal):
            gr.node[subj_str][pred_str] = str(obj)
        elif str(pred) == "http://www.w3.org/1999/02/22-rdf-syntax-ns#type":
            gr.node[subj_str][pred_str] = ns_convert(ns_map, str(obj))
        else:
            gr.add_edge( subj_str, ns_convert(ns_map, str(obj)), attr_dict={"predicate" : pred_str} )
        #print subj, pred, obj
    
    pathway_tools.convert.write_xgmml(gr, sys.stdout)