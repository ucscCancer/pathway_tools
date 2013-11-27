#!/usr/bin/env python

from pathway_tools import convert as network_convert
import sys
import argparse
import os
import re
from glob import glob

def scan_dir(base_dir, name_filter):
    if os.path.isfile(base_dir):
        return [base_dir]
    out = []
    paths = glob(os.path.join(base_dir, "*"))
    for path in paths:
        if os.path.isdir(path):
            out += scan_dir(path)
        else:
            if name_filter is None or re.match(name_filter, os.path.basename(path)):
                out.append(path)
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--xgmml', help="XGMML File", action="store_true", default=False)
    parser.add_argument('--spf', help="SimplePathwayFormat File", action="store_true", default=None)
    parser.add_argument('-p', '--property', help="Property Field to report", default=None)
    parser.add_argument('-f', help="Term file")
    parser.add_argument('-r', action="store_true", help="Recursive Search", default=False)
    parser.add_argument('--filter', help="File name filter", default=None)
    parser.add_argument('-n', action="store_true", help="Report File Name", default=False)    
    parser.add_argument("term", help="Search term")
    parser.add_argument("files", nargs="*", help="Files to search")
    
    args = parser.parse_args()
    
    if len(args.property) == 0:
        args.property.append("pathway")

    start_paths = []
    terms = {}
    if args.f is not None:
        start_paths = [args.term] + args.files
        handle = open(args.f)
        for line in handle:
            terms[line.rstrip("\r\n")] = True
        handle.close()        
    else:
        start_paths = args.files
        terms[ args.term ] = True

    paths = None
    if args.r:
        paths = []
        for a in start_paths:
            paths += list(scan_dir(a, args.filter))
    else:
        paths = start_paths

    

    for file_path in paths:    
        gr = None

        if args.spf:
            handle = open(file_path)
            gr = network_convert.read_spf(handle, strict=False)
            handle.close()
        else:
            handle = open(file_path)
            gr = network_convert.read_xgmml(handle)
            handle.close()


        if gr is not None:
            found = []
            if args.property is None:
                for t in terms:
                    if t in gr.node:
                        found = [ t ]
            else:
                for node in gr.node:
                    if args.property in gr.node[node] and gr.node[node][args.property] in terms:
                        found.append(node)

            if len(found):
                if args.n:
                    print file_path                
                else:
                    print "%s" % ("\t".join(found))

