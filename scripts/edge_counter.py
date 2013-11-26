#!/usr/bin/env python


from pathway_tools import convert
import argparse
import re

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in-xgmml")
    parser.add_argument("--out")
    parser.add_argument("--filter", default=None)

    args = parser.parse_args()

    handle = open(args.in_xgmml)
    gr = convert.read_xgmml(handle)
    handle.close()

    out = open(args.out, "w")

    filter = None
    if args.filter is not None:
        filter = []
        handle = open(args.filter)
        for line in handle:
            tmp = line.rstrip("\r\n").split("\t")
            if len(tmp) == 2:
                filter.append(tmp)
                print "Appending", tmp
        handle.close()

    for node in gr.node:
        count = 0
        for dest in gr.edge[node]:
            for edge in gr.edge[node][dest]:
                add = True
                if filter is not None:
                    for field, value in filter:
                        if field not in gr.edge[node][dest][edge] or not re.match( value, gr.edge[node][dest][edge][field] ):
                            add = False
                if add:
                    count += 1
        out.write("%s\t%d\n" % (node, count))
    out.close()