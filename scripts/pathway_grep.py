#!/usr/bin/env python

import network_convert
import sys
import argparse



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', help="Node Search", action="store_true")
    parser.add_argument('-xgmml', help="XGMML File", default=None)
    parser.add_argument('-property', help="Property Field to report", action="append", default=[])
    
    parser.add_argument("terms", help="Search terms", action="append")
    
    args = parser.parse_args()

    gr = None
    if args.xgmml is not None:
        handle = open(args.xgmml)
        gr = network_convert.read_xgmml(handle)
        handle.close()


    if len(args.property) == 0:
        args.property.append("pathway")
    if gr is not None:
        for query in args.terms:
            if query in gr.node:
                o = []
                for qprop in args.property:
                    if qprop in gr.node[query]:
                        prop = gr.node[query][qprop]
                        o.append(str(prop))
                print "%s: %s" % (query, "\t".join(o))

