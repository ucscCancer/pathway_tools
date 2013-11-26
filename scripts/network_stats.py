#!/usr/bin/env python

import argparse
import networkx
from pathway_tools import convert



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--spf', help="SimplePathwayFormat File")
    parser.add_argument('--xgmml', help="XGMML File")
    parser.add_argument('--top', help="Top Count", type=int, default=10)


    args = parser.parse_args()

    gr = None

    if args.spf:
        handle = open(args.spf)
        gr = convert.read_spf(handle, strict=False)
        handle.close()
    if args.xgmml:
        handle = open(args.xgmml)
        gr = convert.read_xgmml(handle)
        handle.close()


    type_count = {}
    for node in gr.node:
        if 'type' in gr.node[node]:
            node_type = gr.node[node]['type']
            type_count[node_type] = type_count.get(node_type, 0) + 1


    print("Node Count: %d" % (len(gr.nodes())))
    print("Edge Count: %d" % (len(gr.edges())))
    print("Connected Components: %d" % (networkx.number_connected_components(networkx.Graph(gr))))



    print("Top 10 Input Connections:")
    in_degree = gr.in_degree()
    print("\n".join( "\t%s\t%s" % (x[1],x[0]) for x in sorted(in_degree.items(), key=lambda x: x[1], reverse=True )[0:args.top] if (x[1] > 0)))

    print("Top 10 Output Connections:")
    out_degree = gr.out_degree()
    print("\n".join( "\t%s\t%s" % (x[1],x[0]) for x in sorted(out_degree.items(), key=lambda x: x[1], reverse=True )[0:args.top] if (x[1] > 0)))

    for n_type in type_count:
        print("Node Type %s: %d" % (n_type, type_count[n_type]))
