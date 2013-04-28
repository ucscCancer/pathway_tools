#!/usr/bin/env python


import sys
import os
import argparse
from pathway_tools import convert
import networkx

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--xgmml', action="store_true", default=False)
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()

    gr = networkx.MultiDiGraph()

    for path in args.files:
        handle = open(path)
        ngr = convert.read_xgmml(handle)

        for node_name in ngr.node:
            if node_name not in gr.node:
                gr.add_node( node_name, attr_dict=ngr.node[node_name] )
            else:
                pass
                #if gr.node[node_name]['type'] != node_type:
                #    raise Exception("%s failure: Mismatch Node Type: %s :%s --> %s" % (path, node_name, gr.node[node_name]['type'], node_type ))
            if 'pathway' not in gr.node[node_name]:
                gr.node[node_name]['pathway'] = []
            gr.node[node_name]['pathway'].append(os.path.basename(path))


        for src, target, data in ngr.edges(data=True):
            if target in gr.edge[src]:
                has_edge = False
                for i in gr.edge[src][target]:
                    if gr.edge[src][target][i]['interaction'] == data['interaction']:
                        has_edge = True
                if not has_edge:
                    gr.add_edge(src, target, attr_dict=data)

    convert.write_xgmml(gr, sys.stdout)
