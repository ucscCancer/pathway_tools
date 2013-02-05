#!/usr/bin/env python

import network_convert
import sys
import networkx

if __name__ == "__main__":
    handle1 = open(sys.argv[1])
    gr1 = network_convert.read_paradigm_graph(handle1)
    handle1.close()
    
    handle2 = open(sys.argv[2])
    gr2 = network_convert.read_paradigm_graph(handle2)
    handle2.close()
    

    gr = networkx.DiGraph()

    for n in gr1.node:
        if n not in gr2.node:
            sys.stderr.write( "Graph2 Missing Node:" + n + "\n" )
            gr.add_node(n)
            gr.node[n] = gr1.node[n]
    for n in gr1.node:
        for t in gr1.edge[n]:
            if n not in gr2.node or t not in gr2.edge[n]:
                for e in [n,t]:
                    if e not in gr.node:
                        gr.add_node(e)
                        gr.node[e] = gr1.node[e]

                gr.add_edge(n,t)
                gr.edge[n][t] = gr1.edge[n][t]
                sys.stderr.write( "Graph2 Missing Edge:" + str((n,t)) + "\n")
    for n in gr2.node:
        if n not in gr1.node:
            sys.stderr.write( "Graph1 Missing Node:" + n + "\n" )
        else:
            for t in gr2.edge[n]:
                if t not in gr1.edge[n]:
                    sys.stderr.write( "Graph1 Missing Edge:" + str((n,t)) + "\n")
    
    network_convert.write_paradigm_graph(gr, sys.stdout)