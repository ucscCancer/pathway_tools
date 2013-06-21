#!/usr/bin/env python

import sys
import argparse
import networkx
from pathway_tools import convert
from copy import copy


def read_vec_file(path):
    vec = {}
    handle = open(path)
    head = True
    for line in handle:
        if head:
            head = False
        else:
            tmp = line.rstrip().split("\t")
            vec[tmp[0]] = float(tmp[1])
    return vec

class diffuse_heat:
    def __init__(self, delta):
        self.delta = delta

    def calc(self, n):
        key, n_set = n
        x = []
        y = None
        for k, v in n_set.items():
            if k == key:
                y = v
            else:
                x.append(v)
        return ( key, y - (y * float(len(x)) - sum(x)) * self.delta )

def diffuse_net(gr, data, ncycles, delta):

    heat_fun = diffuse_heat(delta)

    current = copy(data)
    for i in range(ncycles):
        vec = []
        for n in gr.node:
            if n in current:
                n_set = {}
                for e in gr.edge[n]:
                    n_set[e] = current[e]
                n_set[n] = current[n]
                vec.append((n, n_set))
        out = map( heat_fun.calc, vec )

        update = {}
        for n in out:
            update[n[0]] = n[1]
        current = update
    return current


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("network", default=None)
    parser.add_argument("values", default=None)
    parser.add_argument("--cycles", type=int, default=10)
    parser.add_argument("--timestep", type=float, default=0.0001)
    parser.add_argument("--header", action="store_true", default=False)
        
    args = parser.parse_args()

    handle = open(args.network)
    gr = networkx.Graph(convert.read_sif(handle))
    handle.close()

    val = read_vec_file(args.values)

    current = {}
    for n in gr.node:
        current[n] = val.get(n, None)

    found = True
    while found:
        update = {}
        found = False
        for n in current:
            if current[n] is not None:
                update[n] = current[n]
            else:
                n_set = []
                for e in gr.edge[n]:
                    if current[e] is not None:
                        n_set.append(current[e])
                if len(n_set):
                    found = True
                    update[n] = sum(n_set) / float(len(n_set))
                else:
                    update[n] = None
        current = update
    for i in current.keys():
        if current[i] is None:
            del current[i]

    handle = sys.stdout
    output = diffuse_net(gr, current, args.cycles, args.timestep)
    if args.header:
        handle.write("symbol\tinput\tstart\tend\tchange\n")
    for n in current:
        handle.write("%s\t%s\t%s\t%s\t%s\n" % (
            n, 
            val.get(n,None), 
            current[n], output[n], 
            output[n] - current[n])
        )
