#!/usr/bin/env python 

import sys
import os
import re
import argparse
import networkx
from glob import glob
from pathway_tools.biopax import BioPaxFile, BioPaxSparql
from pathway_tools import convert

from multiprocessing import Pool

re_clean = re.compile(r'[/\\ \n;]')
re_namesplit = re.compile(r'[/ \#\&]')

class ConvertTask:

    def __init__(self, path_file, outdir=None, pathways=None, singlepath=None):
        self.path_file = path_file
        self.pathways = pathways
        self.outdir = outdir
        self.singlepath = singlepath

    def run(self):
        path_file = self.path_file
        if isinstance(path_file, basestring):
            t = BioPaxFile()
            t.load(path_file)
            path_file = t
        paths = []
        for key,value in path_file.pathways().iteritems():
            add = False
            if self.pathways is None or value in self.pathways:
                add = True
            elif key in self.pathways:
                add = True
            else:
                a_name = re_namesplit.split(key)[-1]
                if a_name in self.pathways:
                    add = True
            if add:
                paths.append(key)

        for subnet in path_file.toNet(paths):                
            gr = networkx.MultiDiGraph()

            for k, v in subnet.meta.items():
                if k != 'type':
                    gr.graph[k] = v
            subnet.to_graph(gr)

            name = re_namesplit.split(gr.graph['url'])[-1]
            if self.outdir is not None:
                handle = open(os.path.join(self.outdir, name + ".xgmml"), "w")
            elif self.singlepath is not None:
                handle = open(self.singlepath, "w")
            else:
                handle = sys.stdout
            
            convert.write_xgmml(gr, handle)
            
            if self.singlepath is not None:
                handle.close()
                return

def runner(x):
    x.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--pathways', action="append")
    parser.add_argument("-s", "--sparql", action="store_true", default=False)
    parser.add_argument("-l", "--list", action="store_true", default=False)
    parser.add_argument("-d", "--out-dir", default=None)
    parser.add_argument("--cpus", type=int, default=1)
   
    
    parser.add_argument("input", nargs="+")
    args = parser.parse_args()

    file_list = []

    if args.sparql:
        file_list = [ BioPaxSparql(args.input) ]
    else:
        for a in args.input:
            if os.path.isdir(a):
                file_list.extend(glob(os.path.join(a, "*")))
            else:
                file_list.append(a)

    if args.list:
        for path_file in file_list:
            if isinstance(path_file, basestring):
                t = BioPaxFile()
                t.load(path_file)
                path_file = t
            for a,b in path_file.pathways().iteritems():
                print a, b.encode('ascii', errors='ignore')
    else:

        if args.cpus > 1:
            runs = []
            for path_file in file_list:
                runs.append( ConvertTask(path_file, args.out_dir, args.pathways) )
            p = Pool(args.cpus)
            p.map(runner, runs)
        else:
            for path_file in file_list:
                ConvertTask(path_file, args.out_dir, args.pathways).run()


