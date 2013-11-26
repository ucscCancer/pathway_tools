#!/usr/bin/env python


import sys
import os
from pathway_tools import convert
import yaml

if __name__ == "__main__":
	handle = open(sys.argv[1])
	gr = convert.read_spf(handle, strict=False)
	handle.close()

	name_tab = {}
	handle = open(sys.argv[2])
	for line in handle:
		tmp = line.rstrip().split("\t")
		name_tab[tmp[0]] = tmp[1:]
	outdir = sys.argv[3]
	
	i = 1
	for pathid in name_tab:
		print pathid
		pathdir = os.path.join(outdir, pathid)
		if not os.path.exists(pathdir):
			os.mkdir(pathdir)
		config = {"DESC" : str(pathid)}
		#if path in name_tab:
		#	config["PID"] = int(name_tab[path][0])
		#	config["SOURCE"] = str(name_tab[path][1])
		
		handle = open(os.path.join(pathdir, "INFO"), "w")
		print pathid
		handle.write(yaml.dump(config, default_flow_style=False))
		handle.close()
		nodes = {}
		for n in gr.node:
			if n in name_tab[pathid]:
				nodes[n] = True
		handle = open(os.path.join(pathdir, "graph"), "w")
		for n in nodes:
			handle.write("%s\t%s\n" % (gr.node[n]['type'], n))
		for s in nodes:
			for t in gr.edge[s]:
				if t in nodes:
					for e in gr.edge[s][t]:
						handle.write("%s\t%s\t%s\n" % (s, t, gr.edge[s][t][e]['interaction']))
		handle.close()
		
	
	
