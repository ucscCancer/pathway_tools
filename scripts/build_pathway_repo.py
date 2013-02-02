#!/usr/bin/env python


import sys
import os
import network_convert

if __name__ == "__main__":
	handle = open(sys.argv[1])
	gr = network_convert.read_xgmml(handle)
	handle.close()

	name_tab = {}
	handle = open(sys.argv[2])
	for line in handle:
		tmp = line.rstrip().split("\t")
		name_tab[tmp[1]] = [tmp[0], tmp[2]]
	outdir = sys.argv[3]
	
	pathways = {}
	for n in gr.node:
		if 'pathway' in gr.node[n]:
			for path in gr.node[n]['pathway']:
				pathways[path] = True
	
	
	for i, path in enumerate(pathways):
		pathid = "TMP%d" % (i+1)
		print pathid
		pathdir = os.path.join(outdir, pathid)
		if not os.path.exists(pathdir):
			os.mkdir(pathdir)
		handle = open(os.path.join(pathdir, "INFO"), "w")
		print path
		handle.write(u"DESC=%s\n" % (unicode(path)))
		if path in name_tab:
			handle.write("PID=%s\n" %(name_tab[path][0]))
			handle.write("SOURCE=%s\n" %(name_tab[path][1]))
		handle.close()
		nodes = {}
		for n in gr.node:
			if 'pathway' in gr.node[n] and path in gr.node[n]['pathway']:
				nodes[n] = True
		handle = open(os.path.join(pathdir, "graph"), "w")
		for n in nodes:
			handle.write("%s\t%s\n" % (gr.node[n]['type'], n))
		for s in nodes:
			for t in gr.edge[s]:
				if t in nodes:
					handle.write("%s\t%s\t%s\n" % (s, t, gr.edge[s][t]))
		handle.close()
		
	
	
