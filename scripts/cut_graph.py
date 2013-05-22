#!/usr/bin/env python

# USAGE: cutGraph -i <input node heats> -n <network .sif format> 
# Desc: Finds graph statistics and the number of connected components for edges remaining after a range of heat thresholds are applied. The parent graph, 
# is supplied, and is treated as undirected.
# Depends: networkx 1.7, python 2.7
# Author: epaull (Evan Paull)
# Date: 4/8/13

import os
import sys
from collections import defaultdict
from optparse import OptionParser
parser = OptionParser()
parser.add_option("-c","--cut_graph", dest="cut_graph",action="store", default=None, help="Just print the graph at supplied cutoff")
parser.add_option("-i","--heats", dest="heats",action="store", default=None, help="Input (diffused) Heats")
parser.add_option("-n","--network",dest="network",action="store", default=None, help="Base Network in UCSC Pathway Format")
parser.add_option("-s","--subdivs",dest="subdivs",action="store", default=100, help="Number of Subdivisions (per heat increment of 1) to test in the Range")
(opts, args) = parser.parse_args()

from pathway_tools import convert
import networkx as nx

def parseHeats(file):
	heats = {}
	for line in open(file, 'r'):
		parts = line.rstrip().split("\t")
		heats[parts[0]] = float(parts[1])

	return heats

def cutGraph(graph, heats, cutoff):

	GC = nx.MultiGraph()
	for edge in graph.edges_iter(data=True):
		source = edge[0]
		target = edge[1]
		interaction = edge[2]['interaction']
		# edge in graph, but not in input heats
		if source not in heats or target not in heats:
			continue

		if (heats[source] >= cutoff) and (heats[target] >= cutoff):
			GC.add_edges_from([(source, target, dict(interaction=interaction))])

	return GC

def maxCC(ccs):

	max = 0
	for cc in ccs:
		if len(cc) > max:
			max = len(cc)

	return max

graph = None
if opts.network.endswith(".sif"):
	graph = convert.read_sif(open(opts.network))
else:
	graph = convert.read_paradigm_graph(open(opts.network))

heats = parseHeats(opts.heats)
if opts.cut_graph:
	cut_val = None
	try:
		cutoff = float(opts.cut_graph)
	except:
		raise Exception("Error: value supplied to cut_graph must be a positive number")

	cutG = cutGraph(graph, heats, cutoff)
	convert.write_sif(cutG, sys.stdout)
	sys.exit(1)

max_heat = max(heats.values())
print "Cutoff\tNum Edges\tNum Connected Components\tLargest Connected ComponentmaxCC(ccs)\tEdge Biggest CC Ratio"
for cutoff in range(0, int(max_heat*opts.subdivs)+1):
	cutoff = cutoff/float(opts.subdivs)
	cutG = cutGraph(graph, heats, cutoff)
	l_ccs = nx.number_connected_components(cutG)
	ccs = nx.connected_components(cutG)
	max_ccs = maxCC(ccs)
	l_edges = len(cutG.edges())
	edge_bigccs_ratio = l_edges/float(l_ccs)
	print "\t".join([str(i) for i in [cutoff, l_ccs, max_ccs, l_edges, edge_bigccs_ratio]])
