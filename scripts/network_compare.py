#!/usr/bin/env python

import network_convert
import sys


if __name__ == "__main__":
	handle1 = open(sys.argv[1])
	gr1 = network_convert.read_paradigm_graph(handle1)
	handle1.close()
	
	handle2 = open(sys.argv[2])
	gr2 = network_convert.read_paradigm_graph(handle2)
	handle2.close()
	
	for n in gr1.node:
		if n not in gr2.node:
			print "Graph2 Missing:", n 
