#!/usr/bin/env python 

import sys
from pathway_tools.biopax import BioPax
from pathway_tools import convert

if __name__ == "__main__":
	b = BioPax()
	b.load(sys.argv[1])
	
	for gr in b.toNet():
		convert.write_paradigm_graph(gr, sys.stdout)