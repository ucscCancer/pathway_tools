#!/usr/bin/env python

import sys
import os

db_order = [
	"pid",
	"Reactome",
	"BioCyc",
	"PANTHER Pathway"
]

db_name = {
	"pid" : "PID",
	"Reactome" : "REACTOME",
	"BioCyc" : "BIOCYC",
	"PANTHER Pathway" : "PANTHER"
}

if __name__ == "__main__":

	pathmap = {}

	handle = open(sys.argv[1])
	for line in handle:
		path, db, db_id = line.rstrip().split("\t")
		if path not in pathmap:
			pathmap[path] = []
		if db in db_order:
			pathmap[path].append( "%s_%s" % (db_name[db], db_id) )
	handle.close()
	for path in pathmap:
		for newname in pathmap[path]:
			print "%s\t%s" % (path, newname)
		
