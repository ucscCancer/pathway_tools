#!/usr/bin/env python

import sys
import urllib
import json
import os

if __name__ == "__main__":
	outdir = sys.argv[1]

	data = json.loads(urllib.urlopen("http://www.pathwaycommons.org/pc2/top_pathways.json").read())
	for i, hit in enumerate(data['searchHit']):
		print "%s of %s : %s" % (i, data['numHits'], hit['uri'])
		urllib.urlretrieve(hit['uri'], os.path.join(outdir, os.path.basename(hit['uri']) + ".owl"))
