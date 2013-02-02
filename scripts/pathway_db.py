#!/usr/bin/env python

import argparse
import sys
import os
import re

reKeyValue = re.compile(r'^(\S+)\s*=\s*(.*)$')

def parse_config(handle):
	config = {}
	for line in handle:
		res = reKeyValue.search(line)
		if res:
			key, value = res.groups()
			config[key] = value
	return config

def check_repo(basedir):
	handle = open(os.path.join(basedir, "INFO"))
	conf = parse_config(handle)
	handle.close()
	
	if 'DESC' not in conf:
		raise Exception("Missing required DESC description field")
		
	if 'SOURCE' not in conf:
		raise Exception("Missing required SOURCE field")
	
	

def main_add(args):
	parser = argparse.ArgumentParser(prog="pathway_db add")
	parser.add_argument("-b", "--base-dir", help="BaseDir", default=None)
	parser.add_argument("input", help="Input Directory")
	
	args = parser.parse_args(args)
	
	if args.base_dir is None or not os.path.exists(args.base_dir):
		sys.stderr.write("Define Pathway REPO\n")
		sys.exit(1)
	
	try:
		check_repo(args.input)
	except Exception, e:
		sys.stderr.write("Pathway Check Error: %s : %s\n" % (args.input, str(e)))
		sys.exit(1)
		
if __name__ == "__main__":
	mode = sys.argv[1]

	if mode == 'add':
		main_add(sys.argv[2:])
