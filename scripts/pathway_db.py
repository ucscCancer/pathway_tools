#!/usr/bin/env python

import argparse
import sys
import os
import re
import yaml
import shutil
import subprocess

reKeyValue = re.compile(r'^(\S+)\s*=\s*(.*)$')

def check_repo(basedir):
	handle = open(os.path.join(basedir, "INFO"))
	conf = yaml.load(handle.read())
	handle.close()
	
	if 'PID' not in conf:
		raise Exception("Missing required PID field")

	if 'DESC' not in conf:
		raise Exception("Missing required DESC description field")
		
	if 'SOURCE' not in conf:
		raise Exception("Missing required SOURCE field")
	
	return conf

def add_repo(input, conf, basedir):
	pathway_name = "PID%s" % (conf['PID'])
	dstdir = os.path.join(basedir, pathway_name)
	if os.path.exists(dstdir):
		raise Exception("Pathway PID%s already exists" % (conf['PID']))
	print "Adding", dstdir
	os.mkdir(dstdir)
	for f in ['INFO', 'graph']:
		shutil.copy(os.path.join(input, f), dstdir)
	subprocess.check_call("cd %s; git add %s; git commit -m 'Adding Pathway %s'" % (basedir, pathway_name, pathway_name), shell=True)


def main_add(args):
	parser = argparse.ArgumentParser(prog="pathway_db add")
	parser.add_argument("-b", "--base-dir", help="BaseDir", default=None)
	parser.add_argument("input", help="Input Directory")
	
	args = parser.parse_args(args)
	
	if args.base_dir is None or not os.path.exists(args.base_dir):
		sys.stderr.write("Define Pathway REPO\n")
		sys.exit(1)
	
	try:
		conf = check_repo(args.input)
		add_repo(args.input, conf, args.base_dir)
	except Exception, e:
		sys.stderr.write("Pathway Check Error: %s : %s\n" % (args.input, str(e)))
		sys.exit(1)
	

if __name__ == "__main__":
	mode = sys.argv[1]

	if mode == 'add':
		main_add(sys.argv[2:])
