#!/usr/bin/env python
"""check-pathway.py: 

Usage:
  check-pathway.py [options] pathwayFile

Options:
  -q            run quietly
  -o str        outputs repaired version if specified
"""
## Written By: Sam Ng
## Last Updated: ##/##/####
import os, os.path, sys, getopt, re
from pathway_tools import mPathway

verbose = True

def usage(code = 0):
    print __doc__
    if code != None: sys.exit(code)

def log(msg, die = False):
    if verbose:
        sys.stderr.write(msg)
    if die:
        sys.exit(1)

def syscmd(cmd):
    log("running:\n\t"+cmd+"\n")
    exitstatus = os.system(cmd)
    if exitstatus != 0:
        print "Failed with exit status %i" % exitstatus
        sys.exit(10)
    log("... done\n")

def main(args):
    ## parse arguments
    try:
        opts, args = getopt.getopt(args, "o:q")
    except getopt.GetoptError, err:
        print str(err)
        usage(2)
    
    if len(args) != 1:
        log("ERROR: incorrect number of arguments", die = True)
    
    inf = args[0]
    
    outf = None
    global verbose
    for o, a in opts:
        if o == "-o":
            outf = a
        elif o == "-q":
            verbose = False
    
    ## execute
    (n, i) = mPathway.rPathway(inf)
    p = mPathway.Pathway(n, i)
    p.selfTest()
    
    if outf != None:
        mPathway.wPathway(outf, p.nodes, p.interactions)

if __name__ == "__main__":
    main(sys.argv[1:])
