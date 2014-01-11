#!/usr/bin/env python

import rdflib
import sys
from glob import glob
import os
from multiprocessing import Pool


def scan_biopax_file(path):
    g = rdflib.Graph()
    g.parse(path)
    out = []
    for a in g.query("""
SELECT ?db ?id WHERE {
    ?path a <http://www.biopax.org/release/biopax-level3.owl#Pathway> .
    ?path <http://www.biopax.org/release/biopax-level3.owl#xref> ?xref .
    ?xref <http://www.biopax.org/release/biopax-level3.owl#db> ?db .
    ?xref <http://www.biopax.org/release/biopax-level3.owl#id> ?id 
}
            """):
        out.append( [os.path.basename(path), a[0], a[1]] )
    for a in g.query("""
SELECT ?dataSource ?name WHERE {
    ?path a <http://www.biopax.org/release/biopax-level3.owl#Pathway> .
    ?path <http://www.biopax.org/release/biopax-level3.owl#dataSource> ?dataSource .
    ?path <http://www.biopax.org/release/biopax-level3.owl#name> ?name .
}
            """):
        out.append( [os.path.basename(path), os.path.basename(a[0]), a[1]] )
    return out

if __name__ == "__main__":
    indir = sys.argv[1]
    fileset = glob(os.path.join(indir, "*.owl"))

    p = Pool(8)
    out = p.map(scan_biopax_file, fileset)
    for scan in out:
        for path, db, dbid in scan:
            print "%s\t%s\t%s" % (path, db, dbid)
