#!/usr/bin/env python

import sys
import os
import re
import rdflib
import argparse
import requests



def sparql_query(url, query):
    data = requests.post(url, data={ "query" : query }, headers={"Accept" : "application/sparql-results+json"}).json()
    for row in data['results']['bindings']:
        yield row



class SparqlServer:

    def __init__(self, url):
        self.url = url

    def typed_node_query(self, type):
        query = """select * where {
            ?a a <%s>;
        }""" % (type)
        
        for row in sparql_query(self.url, query):
            yield { 'type' : row['a']['type'], 'value' : row['a']['value'] }


    def node_query(self, node):
        query = """select * where {
            <%s> ?p ?e ;
        }""" % (node)
        for row in sparql_query(self.url, query):
            yield ({'value' : row['p']['value']}, { 'type' : row['e']['type'], 'value' : row['e']['value'] })



BIOPAX_BASE = "http://www.biopax.org/release/biopax-level3.owl#"
re_namesplit = re.compile(r'[/ \#\&]')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stop-list")
    parser.add_argument("-m", "--max-depth", type=int, default=25)
    parser.add_argument("server")
    parser.add_argument("outdir")

    args = parser.parse_args()

    server = SparqlServer(args.server)

    for pathway_url in server.typed_node_query(BIOPAX_BASE + "Pathway"):
        pathway = pathway_url['value'].encode('ascii', errors='ignore')
        print "Scanning", pathway
        hit_map = { pathway : 0 }
        added = [ pathway ]

        output = rdflib.Graph()

        while len(added):
            next = []
            for n in added:
                for edge, dst in server.node_query(n):
                    if dst['type'] == 'uri':
                        dst_str = dst['value'].encode('ascii', errors='ignore')
                        output.add( (rdflib.URIRef(n), rdflib.URIRef(str(edge['value'])), rdflib.URIRef(dst_str)) )
                        if dst_str not in hit_map:
                            nd = hit_map[n] + 1
                            hit_map[dst_str] = nd
                            if nd < args.max_depth:
                                next.append(dst_str)
                    else:
                        dst_str = dst['value'].encode('ascii', errors='ignore')
                        output.add( (rdflib.URIRef(n), rdflib.URIRef(str(edge['value'])), rdflib.Literal(dst_str)) )
                added = next

        name = re_namesplit.split(pathway)[-1]
        handle = open(os.path.join(args.outdir, name + ".rdf"), "w")
        handle.write(output.serialize(format="pretty-xml"))
        handle.close()
        
    
