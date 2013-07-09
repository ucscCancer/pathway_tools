#!/usr/bin/env python

import sys
import sparql
import rdflib
import argparse

class SparqlServer:

    def __init__(self, url):
        self.url = url
        self.sparql = sparql.Service(url)

    def node_query(self, node):
        query = """select * where {
            <%s> ?p ?e ;
        }""" % (node)
        results = self.sparql.query(query)
        for row in results.fetchall():
            yield row[0], row[1]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--stop-list")
    parser.add_argument("-m", "--max-depth", type=int, default=10)
    parser.add_argument("server")
    parser.add_argument("startpoint")
    args = parser.parse_args()

    start = args.startpoint

    hit_map = { start : 0 }
    added = [ start ]
    server = SparqlServer(args.server)

    output = rdflib.Graph()

    while len(added):
        next = []
        for n in added:
            for edge, dst in server.node_query(n):
                if isinstance(dst, sparql.IRI):
                    dst_str = dst.value.encode('ascii', errors='ignore')
                    output.add( (rdflib.URIRef(n), rdflib.URIRef(str(edge)), rdflib.URIRef(dst_str)) )
                    if dst_str not in hit_map:
                        nd = hit_map[n] + 1
                        hit_map[dst_str] = nd
                        if nd < args.max_depth:
                            next.append(dst_str)
                else:
                    dst_str = dst.value.encode('ascii', errors='ignore')
                    output.add( (rdflib.URIRef(n), rdflib.URIRef(str(edge)), rdflib.Literal(dst_str)) )
            added = next

    sys.stdout.write(output.serialize(format="pretty-xml"))
