#!/usr/bin/env python

"""
Some RDF/XML produced by PathwayCommon fails under certain parsers. This script
removes special characters from RDF IDs
"""

import sys
import xml.sax
import cgi


def clean_id(n):
    return n.replace("+", "_").replace("%", "_")

class BioPaxCleaner(xml.sax.ContentHandler):

    def __init__(self, handle):
        self.handle = handle

    def startDocument(self):
        self.handle.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>")

    def startElement(self, name, attrs):
        items = []
        for k,v in attrs.items():
            if k == "rdf:ID":
                items.append( (k,clean_id(v)) )
            else:
                items.append( (k,v) )
        line = "<%s %s>" % (clean_id(name), " ".join( "%s=\"%s\"" % (k,v) for k,v in items ))
        self.handle.write(line)

    def characters(self, chars):
        self.handle.write(cgi.escape(chars.encode('ascii', 'ignore')))

    def endElement(self, name):
        self.handle.write("</%s>" % (name))

if __name__ == "__main__":
    parser = xml.sax.make_parser()
    parser.setContentHandler(BioPaxCleaner(sys.stdout))
    parser.parse(open(sys.argv[1],"r"))
