#!/usr/bin/env python

import networkx
from xml.dom.minidom import Document
import xml.sax
import argparse
import sys
import zipfile
import os
from glob import glob

def read_paradigm_graph(handle):
    gr = networkx.DiGraph()
    for line in handle:
        tmp = line.rstrip().split("\t")
        if len(tmp) == 2:
            gr.add_node( tmp[1], type=tmp[0] )
        elif len(tmp) == 3:
            gr.add_edge(tmp[0], tmp[1], interaction=tmp[2])
    return gr

def write_paradigm_graph(gr, handle, node_type_field='type', node_type_default='protein', edge_type_field='type', edge_type_default='-a>'):
    for e in sorted(gr.node):
        handle.write("%s\t%s\n" % (gr.node[e].get(node_type_field, node_type_default), e))

    for src in gr.edge:
        for dst in gr.edge[src]:
            handle.write("%s\t%s\t%s\n" % (src, dst, gr.edge[src][dst].get(edge_type_field, edge_type_default)))


def load_paradigm_dir(base_dir):
    pid_names = {}
    handle = open(os.path.join(base_dir, "names.tab"))
    for line in handle:
        tmp = line.rstrip().split("\t")
        pid_names[tmp[0]] = [  
            unicode(tmp[1], 'ascii', errors="ignore"), 
            unicode(tmp[2], 'ascii', errors="ignore")
        ]   
    handle.close()
    gr = networkx.DiGraph()
    for a in glob(os.path.join(base_dir, "pid*.tab")):
        pid_name = os.path.basename(a).split("_")[1]
        handle = open(a)
        for line in handle:
            tmp = line.rstrip().split("\t")
            if len(tmp) == 2:
                if tmp[1] not in gr.node:
                    gr.add_node( tmp[1], { 'type' : tmp[0], 'pathway' : [] } )
                gr.node[tmp[1]]['pathway'].append( pid_names[pid_name][0] )
            elif len(tmp) == 3:
                gr.add_edge(tmp[0], tmp[1], interaction=tmp[2])
        handle.close()
    return gr


def write_xgmml(gr, handle):
    doc = Document()
    graph_node = doc.createElement('graph')
    graph_node.setAttribute("xmlns", "http://www.cs.rpi.edu/XGMML")


    graph_node.setAttribute("xmlns:dc", "http://purl.org/dc/elements/1.1/")
    graph_node.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink" )
    graph_node.setAttribute("xmlns:rdf", "http://www.w3.org/1999/02/22-rdf-syntax-ns#" )
    graph_node.setAttribute("xmlns:cy", "http://www.cytoscape.org" )
    graph_node.setAttribute("directed", "1")
    doc.appendChild(graph_node)

    name_map = {}
    for i, n in enumerate(gr.node):
        name_map[n] = i
        node = doc.createElement('node')
        node.setAttribute('label', n)
        node.setAttribute('id', str(i))
        for key, value in gr.node[n].items():
            att_node = doc.createElement('att')
            att_node.setAttribute('name', key)
            if type(value) == float:
                att_node.setAttribute('value', str(value))
                att_node.setAttribute('type', "real")
            elif type(value) == int:
                att_node.setAttribute('value', str(value))
                att_node.setAttribute('type', "integer")
            elif type(value)== list:
                att_node.setAttribute('type', "list")
                for elm in value:
                    list_node = doc.createElement("att")
                    list_node.setAttribute("name", key)
                    list_node.setAttribute('value', str(elm))
                    if type(elm) == float:
                        list_node.setAttribute("type", "real")
                    elif type(elm) == int:
                        list_node.setAttribute("type", "integer")
                    else:
                        list_node.setAttribute("type", "string")
                    att_node.appendChild(list_node)
            else:
                att_node.setAttribute('value', str(value))
                att_node.setAttribute('type', "string")
            node.appendChild(att_node)
        graph_node.appendChild(node)

    for source in gr.edge:
        for target in gr.edge[source]:
            edge = doc.createElement("edge")
            edge.setAttribute("label", "%s - %s" % (source, target))
            edge.setAttribute("source", str(name_map[source]))
            edge.setAttribute("target", str(name_map[target]))
            for key, value in gr.edge[source][target].items():
                att_node = doc.createElement('att')
                att_node.setAttribute('name', key)
                att_node.setAttribute('value', str(value))
                if type(value) == float:
                    att_node.setAttribute('type', "real")
                elif type(value) == int:
                    att_node.setAttribute('type', "integer")
                else:
                    att_node.setAttribute('type', "string")
                edge.appendChild(att_node)
            graph_node.appendChild(edge)

    doc.writexml(handle, addindent=" ", newl="\n")

class GraphComp:
    def __init__(self, parser, attrs):
        self.graph = networkx.DiGraph()

    def pop(self):
        return self.graph

class NodeComp:
    def __init__(self, parser, parent, attrs):
        self.label = attrs.getValue('label')
        self.parser = parser
        self.parent = parent
        parser.id_map[attrs.getValue('id')] = self.label
        parent.graph.add_node( self.label )

    def add_att(self, name, value):
        self.parent.graph.node[self.label][name] =  value

    def pop(self):
        pass


class EdgeComp:
    def __init__(self, parser, parent, attrs):
        self.parent = parent
        self.src = parser.id_map[attrs.getValue('source')]
        self.target = parser.id_map[attrs.getValue('target')]
        parent.graph.add_edge(self.src,self.target)

    def add_att(self, name, value):
        self.parent.graph.edge[self.src][self.target][name] =  value

    def pop(self):
        pass


class AttComp:
    def __init__(self, parser, parent, attrs):
        self.parent = parent
        self.name = attrs.getValue('name')
        if attrs.has_key('value'):
            parent.add_att(attrs.getValue('name'), attrs.getValue('value'))
        self.att_list = None
        if attrs.getValue("type") == "list":
            self.att_list = []

    def add_att(self, name, value):
        self.att_list.append( value )

    def pop(self):
        if self.att_list is not None:
            self.parent.add_att(self.name, self.att_list)


class XGMMLHandler(xml.sax.ContentHandler):
    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        self.elem_stack = []
        self.id_map = {}
        self.last_value = None

    def startElement(self, name, attrs):
        #print("startElement '" + name + "'")
        if name == "graph":
            self.elem_stack.append(GraphComp(self, attrs))
        elif name == 'node':
            self.elem_stack.append(NodeComp(self, self.elem_stack[-1], attrs))
        elif name == 'edge':
            self.elem_stack.append(EdgeComp(self, self.elem_stack[-1], attrs))
        elif name == 'att':
            self.elem_stack.append(AttComp(self, self.elem_stack[-1], attrs))
        else:
            self.elem_stack.append(None)
            
 
    def endElement(self, name):
        self.last_value = self.elem_stack.pop().pop()

    def result(self):
        return self.last_value



def read_xgmml(handle):
    handler = XGMMLHandler()
    xml.sax.parse(handle, handler)
    return handler.result()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-etd', '--edge-type-default', help="Default Edge Type", default='-a>')
    parser.add_argument('-ntd', '--node-type-default', help='Default Node Type', default='protein')
    parser.add_argument('-ntf', '--node-type-field', help='Node Type Field', default='type')
    parser.add_argument('-etf', '--edge-type-field', help='Edge Type Field', default='interaction')
    parser.add_argument('-pmf', '--pathway-member-field', help="Pathway Membership Field", default="pathway")

    parser.add_argument('-in-xgmml', help='Input XGMML file', default=None)
    parser.add_argument('-in-paradigm', help='Input Paradigm file', default=None)
    parser.add_argument('-in-gmt', help='Input Gene Mapping File', default=None)
    parser.add_argument('-in-paradigm-dir', help='Input Paradigm Directory', default=None)
        
    #parser.add_argument('-cys', help='Input Cytoscape file', default=None)    

    parser.add_argument('-out-paradigm', help='Output Paradigm File', default=None)    
    parser.add_argument('-out-xgmml', help='Output XGMML File', default=None)    
    parser.add_argument('-out-gmt', help="Output Gene Mapping Table", default=None)
    
    args = parser.parse_args()

    gr = None
    if args.in_xgmml is not None:
        handle = open(args.in_xgmml)
        gr = read_xgmml(handle)
        handle.close()

    if args.in_paradigm is not None:
        handle = open(args.in_paradigm)
        gr = read_paradigm_graph(handle)
        handle.close()

    if args.in_gmt is not None:
        handle = open(args.in_gmt)
        for line in handle:
            tmp = line.rstrip().split("\t")
            for elem in tmp[1:]:
                if elem in gr.node:
                    if args.pathway_member_field not in gr.node[elem]:
                        gr.node[elem][args.pathway_member_field] = []
                    gr.node[elem][args.pathway_member_field].append(tmp[0])
        handle.close()

    if args.in_paradigm_dir is not None:
        gr = load_paradigm_dir(args.in_paradigm_dir)

    """
    #not yet complete, because while cytoscape stores an xgmml file internally, it doesn't populate the 
    #att fields until export...
    if args.cys is not None:
        z = zipfile.ZipFile(args.cys)
        for n in z.namelist():
            if n.endswith(".xgmml") and os.path.dirname(n).endswith("networks"):
                handle = z.open(n)
                gr = read_xgmml(handle)
                handle.close()
    """

    if gr is not None:
        if args.out_paradigm is not None:
            if args.out_paradigm == "-":
                ohandle = sys.stdout
            else:
                ohandle = open(args.out_paradigm, "w")
            write_paradigm_graph(gr, ohandle, 
                node_type_field=args.node_type_field, node_type_default=args.node_type_default, 
                edge_type_field=args.edge_type_field, edge_type_default=args.edge_type_default
            )
        if args.out_xgmml is not None:
            if args.out_xgmml == "-":
                ohandle = sys.stdout
            else:
                ohandle = open(args.out_xgmml, "w")
            write_xgmml(gr, ohandle)
        if args.out_gmt:
            if args.out_gmt == "-":
                ohandle = sys.stdout
            else:
                ohandle = open(args.out_gmt, "w")

            pathmap = {}
            for node in gr.node:
                if 'pathway' in gr.node[node]:
                    for p in gr.node[node]['pathway']:
                        if p not in pathmap:
                            pathmap[p] = [node]
                        else:
                            pathmap[p].append(node)
            for path in pathmap:
                ohandle.write("%s\t%s\n" % (path, "\t".join(pathmap[path])))



