#!/usr/bin/env python

import sys
import rdflib
import rdflib.term
import os
from glob import glob
import networkx
from pathway_tools import convert

TYPE_PRED = "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
BIOPAX_BASE = "http://www.biopax.org/release/biopax-level3.owl#"


def log(message):
	sys.stderr.write(message + "\n")

class Link(object):
	def __init__(self, src_type, src, predicate, dst_type, dst):
		self.src_type = src_type
		self.src = src
		self.predicate = predicate
		self.dst_type = dst_type
		self.dst = dst

	def __str__(self):
		return "(%s,%s) %s (%s,%s)" % (self.src_type, self.src, self.predicate, self.dst_type, self.dst)

class ElementBase(object):
	def __init__(self, pax, node):
		self.node = node
		self.name = None
		for a in pax.query(src=node, predicate=BIOPAX_BASE + "displayName"):
			self.name = a.dst

	def __str__(self):
		return self.name


class Protein(ElementBase):
	type = "protein"


class SmallMolecule(ElementBase):
	type = "smallmolecule"
	

class Complex(ElementBase):
	type = "complex"

class PhysicalEntity(ElementBase):
	type = "physicalentity"


class Rna(ElementBase):
	type = "rna"


element_mapping = {
	BIOPAX_BASE + "Protein" : Protein,
	BIOPAX_BASE + "SmallMolecule" : SmallMolecule,
	BIOPAX_BASE + "Complex" : Complex,
	BIOPAX_BASE + "PhysicalEntity" : PhysicalEntity,
	BIOPAX_BASE + "Rna" : Rna
}


class BioPax:
	
	def __init__(self):
		self.graph = {}
	
	def load(self, path):
		g = rdflib.Graph()
		result = g.parse(sys.argv[1])
		for subj, pred, obj in g:
			subj_text = ("%s" % (subj)).encode('ascii', errors='ignore').rstrip()
			pred_text = ("%s" % (pred)).encode('ascii', errors='ignore').rstrip()
			obj_text =  ("%s" % (obj)).encode('ascii', errors='ignore').rstrip()
			if subj_text not in self.graph:
				self.graph[subj_text] = {}
			if pred_text not in self.graph[subj_text]:
				self.graph[subj_text][pred_text] = []
			self.graph[subj_text][pred_text].append( obj_text )
			

	def toNet(self):
		pathway_list = {}
		for a in self.query(src_type=BIOPAX_BASE + "Pathway"):
			pathway_list[a.src] = True

		for pathway in pathway_list:
			gr = networkx.MultiDiGraph()
			log("Pathway: " + pathway)
			name = None
			for a in self.query(src=pathway, predicate= BIOPAX_BASE + "displayName"):
				name = a.dst
			
			component_list = {}
			for component in self.query(src=pathway, predicate=BIOPAX_BASE + "pathwayComponent"):
				component_list[component.dst] = component.dst_type
			
			for component in component_list:
				if component_list[component] == BIOPAX_BASE + "BiochemicalReaction":
					left = []
					for c_info in self.query(src=component, predicate=BIOPAX_BASE + "left"):
						elem = element_mapping[self.get_node_type(c_info.dst)](self, c_info.dst)
						left.append( elem )
						
					right = []
					for c_info in self.query(src=component, predicate=BIOPAX_BASE + "right"):
						elem = element_mapping[self.get_node_type(c_info.dst)](self, c_info.dst)
						right.append( elem )
					
					for l in left:
						if l.name not in gr.node:
							gr.add_node(l.name, type=l.type)

					for r in right:
						if r.name not in gr.node:
							gr.add_node(r.name, type=r.type)
						
					for l in left:
						for r in right:
							interaction = "-a>"
							gr.add_edge( l.name, r.name, None, {'interaction' : interaction} )
					
				elif  component_list[component] == BIOPAX_BASE + "Catalysis":
					pass
				else:
					log( "Unknown Pathway component: " + component_list[component] )
				#for c_info in self.query(src=component):
				#	print c_info
				#print component, component_list[component]
			yield gr
	
	def get_node_type(self, node):
		return self.graph.get(node, {}).get(TYPE_PRED, [None])[0]
	
	def str_compare(self, str1, str2, check_case):
		if check_case:
			return str1 == str2
		print str1, str2
		return str1.lower() == str2.lower()

	def query(self, src_type=None, src=None, predicate=None, dst_type=None, dst=None, check_case=True):
		out = []
		include = [True,True,True,True,True]
		for src_cmp in self.graph:
			include[0] = False
			if src_type is None or self.str_compare(self.get_node_type(src_cmp), src_type, check_case):
				include[0] = True
			include[1] = False
			if src is None or self.str_compare(src_cmp, src, check_case):
				include[1] = True
				
			for pred_cmp in self.graph[src_cmp]:
				include[2] = False
				if predicate is None or self.str_compare(pred_cmp, predicate, check_case):
					include[2] = True
				for dst_cmp in self.graph[src_cmp][pred_cmp]:
					include[3] = False
					if dst_type is None or self.str_compare(self.get_node_type(dst_cmp), dst_type, check_case):
						include[3] = True
					include[4] = False
					if dst is None or self.str_compare(dst, dst_cmp, check_case):
						include[4] = True
					if False not in include:
						out.append( Link(
							self.get_node_type(src_cmp), src_cmp, 
							pred_cmp, 
							self.get_node_type(dst_cmp), dst_cmp) 
						)
		return out			
	
