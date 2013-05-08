#!/usr/bin/env python

import sys
import rdflib
import rdflib.term
import os
from glob import glob
import networkx
from pathway_tools import convert
from StringIO import StringIO

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

class BioPax_ElementBase(object):
	def __init__(self, pax, node):
		self.pax = pax
		self.node = node
		self.name = None
		for a in self.pax.query(src=node, predicate=BIOPAX_BASE + "displayName"):
			self.name = a.dst
		if self.name is None:
			self.name = self.node

	def __str__(self):
		return "%s <%s>" % (self.name, self.type)

class BioPax_Pathway(BioPax_ElementBase):
	type = "pathway"

	def process(self):
		component_list = {}
		for component in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "pathwayComponent"):
			component_list[component.dst] = element_mapping[component.dst_type](self.pax, component.dst)

		gr = networkx.MultiDiGraph()
		re_map = {}
		for c in component_list:
			for r in component_list[c].process():
				if isinstance(r, Output_Edge):
					re_map[r.src.id] = r.src.label
					re_map[r.dst.id] = r.dst.label
					if r.src.id not in gr.node:
						gr.add_node( r.src.id, attr_dict=r.src.data )
					if r.dst.id not in gr.node:
						gr.add_node( r.dst.id, attr_dict=r.dst.data )
					gr.add_edge( r.src.id, r.dst.id, attr_dict=r.data )
				elif isinstance(r, Output_Node):
					if r.id not in gr.node:
						gr.add_node( r.id, attr_dict=r.data )
				else:
					#raise Exception("Unexpected parsing: %s %s" % (c, type(component_list[c])) )
					pass
		gr = networkx.relabel_nodes(gr, re_map)
		yield gr

class BioPax_SmallMolecule(BioPax_ElementBase):
	type = "smallmolecule"

	def process(self):
		node_label = None
		for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName"):
			node_label = name.dst.replace("\n", " ")

		xref_list = []
		for ref in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "entityReference"):
			for xref in self.pax.query(src=ref.dst, predicate=BIOPAX_BASE + "xref"):
				db = None
				db_id = None
				for rel in self.pax.query(src=xref.dst):
					if rel.predicate == BIOPAX_BASE + "db":
						db = rel.dst
					if rel.predicate == BIOPAX_BASE + "id":
						db_id = rel.dst
				db_xref = "%s:%s" % (db, db_id)
				xref_list.append(db_xref)

		yield Output_Node(self.node, node_label, {'db_xref' : xref_list, 'type' : self.type})


class BioPax_Complex(BioPax_ElementBase):
	type = "complex"

	def process(self):
		node_label = None
		#bug, using displayName may cause errors for complexes
		for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName"):
			node_label = name.dst.replace("\n", " ")

		components = []
		for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "component"):
			elem = element_mapping[c_info.dst_type](self.pax, c_info.dst)
			for e in elem.process():
				components.append( e )

		complex_node = Output_Node(self.node, node_label, {"type" : self.type})
		for c in components:
			yield Output_Edge( complex_node, c, {"type" : "component>"} )



class BioPax_PhysicalEntity(BioPax_ElementBase):
	type = "physicalentity"

	def process(self):
		yield None

class BioPax_Rna(BioPax_ElementBase):
	type = "rna"

	def process(self):
		yield None

class BioPax_Protein(BioPax_ElementBase):
	type = "protein"

	def process(self):
		node_label = None
		for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName"):
			node_label = name.dst.replace("\n", " ")

		xref_list = []
		for ref in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "entityReference"):
			for xref in self.pax.query(src=ref.dst, predicate=BIOPAX_BASE + "xref"):
				db = None
				db_id = None
				for rel in self.pax.query(src=xref.dst):
					if rel.predicate == BIOPAX_BASE + "db":
						db = rel.dst
					if rel.predicate == BIOPAX_BASE + "id":
						db_id = rel.dst
				db_xref = "%s:%s" % (db, db_id)
				xref_list.append(db_xref)

		yield Output_Node(self.node, node_label, {'db_xref' : xref_list, 'type' : self.type})


class BioPax_Catalysis(BioPax_ElementBase):
	type = "catalysis"

	def process(self):

		controller = []
		for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "controller"):
			elem = element_mapping[c_info.dst_type](self.pax, c_info.dst)
			controller.append( elem )
			
		controlled = []
		for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "controlled"):
			elem = element_mapping[c_info.dst_type](self.pax, c_info.dst)
			controlled.append( elem )


		for component in self.pax.query(src=self.node):
			#print "NA", component
			#yield "NA", component
			pass
		yield None

class BioPax_BiochemicalReaction(BioPax_ElementBase):
	type = "biochemical_reaction"

	def process(self):
		left = []
		for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "left"):
			elem = element_mapping[c_info.dst_type](self.pax, c_info.dst)
			for e in elem.process():
				left.append( e )
			
		right = []
		for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "right"):
			elem = element_mapping[c_info.dst_type](self.pax, c_info.dst)
			for e in elem.process():
				right.append( e )
		
		for l in left:
			for r in right:
				interaction = "-a>"
				if isinstance(l, Output_Node) and isinstance(r, Output_Node):
					yield Output_Edge( l, r, {'interaction' : interaction} )
				else:
					print "Skipped", l, r


class BioPax_MolecularInteraction(BioPax_ElementBase):
	type = "molecular_interaction"

	def process(self):
		part = []
		for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "participant"):
			elem = element_mapping[c_info.dst_type](self.pax, c_info.dst)
			for e in elem.process():
				part.append( e )

		for p1 in part:
			for p2 in part:
				if p1 != p2:
					if isinstance(p1, Output_Node) and isinstance(p2, Output_Node):
						yield Output_Edge( p1, p2, {'interaction' : "--"} )
					else:
						print "Skipped", p1, p2

class BioPax_Transport(BioPax_ElementBase):
	type = "transport"

	def process(self):
		yield None


element_mapping = {
	BIOPAX_BASE + "Protein" : BioPax_Protein,
	BIOPAX_BASE + "SmallMolecule" : BioPax_SmallMolecule,
	BIOPAX_BASE + "Complex" : BioPax_Complex,
	BIOPAX_BASE + "PhysicalEntity" : BioPax_PhysicalEntity,
	BIOPAX_BASE + "Rna" : BioPax_Rna,
	BIOPAX_BASE + "Catalysis" : BioPax_Catalysis,
	BIOPAX_BASE + "BiochemicalReaction" : BioPax_BiochemicalReaction,
	BIOPAX_BASE + "MolecularInteraction" : BioPax_MolecularInteraction,
	BIOPAX_BASE + "Transport" : BioPax_Transport,
	BIOPAX_BASE + "Pathway" : BioPax_Pathway
}


class Output_Edge:
	def __init__(self, src, dst, data):
		self.src = src
		self.dst = dst
		self.data = data

class Output_Node:
	def __init__(self, id, label, data):
		self.id = id
		self.label = label
		self.data = data

class BioPax:
	
	def __init__(self):
		self.graph = {}
	
	def load(self, path):
		g = rdflib.Graph()
		result = g.parse(path)
		self._graph_parse(g)
	
	def parse(self, text):
		handle = StringIO(text)
		g = rdflib.Graph()
		result = g.parse(handle)
		self._graph_parse(g)
	
	def _graph_parse(self, g):
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
			p_elem = BioPax_Pathway(self, pathway)
			gr = p_elem.process()
			yield gr
			"""
			gr = networkx.MultiDiGraph()
			#log("Pathway: " + pathway)
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
							gr.add_node(l.name.replace("\n", " "), type=l.type)

					for r in right:
						if r.name not in gr.node:
							gr.add_node(r.name.replace("\n", " "), type=r.type)
						
					for l in left:
						for r in right:
							interaction = "-a>"
							gr.add_edge( l.name.replace("\n", " "), r.name.replace("\n", " "), None, {'interaction' : interaction} )
					
				elif  component_list[component] == BIOPAX_BASE + "Catalysis":
					pass
				else:
					log( "Unknown Pathway component: " + component_list[component] )
				#for c_info in self.query(src=component):
				#	print c_info
				#print component, component_list[component]

			for protein in self.query(predicate=TYPE_PRED, dst=BIOPAX_BASE + "Protein"):
				for name in self.query(src=protein.src, predicate=BIOPAX_BASE + "displayName"):
					node_name = name.dst.replace("\n", " ")
					if node_name not in gr.node:
						gr.add_node(node_name, type="protein")
					for ref in self.query(src=protein.src, predicate=BIOPAX_BASE + "entityReference"):
						for xref in self.query(src=ref.dst, predicate=BIOPAX_BASE + "xref"):
							db = None
							db_id = None
							for rel in self.query(src=xref.dst):
								if rel.predicate == BIOPAX_BASE + "db":
									db = rel.dst
								if rel.predicate == BIOPAX_BASE + "id":
									db_id = rel.dst
							db_xref = "%s:%s" % (db, db_id)
							gr.node[node_name]['db_xref'] = db_xref
			yield gr
			"""
	
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
	
