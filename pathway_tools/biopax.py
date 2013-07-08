#!/usr/bin/env python

import sys
import rdflib
import rdflib.term
import sparql
import os
from glob import glob
import networkx
from pathway_tools import convert
from StringIO import StringIO

TYPE_PRED = "http://www.w3.org/1999/02/22-rdf-syntax-ns#type"
BIOPAX_BASE = "http://www.biopax.org/release/biopax-level3.owl#"


DEBUG = True


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
    def __init__(self, pax, node, stack):
        self.pax = pax
        self.node = node
        self.name = None
        self.stack = stack
        for a in self.pax.query(src=node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            self.name = a.dst
        if self.name is None:
            self.name = self.node

    def process_child(self, dst, dst_type):
        self.debug("Entering %s %s" % (dst, dst_type))
        out = element_mapping[dst_type](self.pax, dst, self.stack + [self]).process()
        #self.debug("Exiting " + dst)
        return out

    def debug(self, message):
        if DEBUG:
            sys.stderr.write("DEBUG: %s :%s\n" % (self.stack + [self], message))

    def __repr__(self):
        return "%s (%s) <%s>" % (self.name.encode('ascii', errors='ignore'), self.node, self.type)

class BioPax_Pathway(BioPax_ElementBase):
    type = "Pathway"

    def process(self):
        self.debug("Pathway: %s" % (self.node.encode("ascii", errors='ignore')))

        component_list = {}
        for component in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "pathwayComponent"):
            component_list[component.dst] = self.process_child(component.dst, component.dst_type)

        name = None
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "name", get_type=False):
            name = c_info.dst

        if name is None:
            for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
                name = c_info.dst

        out = Subnet({'name' : name, 'url' : self.node, 'type' : self.type})

        for c in component_list:
            if component_list[c] is not None:
                out.add_node( c, component_list[c] )
        return out

class BioPax_SmallMolecule(BioPax_ElementBase):
    type = "SmallMolecule"

    def process(self):
        node_label = None
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            node_label = name.dst.replace("\n", " ")

        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "name", get_type=False):
            node_label = name.dst.replace("\n", " ")


        xref_list = []
        for ref in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "entityReference", get_type=False):
            for xref in self.pax.query(src=ref.dst, predicate=BIOPAX_BASE + "xref", get_type=False):
                db = None
                db_id = None
                for rel in self.pax.query(src=xref.dst):
                    if rel.predicate == BIOPAX_BASE + "db":
                        db = rel.dst
                    if rel.predicate == BIOPAX_BASE + "id":
                        db_id = rel.dst
                db_xref = "%s:%s" % (db, db_id)
                xref_list.append(db_xref)

            for name in self.pax.query(src=ref.dst, predicate=BIOPAX_BASE + "name", get_type=False):
                node_label = name.dst.replace("\n", " ")



        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})
        data = {'db_xref' : xref_list, 'type' : self.type}
        if node_label is not None:
            data['label'] = node_label
        out.add_node( self.node, data, is_input=True, is_output=True)
        return out


class BioPax_PhysicalEntity(BioPax_ElementBase):
    type = "PhysicalEntity"

    def process(self):
        node_label = None
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "name", get_type=False):
            node_label = name.dst.replace("\n", " ")

        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})
        data = {'type' : self.type}
        if node_label is not None:
            data['label'] = node_label
        out.add_node(self.node, data, is_input=True, is_output=True )
        return out

class BioPax_Rna(BioPax_ElementBase):
    type = "Rna"

    def process(self):
        entity_ref = None
        entity_type = None
        for ref in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "entityReference"):
            entity_ref = ref.dst
            entity_type = ref.dst_type
        if entity_ref is not None:
            return self.process_child(entity_ref, entity_type)

class BioPax_RnaReference(BioPax_ElementBase):
    type = "Rna"

    def process(self):
        
        xref_list = []
        for xref in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "xref", get_type=False):
            db = None
            db_id = None
            for rel in self.pax.query(src=xref.dst):
                if rel.predicate == BIOPAX_BASE + "db":
                    db = rel.dst
                if rel.predicate == BIOPAX_BASE + "id":
                    db_id = rel.dst
            db_xref = "%s:%s" % (db, db_id)
            xref_list.append(db_xref)

        aliases = []
        node_label = None
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            node_label = name.dst.replace("\n", " ")
            aliases.append(node_label)

        for prot_name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "name", get_type=False):
            node_label = prot_name.dst
            aliases.append(node_label)


        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})
        data = {'db_xref' : xref_list, 'type' : 'rna'}
        if len(aliases) > 1:
            data['aliases'] = aliases
        if node_label is not None:
            data['label'] = node_label
        out.add_node(self.node, data, is_input=True, is_output=True )
        return out



class BioPax_Dna(BioPax_ElementBase):
    type = "Dna"

    def process(self):
        entity_ref = None
        entity_type = None
        for ref in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "entityReference"):
            entity_ref = ref.dst
            entity_type = ref.dst_type
        if entity_ref is not None:
            return self.process_child(entity_ref, entity_type)

class BioPax_DnaReference(BioPax_ElementBase):
    type = "DnaReference"

    def process(self):
        
        xref_list = []
        for xref in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "xref", get_type=False):
            db = None
            db_id = None
            for rel in self.pax.query(src=xref.dst):
                if rel.predicate == BIOPAX_BASE + "db":
                    db = rel.dst
                if rel.predicate == BIOPAX_BASE + "id":
                    db_id = rel.dst
            db_xref = "%s:%s" % (db, db_id)
            xref_list.append(db_xref)

        aliases = []
        node_label = None
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            node_label = name.dst.replace("\n", " ")
            aliases.append(node_label)

        for prot_name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "name", get_type=False):
            node_label = prot_name.dst
            aliases.append(node_label)


        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})
        data = {'db_xref' : xref_list, 'type' : 'dna'}
        if len(aliases) > 1:
            data['aliases'] = aliases
        if node_label is not None:
            data['label'] = node_label
        out.add_node(self.node, data, is_input=True, is_output=True )
        return out


class BioPax_Protein(BioPax_ElementBase):
    type = "Protein"

    def process(self):

        entity_ref = None
        entity_type = None
        for ref in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "entityReference"):
            entity_ref = ref.dst
            entity_type = ref.dst_type
        if entity_ref is not None:
            return self.process_child(entity_ref, entity_type)


class BioPax_ProteinReference(BioPax_ElementBase):
    type = "ProteinReference"

    def process(self):

        xref_list = []
        for xref in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "xref", get_type=False):
            db = None
            db_id = None
            for rel in self.pax.query(src=xref.dst):
                if rel.predicate == BIOPAX_BASE + "db":
                    db = rel.dst
                if rel.predicate == BIOPAX_BASE + "id":
                    db_id = rel.dst
            db_xref = "%s:%s" % (db, db_id)
            xref_list.append(db_xref)

        aliases = []
        node_label = None
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            node_label = name.dst.replace("\n", " ")
            aliases.append(node_label)

        for prot_name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "name", get_type=False):
            node_label = prot_name.dst
            aliases.append(node_label)


        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})
        data = {'db_xref' : xref_list, 'type' : 'protein'}
        if len(aliases) > 1:
            data['aliases'] = aliases
        if node_label is not None:
            data['label'] = node_label
        out.add_node(self.node, data, is_input=True, is_output=True )
        return out


class BioPax_Catalysis(BioPax_ElementBase):
    type = "Catalysis"

    def process(self):

        node_label = None
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            node_label = name.dst.replace("\n", " ")

        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})
        controller = []
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "controller"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_input=True)
            controller.append( c_info.dst )
            
        controlled = []
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "controlled"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_output=True)
            controlled.append( c_info.dst )


        interaction = "-a>"
        for control in self.pax.query(src=self.node, predicate=BIOPAX_BASE + 'controlType', get_type=False):
            interaction = control.dst 
        
        for l in controller:
            for r in controlled:
                out.add_edge( l, r, {'interaction' : interaction, 'class' : self.type, 'src_url' : self.node} )
        return out

class BioPax_BiochemicalReaction(BioPax_ElementBase):
    type = "BiochemicalReaction"

    def process(self):
        left = []

        node_label = None
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            node_label = name.dst.replace("\n", " ")

        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "left"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_input=True)
            left.append( c_info.dst )
                
        right = []
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "right"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_output=True)           
            right.append( c_info.dst )

        interaction = "-a>"
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "interactionType", get_type=False):
            interaction = c_info.dst
            for d_info in self.pax.query(src=c_info.dst, predicate=BIOPAX_BASE + "term", get_type=False):
                interaction = d_info.dst

        for l in left:
            for r in right:
                out.add_edge( l, r, {'interaction' : interaction, 'class' : self.type, 'src_url' : self.node} )
        return out


class BioPax_MolecularInteraction(BioPax_ElementBase):
    type = "MolecularInteraction"

    def process(self):

        node_label = None
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            node_label = name.dst.replace("\n", " ")

        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})
        part = []
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "participant"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem)
            part.append( c_info.dst )
        for p1 in part:
            for p2 in part:
                if p1 != p2:
                    out.add_edge(p1, p2, {'interaction' : "--"} ) 
        return out

class BioPax_Transport(BioPax_ElementBase):
    type = "Transport"
    def process(self):
        left = []

        node_label = None
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            node_label = name.dst.replace("\n", " ")

        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "left"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_input=True)
            left.append( c_info.dst )
                
        right = []
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "right"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_output=True)           
            right.append( c_info.dst )

        interaction = "transport"

        for l in left:
            for r in right:
                out.add_edge( l, r, {'interaction' : interaction, 'class' : self.type, 'src_url' : self.node} )
        return out


class BioPax_Degradation(BioPax_ElementBase):
    type = "Degradation"
    def process(self):
        left = []

        node_label = None
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            node_label = name.dst.replace("\n", " ")

        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "left"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_input=True)
            left.append( c_info.dst )

        """     
        right = []
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "right"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_output=True)           
            right.append( c_info.dst )
        interaction = "transport"

        for l in left:
            for r in right:
                out.add_edge( l, r, {'interaction' : interaction, 'class' : self.type, 'src_url' : self.node} )
        """
        return out

class BioPax_Control(BioPax_ElementBase):
    type = "Control"
    def process(self):
        controller = []

        node_label = None
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            node_label = name.dst.replace("\n", " ")

        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "controller"):
            self.debug("call %s %s" % (c_info.dst, c_info.dst_type))
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_input=True)
            controller.append( c_info.dst )
            
        controlled = []
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "controlled"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_output=True)
            controlled.append( c_info.dst )

        interaction = "-a>"
        for control in self.pax.query(src=self.node, predicate=BIOPAX_BASE + 'controlType', get_type=False):
            interaction = control.dst 
        
        for l in controller:
            for r in controlled:
                out.add_edge( l, r, {'interaction' : interaction, 'class' : self.type, 'src_url' : self.node} )

        return out

class BioPax_TemplateReactionRegulation(BioPax_ElementBase):
    type = "TemplateReactionRegulation"
    def process(self):
        #raise Exception("Missing Code: %s" % (self.type))
        return None

class BioPax_TemplateReaction(BioPax_ElementBase):
    type = "TemplateReaction"
    def process(self):
        #raise Exception("Missing Code: %s" % (self.type))
        return None

class BioPax_TransportWithBiochemicalReaction(BioPax_ElementBase):
    type = "TransportWithBiochemicalReaction"
    def process(self):

        node_label = None
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            node_label = name.dst.replace("\n", " ")
        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})

        left = []
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "left"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_input=True)
            left.append( c_info.dst )
                
        right = []
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "right"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_output=True)           
            right.append( c_info.dst )

        interaction = "transport"

        for l in left:
            for r in right:
                out.add_edge( l, r, {'interaction' : interaction, 'class' : self.type, 'src_url' : self.node} )
        return out

class BioPax_Interaction(BioPax_ElementBase):
    type = "Interaction"
    def process(self):
        #raise Exception("Missing Code: %s" % (self.type))
        return None

class BioPax_Complex(BioPax_ElementBase):
    type = "Complex"

    def process(self):
        node_label = None
        #bug, using displayName may cause errors for complexes
        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "displayName", get_type=False):
            node_label = name.dst.replace("\n", " ")

        for name in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "name", get_type=False):
            node_label = name.dst.replace("\n", " ")

        out = Subnet({'name' : node_label, 'url' : self.node, 'type' : self.type})
        data = {"type" : "complex"}
        if node_label is not None:
            data['label'] = node_label

        out.add_node(self.node, data, is_output=True, is_input=True)

        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "component"):
            c = self.process_child(c_info.dst, c_info.dst_type)
            if c is not None:
                self.debug("ComponentLink %s %s" % (c_info.dst, c))
                out.add_node(c_info.dst, c)
                out.add_edge(c_info.dst, self.node, {"interaction" : "component>", 'class' : self.type, 'src_url' : self.node} )
        return out


class BioPax_ComplexAssembly(BioPax_ElementBase):
    type = "ComplexAssembly"

    def process(self):
        
        """
        out = Subnet()
        left = []
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "left"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem)
            left.append( c_info.dst )
        
        right = []
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "right"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            out.add_node(c_info.dst, elem, is_input=True, is_output=True)
            right.append( c_info.dst )

        if len(right) != 1:
            self.debug("Weird product set")
        
        for l in left:
            for r in right:
                out.add_edge(l, r, {"interaction" : "component>", 'class' : self.type, 'src_url' : self.node} )        
        return out
        """
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "right"):
            elem = self.process_child(c_info.dst, c_info.dst_type)
            return elem



class BioPax_Blank(BioPax_ElementBase):
    type = "Blank"

    def process(self):
        out = Subnet()
        name = None
        for c_info in self.pax.query(src=self.node, predicate=BIOPAX_BASE + "name", get_type=False):
            name = c_info.dst
        out.add_node( self.node, {'name' : name, 'url' : self.node}, is_input=True, is_output=True )
        return out

element_mapping = {
    BIOPAX_BASE + "Protein" : BioPax_Protein,
    BIOPAX_BASE + "ProteinReference" : BioPax_ProteinReference,
    BIOPAX_BASE + "SmallMolecule" : BioPax_SmallMolecule,
    BIOPAX_BASE + "Complex" : BioPax_Complex,
    BIOPAX_BASE + "ComplexAssembly" : BioPax_ComplexAssembly,
    BIOPAX_BASE + "PhysicalEntity" : BioPax_PhysicalEntity,
    BIOPAX_BASE + "Rna" : BioPax_Rna,
    BIOPAX_BASE + "RnaReference" : BioPax_RnaReference,
    BIOPAX_BASE + "Dna" : BioPax_Dna,
    BIOPAX_BASE + "DnaReference" : BioPax_DnaReference,
    BIOPAX_BASE + "Catalysis" : BioPax_Catalysis,
    BIOPAX_BASE + "BiochemicalReaction" : BioPax_BiochemicalReaction,
    BIOPAX_BASE + "MolecularInteraction" : BioPax_MolecularInteraction,
    BIOPAX_BASE + "Transport" : BioPax_Transport,
    BIOPAX_BASE + "Degradation" : BioPax_Degradation,
    BIOPAX_BASE + "Pathway" : BioPax_Pathway,
    BIOPAX_BASE + "Control" : BioPax_Control,
    BIOPAX_BASE + "TemplateReactionRegulation" : BioPax_TemplateReactionRegulation,
    BIOPAX_BASE + "TemplateReaction" : BioPax_TemplateReaction,
    BIOPAX_BASE + "TransportWithBiochemicalReaction" : BioPax_TransportWithBiochemicalReaction,
    BIOPAX_BASE + "Interaction" : BioPax_Interaction,
    None : BioPax_Blank
}

class Subnet:
    def __init__(self, meta=None):
        self.nodes = {}
        self.edges = []
        self.output_node = None
        self.input_node = None
        self.meta = meta

    def add_node(self, name, data, is_input=False, is_output=False):
        self.nodes[name] = data
        if is_input:
            self.input_node = name
        if is_output:
            self.output_node = name

    def add_edge(self, src, dst, data):
        self.edges.append( (src, dst, data) )

    def to_graph(self, graph, visited={}):
        out = ""
        for a in self.nodes:
            if isinstance(self.nodes[a], Subnet):
                if a not in visited:
                    self.nodes[a].to_graph(graph, visited)
                #visited[a] = True
            else:
                if isinstance(self.nodes[a], dict):
                    data = dict(self.nodes[a])
                else:
                    data = {}
                data['src_url'] = a
                graph.add_node(a, attr_dict=data)
        for a in self.edges:
            if isinstance(self.nodes[a[0]], Subnet):
                src = self.nodes[a[0]].get_input()
            else:
                src = a[0]
            if isinstance(self.nodes[a[1]], Subnet):
                dst = self.nodes[a[1]].get_output()
            else:
                dst = a[1]
            if src in graph.edge and dst not in graph.edge[src]:
                graph.add_edge( src, dst, attr_dict=a[2] )
    
    def get_input(self):
        if self.input_node is None:
            return None
        if isinstance(self.nodes[self.input_node], Subnet):
            return self.nodes[self.input_node].get_input()
        return self.input_node

    def get_output(self):
        if self.output_node is None:
            return None
        if isinstance(self.nodes[self.output_node], Subnet):
            return self.nodes[self.output_node].get_output()
        return self.output_node

    def __str__(self):
        return "{%s}" % (",".join( [ str(self.nodes[a]) for a in self.nodes ] ))


class BioPax:

    def __init__(self):
        pass

    def pathways(self):
        out = {}
        for a in self.query(src_type=BIOPAX_BASE + "Pathway", predicate=BIOPAX_BASE + "name"):
            out[a.dst] = True
        for a in self.query(src_type=BIOPAX_BASE + "Pathway", predicate=BIOPAX_BASE + "displayName"):
            out[a.dst] = True

        return out.keys()

    def toNet(self, pathways=None):
        start_set = {}
        for a in self.query(src_type=BIOPAX_BASE + "Pathway", predicate=BIOPAX_BASE + "name"):
            if pathways is None or a.dst in pathways:
                start_set[a.src] = True

        for a in self.query(src_type=BIOPAX_BASE + "Pathway", predicate=BIOPAX_BASE + "displayName"):
            if pathways is None or a.dst in pathways:
                start_set[a.src] = True

        for pathway in start_set:
            p_elem = BioPax_Pathway(self, pathway, [])
            gr = p_elem.process()
            yield gr
    


class BioPaxFile(BioPax):
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

    def get_node_type(self, node):
        return self.graph.get(node, {}).get(TYPE_PRED, [None])[0]
    
    def str_compare(self, str1, str2, check_case):
        if check_case:
            return str1 == str2
        #print str1, str2
        return str1.lower() == str2.lower()

    def query(self, src_type=None, src=None, predicate=None, dst_type=None, dst=None, check_case=True, get_type=True):
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


class BioPaxSparql(BioPax):

    def __init__(self, url):
        self.url = url
        self.sparql = sparql.Service(url)

    def get_node_type(self, node):
        query = "select ?t { <%s> a ?t }" % (node)
        results = self.sparql.query(query)
        out = list(results.fetchone())
        if len(out):
            return unicode(out[0][0])
        return None

    def query(self, src_type=None, src=None, predicate=None, dst_type=None, dst=None, check_case=True, get_type=True):
        if src:
            src_label = "<%s>" % (src)
        else:
            src_label = "?s"

        if dst:
            dst_label = "<%s>" % (dst)
        else:
            dst_label = "?d"

        if predicate:
            predicate_label = "<%s>" % (predicate)
        else:
            predicate_label = "?p"

        if src_type:
            src_type_select = "%s a <%s> ." % (src_label, src_type)
        else:
            src_type_select = ""
        if dst_type:
            dst_type_select = "%s a <%s> ." % (dst_label, dst_type)
        else:
            dst_type_select = ""



        link_select = "%s %s %s ." % (src_label, predicate_label, dst_label)

        query = """SELECT *
        WHERE {
            %s
            %s
            %s
        }""" % (src_type_select, dst_type_select, link_select)
        try:
            results = self.sparql.query(query)
            for row in results.fetchall():
                rsrc = None
                rsrc_type = None
                rdst = None
                rdst_type = None
                rpred = None
                for i, v in enumerate(results.variables):
                    if v == 's':
                        rsrc = unicode(row[i])
                        if isinstance(row[i], sparql.IRI) and get_type:
                            rsrc_type = self.get_node_type(rsrc)
                    if v == 'd':
                        rdst = unicode(row[i])
                        if isinstance(row[i], sparql.IRI) and get_type:
                            rdst_type = self.get_node_type(rdst)
                    if v == 'p':
                        rpred = unicode(row[i])
                link = Link(rsrc_type, rsrc, rpred, rdst_type, rdst)
                yield link
        except Exception, e:
            raise Exception("Query Fail: %s : %s " % (query, str(e)) )



def write_biopax(gr, handle):

    gene = rdflib.Namespace("http://ucsc.edu/gene#")
    xsd=rdflib.Namespace("http://www.w3.org/2001/XMLSchema#")
    owl=rdflib.Namespace("http://www.w3.org/2002/07/owl#")
    rdf=rdflib.Namespace("http://www.w3.org/1999/02/22-rdf-syntax-ns#")
    bp=rdflib.Namespace("http://www.biopax.org/release/biopax-level3.owl#")

    pax = rdflib.Graph()

    node_map = {}
    for n in gr.node:
        if gr.node[n]['type'] == 'family':
            pnode = gene[n] #rdflib.BNode()            
            node_map[n] = pnode
            pax.add( (pnode, rdf.type, bp.Protein) )
            pax.add( (pnode, bp.displayName, rdflib.Literal(n)) )

            for src, dst, data in gr.edges( (None, n), data=True):
                prefnode = gene[src + "_ref"]
                pax.add( (pnode, bp.memberEntityReference, prefnode) )

    for n in gr.node:
        if gr.node[n]['type'] == 'protein':
            pnode = gene[n] #rdflib.BNode()
            node_map[n] = pnode
            pax.add( (pnode, rdf.type, bp.Protein) )
            pax.add( (pnode, bp.displayName, rdflib.Literal(n)) )
            prefnode = gene[n + "_ref"]
            pax.add( (prefnode, rdf.type, bp.ProteinReference) )
            pax.add( (pnode, bp.entityReference, prefnode) )
        elif gr.node[n]['type'] == 'complex':
            pnode = gene[n] #rdflib.BNode()
            node_map[n] = pnode
            pax.add( (pnode, rdf.type, bp.Complex) )
            pax.add( (pnode, bp.displayName, rdflib.Literal(n)) )
            prefnode = gene[n + "_ref"]
            pax.add( (prefnode, rdf.type, bp.ProteinReference) )
            pax.add( (pnode, bp.entityReference, prefnode) )
        else:
            pnode = gene[n] #rdflib.BNode()
            node_map[n] = pnode
            pax.add( (pnode, rdf.type, bp.Protein) )
            pax.add( (pnode, bp.displayName, rdflib.Literal(n)) )
        

    i = 0
    reaction_list = {}
    for src, dst, data in gr.edges(data=True):
        reaction = rdflib.URIRef("reaction_%d" % (i) ) # rdflib.BNode()
        i += 1
        pax.add( (reaction, rdf.type, bp.BiochemicalReaction ) )
        pax.add( (reaction, rdf.type, bp.BiochemicalReaction ) )
        pax.add( (reaction, bp.left, node_map[src] ) )
        pax.add( (reaction, bp.right, node_map[dst] ) )
        pax.add( (reaction, bp.displayName, rdflib.Literal("%s - %s" % (src, dst)) ))
        for spid in gr.node[src].get('pid', [None]):
            for dpid in gr.node[dst].get('pid', [None]):
                if spid == dpid:
                    if spid not in reaction_list:
                        reaction_list[spid] = []
                    reaction_list[spid].append(reaction)

    """
    pathways = {}
    for n in gr.node:
        if 'pid' not in gr.node[n]:
            if None not in pathways:
                pathways[None] = rdflib.BNode()
                pax.add( (pathways[None], rdf.type, bp.Pathway) )
        else:
            for pid in gr.node[n]['pid']:
                if pid not in pathways:
                    pathways[pid] = rdflib.BNode()
                    pax.add( (pathways[pid], rdf.type, bp.Pathway) )
    """
    i = 0
    for pid in reaction_list:
        pathnode = rdflib.URIRef("path_%d" % (i))
        i += 1
        for reaction in reaction_list[pid]:
            pax.add( (pnode, bp.pathwayComponent, reaction) )


    handle.write(pax.serialize(format="pretty-xml"))

