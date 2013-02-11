#!/usr/local/bin/python2.6

import mPathway, sys, re, collections, itertools

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p","--pathway_file",type="string",dest="pathway_file", action="store", help="UCSC Pathway File")
parser.add_option("--c_output",type="string",dest="component_output", action="store", help="Output file for Component Map")
parser.add_option("--flattened",type="string",dest="flattened", action="store", help="Join genes with all pathway links for complexes and families they belong to. Print that network with just proteins")
(options, args) = parser.parse_args()


# nodes:
#	name -> type
# interactions:
#	name -> interacting nodes
nodes, Interactions, Proteins = mPathway.rPathway(options.pathway_file, reverse = False, retProteins = True)

rev_nodes, revInteractions = mPathway.rPathway(options.pathway_file, reverse = True, retProteins = False)
# maps complex strings to the components in each
componentMap = mPathway.getComponentMap(rev_nodes, revInteractions)

complexRE = re.compile(".*\((complex|family)\).*")
abstractRE = re.compile(".*\(abstract\).*")

# get the constituents of a complex
def getGenes(complex_str, componentMap, depth):
	all_genes = []

	# give up
	if complex_str not in componentMap:
		return []

	for component in componentMap[complex_str]:
		if complexRE.match(component):
			if component in componentMap and depth < 4:
				genes = getGenes(component, componentMap, depth+1)
				for gene in genes:
					all_genes.append(gene)
		else:
			all_genes.append(component)

	return all_genes

def addEdges(abs_source, abs_target, interaction, map, componentMap):
	"""
	Add the source/target interaction, break them up into
	component parts with the component map and add them to the 
	edge map
	map[ (sourceGene, interactionType) ] = targetGene
	Returns: An updated edge mapping
	"""

	# use to component map to find all the source genes
	sourceGenes = []
	if complexRE.match(abs_source):
		sourceGenes = getGenes(abs_source, componentMap, 1)
	else:
		sourceGenes.append(abs_source)

	# find target genes with the component map
	targetGenes = []
	if complexRE.match(abs_target):
		for gene in getGenes(abs_target, componentMap, 1):
			targetGenes.append(gene)
	else:
		targetGenes.append(abs_target)

	# connect these sink genes to the source gene
	for sr in sourceGenes:
		if (sr,interaction) not in map:
			map[(sr,interaction)] = []
		for target in targetGenes:
			# add target interactions: don't duplicate
			if target not in map[(sr,interaction)]:
				map[(sr,interaction)].append(target)

	return map

# print out a 2-column interactions file of simple PPIs
# protein -> protein
ppi_edges = {}
tf_edges = {}
all_edges = {}

# for each node, find all interactions, and use the component map to refine them to genes
# if the interaction is a -t*, add it to the TF edge list, otherwise if -p* add it to the ppi edge list
for node in Interactions:

	# complex, -t> or -t| interaction type
	type = nodes[node]
		
	if node not in Interactions:
		continue

	for inode in Interactions[node]:

		# type is a ';' separated list of the interactions
		# for these node pairs -- may be member, component
		# -t or -a or any mix of those

		type = Interactions[node][inode]
		all_interactions = type.split(";")
		for int in all_interactions:
			# add to the all-edge map
			all_edges = addEdges(node, inode, int, all_edges, componentMap) 

# test: try writing all edges out
all_out = open(options.flattened, 'w')
for (source, interaction) in all_edges:

	for sink in all_edges[(source, interaction)]:
		# remove self links:
		if source == sink:
			continue
		all_out.write(source+"\t"+sink+"\t"+interaction+"\n")

all_out.close()

c_out = open(options.component_output, 'w')
for component in componentMap:
	c_out.write(component+"\t"+" ".join(componentMap[component])+"\n")
c_out.close()
