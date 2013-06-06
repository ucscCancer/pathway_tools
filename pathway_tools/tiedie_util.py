## Utility function file for hotLink program

import re, math, os, sys, operator, random

def parseHeats(file):
	heats = {}
	signs = {}
	for line in open(file, 'r'):
		parts = line.rstrip().split("\t")
		if len(parts) > 2:
			prot, heat, sign = line.rstrip().split("\t")
			heats[prot] = float(heat)
			signs[prot] = sign
		else:
			heats[parts[0]] = float(parts[1])

	return (heats, signs)

def edgelist2nodes(list):

	nodes = set()
	for (source,i, target) in list:
		nodes.add(source)
		nodes.add(target)

	return nodes
	
def classifyInteraction(i):
	componentRE = re.compile("^-?component>$")
	activatingRE = re.compile("^-?(\S)>$")
	inactivatingRE = re.compile("^-?(\S)\|$")
	rewiredAC = re.compile("^-?REWIRED>$")
	rewiredIN = re.compile("^-?REWIRED\|$")

	if componentRE.match(i):
		return (0, "component")
	elif activatingRE.match(i):
		type = activatingRE.match(i)
		return (1, type.group(1))
	elif inactivatingRE.match(i):
		type = inactivatingRE.match(i)
		return (2, type.group(1))
	elif rewiredAC.match(i):
		type = "REWIRED"
		return (1, type.group(1))
	elif rewiredIN.match(i):
		type = "REWIRED"
		return (2, type.group(1))
	else:
		# default to activating links for HPRD
		return (1, "INTERACTS")

def getOutDegrees(network):

	outDegrees = {}
	for s in network:
		outDegrees[s] = len(network[s]) 
		for (i, t) in network[s]:
			if t not in outDegrees:
				outDegrees[t] = 0

	return outDegrees
	
def edges2degrees(edges):

	nodes = {}
	for (s,t) in edges:
		if s not in nodes:
			nodes[s] = {}
		if t not in nodes:
			nodes[t] = {}

		nodes[s][t] = 1
		nodes[t][s] = 1

	sizes = {}	
	for n in nodes:
		sizes[n] = len(nodes[n])

	return sizes

def isRewired(i):
	rewiredRE = re.compile(".*REWIRED.*")
	rewiredComponentRE = re.compile(".*\-component.*")

	if rewiredRE.match(i):
		return True
	elif rewiredComponentRE.match(i):
		return True

	return False

# do a depth-first search by following directional links 
# until we hit another source
# find edges 
def searchDFS(source, action, discovered, linker_nodes, target_set, net, gene_states, transcriptional_signs, depth, truePaths, falsePaths, falsePathStatus):

	if depth == 0:
		return

	if source not in net:
		return

	for (interaction, target) in net[source]:

		# if we arrived here through a bad link, any continued path is counted as false
		pathStatus_ThisTarget = falsePathStatus
		(i_type, post_t_type) = classifyInteraction(interaction)
		if isRewired(interaction):
			pathStatus_ThisTarget = True

		# don't follow component links
		if i_type == 0:
			continue

		# activating nodes keep the signal, inactivating nodes reverse the signal
		action_this_target = None
		if i_type == 1:
			action_this_target = action
		elif i_type == 2:
			action_this_target = -action

		# for transcriptional states: the expression activity is what we want to measure 
		this_state = None
		if target in gene_states:	
			this_state = gene_states[target]
		#if post_t_type == "t":
		# this depends on weather we monitor the activities of downstream genes, or just the transcription
		# leave commented out for the former
		#	if target not in transcriptional_signs:
		#		continue
		#	this_state = transcriptional_signs[target] 

		# we hit a target that has a matching action/signal from the original source
		if (target in gene_states) and (target in target_set) \
			and (action_this_target == this_state):

			for (s,i,t) in linker_nodes:
				discovered.add((s,i,t))
			discovered.add((source, interaction, target))
			linker_nodes = set()
			new_linkers = set()
			# and keep going

			# add this to our TP or FP score, depending on the path
			if pathStatus_ThisTarget == True:
				falsePaths.append(target)
			else:
				truePaths.append(target)

		# search the target, but with any previous linkers	
		else:
			new_linkers = set()
			new_linkers.add((source, interaction, target))
			new_linkers = new_linkers.union(linker_nodes)	

		# if we come from a transcriptionally activating link, this cuts the cycle. Gene must
		# be upregulated 
		if post_t_type == "t":
			continue

		# add this link and keep searching from the target
		searchDFS(target, action_this_target, discovered, new_linkers, target_set, net, gene_states, transcriptional_signs, depth-1, truePaths, falsePaths, pathStatus_ThisTarget)

# do a depth-first search by following directional links 
# until we hit another source
# find edges 
def searchDFS_NoAction(source, discovered, linker_nodes, target_set, net, gene_states, depth, truePaths, falsePaths, falsePathStatus):

	if depth == 0:
		return

	if source not in net:
		return

	for (interaction, target) in net[source]:

		# if we arrived here through a bad link, any continued path is counted as false
		pathStatus_ThisTarget = falsePathStatus
		(i_type, post_t_type) = classifyInteraction(interaction)
		if isRewired(interaction):
			pathStatus_ThisTarget = True


		# we hit a target that has a matching action/signal from the original source
		if (target in gene_states) and (target in target_set):

			for (s,i,t) in linker_nodes:
				discovered.add((s,i,t))
			discovered.add((source, interaction, target))
			linker_nodes = set()
			new_linkers = set()
			# and keep going

			# add this to our TP or FP score, depending on the path
			if pathStatus_ThisTarget == True:
				falsePaths.append(target)
			else:
				truePaths.append(target)

		# search the target, but with any previous linkers	
		else:
			new_linkers = set()
			new_linkers.add((source, interaction, target))
			new_linkers = new_linkers.union(linker_nodes)	

		# if we come from a transcriptionally activating link, this cuts the cycle. Gene must
		# be upregulated 
		if post_t_type == "t":
			continue

		# add this link and keep searching from the target
		searchDFS_NoAction(target, discovered, new_linkers, target_set, net, gene_states, depth-1, truePaths, falsePaths, pathStatus_ThisTarget)

# encapsulates the guesses of the biological effects of
# each events. 
def classifyState(up_signs, down_signs):

	c = {}
	t_states = {}
	# The order matters here: 
	for gene in down_signs:
		if down_signs[gene] == "+":
			c[gene] = 1
			t_states[gene] = 1
		else:
			c[gene] = -1
			t_states[gene] = -1

	# The order matters here: 
	for gene in up_signs:
		if up_signs[gene] == "+":
			c[gene] = 1
		else:
			c[gene] = -1

	return (c, t_states)


def getNetworkNodes(network):
	nodes = set()
	for s in network:
		nodes.add(s)
		for t in network[s]:
			nodes.add(t)
	return nodes

# build an index, source to targets fro the directed graph
def parseNet(network):

	net = {}
	fh = open(network, 'r')
	for line in fh:

		parts = line.rstrip().split("\t")
		source = parts[0]
		interaction = parts[1]
		target = parts[2]

		(i,t) = classifyInteraction(interaction)
		# skip component links
		#if i == 0:
		#	continue

		if source not in net:
			net[source] = set()

		net[source].add((interaction, target))
		if t == "INTERACTS":
			# these are bi-directional links: we need to add duplicate edges
			if target not in net:
				net[target] = set()
			net[target].add((interaction, source))

	fh.close()
	return net

def findLinkerCutoff(source_set, target_set, up_heat_diffused, down_heat_diffused, size):

	if down_heat_diffused is None:
		cutoff, score = findLinkerCutoffSingle(source_set, up_heat_diffused, size)		
	else:
		try:
			cutoff, score = findLinkerCutoffMulti(source_set, target_set, up_heat_diffused, down_heat_diffused, size)		
		except:
			return (0,0)

	return (cutoff, score)

def findLinkerCutoffSingle(source_set, up_heat_diffused, size):

	source_set = set(source_set)
	up_sorted = sorted(up_heat_diffused, key=up_heat_diffused.get, reverse=True)

	EPSILON = 0.0001

	# we want to find this many diffused genes in total
	target_size = size*len(source_set)

	i = 1
	cutoff = None
	for (gene, heat) in sorted(up_heat_diffused.iteritems(), key=operator.itemgetter(1), reverse=True):

		cutoff = heat+EPSILON
		diffused_set = set(up_sorted[0:i])
		if len(diffused_set.difference(source_set)) > target_size:
			break
		i += 1

	return (cutoff, 0)	


'''			
def findLinkerCutoffMulti(source_set, target_set, up_heat_diffused, down_heat_diffused, size):

	if size == 0:
		return (1000000, 0)

	target_set = set(target_set)
	source_set = set(source_set)
	up_sorted = sorted(up_heat_diffused, key=up_heat_diffused.get, reverse=True)
	down_sorted = sorted(down_heat_diffused, key=down_heat_diffused.get, reverse=True)
	scores = {}
	best = None
	last = 0

	foundPos = False

	# rank the min heats, and decrement the cutoff, adding an additional gene at each step
	EPSILON = 0.0001

	f, min_heats = filterLinkers(up_heat_diffused,down_heat_diffused,1)
	# at a cutoff just below, we'll get the gene 
	for cutoff in [h-EPSILON for (l,h) in sorted(min_heats.iteritems(), key=operator.itemgetter(1), reverse=True)]:
		score, size_corr = scoreLinkers(up_heat_diffused, up_sorted, down_heat_diffused, down_sorted, source_set, target_set, cutoff, size)
		abs_score = score - abs(size_corr)
		scores[cutoff] = abs_score

		if abs_score > 0:
			foundPos = True

		# find the first local max
		if abs_score >= last:
			last = abs_score
		elif best is None:
			best = (cutoff, last)
		elif abs_score < 0 and foundPos:
			break

	best_cutoff, best_score = sorted(scores.iteritems(), key=operator.itemgetter(1),reverse=True)[0]

	return (best_cutoff, best_score)
'''

def findLinkerCutoffMulti(source_set, target_set, up_heat_diffused, down_heat_diffused, size):

	if size == 0:
		return (1000000, 0)

	target_set = set(target_set)
	source_set = set(source_set)
	up_sorted = sorted(up_heat_diffused, key=up_heat_diffused.get, reverse=True)
	down_sorted = sorted(down_heat_diffused, key=down_heat_diffused.get, reverse=True)
	scores = {}
	best = None
	last = 0

	foundPos = False

	# rank the min heats, and decrement the cutoff, adding an additional gene at each step
	EPSILON = 0.0001

	f, min_heats = filterLinkers(up_heat_diffused,down_heat_diffused,1)
	# at a cutoff just below, we'll get the gene 
	for cutoff in [h-EPSILON for (l,h) in sorted(min_heats.iteritems(), key=operator.itemgetter(1), reverse=True)]:
		score, size_frac = scoreLinkers(up_heat_diffused, up_sorted, down_heat_diffused, down_sorted, source_set, target_set, cutoff, size)

		if size_frac > 1:
			return (cutoff, score)

'''	
def scoreLinkers(heats1, sorted1, heats2, sorted2, sourceSet, targetSet, cutoff, size):
	"""
		Get linkers greater than this cutoff according to reverse-sorted list..
	""" 
	filtered_h1 = {}
	for l in sorted1:
		s = heats1[l]
		if s < cutoff:
			break

		filtered_h1[l] = s

	filtered_h2 = {}
	for l in sorted2:
		s = heats2[l]
		if s < cutoff:
			break

		filtered_h2[l] = s

	f1 = set(filtered_h1)
	f2 = set(filtered_h2)

	union = f1.union(f2)
	intersection = f1.intersection(f2)
	connecting = intersection.difference(sourceSet).difference(targetSet)
	score = len(connecting)/float(len(union))
	size_corr = None
	try:
		size_corr = (len(connecting)/size)/float(len(sourceSet.union(targetSet)))
	except:
		return (score, 1)
	
	
	return (score, (1-abs(size_corr)))
'''

def scoreLinkers(heats1, sorted1, heats2, sorted2, sourceSet, targetSet, cutoff, size):
	"""
		Get linkers greater than this cutoff according to reverse-sorted list..
	""" 
	filtered_h1 = {}
	for l in sorted1:
		s = heats1[l]
		if s < cutoff:
			break

		filtered_h1[l] = s

	filtered_h2 = {}
	for l in sorted2:
		s = heats2[l]
		if s < cutoff:
			break

		filtered_h2[l] = s

	f1 = set(filtered_h1)
	f2 = set(filtered_h2)

	union = f1.union(f2)
	intersection = f1.intersection(f2)
	connecting = intersection.difference(sourceSet).difference(targetSet)
	score = len(connecting)/float(len(union))
	size_frac = (len(connecting)/float(len(sourceSet.union(targetSet))))/float(size)
	
	return (score, size_frac)

def scoreLinkersTest(input_heats, min_heats, cutoff, size):
	"""
		Get linkers greater than this cutoff according to reverse-sorted list..
	""" 

	# get the full input set
	all_inputs = set()
	for name in input_heats:
		all_inputs = all_inputs.union(input_heats[name].keys())

	# get the full set of linkers
	all_linkers = set()
	for (gene, heat) in sorted(min_heats.iteritems(), key=operator.itemgetter(1), reverse=True):
		if heat < cutoff:
			break
		all_linkers.add(gene)

	all_genes = all_inputs.union(all_linkers)
	connecting = all_linkers.difference(all_inputs)
	score = len(connecting)/float(len(all_linkers))
	size_frac = (len(connecting)/float(len(all_inputs)))/float(size)
	
	return (score, size_frac)

def getMinHeats(diffused):

	mins = {}
	for file in diffused:
		for (gene, heat) in diffused[file].iteritems():
			if gene in mins:
				if mins[gene] > heat:
					mins[gene] = heat
			else:
				mins[gene] = heat

	return mins


def getMaxHeats(diffused):

	max = {}
	for file in diffused:
		for (gene, heat) in diffused[file].iteritems():
			if gene in max:
				if max[gene] < heat:
					max[gene] = heat
			else:
				max[gene] = heat

	return max

def getLinkersMulti(input_heats, diffused_heats, input_order, size):
	'''
		Get linkers based on the size cutoff and supplied input heats/diffused heats
	'''

	# construct min-heats list for all genes here

	# use the order to construct minimum 
	set1 = [input_order[0], input_order[1]]
	set2 = [input_order[1], input_order[2]]
	set1_input = {}	
	set1_input[set1[0]] = diffused_heats[set1[0]]
	set1_input[set1[1]] = diffused_heats[set1[1]]
	set2_input = {}	
	set2_input[set1[0]] = diffused_heats[set2[0]]
	set2_input[set1[1]] = diffused_heats[set2[1]]

	min1 = getMinHeats(set1_input)
	min2 = getMinHeats(set2_input)
	meta_comparison = {}
	meta_comparison['s1'] = min1
	meta_comparison['s2'] = min2
	# final linker heats here: get MAX of either
	linker_heats = getMaxHeats(meta_comparison)
	
	EPSILON = 0.0001
	cutoff = None
	score = None
	linkers = set()
	linker_scores = {}
	for (l,h) in sorted(linker_heats.iteritems(), key=operator.itemgetter(1), reverse=True):
		c = h-EPSILON
		score, size_frac = scoreLinkersTest(input_heats, linker_heats, c, size)
		cutoff = c
		linkers.add(l)	
		linker_scores[l] = h
		if size_frac > 1:
			break

	return (cutoff, linkers, linker_scores)

def getLinkers(input_heats, diffused_heats, size, maxHeats):
	'''
		Get linkers based on the size cutoff and supplied input heats/diffused heats
	'''

	# construct min-heats list for all genes here

	linker_heats = None	
	if not maxHeats:
		# this is the default
		linker_heats = getMinHeats(diffused_heats)
	else:
		linker_heats = getMaxHeats(diffused_heats)
	
	EPSILON = 0.0001
	cutoff = None
	score = None
	linkers = set()
	linker_scores = {}
	for (l,h) in sorted(linker_heats.iteritems(), key=operator.itemgetter(1), reverse=True):
		c = h-EPSILON
		score, size_frac = scoreLinkersTest(input_heats, linker_heats, c, size)
		cutoff = c
		linkers.add(l)	
		linker_scores[l] = h
		if size_frac > 1:
			break

	return (cutoff, linkers, linker_scores)


def filterLinkers(up_heats_diffused, down_heats_diffused, cutoff):

	linkers = {}
	filtered = []
	for node in up_heats_diffused:
		if down_heats_diffused is None:
			min_heat = up_heats_diffused[node]	
		else:
			if node not in down_heats_diffused:
				continue
			min_heat = min(up_heats_diffused[node], down_heats_diffused[node])
		linkers[node] = min_heat
		if min_heat > cutoff:
			filtered.append(node)

	return (filtered, linkers)

def mapUGraphToNetwork(edge_list, network):
	"""
		Map undirected edges to the network to form a subnetwork
	"""

	subnetwork = {}
	
	for (s,t) in edge_list:

		# edges: 
		if s in network:
			for (i, nt) in network[s]:
				if nt == t:
					if s not in subnetwork:
						subnetwork[s] = set()
					subnetwork[s].add((i,t))	

		if t in network:
			for (i, nt) in network[t]:
				if nt == s:
					if t not in subnetwork:
						subnetwork[t] = set()
					subnetwork[t].add((i,s))	
	
	return subnetwork	


def connectedSubnets(network, subnet_nodes):

	edgelist = set()

	for s in network:
		for (i,t) in network[s]:
			if s in subnet_nodes and t in subnet_nodes:
				edgelist.add((s,t))

	return edgelist

def connectedNodes(network, hot_nodes):

	nodes = set()
	for (s, t) in connectedSubnets(network, hot_nodes):
		nodes.add(s)
		nodes.add(t)
	return nodes

def runPCST(up_heats, down_heats, linker_genes, network_file):
	"""
		Convert input to format used for PCST program.
	"""
	
	# convert up/down heats to p-values	
	# find the maximum heat for any value
	s_up = sorted([v for k, v in up_heats.iteritems()], reverse=True)
	s_down = sorted([v for k, v in down_heats.iteritems()], reverse=True)

	if len(up_heats) > 0:	
		max_heat = s_up[0]
		min_heat = s_up[-1]
	
		if len(s_down) > 0:
			if s_down[0] > max_heat:
				max_heat = s_down[0]
				min_heat = s_up[-1]
			if s_down[-1] > min_heat:
				min_heat = s_down[-1]
	else:
		max_heat = 1
		min_heat = 1

	# take the sqrt of the fold difference over the min
	normalized_max = math.sqrt(max_heat/min_heat)
	scores = {}
	# the order is important here: there may be overlap between the source, target
	# and linker sets. The linkers are the highest priority, over the source/target.
	for node in down_heats:
		heat = down_heats[node]
		normalized_heat = math.sqrt(heat/min_heat)
		pval = math.exp( normalized_heat*math.log(float("1e-10"))/normalized_max )
		scores[node] = str(pval)
	for node in up_heats:
		heat = up_heats[node]
		normalized_heat = math.sqrt(heat/min_heat)
		pval = math.exp( normalized_heat*math.log(float("1e-10"))/normalized_max )
		scores[node] = str(pval)
	for node in linker_genes:
		scores[node] = "1e-10"	


	pid = str(os.geteuid())

	tmp_act = open("/tmp/tmp_act_"+pid+".tab",'w')
	for node in scores:
		tmp_act.write(node+"\t"+scores[node]+"\n")
	tmp_act.close()	

	os.system(sys.path[0]+"/span.R --activities /tmp/tmp_act_"+pid+".tab --network "+network_file+" > /tmp/pcst_"+pid+".tab 2>/dev/null")

	pcst_network = []
	pcst_line = re.compile("\[\d+\]\s+(\S+)\s+\-\-\s+(\S+)\s+")
	pcst = open("/tmp/pcst_"+pid+".tab",'r')	
	for line in pcst:
		m = pcst_line.match(line)	
		if m:
			pcst_network.append((m.group(1),m.group(2)))	
	pcst.close()

	return pcst_network	


def writeNetwork(net, out_file):

	out = open(out_file, 'w')
	for source in net:
		for (int, target) in net[source]:
			out.write("\t".join([source, int, target])+"\n")

	out.close()

def randomSubnet(network, num_sources):

	sub = {}
	for source in random.sample(network, num_sources):
		sub[source] = network[source]
			
	return sub

def writeEL(el, so, down_set, out_file):

	out = open(out_file, 'w')
	for (source,int,target) in el:
		out.write("\t".join([source, int, target])+"\n")
	out.close()

	if so is None or down_set is None:
		return

	out = open(out_file+".txt", 'w')
	set2 = set()
	for (source,int,target) in el:
		if target in down_set:
			set2.add(target)
	out.write(so+"\t"+"\t".join(set2)+"\n")
	out.close()

def sampleHeats(heats):

	ss = int(len(heats)*0.8)
	keys = random.sample(heats, ss)
	subset = {}
	for k in keys:
		subset[k] = heats[k] 

	return subset


