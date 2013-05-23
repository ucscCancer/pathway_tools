#!/usr/bin/env	python2.7

import sys

class DataMatrix:

	def __init__(self, matrix_file=None):


		dfh = None
		try:
			dhf = open(matrix_file,'r')	
		except:
			raise Exception("Error: can't open data file:"+matrix_file)

		features, samples = self.parseMatrix(dhf)

		self.features = features
		self.samples = samples
	

	def getFeatureValues(self, feature, restrict_samples=None):

		values = []
		samples = []
		if feature in self.features:
			for sample in self.samples.keys():
				if restrict_samples is None or sample in restrict_samples:
					samples.append(samples)
					values.append(self.features[feature][sample])
				else:
					continue
	
		return (values)

	def getFeatureList(self):
		return self.features.keys()

	def getSampleValues(self, sample):

		# build a hash indexed by 
		values = {}
		for feature in self.features:
			values[feature] = self.features[feature][sample]
	
		return (values)
	
	def getValue(self, idx1, idx2):
		'''
			Either index 1 or 2 could be a feature, the other should be a sample.
			Find the value associated with it in the matrix, and return.
		'''
		if idx1 in self.features:
			if idx2 not in self.features[idx1]:
				return None
			else:
				return self.features[idx1][idx2]
		elif idx1 in self.samples:
			if idx2 not in self.samples[idx1]:
				return None
			else:
				return self.samples[idx1][idx2]
		else:
			return None	
			
		
	def getSampleList(self):
		return self.samples.keys()
	
	def parseMatrix(self, dfh, convert_ids=False):
		'''
			Columns are TCGA Sample IDS/Patient IDS
			Rows are feature names	
		'''

		first = True
		header = None

		samples = {}
		features = {}

		lineno = 1

		for line in dfh:
			parts = line.rstrip().split("\t")
			feature = parts[0]
			vals = parts[1:]
			if first:
				first = False
				header = vals
				if convert_ids:
					fixedNames = []
					for i in range(0, len(header)):
						patient = header[i][0:12]
						fixedNames.append(patient)
					header = fixedNames

				# initialize a defaultdict for each sample
				for id in header:
	
					# check for duplicate samples
					if id in samples:
						raise Exception("Error: duplicate sample/patient id:"+id)

					samples[id] = {}

				continue

			lineno += 1	

			# check for duplicate features
			if feature in features:
				raise Exception("Error: duplicate feature:"+feature+" lineno: "+str(lineno))

			features[feature] = {}

			for i in range(0,len(vals)):

				val = None
				try:
					val = float(vals[i])
				except:
					sys.stderr.write("Error: cannot cast value to float: line"+str(lineno)+" column: "+str(i)+" value: "+str(val)+"\n")
					val = None

				sample_id = header[i]		

				# index by both features and samples
				features[feature][sample_id] = val
				samples[sample_id][feature] = val

		feature_names = set(features.keys())
		sample_names = set(samples.keys())
		if len(feature_names.intersection(sample_names)) > 0:
			raise Exception("Error: feature and sample names are not exclusive!")

		return (features, samples)	
			
	
