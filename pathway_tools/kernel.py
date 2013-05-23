#!/usr/bin/env  python2.7

from numpy import genfromtxt, dot
import sys
import math

class Kernel:

    def __init__(self, kernel_files):
        """ 
            Input:
                kernel_file - a tab-delimited matrix file with both a header
                and first-row labels, in the same order. 
        """

        # might have multiple kernels here
        self.kernels = {}
        self.labels = {}
        self.ncols = {}
        self.nrows = {}
        for kernel in kernel_files.split(":"):
            self.kernels[kernel] = genfromtxt(kernel,delimiter="\t")[1:,1:]
            self.labels[kernel] = None
            fh = open(kernel,'r')
            # get the header line
            for line in fh:
                self.labels[kernel] = line.rstrip().split("\t")[1:]
                break
            fh.close()

            self.ncols[kernel] = self.kernels[kernel].shape[1]-1
            self.nrows[kernel] = self.kernels[kernel].shape[0]-1

    def kernelMultiplyOne(self, kernel, vector):
        """
            Input:
                vector: A hash mapping gene labels to floating point values 
        """
        array = []
        for label in self.labels[kernel]:
            if label in vector:
                array.append(vector[label])
            else:
                array.append(0)

        value = dot(self.kernels[kernel], array)

        return_vec = {}
        idx = 0
        for label in self.labels[kernel]:
            return_vec[label] = float(value[idx])
            idx += 1

        return return_vec

    def addVectors(self, vector_list):
        sum = {}

        for vec in vector_list:
            for key in vec:
                val = vec[key]
                if key not in sum:
                    sum[key] = val
                else:
                    sum[key] += val

        return sum

    @staticmethod
    def getAngle(v1, v2):

        arry1 = []
        arry2 = []
        for key in v1:
            arry1.append(float(v1[key]))
            arry2.append(float(v2[key]))

        mag_1 = math.sqrt(dot(arry1,arry1))
        mag_2 = math.sqrt(dot(arry2,arry2))

        cos_theta = dot(arry1,arry2)/(mag_1*mag_2)

        return math.acos(cos_theta)

    def diffuse(self, vector, reverse=False):

        # reverse is not used: heat diffusion is undirected
        return_vectors = []
        for kernel in self.kernels:
            diffused_vector = self.kernelMultiplyOne(kernel, vector)
            return_vectors.append(diffused_vector)

        return self.addVectors(return_vectors)

class MultiDiffuser:

    def __init__(self, absolute_diffuser):

        self.diffuser_engine = absolute_diffuser


    def diffuse(self, vector, reverse=False):
        '''
            Split the input vector into positive and negative components, then merge back after diffusing 
            each
        '''
        pos_vector = {}
        neg_vector = {}
        for feature in vector:
            val = vector[feature]
            if float(val) == val and val < 0:
                neg_vector[feature] = abs(val)
                pos_vector[feature] = 0
            else:
                pos_vector[feature] = val
                neg_vector[feature] = 0

        pos_diffused = self.diffuser_engine.diffuse(pos_vector, reverse)
        neg_diffused = self.diffuser_engine.diffuse(neg_vector, reverse)
        # merge back
        merged_diffused = {}
        for feature in vector:
            v = pos_diffused[feature]-neg_diffused[feature]
            merged_diffused[feature] = v

        return merged_diffused    
