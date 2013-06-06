#!/usr/bin/env python

import sys
from optparse import OptionParser




# local imports
from pathway_tools.kernel import Kernel, MultiDiffuser
from pathway_tools.data_matrix import DataMatrix
#from tiedie_util import *

# Program Constants
SCORE_MU = 0.1


def main():
    parser = OptionParser()
    parser.add_option("-i","--input_heats",dest="input_heats",action="store",default=None,help="\
    File with columns as samples, rows as pathway features")
    parser.add_option("-k","--kernel",dest="kernel",action="store",default=None,help="\
    Pre-computed heat diffusion kernel in tab-delimited form. Should have both a header and row labels")
    parser.add_option("-o","--output",dest="output",action="store",default="diffused.tab")
    parser.add_option("-n","--network",dest="network",action="store",default=None,help="\
    .sif network file for the curated pathway to search. <source>   <(-a>,-a|,-t>,-t|,-component>)> <target>")
    parser.add_option("--pagerank",dest="pagerank",action="store_true",default=False,help="Use Personalized PageRank to Diffuse")
    (opts, args) = parser.parse_args()

    if opts.input_heats is None:
        print "Please define input_heat"
        return

    input_heats = DataMatrix(opts.input_heats)

    #
    # Diffusion Step:
    #   Load the heat diffusion kernel and perform a kernel-multiply, or alternatively use supplied
    #   page-rank diffused vectors
    #

    print "Parsing Network File..."
    network = parseNet(opts.network)
    network_nodes = getNetworkNodes(network) 
    diffuser = None
    if opts.pagerank:
        from ppr import PPrDiffuser
        diffuser = PPrDiffuser(network)
    else:
        print "Loading Heat Diffusion Kernel.."
        if not opts.kernel:
            # requires SCIPY installation
            from kernel_scipy import SciPYKernel
            diffuser = SciPYKernel(opts.network)
        else:   
            diffuser = Kernel(opts.kernel)

    # remove any feature not in the network
    feature_list = input_heats.getFeatureList()
    restricted = []
    for f in set(feature_list).intersection(network_nodes):
        restricted.append(f)
    feature_list = restricted

    # negative values require dual-diffusion and a merge/sum step
    m_diffuser = MultiDiffuser(diffuser)

    out = open(opts.output, 'w')
    out.write("Sample\t"+"\t".join(feature_list)+"\n")
    for sample in input_heats.getSampleList():
        in_heat_vec = input_heats.getSampleValues(sample)   
        diffused_heat_vec = m_diffuser.diffuse(in_heat_vec, reverse=False)
        # print diffused heats in order, to line up with feature columns
        out.write(sample+"\t"+"\t".join([str(diffused_heat_vec[feature]) for feature in feature_list])+"\n")
        
if __name__ == "__main__":
    main()