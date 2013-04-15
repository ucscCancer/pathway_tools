#!/usr/bin/env python

from pathway_tools.infer import pathway_inference

from argparse import ArgumentParser

                
if __name__ == "__main__":
    parser = ArgumentParser()    
    parser.add_argument("-d", "--dogma", 
        dest="dogma", help="DogmaFile", default=None)
    parser.add_argument("-p", "--pathway", 
        dest="pathway", help="PathwayFile", default=None)
    parser.add_argument("-e", "--evidence", nargs=2, action="append", default=[])
    parser.add_argument("--report", action="store_true", help="Print out factor graph and quit", default=False)
    parser.add_argument("--dai-file", help="Build factor graph in DAI, then print out to file and quit", default=None)
    parser.add_argument("-s", "--sample", help="Sample Name", default=None)
    parser.add_argument("--em-mode", action="store_true", default=False)
    
    args = parser.parse_args()
    pathway_inference(args)