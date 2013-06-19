#!/usr/bin/env python

import sys
import argparse
from pathway_tools import convert, biopax


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', "--output", default=None)
    parser.add_argument("input")
    args = parser.parse_args()

    handle = open(args.input)
    gr = convert.read_xgmml(handle)

    biopax.write_biopax(gr, sys.stdout)
