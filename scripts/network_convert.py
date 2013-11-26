#!/usr/bin/env python

import argparse
import zipfile
import sys

from pathway_tools import convert

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-etd', '--edge-type-default', help="Default Edge Type", default='-a>')
    parser.add_argument('-ntd', '--node-type-default', help='Default Node Type', default='protein')
    parser.add_argument('-ntf', '--node-type-field', help='Node Type Field', default='type')
    parser.add_argument('-etf', '--edge-type-field', help='Edge Type Field', default='interaction')
    parser.add_argument('-pmf', '--pathway-member-field', help="Pathway Membership Field", default="pathway")

    parser.add_argument('--in-xgmml', help='Input XGMML file', default=None)
    parser.add_argument('--in-spf', help='Input Simple Pathway Format', default=None)
    parser.add_argument('--in-gmt', help='Input Gene Mapping File', default=None)
    parser.add_argument('--in-spf-dir', help='Input SimplePathwayFormat Directory', default=None)

    parser.add_argument('--strict', help='Strict SimplePathwayFormat input', action="store_true", default=False)
    
        
    #parser.add_argument('-cys', help='Input Cytoscape file', default=None)    

    parser.add_argument('--out-spf', help='Output SimplePathwayFormat File', default=None)
    parser.add_argument('--out-sif', help='Output SimpleInteractionFormat File', default=None)    
    parser.add_argument('--out-xgmml', help='Output XGMML File', default=None)    
    parser.add_argument('--out-gmt', help="Output Gene Mapping Table", default=None)
    
    args = parser.parse_args()

    gr = None
    if args.in_xgmml is not None:
        if args.in_xgmml == "-":
            handle = sys.stdin
        else:
            handle = open(args.in_xgmml)
        gr = convert.read_xgmml(handle)
        handle.close()

    if args.in_spf is not None:
        handle = open(args.in_spf)
        gr = convert.read_spf_graph(handle, strict=args.strict)
        handle.close()

    if args.in_gmt is not None:
        handle = open(args.in_gmt)
        for line in handle:
            tmp = line.rstrip().split("\t")
            for elem in tmp[1:]:
                if elem in gr.node:
                    if args.pathway_member_field not in gr.node[elem]:
                        gr.node[elem][args.pathway_member_field] = []
                    gr.node[elem][args.pathway_member_field].append(tmp[0])
        handle.close()

    if args.in_spf_dir is not None:
        gr = load_spf_dir(args.in_spf_dir)

    """
    #not yet complete, because while cytoscape stores an xgmml file internally, it doesn't populate the 
    #att fields until export...
    if args.cys is not None:
        z = zipfile.ZipFile(args.cys)
        for n in z.namelist():
            if n.endswith(".xgmml") and os.path.dirname(n).endswith("networks"):
                handle = z.open(n)
                gr = convert.read_xgmml(handle)
                handle.close()
    """

    if gr is not None:
        if args.out_spf is not None:
            if args.out_spf == "-":
                ohandle = sys.stdout
            else:
                ohandle = open(args.out_spf, "w")
            convert.write_spf(gr, ohandle, 
                node_type_field=args.node_type_field, node_type_default=args.node_type_default, 
                edge_type_field=args.edge_type_field, edge_type_default=args.edge_type_default
            )
        if args.out_xgmml is not None:
            if args.out_xgmml == "-":
                ohandle = sys.stdout
            else:
                ohandle = open(args.out_xgmml, "w")
            convert.write_xgmml(gr, ohandle)
        if args.out_gmt:
            if args.out_gmt == "-":
                ohandle = sys.stdout
            else:
                ohandle = open(args.out_gmt, "w")

            pathmap = {}
            for node in gr.node:
                if 'pathway' in gr.node[node]:
                    for p in gr.node[node]['pathway']:
                        if p not in pathmap:
                            pathmap[p] = [node]
                        else:
                            pathmap[p].append(node)
            for path in pathmap:
                ohandle.write("%s\t%s\n" % (path, "\t".join(pathmap[path])))



