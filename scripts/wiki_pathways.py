#!/usr/bin/env python

# Load SOAPpy and dependent modules (fpconst) and access the remote
# SOAP server through a proxy class, SOAPProxy - see:
# (http://diveintopython.org/soap_web_services/first_steps.html)

import sys
import os
import csv
import base64
from SOAPpy import SOAPProxy      
url = 'http://www.wikipathways.org/wpi/webservice/webservice.php'
namespace = 'http://www.wikipathways.org/webservice'

import argparse

from pathway_tools import biopax
from pathway_tools import convert
import networkx as nx

# Print out function for query results (see code below function first)
def printOutput(ws_output):
    #Loops through a list of dictionary items
    index=1
    for object in ws_output:
        #calls select dictionary keys to print out values
        print_output = '   '+str(index)+')'+'species:'+object['species']+'\t '
        print_output+= 'id:'+object['id']+'\t '+'name:'+object['name']
        print print_output
        index+=1


#mapping of databases referenced by WikiPathways, and what they can be translated to
db_map = {
    "Ensembl" : "Ensembl ID",
    "Ensembl Human" : "Ensembl ID",
    "Uniprot" : "UniProt ID",
    "UniProt" : "UniProt ID",
    'Entrez Gene' : "Entrez Gene ID",
    "Enzyme Nomenclature" : None,
    "SwissProt" : None,
    'Affy' : None,
    "ec-code" : None,
    "PubChem" : None,
    "Wikipedia" : None,
    "Other" : None,
    "EMBL" : None,
    "miRBase Sequence" : None,
    "HGNC" : None,
    "GenBank" : None,
    "UniGene" : None,
    "PubChem-compound" : None,
    "Reactome" : None,
    "COMPOUND" : None,
    "ChEBI" : None,
    "GenBank" : "Accession Numbers",
    "GLYCAN" : None,
    "CAS" : None,
    "Kegg ortholog" : None,
    "KEGG Orthology" : None,
    "KEGG Genes" : None,
    "OMIM" : "OMIM ID",
    "RefSeq" : "RefSeq IDs",
    "enzyme" : "Enzyme IDs",
    "CTD Gene" : None,
    "GeneDB" : None,
    "Pfam" : None,
    "WikiPathways" : None
}

def find_translation(label):
    out = None
    try:
        db, db_id = gr.node[n]['db_xref'].split(":")
        found = False
        for row in translate_table:
            if db not in db_map:
                print "WARNING: Database not found: %s" % (db)
            else:
                if db_map[db] in row:
                    if row[db_map[db]] == db_id:
                        if args.verbose:
                            print db_id, row['Approved Symbol']
                        found = True
                        out = row['Approved Symbol']
        #print found
    except ValueError:
        pass
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', "--outdir")
    parser.add_argument('-t', "--translate")
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-f", "--file-mode", action="store_true", default=False)
    parser.add_argument("pathways", nargs="*")
    args = parser.parse_args()

    server = SOAPProxy(url, namespace)
    
    if len(args.pathways) == 0:    
        pathways = []
        for row in server.listPathways():
            if row['species'] == "Homo sapiens":
                pathways.append(row['id'])
    else:
        pathways = args.pathways

    translate_table = None
    if args.translate:
        handle = open(args.translate)
        reader = csv.DictReader(handle, delimiter="\t")
        translate_table = []
        for row in reader:
            translate_table.append(row)
        handle.close()

    if args.file_mode:
        for path in pathways:
            handle = open(path)
            gr = convert.read_gpml(handle)
            handle.close()
            if translate_table:
                re_map = {}
                for n in gr.node:
                    if 'db_xref' in gr.node[n]:
                        relabel = find_translation(gr.node[n]['db_xref'])
                        if relabel:
                            re_map[n] = relabel
                gr = nx.relabel_nodes(gr, re_map)
            convert.write_xgmml(gr, sys.stdout)            

    else:
        for path in pathways:
            print "Getting", path
            pathdata_str = server.getPathwayAs(fileType="owl", pwId=path)
            pathdata_xml = base64.b64decode(pathdata_str)
            handle = open( os.path.join(args.outdir, path + ".owl"), "w" )
            handle.write(pathdata_xml)
            handle.close()
            b = biopax.BioPax()
            b.parse(pathdata_xml)
            for gr in b.toNet():
                if translate_table:
                    re_map = {}
                    for n in gr.node:
                        if 'db_xref' in gr.node[n]:
                            relabel = find_translation(gr.node[n]['db_xref'])
                            if relabel:
                                re_map[n] = relabel
                    gr = nx.relabel_nodes(gr, re_map)
                handle = open( os.path.join(args.outdir, path + ".xgmml"), "w" )
                #convert.write_paradigm_graph(gr, handle)            
                convert.write_xgmml(gr, handle)            
                handle.close()
            

    """
    # listOrganisms (no args; returns list of strings)
    wp_organisms = server.listOrganisms()
    print '\nSupported organisms at WikiPathways'
    index=1
    for organism in wp_organisms:
        print '   '+str(index)+')',organism; index+=1

    # getPathwayInfo (one arg; returns dictionary reference)
    pathway = 'WP274'
    pathway_info = server.getPathwayInfo(pwId = pathway)
    print '\nPathway information for %s' % pathway
    # Access data via dictionary reference
    printOutput([pathway_info])

    # findPathwaysByText (multiple args; returns list of dictionary references)
    pathway = 'apoptosis'    
    apoptosis_containing = server.findPathwaysByText(query= pathway, species = "")
    print '\nPathways containing the term %s' % pathway 
    printOutput(apoptosis_containing)

    # Define the order of args: needed for this service
    server.config.argsOrdering = {'findPathwaysByXref': ('ids', 'codes') }

    # findPathwaysByXref (multiple args; returns list of dictionary references)
    sc = 'X'; gi = '201746_at'
    probeset_containing = server.findPathwaysByXref(codes=sc, ids=gi )
    print '\nPathways containing the gene ID "%s" for system code "%s"' % (gi,sc) 
    printOutput(probeset_containing)
    """
