#!/usr/bin/env python

# Load SOAPpy and dependent modules (fpconst) and access the remote
# SOAP server through a proxy class, SOAPProxy - see:
# (http://diveintopython.org/soap_web_services/first_steps.html)

import sys
import os
import base64
from SOAPpy import SOAPProxy      
url = 'http://www.wikipathways.org/wpi/webservice/webservice.php'
namespace = 'http://www.wikipathways.org/webservice'

from pathway_tools import biopax
from pathway_tools import convert

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


if __name__ == "__main__":
    server = SOAPProxy(url, namespace)
    outdir = sys.argv[1]
    
    for row in server.listPathways():
        if row['species'] == "Homo sapiens":
            print row
            pathdata_str = server.getPathwayAs(fileType="owl", pwId=row['id'])
            pathdata_xml = base64.b64decode(pathdata_str)
            b = biopax.BioPax()
            b.parse(pathdata_xml)
            for gr in b.toNet():
                handle = open( os.path.join(outdir, row['id']), "w" )
                convert.write_paradigm_graph(gr, handle)            
                handle.close()
            
    #print server.getPathway(pwId="WP98")

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
