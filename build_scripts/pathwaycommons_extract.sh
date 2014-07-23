#!/bin/bash

mkdir -p work/Pathway_Commons.4.All.BIOPAX.split
curl -o work/Pathway_Commons.4.All.BIOPAX.owl.gz http://www.pathwaycommons.org/pc2/downloads/Pathway%20Commons.4.All.BIOPAX.owl.gz
./sbt/sbt "run work/Pathway_Commons.4.All.BIOPAX.owl.gz work/Pathway_Commons.4.All.BIOPAX.split"
./scripts/pathway_commons_rename.py work/Pathway_Commons.4.All.BIOPAX.split/ > work/Pathway_Commons.4.All.BIOPAX.split.names
./scripts/pathway_db.py biopax-format work/Pathway_Commons.4.All.BIOPAX.library biopax-dir work/Pathway_Commons.4.All.BIOPAX.split --rename work/Pathway_Commons.4.All.BIOPAX.split.names 
./scripts/pathway_db.py biopax-convert --cpus=8 work/Pathway_Commons.4.All.BIOPAX.library
