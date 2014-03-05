#!/bin/bash

./scripts/pathway_db.py build --dogma data/simple_dogma.yaml --exclude data/exclude.list --merge-file data/merge.gmt --all library Pathway_Commons.4.All.BIOPAX.library --spf --min-subgraph 2 --exclude-type chemical > Pathway_Commons.4.All.BIOPAX.library.spf
./scripts/pathway_db.py library-compile --dogma data/simple_dogma.yaml --exclude data/exclude.list --merge-file data/merge.gmt --all  Pathway_Commons.4.All.BIOPAX.library --exclude-type chemical
./scripts/pathway_db.py library-gmt Pathway_Commons.4.All.BIOPAX.library > Pathway_Commons.4.All.BIOPAX.library.gmt
./scripts/pathway_db.py library-gmt Pathway_Commons.4.All.BIOPAX.library --filter hugo.list > Pathway_Commons.4.All.BIOPAX.library.gmt_hugo 
