#!/bin/bash

DATADIR=work

./scripts/pathway_db.py library-compile --dogma data/simple_dogma.yaml --exclude data/exclude.list --merge-file data/merge.gmt --all $DATADIR/Pathway_Commons.4.All.BIOPAX.library --exclude-type chemical

./scripts/pathway_db.py build --organism 9606 --dogma data/simple_dogma.yaml --exclude data/exclude.list --merge-file data/merge.gmt --all library $DATADIR/Pathway_Commons.4.All.BIOPAX.library --spf --min-subgraph 2 --exclude-type chemical > $DATADIR/Pathway_Commons.4.All.BIOPAX.library.9606.spf
./scripts/pathway_db.py library-gmt --skip-empty $DATADIR/Pathway_Commons.4.All.BIOPAX.library > $DATADIR/Pathway_Commons.4.All.BIOPAX.library.9606.gmt
cat data/hugo.tab  | awk '{print $2}' > data/hugo.list
./scripts/pathway_db.py library-gmt $DATADIR/Pathway_Commons.4.All.BIOPAX.library --skip-empty --filter data/hugo.list > $DATADIR/Pathway_Commons.4.All.BIOPAX.library.9606.gmt_hugo
./scripts/pathway_db.py library-table $DATADIR/Pathway_Commons.4.All.BIOPAX.library > $DATADIR/Pathway_Commons.4.All.BIOPAX.library.table
./scripts/pathway_db.py library-copy --dogma data/simple_dogma.yaml --exclude data/exclude.list --merge-file data/merge.gmt --all $DATADIR/Pathway_Commons.4.All.BIOPAX.library $DATADIR/Pathway_Commons.4.All.BIOPAX.spf.zip
