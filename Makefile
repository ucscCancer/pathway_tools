

galaxy : galaxy_pathway_tools
	cp scripts/edge_counter.py galaxy_pathway_tools/
	cp galaxy/edge_counter.xml galaxy_pathway_tools/
	cp scripts/network_convert.py galaxy_pathway_tools/
	cp galaxy/network_convert.xml galaxy_pathway_tools/
	cp -r pathway_tools galaxy_pathway_tools/

galaxy_pathway_tools : 
	mkdir galaxy_pathway_tools

clean :
	rm -rf galaxy_pathway_tools