#! /bin/env bash

# this assumes scripts live in a dir "GH" and data in a dir "data"

# parse Tas data to phylip file
python GH/TasParse.py data/Tas/150228_Clone_3-8.fasta > data/Tas/150228_Clone_3-8.phylip

# run phylip's dna parsimony program and rename its outputs
dnapars <<- 	STDIN
		/fh/fast/matsen_e/wdewitt/gctree/data/Tas/150228_Clone_3-8.phylip
                O
                1
		J
		13
		10
                4
                5
                .
		Y

		STDIN
mv outfile data/Tas/150228_Clone_3-8.dnapars.outfile
mv outtree data/Tas/150228_Clone_3-8.dnapars.outtree

# branching process likelihoods of parsimony results
python GH/recurse.py --outtree data/Tas/150228_Clone_3-8.dnapars.outtree --outfile data/Tas/150228_Clone_3-8.dnapars.outfile --plot_file data/Tas/150228_Clone_3-8.collapsed_tree --germline GL
