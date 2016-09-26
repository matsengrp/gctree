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
                3
                4
                5
		Y

		STDIN
mv outfile data/Tas/150228_Clone_3-8.dnapars.outfile
mv outtree data/Tas/150228_Clone_3-8.dnapars.outtree

# branching process likelihoods of parsimony results
python GH/recurse.py data/Tas/150228_Clone_3-8.dnapars.outtree --plot_file data/Tas/150228_Clone_3-8.collapsed_tree.pdf --germline GL
