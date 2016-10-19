#! /bin/env bash

# this assumes scripts live in a dir "GH" and data in a dir "data"

fasta=$1
gctree_outdir=$2
germline=$3

# parse fasta data to phylip file
phylip=${fasta}.phylip
python gctree/TasParse.py ${fasta} > ${phylip}

rm outfile
rm outtree
rm ${gctree_outdir}/*png

# run phylip's dna parsimony program and rename its outputs
dnapars <<STDIN
`pwd`/${phylip}
O
1
V
100000
J
13
10
4
5
.
Y

STDIN

outfile=${phylip}.dnapars.outfile
outtree=${phylip}.dnapars.outtree

mv outfile ${outfile}
mv outtree ${outtree}

mkdir -p ${gctree_outdir}

# branching process likelihoods of parsimony results
python gctree/gctree.py --phylipfile ${outfile} --plot_file ${gctree_outdir}/gctree --germline ${germline}
