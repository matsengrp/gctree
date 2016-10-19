#! /bin/env bash

# this assumes scripts live in a dir "GH" and data in a dir "data"

fasta=$1
gctree_outdir=$2
germline=$3

# parse fasta data to phylip file
phylip=${fasta}.phylip
python gctree/TasParse.py ${fasta} ${germline} > ${phylip}

rm -f outfile
rm -f outtree

echo -n "computing parsimony trees... "
# run phylip's dna parsimony program and rename its outputs
dnapars <<STDIN > /dev/null
`pwd`/${phylip}
O
1
V
10000
J
13
10
4
5
.
Y

STDIN
echo "done"

outfile=${phylip}.dnapars.outfile
outtree=${phylip}.dnapars.outtree

mv outfile ${outfile}
mv outtree ${outtree}

rm -f ${gctree_outdir}/*
mkdir -p ${gctree_outdir}

echo "branching process likelihood ranking of parsimony trees:"
python gctree/gctree.py --phylipfile ${outfile} --plot_file ${gctree_outdir}/gctree --germline ${germline}
