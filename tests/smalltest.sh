#!/bin/bash
set -eu

mut_file=/fh/fast/matsen_e/wdumm/gctree_validation/MK_RS5NF_mutability.csv
sub_file=/fh/fast/matsen_e/wdumm/gctree_validation/MK_RS5NF_substitution.csv
export QT_QPA_PLATFORM=offscreen
export XDG_RUNTIME_DIR=/tmp/runtime-runner
export MPLBACKEND=agg
mkdir -p tests/smalltest_output
wget -O HS5F_Mutability.csv https://bitbucket.org/kleinstein/shazam/raw/ba4b30fc6791e2cfd5712e9024803c53b136e664/data-raw/HS5F_Mutability.csv
wget -O HS5F_Substitution.csv https://bitbucket.org/kleinstein/shazam/raw/ba4b30fc6791e2cfd5712e9024803c53b136e664/data-raw/HS5F_Substitution.csv
gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer  --root GL --frame 1 --verbose --idlabel
gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt
gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer --root GL --frame 1 --verbose --idlabel --mutability $mut_file --substitution $sub_file
gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability $mut_file --substitution $sub_file
