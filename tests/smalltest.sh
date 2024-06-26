#!/bin/bash
set -eu

export QT_QPA_PLATFORM=offscreen
export XDG_RUNTIME_DIR=/tmp/runtime-runner
export MPLBACKEND=agg
mkdir -p tests/smalltest_output
wget -O HS5F_Mutability.csv https://bitbucket.org/kleinstein/shazam/raw/ba4b30fc6791e2cfd5712e9024803c53b136e664/data-raw/HS5F_Mutability.csv
wget -O HS5F_Substitution.csv https://bitbucket.org/kleinstein/shazam/raw/ba4b30fc6791e2cfd5712e9024803c53b136e664/data-raw/HS5F_Substitution.csv

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer  --root GL --frame 1 --verbose --idlabel

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability HS5F_Mutability.csv --substitution HS5F_Substitution.csv --ranking_coeffs 0 1 0 --use_old_mut_parsimony --branching_process_ranking_coeff 0 

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability HS5F_Mutability.csv --substitution HS5F_Substitution.csv --ranking_coeffs 1 1 0 --use_old_mut_parsimony --branching_process_ranking_coeff 0 

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability HS5F_Mutability.csv --substitution HS5F_Substitution.csv --ranking_coeffs .01 -1 0 --branching_process_ranking_coeff -1 --summarize_forest --tree_stats


gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer --root GL --frame 1 --verbose --idlabel --mutability HS5F_Mutability.csv --substitution HS5F_Substitution.csv

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability HS5F_Mutability.csv --substitution HS5F_Substitution.csv
