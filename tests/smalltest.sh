#!/bin/bash
set -eu

export QT_QPA_PLATFORM=offscreen
export XDG_RUNTIME_DIR=/tmp/runtime-runner
export MPLBACKEND=agg
mkdir -p tests/smalltest_output

mutabilities=HS5F_Mutability.csv
substitutions=HS5F_Substitution.csv

wget -O $mutabilities https://bitbucket.org/kleinstein/shazam/raw/ba4b30fc6791e2cfd5712e9024803c53b136e664/data-raw/HS5F_Mutability.csv
wget -O $substitutions https://bitbucket.org/kleinstein/shazam/raw/ba4b30fc6791e2cfd5712e9024803c53b136e664/data-raw/HS5F_Substitution.csv

# testing backward compatibility:

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/nonexistant_dir/gctree.infer --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability HS5F_Mutability.csv --substitution HS5F_Substitution.csv --ranking_coeffs 1 1 0 --use_old_mut_parsimony --branching_process_ranking_coeff 0 

gctree infer tests/smalltest_output/gctree.infer.inference.parsimony_forest.p --outbase tests/smalltest_output/gctree.infer --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability HS5F_Mutability.csv --substitution HS5F_Substitution.csv --ranking_coeffs 1 -1 0 --use_old_mut_parsimony

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability HS5F_Mutability.csv --substitution HS5F_Substitution.csv --use_old_mut_parsimony

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability HS5F_Mutability.csv --substitution HS5F_Substitution.csv --ranking_coeffs -1 1 0 --use_old_mut_parsimony --branching_process_ranking_coeff 0 --ranking_strategy "M"

# End testing backward compatibility ^^

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer1  --root GL --frame 1 --verbose --idlabel

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer2 --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability $mutabilities --substitution $substitutions --ranking_strategy "M"

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer3 --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability $mutabilities --substitution $substitutions --ranking_strategy "I+M"

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer4 --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability $mutabilities --substitution $substitutions --ranking_strategy="-B+0.01I+C" --summarize_forest --tree_stats


gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer5 --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer6 --root GL --frame 1 --verbose --idlabel --mutability $mutabilities --substitution $substitutions

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer7 --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability $mutabilities --substitution $substitutions

gctree infer tests/small_outfile tests/abundances.csv --outbase tests/smalltest_output/gctree.infer8 --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability $mutabilities --substitution $substitutions --ranking_strategy "R,-1.0B"

gctree infer tests/smalltest_output/gctree.infer.inference.parsimony_forest.p --outbase tests/smalltest_output/gctree.infer9 --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability $mutabilities --substitution $substitutions --ranking_strategy "A, 0B, 0R" --tree_stats --summarize_forest
