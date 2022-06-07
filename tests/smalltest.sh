#!/bin/bash

export QT_QPA_PLATFORM=offscreen
export XDG_RUNTIME_DIR=/tmp/runtime-runner
export MPLBACKEND=agg
gctree infer tests/small_outfile tests/abundances.csv --root GL --frame 1 --verbose --idlabel 
gctree infer tests/small_outfile tests/abundances.csv --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt
gctree infer tests/small_outfile tests/abundances.csv --root GL --frame 1 --verbose --idlabel --mutability S5F/Mutability.csv --substitution S5F/Substitution.csv
gctree infer tests/small_outfile tests/abundances.csv --root GL --frame 1 --verbose --idlabel --idmapfile tests/idmap.txt --isotype_mapfile tests/isotypemap.txt --mutability S5F/Mutability.csv --substitution S5F/Substitution.csv
