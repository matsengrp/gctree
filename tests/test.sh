#!/bin/bash

deduplicate example/150228_Clone_3-8.fasta --root GL --abundance_file abundances.csv --idmapfile idmap.txt > deduplicated.phylip
mkconfig deduplicated.phylip dnapars > dnapars.cfg
dnapars < dnapars.cfg > dnapars.log
export QT_QPA_PLATFORM=offscreen
export XDG_RUNTIME_DIR=/tmp/runtime-runner
export MPLBACKEND=agg
xvfb-run -a gctree infer outfile abundances.csv --root GL --frame 1 --verbose
