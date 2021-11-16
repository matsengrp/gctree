#! /bin/bash
ml libGLU
ml Xvfb

TMPDIR=/tmp xvfb-run -a \
    isotype --parsimony_forest littleforest.p \
    --inference_log gctree.inference.log --isotype_mapfile example.isotypemap \
    --idmapfile example.idmap --out_directory isotypes

