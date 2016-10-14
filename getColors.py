#!/bin/env python

import sys, scipy
from Bio import AlignIO

aln = AlignIO.read(sys.argv[1], 'fasta')

red = scipy.array([1, 0, 0])
yellow = scipy.array([1, 1, 0])
cyan = scipy.array([0, 1, 1])

colors = {}
for seq in aln:
    if seq.id[-4:] == 'IgGR':
        color = red
    elif seq.id[-4:] == 'IgGY':
        color = yellow
    elif seq.id[-4:] == 'IgGC':
        color = cyan
    else:
        continue
    if seq.seq in colors:
        colors[str(seq.seq)] = scipy.append(colors[str(seq.seq)], [color], axis=0)
    else:
        colors[str(seq.seq)] = scipy.array([color])

# convert to hex color mixes
for seq in colors:
    print seq + '\t#%02x%02x%02x' % tuple(255*colors[seq].mean(axis=0))
