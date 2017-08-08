#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''compute GCtree for a phylip consensus tree'''

from __future__ import division, print_function
from gctree import *
import argparse



def main():
    parser = argparse.ArgumentParser(description='consensus tree support')
    parser.add_argument('consense', type=str, help='dnapars outfile for fixed consense tree')
    parser.add_argument('countfile', type=str, help='File containing allele frequencies (sequence counts) in the format: "SeqID,Nobs"')
    parser.add_argument('bootstrap', type=str, help='dnapars outfile from seqboot (multiple data sets)')
    parser.add_argument('--naive', type=str, default=None, help='name of naive sequence (outgroup root)')
    parser.add_argument('--frame', type=int, default=None, choices=(1, 2, 3), help='codon frame')
    parser.add_argument('--outbase', type=str, default='gctree.out', help='output file base name')
    args = parser.parse_args()

    consense = phylip_parse.parse_outfile(args.consense, args.countfile, args.naive)
    assert len(consense) == 1
    consensus_gctree = CollapsedTree(tree=consense[0], frame=args.frame)

    bootstrap_gctrees = []
    # fractional trees for parsimony degeneracy
    weights = []

    for bootstrap in phylip_parse.parse_outfile(args.bootstrap, args.countfile, args.naive):
        new_gctrees = []
        for tree in bootstrap:
            gctree = CollapsedTree(tree=tree, frame=args.frame, allow_repeats=True)
            if sum(gctree.compare(gctree2, method='identity') for gctree2 in new_gctrees) == 0:
                new_gctrees.append(gctree)
        bootstrap_gctrees.extend(new_gctrees)
        weights.extend([1/len(new_gctrees) for _ in range(len(new_gctrees))])

        assert len(weights) == len(bootstrap_gctrees)

    consensus_gctree.support(bootstrap_gctrees, weights=weights)
    consensus_gctree.render(args.outbase+'.bootstrap_support.svg', show_support=True)
    consensus_gctree.support(bootstrap_gctrees, weights=weights, compatibility=True)
    consensus_gctree.render(args.outbase+'.bootstrap_compatibility.svg', show_support=True)

if __name__ == "__main__":
    main()
