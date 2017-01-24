#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
aggregation plots of validation output from several simulation/validation runs
'''

from __future__ import division, print_function
import scipy, matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rc, ticker
import pandas as pd
import seaborn as sns
import argparse

parser = argparse.ArgumentParser(description='germinal center tree inference and simulation')
parser.add_argument('input', type=str, nargs='+', help='gctree.validation.tsv files')
parser.add_argument('--outbase', type=str, help='output file base name')
args = parser.parse_args()

aggdat = pd.DataFrame(columns=('forest size',
                               'max l', 'max l distance', 'max l distance rank',
                               'min l', 'min l distance', 'min l distance rank'))

for i, fname in enumerate(args.input):
    df = pd.read_csv(fname, sep='\t')
    maxl_idx = df['log-likelihood'].idxmax()
    minl_idx = df['log-likelihood'].idxmin()
    d_ranks = df['distance'].rank(method='min')
    forest_size = len(df.index)
    l_max = df['log-likelihood'][maxl_idx]
    l_max_distance = df['distance'][maxl_idx]
    l_max_distance_rank = d_ranks[maxl_idx]
    l_min = df['log-likelihood'][minl_idx]
    l_min_distance = df['distance'][minl_idx]
    l_min_distance_rank = d_ranks[minl_idx]

    aggdat.loc[i] = (forest_size,
                     l_max, l_max_distance, l_max_distance_rank,
                     l_min, l_min_distance, l_min_distance_rank)

aggdat['delta l'] = aggdat['max l'] - aggdat['min l']
aggdat['delta r'] = (aggdat['min l distance rank'] - aggdat['max l distance rank'])/aggdat['forest size']

aggdat.to_csv(args.outbase+'.tsv', sep='\t', index=False)
sns.jointplot(x='forest size', y='delta l', data=aggdat, kind="kde").savefig(args.outbase+'.delta_l.pdf')
sns.jointplot(x='forest size', y='delta r', data=aggdat, kind="kde").savefig(args.outbase+'.delta_r.pdf')
