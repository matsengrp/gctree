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
                               'max l', 'max l RF', 'max l RF rank', 'max l MRCA', 'max l MRCA rank',
                               'min l', 'min l RF', 'min l RF rank', 'min l MRCA', 'min l MRCA rank'))

for i, fname in enumerate(args.input):
    df = pd.read_csv(fname, sep='\t')
    maxl_idx = df['log-likelihood'].idxmax()
    minl_idx = df['log-likelihood'].idxmin()
    d_ranks = df['RF'].rank(method='min')
    MRCA_ranks = df['MRCA'].rank(method='min')
    forest_size = len(df.index)

    l_max = df['log-likelihood'][maxl_idx]
    l_max_RF = df['RF'][maxl_idx]
    l_max_RF_rank = d_ranks[maxl_idx]
    l_max_MRCA = df['MRCA'][maxl_idx]
    l_max_MRCA_rank = MRCA_ranks[maxl_idx]

    l_min = df['log-likelihood'][minl_idx]
    l_min_RF = df['RF'][minl_idx]
    l_min_RF_rank = d_ranks[minl_idx]
    l_min_MRCA = df['MRCA'][minl_idx]
    l_min_MRCA_rank = MRCA_ranks[minl_idx]

    aggdat.loc[i] = (forest_size,
                     l_max, l_max_RF, l_max_RF_rank, l_max_MRCA, l_max_MRCA_rank,
                     l_min, l_min_RF, l_min_RF_rank, l_min_MRCA, l_min_MRCA_rank)

aggdat['delta l'] = aggdat['max l'] - aggdat['min l']
aggdat['delta RF rank'] = (aggdat['min l RF rank'] - aggdat['max l RF rank'])/aggdat['forest size']
aggdat['delta MRCA rank'] = (aggdat['min l MRCA rank'] - aggdat['max l MRCA rank'])/aggdat['forest size']

aggdat.to_csv(args.outbase+'.tsv', sep='\t', index=False)
g = sns.pairplot(aggdat, kind='reg', x_vars='forest size', y_vars=('delta l', 'delta RF rank', 'delta MRCA rank'),
                 aspect=1.5, plot_kws={'fit_reg':None})
g.set(xlim=(0, None))
bound = scipy.ceil(max(aggdat['delta l']))
g.axes[0, 0].set_ylim(-bound, bound)
g.axes[1, 0].set_ylim(-1, 1)
g.axes[2, 0].set_ylim(-1, 1)
g.savefig(args.outbase+'.pdf')
# sns.jointplot(x='forest size', y='delta l', data=aggdat, kind="kde").savefig(args.outbase+'.delta_l.pdf')
# sns.jointplot(x='forest size', y='delta RF r', data=aggdat, kind="kde").savefig(args.outbase+'.delta_RF_rank.pdf')
# sns.jointplot(x='forest size', y='delta MRCA r', data=aggdat, kind="kde").savefig(args.outbase+'.delta_RF_rank.pdf')
