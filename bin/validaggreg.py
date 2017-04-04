#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
aggregation plots of validation output from several simulation/validation runs
'''

from __future__ import division, print_function
import scipy, matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
from matplotlib import rc, ticker
import pandas as pd
import argparse
import seaborn as sns; sns.set(style="white", color_codes=True)

parser = argparse.ArgumentParser(description='aggregate validation of repeated runs with same parameters')
parser.add_argument('input', type=str, nargs='+', help='gctree.validation.tsv files')
parser.add_argument('--outbase', type=str, help='output file base name')
args = parser.parse_args()

# aggdat = pd.DataFrame(columns=('parsimony forest size', 'mean allele frequency', 'mean branch length', 'MRCA distance to true tree', 'RF distance to true tree', 'trees with MRCA less than or equal to optimal tree'))

aggdat = pd.DataFrame(columns=('simulation', 'parsimony forest size', 'MRCA distance to true tree', 'RF distance to true tree', 'log-likelihood', 'ismle'))

rowct = 0
frames = []
for i, fname in enumerate(args.input):
    frames.append(pd.read_csv(fname, sep='\t', usecols=('MRCA', 'RF', 'log-likelihood')))
frames.sort(key=lambda df: len(df.index), reverse=True)
for i, df in enumerate(frames):
    mle = df['log-likelihood'].max()
    for _, row in df.iterrows():
        aggdat.loc[rowct] = (i, len(df.index), row['MRCA'], row['RF'], row['log-likelihood'], row['log-likelihood'] == mle)
        rowct += 1

aggdat.to_csv(args.outbase+'.tsv', sep='\t', index=False)

plt.figure(figsize=(20 ,6))
for i, metric in enumerate(('RF distance to true tree', 'MRCA distance to true tree'), 1):
    plt.subplot(2, 1, i)
    ax = sns.stripplot(x='simulation', y=metric, hue='ismle', data=aggdat[aggdat['parsimony forest size'] > 1], palette={False:'black', True:'red'}, jitter=True, alpha=.5, hue_order=[False, True])
    ax.legend_.remove()
    plt.ylim([-.5, None])
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off

# sns.regplot(x='parsimony forest size', y='trees with MRCA less than or equal to optimal tree',
#             data=aggdat, fit_reg=False, scatter_kws={'alpha':0.4})
# plt.xlim(.5, None)
# plt.ylim(.5, plt.xlim()[1]/2)
# linex = scipy.arange(1, int(aggdat['parsimony forest size'].max())+1)
# liney = (linex+1)/2
# plt.plot(linex, liney, '--k', lw=1)
# plt.xscale('log')
# plt.yscale('log')
plt.savefig(args.outbase+'.pdf')

# g = sns.pairplot(aggdat, kind='reg',
#                  aspect=1.5, plot_kws={'fit_reg':None})
# g.set(xlim=(0, None))
# g.axes[1, 0].set_ylim(-1, 101)
# g.axes[1, 0].set_ylim(-1, 101)
# g.axes[2, 0].set_ylim(-1, 101)
# g.savefig(args.outbase+'.pdf')
