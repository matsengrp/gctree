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
import argparse
from scipy.stats import percentileofscore
import seaborn as sns

parser = argparse.ArgumentParser(description='germinal center tree inference and simulation')
parser.add_argument('input', type=str, nargs='+', help='gctree.validation.tsv files')
parser.add_argument('--outbase', type=str, help='output file base name')
args = parser.parse_args()

aggdat = pd.DataFrame(columns=('parsimony forest size', 'trees with MRCA less than or equal to optimal tree'))

for i, fname in enumerate(args.input):
    df = pd.read_csv(fname, sep='\t')
    forest_size = len(df.index)
    aggdat.loc[i] = (forest_size, sum(x <= df['MRCA'][0] for x in df['MRCA']))
    #percentileofscore(df['MRCA'], df['MRCA'][0], kind='mean'))

aggdat.to_csv(args.outbase+'.tsv', sep='\t', index=False)

plt.figure()
sns.regplot(x='parsimony forest size', y='trees with MRCA less than or equal to optimal tree',
            data=aggdat, fit_reg=False, scatter_kws={'alpha':0.4})
plt.xlim(.5, None)
plt.ylim(.5, plt.xlim()[1]/2)
linex = scipy.arange(1, int(aggdat['parsimony forest size'].max())+1)
liney = (linex+1)/2
plt.plot(linex, liney, '--k', lw=1)
plt.xscale('log')
plt.yscale('log')
plt.savefig(args.outbase+'.pdf')

# g = sns.pairplot(aggdat, kind='reg',
#                  aspect=1.5, plot_kws={'fit_reg':None})
# g.set(xlim=(0, None))
# g.axes[1, 0].set_ylim(-1, 101)
# g.axes[1, 0].set_ylim(-1, 101)
# g.axes[2, 0].set_ylim(-1, 101)
# g.savefig(args.outbase+'.pdf')
