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
from scipy.stats import pearsonr, mannwhitneyu

parser = argparse.ArgumentParser(description='aggregate validation of repeated runs with same parameters')
parser.add_argument('input', type=str, nargs='+', help='gctree.validation.tsv files')
parser.add_argument('--outbase', type=str, help='output file base name')
args = parser.parse_args()

# aggdat = pd.DataFrame(columns=('parsimony forest size', 'mean allele frequency', 'mean branch length', 'MRCA distance to true tree', 'RF distance to true tree', 'trees with MRCA less than or equal to optimal tree'))

aggdat = pd.DataFrame(columns=('simulations ranked by parsimony degeneracy', 'parsimony forest size', 'MRCA distance to true tree', 'RF distance to true tree', 'log-likelihood', 'ismle'))

rowct = 0
frames = []
for i, fname in enumerate(args.input):
    df = pd.read_csv(fname, sep='\t', usecols=('MRCA', 'RF', 'log-likelihood'))
    frames.append(df)
frames.sort(key=lambda df: len(df.index), reverse=True)
for i, df in enumerate(frames):
    mle = df['log-likelihood'].max()
    for _, row in df.iterrows():
        aggdat.loc[rowct] = (i, len(df.index), row['MRCA'], row['RF'], row['log-likelihood'], row['log-likelihood'] == mle)
        rowct += 1

aggdat.to_csv(args.outbase+'.tsv', sep='\t', index=False)

aggdat = aggdat[(aggdat['parsimony forest size'] >= 2)]

plt.figure(figsize=(3, 6))
for metric in ('RF', 'MRCA'):
    x, y = zip(*[pearsonr(aggdat[aggdat['simulations ranked by parsimony degeneracy']==ct]['log-likelihood'], aggdat[aggdat['simulations ranked by parsimony degeneracy']==ct][metric+' distance to true tree']) for ct in set(aggdat['simulations ranked by parsimony degeneracy'])])
    plt.plot(x, -scipy.log10(y), 'o', alpha=.75, label=metric, clip_on=False)
#sns.regplot('Pearson r', '-log10 P', aggdat, fit_reg=False, color='black', scatter_kws={'clip_on':False})
plt.xlabel('correlation of log-likelihood with\ndistance from true tree (Pearson r)')
plt.ylabel('-log10 p-value')
plt.xlim([-1, 1])
plt.legend(frameon=True)
#plt.axvline(linewidth=1, ls='--', color='gray')
plt.gca().spines['left'].set_position('center')
plt.gca().spines['right'].set_color('none')
plt.gca().spines['top'].set_color('none')
plt.savefig(args.outbase+'.volcano.pdf')

for metric in ('RF', 'MRCA'):
    plt.figure()
    g = sns.FacetGrid(aggdat, col='simulations ranked by parsimony degeneracy', size=3, ylim=[0, aggdat[metric+' distance to true tree'].max()], sharex=False)#, sharey=False)
    g.map(sns.regplot, 'log-likelihood', metric+' distance to true tree', scatter_kws={'color':'black', 'alpha':.5}, line_kws={'color':'grey', 'alpha':.5})
    g.set(clip_on=False, title='')
    g.fig.subplots_adjust(wspace=.15, hspace=.15)
    #for ax in g.axes[0]:
    #    ax.set_ylim([0, None])
    plt.savefig(args.outbase+'.'+metric+'.pdf')

maxy = {}

plt.figure(figsize=(7, 6))
for i, metric in enumerate(('RF distance to true tree', 'MRCA distance to true tree'), 1):
    maxy[metric] = 1.1*aggdat[metric].max()
    plt.subplot(2, 1, i)
    sns.boxplot(x='simulations ranked by parsimony degeneracy', y=metric, data=aggdat[aggdat['ismle']==False], color='gray')
    sns.stripplot(x='simulations ranked by parsimony degeneracy', y=metric, color='red', data=aggdat[aggdat['ismle']])
    plt.ylim(-.5, maxy[metric])
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off
    # plt.gca().spines['left'].set_color('none')
    # plt.gca().yaxis.tick_right()
    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.tight_layout()
    plt.savefig(args.outbase+'.boxplot.pdf')

plt.figure(figsize=(2, 6))
#bins = {}
for i, metric in enumerate(('RF distance to true tree', 'MRCA distance to true tree'), 1):
    plt.subplot(2, 1, i)
    sns.boxplot(x='ismle', y=metric, data=aggdat, palette={False:'gray', True:'red'}, order=(True, False))
    #binmax = int(aggdat[metric].max()+1)
    #n, bins[metric], patches = plt.hist([aggdat[metric][aggdat['ismle']],
    #                                     aggdat[metric][aggdat['ismle']==False]],
    #                                    bins=range(binmax) if binmax <= 20 else 20,
    #                                    color=['red', 'gray'],
    #                                    label=['inferred GCtrees', 'other parsimony trees'],
    #                                    histtype='barstacked',
    #                                    align='left',
    #                                    orientation='horizontal')
    # plt.legend(frameon=True)
    #plt.ylabel(metric)
    #plt.ylim(bins[metric][0] - .5*(bins[metric][1] - bins[metric][0]), bins[metric][-1])
    plt.ylim(-.5, maxy[metric])
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off
    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.xlabel('aggregated')
    plt.tight_layout()
plt.savefig(args.outbase+'.histogram.pdf')

for metric in ('RF distance to true tree', 'MRCA distance to true tree'):
    p = mannwhitneyu(aggdat[metric][aggdat['ismle']==True], aggdat[metric][aggdat['ismle']==False], alternative='less')[1]
    print('U test significance ({}): {}'.format(metric, p))


# sns.regplot(x='parsimony forest size', y='trees with MRCA less than or equal to optimal tree',
#             data=aggdat, fit_reg=False, scatter_kws={'alpha':0.4})
# plt.xlim(.5, None)
# plt.ylim(.5, plt.xlim()[1]/2)
# linex = scipy.arange(1, int(aggdat['parsimony forest size'].max())+1)
# liney = (linex+1)/2
# plt.plot(linex, liney, '--k', lw=1)
# plt.xscale('log')
# plt.yscale('log')
#plt.savefig(args.outbase+'.pdf')

# g = sns.pairplot(aggdat, kind='reg',
#                  aspect=1.5, plot_kws={'fit_reg':None})
# g.set(xlim=(0, None))
# g.axes[1, 0].set_ylim(-1, 101)
# g.axes[1, 0].set_ylim(-1, 101)
# g.axes[2, 0].set_ylim(-1, 101)
# g.savefig(args.outbase+'.pdf')
