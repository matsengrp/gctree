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
from itertools import cycle


parser = argparse.ArgumentParser(description='aggregate validation of repeated runs with same parameters')
parser.add_argument('input', type=str, nargs='+', help='gctree.validation.tsv files')
parser.add_argument('--outbase', type=str, help='output file base name')
parser.add_argument('--allmetrics', action='store_true', default=False, help='use validation.tsv files containing other metrics')
args = parser.parse_args()

# aggdat = pd.DataFrame(columns=('parsimony forest size', 'mean allele frequency', 'mean branch length', 'MRCA distance to true tree', 'RF distance to true tree', 'trees with MRCA less than or equal to optimal tree'))

aggdat = pd.DataFrame(columns=('simulations ranked by parsimony degeneracy', 'parsimony forest size', 'MRCA distance to true tree', 'RF distance to true tree', 'COAR distance to true tree', 'COAR_fw distance to true tree', 'log-likelihood', 'ismle'))

rowct = 0
frames = []
for fname in args.input:
    df = pd.read_csv(fname, sep='\t', usecols=('COAR', 'COAR_fw', 'MRCA', 'RF', 'log-likelihood'))
    frames.append(df)
sorted_indices = sorted(range(len(frames)), key=lambda i: len(frames[i].index), reverse=True)
frames = [frames[i] for i in sorted_indices]

for i, df in enumerate(frames):
    mle = df['log-likelihood'].max()
    for _, row in df.iterrows():
        aggdat.loc[rowct] = (i, len(df.index), row['MRCA'], row['RF'], row['COAR'], row['COAR_fw'], row['log-likelihood'], row['log-likelihood'] == mle)
        rowct += 1

aggdat.to_csv(args.outbase+'.tsv', sep='\t', index=False)

if args.allmetrics:

    aggdat['method'] = 'parsimony'

    aggdat_other = pd.DataFrame(columns=('simulations ranked by parsimony degeneracy', 'method', 'MRCA distance to true tree', 'RF distance to true tree', 'COAR distance to true tree', 'COAR_fw distance to true tree'))

    frames_other = []
    for fname in [f.replace('.gctree', '') for f in args.input]:
        df = pd.read_csv(fname, sep='\t', usecols=('method', 'COAR', 'COAR_fw', 'MRCA', 'RF'))
        frames_other.append(df)
    frames_other = [frames_other[i] for i in sorted_indices]

    for i, df in enumerate(frames_other):
        for _, row in df.iterrows():
            aggdat_other.loc[rowct] = (i, row['method'], row['MRCA'], row['RF'], row['COAR'], row['COAR_fw'])
            rowct += 1

    aggdat_other.to_csv(args.outbase+'.other_metrics.tsv', sep='\t', index=False)

if not args.allmetrics:
    aggdat = aggdat[(aggdat['parsimony forest size'] >= 2)]
if aggdat.shape[0] == 0:
    print('no simulations with parsimony degeneracy!')
else:
    # plt.figure(figsize=(3, 12))
    # for metric in ('RF', 'MRCA'):
    #     x, y = zip(*[pearsonr(aggdat[aggdat['simulations ranked by parsimony degeneracy']==ct]['log-likelihood'], aggdat[aggdat['simulations ranked by parsimony degeneracy']==ct][metric+' distance to true tree']) for ct in set(aggdat['simulations ranked by parsimony degeneracy'])])
    #     plt.plot(x, -scipy.log10(y), 'o', alpha=.75, label=metric, clip_on=False)
    # #sns.regplot('Pearson r', '-log10 P', aggdat, fit_reg=False, color='black', scatter_kws={'clip_on':False})
    # plt.xlabel('correlation of log-likelihood with\ndistance from true tree (Pearson r)')
    # plt.ylabel('-log10 p-value')
    # plt.xlim([-1, 1])
    # plt.legend(frameon=True)
    # #plt.axvline(linewidth=1, ls='--', color='gray')
    # plt.gca().spines['left'].set_position('center')
    # plt.gca().spines['right'].set_color('none')
    # plt.gca().spines['top'].set_color('none')
    # plt.savefig(args.outbase+'.volcano.pdf')

    # for metric in ('RF', 'MRCA', 'COAR', 'COAR_fw'):
    #     plt.figure()
    #     g = sns.FacetGrid(aggdat, col='simulations ranked by parsimony degeneracy', size=3, ylim=[0, aggdat[metric+' distance to true tree'].max()], sharex=False)#, sharey=False)
    #     g.map(sns.regplot, 'log-likelihood', metric+' distance to true tree', scatter_kws={'color':'black', 'alpha':.5}, line_kws={'color':'grey', 'alpha':.5})
    #     g.set(clip_on=False, title='')
    #     g.fig.subplots_adjust(wspace=.15, hspace=.15)
    #     #for ax in g.axes[0]:
    #     #    ax.set_ylim([0, None])
    #     plt.savefig(args.outbase+'.'+metric+'.pdf')

    maxy = {}

    plt.figure(figsize=(7, 12))
    for i, metric in enumerate(('RF distance to true tree', 'MRCA distance to true tree', 'COAR distance to true tree', 'COAR_fw distance to true tree'), 1):
        maxy[metric] = 1.1*aggdat[metric].max()
        plt.subplot(4, 1, i)
        if not args.allmetrics:
            sns.boxplot(x='simulations ranked by parsimony degeneracy', y=metric, data=aggdat[aggdat['ismle']==False], color='gray')
            sns.stripplot(x='simulations ranked by parsimony degeneracy', y=metric, color='red', data=aggdat[aggdat['ismle']])
        else:
            # sns.boxplot(x='simulations ranked by parsimony degeneracy', y=metric, data=aggdat)#, color='gray')
            #sns.swarmplot(x='simulations ranked by parsimony degeneracy', y=metric, hue='method', size=2, data=pd.concat([aggdat, aggdat_other]))
            palette = cycle(list(sns.color_palette()))
            sns.boxplot(x='simulations ranked by parsimony degeneracy', y=metric, data=aggdat, color=next(palette))
            sns.stripplot(x='simulations ranked by parsimony degeneracy', y=metric, size=4, alpha=.7, jitter=2, color=next(palette), data=aggdat_other[aggdat_other['method'] == 'gctree'])
            sns.stripplot(x='simulations ranked by parsimony degeneracy', y=metric, size=4, alpha=.7, jitter=2, color=next(palette), data=aggdat_other[aggdat_other['method'] == 'dnaml'])
            sns.stripplot(x='simulations ranked by parsimony degeneracy', y=metric, size=4, alpha=.7, jitter=2, color=next(palette), data=aggdat_other[aggdat_other['method'] == 'igphyml'])
        if 'COAR' in metric:
            plt.ylim(0, maxy[metric])
        else:
            plt.ylim(0, maxy[metric])
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

    plt.figure(figsize=(2, 12))
    #bins = {}
    for i, metric in enumerate(('RF distance to true tree', 'MRCA distance to true tree', 'COAR distance to true tree', 'COAR_fw distance to true tree'), 1):
        plt.subplot(4, 1, i)
        if not args.allmetrics:
            sns.boxplot(x='ismle', y=metric, data=aggdat, palette={False:'gray', True:'red'}, order=(True, False))
        else:
            ax = sns.boxplot(x='method', y=metric, data=pd.concat([aggdat, aggdat_other]), order=('parsimony', 'gctree', 'dnaml', 'igphyml'))#, palette={False:'gray', True:'red'}, order=(True, False))
            #ax = sns.boxplot(x='method', y=metric, data=aggdat, order=('parsimony', 'gctree', 'dnaml', 'igphyml'))#, palette={False:'gray', True:'red'}, order=(True, False))
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
        if 'COAR' in metric:
            plt.ylim(0, maxy[metric])
        else:
            plt.ylim(0, maxy[metric])
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
    plt.savefig(args.outbase+'.boxplot_all.pdf')

    if args.allmetrics:
        methods = list(set(aggdat_other['method']))
        nmethods = len(methods)
        for metric in ('RF distance to true tree', 'MRCA distance to true tree', 'COAR distance to true tree', 'COAR_fw distance to true tree'):
            for i in range(nmethods):
                for j in range(i+1, nmethods):
                    p = mannwhitneyu(aggdat_other[metric][aggdat_other['method']==methods[i]], aggdat_other[metric][aggdat_other['method']==methods[j]], alternative='two-sided')[1]
                    print('U test significance ({}, {} vs {}): {}'.format(metric, methods[i], methods[j], p))
                p = mannwhitneyu(aggdat_other[metric][aggdat_other['method']==methods[i]], aggdat[metric], alternative='two-sided')[1]
                print('U test significance ({}, {} vs parsimony): {}'.format(metric, methods[i], p))

    else:
        for metric in ('RF distance to true tree', 'MRCA distance to true tree', 'COAR distance to true tree', 'COAR_fw distance to true tree'):
            if len(set(aggdat[metric])) > 1:
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
