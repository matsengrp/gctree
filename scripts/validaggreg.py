#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
aggregation plots of validation output from several simulation/validation runs
'''

from matplotlib import pyplot as plt
import pandas as pd
import argparse
import seaborn as sns
sns.set(style="white", color_codes=True)
from scipy.stats import mannwhitneyu
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
    maxy = {}

    plt.figure(figsize=(7, 12))
    for i, metric in enumerate(('RF distance to true tree', 'MRCA distance to true tree', 'COAR distance to true tree', 'COAR_fw distance to true tree'), 1):
        maxy[metric] = 1.1*aggdat[metric].max()
        plt.subplot(4, 1, i)
        if not args.allmetrics:
            sns.boxplot(x='simulations ranked by parsimony degeneracy', y=metric, data=aggdat[aggdat['ismle']==False], color='gray')
            sns.stripplot(x='simulations ranked by parsimony degeneracy', y=metric, color='red', data=aggdat[aggdat['ismle']])
        else:
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
        plt.gca().spines['right'].set_color('none')
        plt.gca().spines['top'].set_color('none')
        plt.tight_layout()
        plt.savefig(args.outbase+'.boxplot.pdf')

    plt.figure(figsize=(2, 12))
    for i, metric in enumerate(('RF distance to true tree', 'MRCA distance to true tree', 'COAR distance to true tree', 'COAR_fw distance to true tree'), 1):
        plt.subplot(4, 1, i)
        if not args.allmetrics:
            sns.boxplot(x='ismle', y=metric, data=aggdat, palette={False:'gray', True:'red'}, order=(True, False))
        else:
            ax = sns.boxplot(x='method', y=metric, data=pd.concat([aggdat, aggdat_other]), order=('parsimony', 'gctree', 'dnaml', 'igphyml'))
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
