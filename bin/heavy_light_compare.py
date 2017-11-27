#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
comparison of heavy and light trees
'''

from __future__ import division, print_function
from gctree import CollapsedTree, CollapsedForest
import pandas as pd
import scipy
#import matplotlib
#matplotlib.use('PDF')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='white', color_codes=True)
import pickle
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from collections import defaultdict


def main():

    import argparse

    parser = argparse.ArgumentParser(description='compare heavy and light chain trees',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('heavy_forest', type=str, help='.p forest file')
    parser.add_argument('light_forest', type=str, help='.p forest file')
    parser.add_argument('heavy_fasta', type=str, help='fasta')
    parser.add_argument('light_fasta', type=str, help='fasta')
    parser.add_argument('heavy_log', type=str, help='log')
    parser.add_argument('light_log', type=str, help='log')
    parser.add_argument('outbase', type=str, help='output file base name')
    parser.add_argument('--naiveid', type=str, default='GL', help='fasta')

    # parser.add_argument('--outbase', type=str, required=True, help='output file base name')
    args = parser.parse_args()

    with open(args.heavy_forest, 'rb') as f:
        heavy_forest = pickle.load(f)
    with open(args.light_forest, 'rb') as f:
        light_forest = pickle.load(f)

    heavy_map = defaultdict(list)
    for seq in AlignIO.read(args.heavy_fasta, 'fasta'):
        if seq.id != args.naiveid:
            heavy_map[str(seq.seq).lower()].append(seq.id)
    light_map = defaultdict(list)
    for seq in AlignIO.read(args.light_fasta, 'fasta'):
        if seq.id != args.naiveid:
            light_map[str(seq.seq).lower()].append(seq.id)

    for heavy_tree in heavy_forest.forest:
        for node in list(heavy_tree.tree.traverse()):
            assert node.frequency == len(heavy_map[node.sequence.lower()])
            for cell_name in heavy_map[node.sequence.lower()]:
                node.add_child(name=cell_name)
    for light_tree in light_forest.forest:
        for node in list(light_tree.tree.traverse()):
            assert node.frequency == len(light_map[node.sequence.lower()])
            for cell_name in light_map[node.sequence.lower()]:
                node.add_child(name=cell_name)

    heavy_log = pd.read_csv(args.heavy_log, sep='\t', usecols=('logLikelihood',), squeeze=True)
    light_log = pd.read_csv(args.light_log, sep='\t', usecols=('logLikelihood',), squeeze=True)

    dat = pd.DataFrame(columns=('l_heavy', 'l_light', 'l_heavy+l_light', 'RF', 'mle'))
    # dat = scipy.zeros((len(heavy_forest.forest), len(light_forest.forest)), dtype=int)
    ct = 0
    for i, heavy_tree in enumerate(heavy_forest.forest):
        for j, light_tree in enumerate(light_forest.forest):
            # dat[i, j] = heavy_tree.tree.robinson_foulds(light_tree.tree, attr_t1='name', attr_t2='name', unrooted_trees=True)[0]
            dat.loc[ct] = [heavy_log[i],
                           light_log[j],
                           heavy_log[i]+light_log[j],
                           heavy_tree.tree.robinson_foulds(light_tree.tree, attr_t1='name', attr_t2='name', unrooted_trees=True)[0],
                           'both' if i == j == 0 else 'heavy' if i == 0 else 'light' if j == 0 else 'neither']
            ct += 1

    dat2 = pd.DataFrame(columns=('logLikelihood', 'chain', 'mle'))
    # dat = scipy.zeros((len(heavy_forest.forest), len(light_forest.forest)), dtype=int)
    ct = 0
    for i, heavy_tree in enumerate(heavy_log):
        dat2.loc[i] = [heavy_log[i], 'IgH\n(n={})'.format(len(heavy_log)), True if i == 0 else False]
    for j, light_tree in enumerate(light_log):
        dat2.loc[i + 1 + j] = [light_log[j], 'IgL\n(n={})'.format(len(light_log)), True if j == 0 else False]

    plt.figure(figsize=(1.75, 3))
    sns.boxplot(x='chain', y='logLikelihood', data=dat2[dat2['mle']==False], color='gray')
    ax = sns.stripplot(x='chain', y='logLikelihood', data=dat2[dat2['mle']==True], color='red')
    plt.xlabel('parsimony trees')
    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.tight_layout()
    plt.savefig(args.outbase+'.2.pdf')

    plt.figure(figsize=(4, 2.5))
    # sns.barplot(x=scipy.argsort(-dat['RF']), y='RF', hue='mle', data=dat, palette=('red', 'orange', 'yellow', 'gray'))
    plt.hist((dat['RF'][dat['mle']=='both'], dat['RF'][dat['mle']=='heavy'], dat['RF'][dat['mle']=='light'], dat['RF'][dat['mle']=='neither']),
             bins=scipy.arange(-.5, max(dat['RF'])+1.5, 1.),
             histtype='barstacked',
             label=('IgH and IgL MLE', 'IgH MLE', 'IgL MLE'),
             color=('red', 'indianred', 'rosybrown', 'gray'),
             edgecolor='black')
    plt.legend(loc=(0,.68))
    plt.xlim([min(dat['RF'])-.5, max(dat['RF'])+.5])
    plt.xlabel('RF distance between IgH tree and IgL tree\nlabelled by cell identity')
    plt.ylabel('IgH/IgL parsimony tree pairs\n(n={})'.format(len(dat)))
    plt.gca().spines['right'].set_color('none')
    plt.gca().spines['top'].set_color('none')
    plt.tight_layout()
    plt.savefig(args.outbase+'.pdf')

    # return

    # dat = pd.DataFrame(columns=('l_heavy', 'l_light', 'RF'))
    dat = scipy.zeros((len(light_forest.forest), len(heavy_forest.forest)), dtype=int)
    for j, heavy_tree in enumerate(heavy_forest.forest):
        for i, light_tree in enumerate(light_forest.forest):
            dat[i, j] = heavy_tree.tree.robinson_foulds(light_tree.tree, attr_t1='name', attr_t2='name', unrooted_trees=True)[0]
            # dat.loc[-1] = [heavy_log[i], light_log[j], heavy_tree.tree.robinson_foulds(light_tree.tree, attr_t1='name', attr_t2='name', unrooted_trees=True)[0]]
            # dat.index += 1

    plt.figure(figsize=(.3*len(heavy_forest.forest), .3*len(light_forest.forest)))
    # sns.lmplot(x='l_heavy', y='l_light', hue='RF', data=dat, fit_reg=False, scatter_kws={'alpha':.5}, legend_out=False)
    # for x in heavy_log:
    #     plt.axvline(x, ls='--', color='black', lw=.1, zorder=0)
    # for x in light_log:
    #     plt.axhline(x, ls='--', color='black', lw=.1, zorder=0)
    # plt.subplot(2,2,3)
    ax = sns.heatmap(dat, annot=True, fmt='d', cmap='Reds', cbar=False, linewidths=1)
    ax.invert_yaxis()
    # plt.subplot(2,2,4)
    # plt.plot(list(reversed(heavy_log)), list(reversed(range(len(heavy_log)))), ls='none', marker='o', clip_on=False, c='k')
    # plt.xlabel('GCtree log likelihood')
    # plt.subplot(2,2,1)
    # plt.plot(range(len(light_log)), light_log, ls='none', marker='o', clip_on=False, c='k')
    # plt.ylabel('GCtree log likelihood')
    plt.tight_layout()
    sns.despine()
    plt.savefig(args.outbase+'.3.pdf')

if __name__ == '__main__':
    main()
