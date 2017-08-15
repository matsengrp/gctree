#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
comparison of inference and simulated trees
'''

from __future__ import division, print_function
from gctree import CollapsedTree, CollapsedForest
import gctree
from utils import hamming_distance
from random import randint
try:
    import cPickle as pickle
except:
    import pickle
import pandas as pd
import scipy
import matplotlib
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='white', color_codes=True)
import os, sys
import numpy as np
try:
    import jellyfish  # Last time I checked this module was the fastest kid on the block, however 2x slower than a simple Cython function
    def fast_hamming_dist(s1, s2):
        if s1 == s2:
            return 0
        else:
            return jellyfish.hamming_distance(unicode(s1), unicode(s2))
except:
    print('Couldn\'t find the python module "jellyfish" which is used for fast string comparison. Falling back to pure python function.')
    fast_hamming_dist = hamming_distance


def reconstruct_lineage(tree, node):
    lineage = list()
    while True:
        lineage.append(node.sequence)
        if node.up is None:
            return lineage
        node = node.up


def find_node_by_seq(tree, sequence):
    nodes = [node for node in tree.traverse() if node.sequence == sequence and node.frequency > 0]
    if len(nodes) > 1:
        nodes = [nodes[randint(0, len(nodes)-1)]]
    try:
        assert(len(nodes) == 1)
    except Exception as e:
        print('Nodes list:')
        print(nodes)
        print(sequence)
        print(tree)
        print([(node.frequency, node.name, node.sequence) for node in tree.traverse()])
        print(nodes[0])
        print(nodes[1])

        raise e
    return nodes[0]


def align_lineages(seq, tree_t, tree_i, gap_penalty_pct=0, known_root=True, allow_double_gap=False):
    '''
    Standard implementation of a Needleman-Wunsch algorithm as described here:
    http://telliott99.blogspot.com/2009/08/alignment-needleman-wunsch.html
    https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
    And implemented here:
    https://github.com/alevchuk/pairwise-alignment-in-python/blob/master/alignment.py

    gap_penalty_pct is the gap penalty relative to the sequence length of the sequences on the tree.
    '''
    nt = find_node_by_seq(tree_t, seq)
    lt = reconstruct_lineage(tree_t, nt)
    ni = find_node_by_seq(tree_i, seq)
    li = reconstruct_lineage(tree_i, ni)
    # One lineages must be longer than just the root and the terminal node
    if len(lt) <= 2 and len(li) <= 2:
        return False

    # Gap penalty chosen not too large:
    gap_penalty = -1 * int((len(seq) / 100.0) * gap_penalty_pct)
    assert(gap_penalty <= 0)  # Penalties must be negative
    if gap_penalty == 0:  # If gap penalty is zero only gaps in the shortes sequence will be allowed
        assert(allow_double_gap is False)

    # Generate a score matrix matrix:
    kt = len(lt)
    ki = len(li)
    # Disallow gaps in the longest list:
    if allow_double_gap is False and kt > ki:
        # If true is longer than inferred allow gap only in inferred:
        gap_penalty_i = gap_penalty
        gap_penalty_j = -1 * float('inf')
    elif allow_double_gap is False and kt < ki:
        # If inferred is longer than true allow gap only in true:
        gap_penalty_i = -1 * float('inf')
        gap_penalty_j = gap_penalty
    elif allow_double_gap is False and kt == ki:
        # If lists are equally long no gaps are allowed:
        gap_penalty_i = -1 * float('inf')
        gap_penalty_j = -1 * float('inf')
    else:
        gap_penalty_i = gap_penalty
        gap_penalty_j = gap_penalty

    sc_mat = np.zeros((kt, ki), dtype=np.float64)
    for i in range(kt):
        for j in range(ki):
            # Notice the score is defined by number of mismatches:
            #sc_mat[i, j] = len(lt[i]) - fast_hamming_dist(lt[i], li[j])
            sc_mat[i, j] = -1 * fast_hamming_dist(lt[i], li[j])

###    print(sc_mat)
    # Calculate the alignment scores:
    aln_sc = np.zeros((kt+1, ki+1), dtype=np.float64)
    for i in range(0, kt+1):
        if known_root is True:
            aln_sc[i][0] = -1 * float('inf')
        else:
            aln_sc[i][0] = gap_penalty_i * i
    for j in range(0, ki+1):
        if known_root is True:
            aln_sc[0][j] = -1 * float('inf')
        else:
            aln_sc[0][j] = gap_penalty_j * j
    aln_sc[0][0] = 0  # The top left is fixed to zero
###    print(aln_sc)
    for i in range(1, kt+1):
        for j in range(1, ki+1):
            match = aln_sc[i-1][j-1] + sc_mat[i-1, j-1]
            gap_in_inferred = aln_sc[i-1][j] + gap_penalty_i
            gap_in_true = aln_sc[i][j-1] + gap_penalty_j
            aln_sc[i][j] = max(match, gap_in_inferred, gap_in_true)
###    print(aln_sc)
    # Traceback to compute the alignment:
    align_t, align_i, asr_align = list(), list(), list()
    i, j = kt, ki
    alignment_score = aln_sc[i][j]
    while i > 0 and j > 0:
        sc_current = aln_sc[i][j]
        sc_diagonal = aln_sc[i-1][j-1]
        sc_up = aln_sc[i][j-1]
        sc_left = aln_sc[i-1][j]

        if sc_current == (sc_diagonal + sc_mat[i-1, j-1]):
            align_t.append(lt[i-1])
            align_i.append(li[j-1])
            i -= 1
            j -= 1
        elif sc_current == (sc_left + gap_penalty_i):
            align_t.append(lt[i-1])
            align_i.append('-')
            i -= 1
        elif sc_current == (sc_up + gap_penalty_j):
            align_t.append('-')
            align_i.append(li[j-1])
            j -= 1

    # If space left fill it with gaps:
    while i > 0:
        asr_align.append(gap_penalty_i)
        align_t.append(lt[i-1])
        align_i.append('-')
        i -= 1
    while j > 0:
        asr_align.append(gap_penalty_j)
        align_t.append('-')
        align_i.append(li[j-1])
        j -= 1

    max_penalty = 0
    for a, b in zip(align_t, align_i):
        if a == '-' or b == '-':
            max_penalty += gap_penalty
        else:
            max_penalty += -len(a)
    # Notice that the root and the terminal node is excluded from this comparison.
    # by adding their length to the max_penalty:
    if known_root is True:
        max_penalty += 2 * len(lt[0])
    else:  # Or in the case of an unknown root, just add the terminal node
        max_penalty += len(lt[0])

    return [align_t, align_i, alignment_score, max_penalty]


def lineage_dist(true_tree, inferred_tree, freq_weigthing=False, known_root=True, allow_double_gap=False):
    norm_lineage_dist = list()
    nlineages = 0
#    for node in true_tree.tree.traverse():
    for node in true_tree.tree.iter_leaves():  # Iterate only through the leaves
        if not node.frequency > 0:
            continue

        aln_res = align_lineages(node.sequence, true_tree.tree, inferred_tree.tree, known_root=known_root, allow_double_gap=allow_double_gap)
        if aln_res is  False:  # Skip lineages less than three members long
            continue
        align_t, align_i, final_score, max_penalty = aln_res
        if freq_weigthing is True:
            total_max_penalty = max_penalty * node.frequency
            total_lineage_dist = final_score * node.frequency
            # Normalize with the max penalty:
            if total_max_penalty < 0:
                norm_lineage_dist.append(total_lineage_dist/total_max_penalty)
        else:
            total_max_penalty = max_penalty
            total_lineage_dist = final_score
            # Normalize with the max penalty:
            if total_max_penalty < 0:
                norm_lineage_dist.append(total_lineage_dist/total_max_penalty)

    if len(norm_lineage_dist) == 0:  # There can be total_max_penalty == 0 when all lineages have less than three members
        return 0
    # Take the mean of the distances:
    mean_norm_lineage_dist = sum(norm_lineage_dist) / len(norm_lineage_dist)
    return mean_norm_lineage_dist


def validate(true_tree, inferences, true_tree_colormap, outbase):
    '''
    inferences is a dict mapping inference name, like "gctree" to pickle files of
    CollapsedForest
    '''

    # With/without frequency weighting:
    all_lineage_dist = lambda x, y: [lineage_dist(x, y, freq_weigthing=fw) for fw in [False, True]]

    # if gctre is among the inferences, let's evaluate the likelihood ranking
    # among the parsimony trees
    if 'gctree' in inferences:
        n_trees = len(inferences['gctree'].forest)
        # note: the unrooted_trees flag is needed because, for some reason, the RF
        #       function sometimes thinks the collapsed trees are unrooted and barfs
        distances, likelihoods = zip(*[(true_tree.compare(tree, method='RF'),
                                        tree.l(inferences['gctree'].params)[0]) for tree in inferences['gctree'].forest])
        MRCAs = [true_tree.compare(tree, method='MRCA') for tree in inferences['gctree'].forest]
        lineage_distances = [all_lineage_dist(true_tree, tree) for tree in inferences['gctree'].forest]
        lineage_distances = zip(*lineage_distances)  # Unzip the forest tuple to get lineage_distances[ld0-3][tree_n]
        mean_frequencies = [scipy.mean([node.frequency for node in tree.tree.traverse()]) for tree in inferences['gctree'].forest]
        mean_branch_lengths = [scipy.mean([node.dist for node in tree.tree.iter_descendants()]) for tree in inferences['gctree'].forest]
        df = pd.DataFrame({'log-likelihood':likelihoods,
                           'RF':distances,
                           'MRCA':MRCAs,
                           'COAR':lineage_distances[0],
                           'COAR_fw':lineage_distances[1],
                           'mean_frequency':mean_frequencies,
                           'mean_branch_length':mean_branch_lengths})

        if n_trees > 1:
            # plots
            maxll = df['log-likelihood'].max()
            if len(df[df['log-likelihood']!=maxll]) >= 2:
                fit_reg = True
            else:
                fit_reg = False
            plt.figure(figsize=(10, 10))
            for i, metric in enumerate(('RF', 'MRCA', 'COAR', 'COAR_fw'), 1):
                plt.subplot(2, 2, i)
                ax = sns.regplot('log-likelihood', metric, data=df[df['log-likelihood']!=maxll], fit_reg=fit_reg, color='black', scatter_kws={'alpha':.8, 'clip_on':False})
                sns.regplot('log-likelihood', metric, data=df[df['log-likelihood']==maxll], fit_reg=False, color='red', scatter_kws={'alpha':.8, 'clip_on':False}, ax=ax)
                plt.ylim(0, 1.1*df[metric].max())
                plt.xlim(df['log-likelihood'].min(), df['log-likelihood'].max())
                plt.tight_layout()
            plt.savefig(outbase+'.gctree.pdf')

        df.to_csv(outbase+'.gctree.tsv', sep='\t', index=False)

        plt.figure(figsize=(14, 14))
        sns.pairplot(df, kind="reg")
        plt.savefig(outbase+'_pairplot.pdf')

        for i, tree in enumerate(inferences['gctree'].forest, 1):
            colormap = {}
            for node in tree.tree.traverse():
                if node.sequence in true_tree_colormap:
                    colormap[node.name] = true_tree_colormap[node.sequence]
                else:
                    assert node.frequency == 0
                    colormap[node.name] = 'lightgray'
            tree.render(outbase+'.gctree.colored_tree.{}.svg'.format(i), colormap=colormap)




    # compare the inference methods
    # assume the first tree in the forest is the inferred tree

    methods, n_taxa, distances, MRCAs, lineage_distances = zip(
        *[(method,
           len(list(true_tree.tree.traverse())),  # Get all taxa in the tree
           true_tree.compare(inferences[method].forest[0], method='RF'),
           true_tree.compare(inferences[method].forest[0], method='MRCA'),
           all_lineage_dist(true_tree, inferences[method].forest[0])) for method in inferences])
    lineage_distances = zip(*lineage_distances)  # Unzip the methods tuple to get lineage_distances[ld0-3][method]
    df = pd.DataFrame({'method':methods, 'N_taxa':n_taxa, 'RF':distances, 'MRCA':MRCAs, 'COAR':lineage_distances[0], 'COAR_fw':lineage_distances[1]},
                      columns=('method', 'N_taxa', 'RF', 'MRCA', 'COAR', 'COAR_fw'))
    df.to_csv(outbase+'.tsv', sep='\t', index=False)


def main():

    import argparse

    parser = argparse.ArgumentParser(description='validate results of inference on simulation data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('true_tree', type=str, help='.p file containing true tree')
    parser.add_argument('true_tree_colormap', type=str, help='.tsv colormap file for true tree')
    parser.add_argument('forest_files', type=str, nargs='*', help='.p files containing forests from each inference method')
    # Uncomment and insert into function calls to allow for these settings:
    # parser.add_argument('--known_root', type=bool, default=True, help='This is root sequence known in the inferred phylogeny? If yes, the alignment is forced to end in the top left corner of the alignment grid.')
    # parser.add_argument('--allow_double_gap', type=bool, default=False, help='Allow the Needleman-Wunsch algorithm to move freely and possibly introducing gaps in both true and inferred list under the same comparison.')
    parser.add_argument('--outbase', type=str, required=True, help='output file base name')
    global args
    args = parser.parse_args()

    with open(args.true_tree, 'rb') as f:
        true_tree = pickle.load(f)
    inferences = {}
    for forest_file in args.forest_files:
        with open(forest_file, 'rb') as f:
            inferences[os.path.basename(forest_file).split('.')[0]] = pickle.load(f)
    # now we rerender the inferred trees, but using colormap from true tree, makes visual comaprison easier
    true_tree_colormap = {} # map for tree sequences
    with open(args.true_tree_colormap, 'r') as f:
        for line in f:
            name, color = line.rstrip().split('\t')
            if ',' not in name:
                true_tree_colormap[true_tree.tree.search_nodes(name=name)[0].sequence] = color
            else:
                search = true_tree.tree.search_nodes(name=tuple(name.split(',')))
                if search:
                    true_tree_colormap[search[0].sequence] = color

    validate(true_tree, inferences, true_tree_colormap, args.outbase)
    print('Done')  # Print something to the log to make the wait_func run smooth


if __name__ == '__main__':
    main()
