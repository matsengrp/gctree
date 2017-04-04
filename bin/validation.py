#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
comparison of inference and simulated trees
'''

from __future__ import division, print_function
from gctree import CollapsedTree, CollapsedForest, hamming_distance
try:
    import cPickle as pickle
except:
    import pickle
import pandas as pd
import scipy
import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='whitegrid', color_codes=True)
import os
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
    node = [node for node in tree.traverse() if node.sequence == sequence and node.frequency > 0]
    try:
        assert(len(node) == 1)
    except Exception as e:
        print('Nodes list:')
        print(node)

        print(tree)
        print(node[0])
        print(node[1])

        raise e
    return node[0]


def align_lineages(seq, tree_t, tree_i, penalty_cap=None):
    assert(penalty_cap < 0)  # Penalties must be negative
    nt = find_node_by_seq(tree_t, seq)
    lt = reconstruct_lineage(tree_t, nt)
    ni = find_node_by_seq(tree_i, seq)
    li = reconstruct_lineage(tree_i, ni)
    # One lineages must be longer than just the root and the terminal node
    if len(lt) <= 2 and len(li) <= 2:
        return False

    # Gap penalty chosen not too large:
    gap_penalty = -10
    assert(gap_penalty < 0)  # Penalties must be negative
    if penalty_cap is not None and gap_penalty < penalty_cap:
        gap_penalty = penalty_cap

    # Regerate a score matrix matrix:
    kt = len(lt)
    ki = len(li)
    sc_mat = np.zeros((kt, ki), dtype=np.int64)
    for i in range(kt):
        for j in range(ki):
            # Notice the score is defined by number of mismatches:
            #sc_mat[i, j] = len(lt[i]) - fast_hamming_dist(lt[i], li[j])
            sc_mat[i, j] =  -1 * fast_hamming_dist(lt[i], li[j])

    # Calculate the alignment scores:
    aln_sc = np.zeros((kt+1, ki+1), dtype=np.int64)
    for i in range(0, kt+1):
        aln_sc[i][0] = gap_penalty * i
    for j in range(0, ki+1):
        aln_sc[0][j] = gap_penalty * j
    for i in range(1, kt+1):
        for j in range(1, ki+1):
            match = aln_sc[i-1][j-1] + sc_mat[i-1, j-1]
            delete = aln_sc[i-1][j] + gap_penalty
            insert = aln_sc[i][j-1] + gap_penalty
            aln_sc[i][j] = max(match, delete, insert)

    # Traceback to compute the alignment:
    align_t, align_i, asr_align = list(), list(), list()
    i, j = kt, ki
    while i > 0 and j > 0:
        sc_current = aln_sc[i][j]
        sc_diagonal = aln_sc[i-1][j-1]
        sc_up = aln_sc[i][j-1]
        sc_left = aln_sc[i-1][j]

        if sc_current == (sc_diagonal + sc_mat[i-1, j-1]):
            if penalty_cap is not None and sc_mat[i-1, j-1] > penalty_cap:
                asr_align.append(penalty_cap)
            else:
                asr_align.append(sc_mat[i-1, j-1])
            align_t.append(lt[i-1])
            align_i.append(li[j-1])
            i -= 1
            j -= 1
        elif sc_current == (sc_left + gap_penalty):
            asr_align.append(gap_penalty)
            align_t.append(lt[i-1])
            align_i.append('-')
            i -= 1
        elif sc_current == (sc_up + gap_penalty):
            asr_align.append(gap_penalty)
            align_t.append('-')
            align_i.append(li[j-1])
            j -= 1

    # If space left fill it with gaps:
    while i > 0:
        asr_align.append(gap_penalty)
        align_t.append(lt[i-1])
        align_i.append('-')
        i -= 1
    while j > 0:
        asr_align.append(gap_penalty)
        align_t.append('-')
        align_i.append(li[j-1])
        j -= 1

    max_penalty = 0
    penalty = sum(asr_align)
    for a in align_t:
        if a == '_':
            max_penalty += gap_penalty
        else:
            if penalty_cap is not None:
                max_penalty += penalty_cap
            else:
                max_penalty += -len(a)
    # Notice that the root and the terminal node is excluded from this comparison.
    # by adding their length to the max_penalty:
    max_penalty += 2 * len(lt[0])
    return [align_t, align_i, asr_align, penalty, max_penalty]


def lineage_dist(true_tree, inferred_tree, penalty_cap=None, freq_weigthing=False):
    total_lineage_dist = 0
    total_max_penalty = 0
    nlineages = 0
    for node in true_tree.tree.traverse():
        if not node.frequency > 0:
            continue

        aln_res = align_lineages(node.sequence, true_tree.tree, inferred_tree.tree, penalty_cap=None)
        if aln_res is  False:  # Skip lineages less than three members long
            continue
        align_t, align_i, asr_align, final_score, max_penalty = aln_res
        assert(sum(asr_align) == final_score)
        if freq_weigthing is True:
            total_max_penalty += max_penalty * node.frequency
            total_lineage_dist += final_score * node.frequency
        else:
            total_max_penalty += max_penalty
            total_lineage_dist += final_score

    if total_max_penalty == 0:  # There can be total_max_penalty == 0 when all lineages have less than three members
        return 0
    norm_lineage_dist = total_lineage_dist / total_max_penalty  # Normalize with max penalty to get a number between 0 and 1
    return norm_lineage_dist


def validate(true_tree, inferences, true_tree_colormap, outbase):
    '''
    inferences is a dict mapping infernce name, like "gctree" to pickle files of
    CollapsedForest
    '''

    # [(None, False), (1, False), (None, True), (1, True)]
    all_lineage_dist = lambda x, y: [lineage_dist(x, y, i1, i2) for i2 in [False, True] for i1 in [None, 1]]

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
                           'ld1':lineage_distances[0],
                           'ld2':lineage_distances[1],
                           'ld3':lineage_distances[2],
                           'ld4':lineage_distances[3],
                           'mean_frequency':mean_frequencies,
                           'mean_branch_length':mean_branch_lengths})

        if n_trees > 1:
            # plots
            plt.figure()#figsize=plt.figaspect(2))
            g = sns.pairplot(df, kind='reg', x_vars='log-likelihood', y_vars=('MRCA', 'RF'), plot_kws={'color':'black', 'scatter_kws':{'clip_on':False}, 'line_kws':{'alpha':.5}}, size=5)
            axes = g.axes
            axes[0,0].set_ylim(0,)
            axes[1,0].set_ylim(0,)
            plt.savefig(outbase+'.gctree.pdf')

        df.to_csv(outbase+'.gctree.tsv', sep='\t', index=False)

        for i, tree in enumerate(inferences['gctree'].forest, 1):
            colormap = dict((node.name, true_tree_colormap[node.sequence]) if node.sequence in true_tree_colormap else (node.name, 'lightgray') for node in tree.tree.traverse())
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
    df = pd.DataFrame({'method':methods, 'N_taxa':n_taxa, 'RF':distances, 'MRCA':MRCAs, 'ld1':lineage_distances[0], 'ld2':lineage_distances[1], 'ld3':lineage_distances[2], 'ld4':lineage_distances[3]},
                      columns=('method', 'N_taxa', 'RF', 'MRCA', 'ld1', 'ld2', 'ld3', 'ld4'))
    df.to_csv(outbase+'.tsv', sep='\t', index=False)


def main():

    import argparse

    parser = argparse.ArgumentParser(description='validate results of inference on simulation data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('true_tree', type=str, help='.p file containing true tree')
    parser.add_argument('true_tree_colormap', type=str, help='.tsv colormap file for true tree')
    parser.add_argument('forest_files', type=str, nargs='*', help='.p files containing forests from each inference method')
    parser.add_argument('--outbase', type=str, required=True, help='output file base name')
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
            true_tree_colormap[true_tree.tree.search_nodes(name=name)[0].sequence] = color

    validate(true_tree, inferences, true_tree_colormap, args.outbase)
    print('Done')  # Print something to the log to make the wait_func run smooth


if __name__ == '__main__':
    main()
