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
import os
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
    node = [node for node in tree.traverse() if node.sequence == sequence]
    assert(len(node) == 1)
    return node[0]


def align_lineages(seq, tree_t, tree_i):
    nt = find_node_by_seq(tree_t, seq)
    lt = reconstruct_lineage(tree_t, nt)
    ni = find_node_by_seq(tree_i, seq)
    li = reconstruct_lineage(tree_i, ni)
    # SW parameters:
    gap_penalty = -10

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
    for a in align_t:
        if a == '_':
            max_penalty += gap_penalty
        else:
            max_penalty += -len(a)
    return align_t, align_i, asr_align, aln_sc[kt, ki], max_penalty


def lineage_dist(true_tree, inferred_tree):
    total_lineage_dist = 0
    leaf_seqs = [node.sequence for node in true_tree.tree.traverse() if node.is_leaf()]
    for leaf_seq in leaf_seqs:
        align_t, align_i, asr_align, final_score, max_penalty = align_lineages(leaf_seq, true_tree.tree, inferred_tree.tree)
        total_lineage_dist += final_score / max_penalty  # Normalize with max penalty to get a number between 0 and 1
        assert(sum(asr_align) == final_score)
    norm_lineage_dist = total_lineage_dist / len(leaf_seqs)
    return total_lineage_dist


def MRCA_distance(true_tree, tree):
    '''
    here's Erick's idea of matrix of hamming distance of common ancestors of taxa
    takes a true and inferred tree as CollapsedTree objects
    '''
    taxa = [node.sequence for node in true_tree.tree.traverse() if node.frequency]
    n_taxa = len(taxa)
    d = scipy.zeros(shape=(n_taxa, n_taxa))
    for i in range(n_taxa):
        nodei_true = true_tree.tree.iter_search_nodes(sequence=taxa[i]).next()
        nodei      =      tree.tree.iter_search_nodes(sequence=taxa[i]).next()
        for j in range(i + 1, n_taxa):
            nodej_true = true_tree.tree.iter_search_nodes(sequence=taxa[j]).next()
            nodej      =      tree.tree.iter_search_nodes(sequence=taxa[j]).next()
            MRCA_true = true_tree.tree.get_common_ancestor((nodei_true, nodej_true)).sequence
            MRCA =           tree.tree.get_common_ancestor((nodei, nodej)).sequence
            d[i, j] = nodei.frequency*nodej.frequency*hamming_distance(MRCA_true, MRCA)
    return d


def validate(true_tree, inferences, outbase):
    '''
    inferences is a dict mapping infernce name, like "gctree" to pickle files of
    CollapsedForest
    '''

    # if gctre is among the inferences, let's evaluate the likelihood ranking
    # among the parsimony trees
    if 'gctree' in inferences:
        n_trees = len(inferences['gctree'].forest)
        # NOTE: the unrooted_trees flag is needed because, for some reason, the RF
        #       function sometimes thinks the collapsed trees are unrooted and barfs
        distances, likelihoods = zip(*[(true_tree.tree.robinson_foulds(tree.tree, attr_t1='sequence', attr_t2='sequence', unrooted_trees=True)[0],
                                        tree.l(inferences['gctree'].params)[0]) for tree in inferences['gctree'].forest])
        MRCAs = [MRCA_distance(true_tree, tree).sum() for tree in inferences['gctree'].forest]
        lineage_distance = [lineage_dist(true_tree, tree) for tree in inferences['gctree'].forest]

        mean_frequencies = [scipy.mean([node.frequency for node in tree.tree.traverse()]) for tree in inferences['gctree'].forest]
        mean_branch_lengths = [scipy.mean([node.dist for node in tree.tree.iter_descendants()]) for tree in inferences['gctree'].forest]
        df = pd.DataFrame({'log-likelihood':likelihoods,
                           'RF':distances,
                           'MRCA':MRCAs,
                           'lineage_distance':lineage_distance,
                           'mean_frequency':mean_frequencies,
                           'mean_branch_length':mean_branch_lengths})

        if n_trees > 1:
            # plots
            plt.figure()
            sns.pairplot(df, kind='reg', x_vars='log-likelihood', y_vars=('MRCA', 'RF'), aspect=1.5)
            plt.savefig(outbase+'.gctree.pdf')

        df.to_csv(outbase+'.gctree.tsv', sep='\t', index=False)

    # compare the inference methods
    # assume the first tree in the forest is the inferred tree
    methods, n_taxa, distances, MRCAs, lineage_distance = zip(
        *[(method,
           len(list(true_tree.tree.traverse())),  # Get all taxa in the tree
           true_tree.tree.robinson_foulds(inferences[method].forest[0].tree,
                                          attr_t1='sequence',
                                          attr_t2='sequence',
                                          unrooted_trees=True)[0],
           MRCA_distance(true_tree, inferences[method].forest[0]).sum(),
           lineage_dist(true_tree, inferences[method].forest[0]))
           for method in inferences])
    df = pd.DataFrame({'method':methods, 'N_taxa':n_taxa, 'RF':distances, 'MRCA':MRCAs, 'lineage_distance':lineage_distance},
                      columns=('method', 'N_taxa', 'RF', 'MRCA', 'lineage_distance'))
    df.to_csv(outbase+'.tsv', sep='\t', index=False)


def main():

    import argparse

    parser = argparse.ArgumentParser(description='validate results of inference on simulation data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('true_tree', type=str, help='.p file containing true tree')
    parser.add_argument('forest_files', type=str, nargs='*', help='.p files containing forests from each inference method')
    parser.add_argument('--outbase', type=str, required=True, help='output file base name')
    args = parser.parse_args()

    with open(args.true_tree, 'rb') as f:
        true_tree = pickle.load(f)
    inferences = {}
    for forest_file in args.forest_files:
        with open(forest_file, 'rb') as f:
            inferences[os.path.basename(forest_file).split('.')[0]] = pickle.load(f)
    validate(true_tree, inferences, args.outbase)
    print('Done')  # Print something to the log to make the wait_func run smooth


if __name__ == '__main__':
    main()
