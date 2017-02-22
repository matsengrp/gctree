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

        mean_frequencies = [scipy.mean([node.frequency for node in tree.tree.traverse()]) for tree in inferences['gctree'].forest]
        mean_branch_lengths = [scipy.mean([node.dist for node in tree.tree.iter_descendants()]) for tree in inferences['gctree'].forest]
        df = pd.DataFrame({'log-likelihood':likelihoods,
                           'RF':distances,
                           'MRCA':MRCAs,
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
    methods, n_taxa, distances, MRCAs = zip(*[(method,
                                       len(list(true_tree.tree.traverse())),  # Get all taxa in the tree
                                       true_tree.tree.robinson_foulds(inferences[method].forest[0].tree,
                                                                      attr_t1='sequence',
                                                                      attr_t2='sequence',
                                                                      unrooted_trees=True)[0],
                                       MRCA_distance(true_tree, inferences[method].forest[0]).sum())
                                       for method in inferences])
    df = pd.DataFrame({'method':methods, 'N_taxa':n_taxa, 'RF':distances, 'MRCA':MRCAs}, columns=('method', 'N_taxa', 'RF', 'MRCA'))
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
