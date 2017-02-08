#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
comparison of inference and simulated trees
'''

from __future__ import division, print_function
import gctree

def validate(args):
    with open(args.truetree, 'rb') as f:
        true_tree = pickle.load(f)
    with open(args.parfor, 'rb') as f:
        parsimony_forest = pickle.load(f)

    # NOTE: the unrooted_trees flag is needed because, for some reason, the RF
    #       function sometimes thinks the collapsed trees are unrooted and barfs
    distances, likelihoods = zip(*[(true_tree.tree.robinson_foulds(tree.tree, attr_t1='sequence', attr_t2='sequence', unrooted_trees=True)[0],
                                    tree.l(parsimony_forest.params)[0]) for tree in parsimony_forest.forest])

    df = pd.DataFrame({'RF':distances, 'log-likelihood':likelihoods})

    # here's Erick's idea of matrix of hamming distance of common ancestors of taxa
    taxa = [node.sequence for node in true_tree.tree.traverse() if node.frequency]
    sequence_length = len(taxa[0])
    n_taxa = len(taxa)
    MRCA_sum_metric = []
    for ct, tree in enumerate(parsimony_forest.forest, 1):
        d = scipy.zeros(shape=(n_taxa, n_taxa))
        for i in range(n_taxa):
            nodei_true = true_tree.tree.iter_search_nodes(sequence=taxa[i]).next()
            nodei      =      tree.tree.iter_search_nodes(sequence=taxa[i]).next()
            for j in range(i + 1, n_taxa):
                nodej_true = true_tree.tree.iter_search_nodes(sequence=taxa[j]).next()
                nodej      =      tree.tree.iter_search_nodes(sequence=taxa[j]).next()
                MRCA_true = true_tree.tree.get_common_ancestor((nodei_true, nodej_true)).sequence
                MRCA =           tree.tree.get_common_ancestor((nodei, nodej)).sequence
                d[i, j] = nodei.frequency*nodej.frequency*hamming_distance(MRCA_true, MRCA)#/sequence_length
        sns.heatmap(d, vmin=0)
        plt.savefig(args.outbase+'.validation.ancestor.{}.pdf'.format(ct))
        plt.clf()
        MRCA_sum_metric.append(d.sum())
    df['MRCA'] = MRCA_sum_metric

    # plots
    sns.pairplot(df, kind='reg', x_vars='log-likelihood', y_vars=('MRCA', 'RF'), aspect=1.5).savefig(args.outbase+'.validation.pdf')
    # sns.regplot(x='log-likelihood', y='distance', data=df).get_figure().savefig(args.outbase+'.validation.RF.pdf')
    # sns.regplot(x='log-likelihood', y='MRCA', data=df).get_figure().savefig(args.outbase+'.validation.MRCA.pdf')

    df.to_csv(args.outbase+'.validation.tsv', sep='\t', index=False)


def main():

    import argparse

    parser = argparse.ArgumentParser(description='validate results of inference on simulation data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('truetree', type=str, help='.p file containing true tree')
    parser.add_argument('parfor', type=str, help='.p file containing parsimony forest from inference')
    parser.add_argument('--outbase', type=str, default='gctree.out', help='output file base name')
    parser.set_defaults(func=validate)


if __name__ == '__main__':
    main()
