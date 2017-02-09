#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
Parse the output ancestral state reconstruction from IgPhyML together with the topology
to create an ete3 tree with the ancestral sequences.
'''

from __future__ import division, print_function
import sys
import os
import argparse
try:
    import cPickle as pickle
except:
    import pickle
from ete3 import Tree, TreeNode, NodeStyle, TreeStyle, TextFace, add_face_to_node, CircleFace, faces, AttrFace
from Bio import SeqIO
sys.path.append(os.path.abspath("/fh/fast/matsen_e/kdavidse/gctree/bin"))
from gctree import CollapsedForest, CollapsedTree
import warnings


def find_node(tree, pattern):
    regex = re.compile(pattern).search
    nodes =  [ node for node in tree.traverse() for m in [regex(node.name)] if m]
    if not nodes:
        warn("Cannot find matching node; looking for name matching '{}'".format(pattern))
        return
    else:
        if len(nodes) > 1:
            warn("multiple nodes found; using first one.\nfound: {}".format([n.name for n in nodes]))
        return nodes[0]


# reroot the tree on node matching regex pattern.
# Usually this is used to root on the naive germline sequence with a name matching '.*naive.*'
def reroot_tree(tree, pattern='.*naive.*', outgroup=0):
    # find all nodes matching pattern
    node = find_node(tree, pattern)
    tree.set_outgroup(node)
    if tree != node and not outgroup:
        tree.remove_child(node)
        node.add_child(tree)
        tree.dist = node.dist
        node.dist = 0
        tree = node
    return tree


class TreeFileParsingError(Exception):
    '''When ete3 fails to read the input tree.'''


def map_asr_to_tree(asr_seq, tree):
    '''Takes a IgPhyML fasta header and returns the matching ete3 tree node.'''
    for record in SeqIO.parse(asr_seq, "fasta"):
        descendants = record.id.split(';')[1].split(',')
        # Fix IgPhyML ASR turning fasta headers to upper case:
        descendants = [d.lower() for d in descendants]
        assert(descendants)
        if len(descendants) > 1:
            ancestor = tree.get_common_ancestor(descendants)
            frequency = 0
        elif descendants[0] == 'naive':
            ancestor = tree.get_leaves_by_name(descendants[0])
            ancestor = ancestor[0]
            frequency = 0
        else:
            ancestor = tree.get_leaves_by_name(descendants[0])
            ancestor = ancestor[0]
            frequency = int(descendants[0].split('_')[-1])

        assert(ancestor)
        ancestor.add_feature('frequency', frequency)
        ancestor.add_feature('sequence', record.seq)

    return tree


def main():
    parser = argparse.ArgumentParser(description='Reroot a tree based on node containing a keyword.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--tree',
        required=True,
        metavar='NEWICK TREE',
        help='Input tree used for topology.')
    parser.add_argument(
        '--asr_seq',
        required=True,
        help='Input ancestral sequences.')
    parser.add_argument(
        '--outbase',
        required=True,
        metavar='FILENAME',
        help='Filename for the output ASR tree.')
    args = parser.parse_args()

    try:
        tree = Tree(args.tree)
    except:
        raise TreeFileParsingError('Could not read the input tree. Is this really newick format?')

    tree.add_feature('frequency', 0)
    tree = map_asr_to_tree(args.asr_seq, tree)

    igphyml_tree = CollapsedTree(tree=tree)
    igphyml_tree.render(args.outbase + '.svg')
    igphyml_forest = CollapsedForest(forest=[igphyml_tree])
    print('number of trees with integer branch lengths:', igphyml_forest.n_trees)

    # check for unifurcations at root
    unifurcations = sum(tree.tree.frequency == 0 and len(tree.tree.children) == 1 for tree in igphyml_forest.forest)
    if unifurcations:
        print('WARNING: {} trees exhibit unifurcation from root, which is not possible under current model. Such nodes will be ommitted from likelihood calculation'.format(unifurcations))

    with open(args.outbase + '.p', 'wb') as f:
        pickle.dump(igphyml_forest, f)


if __name__ == "__main__":
    main()



# Reroot to make the naive sequence the real root instead of just an outgroup???
    # rerotted_tree = reroot_tree(tree, pattern=args.pattern, outgroup=args.outgroup)
    # rerotted_tree.write(outfile=args.reroot_tree)
