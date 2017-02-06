from __future__ import print_function
import argparse
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, add_face_to_node
import os
import re


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


def main():
    parser = argparse.ArgumentParser(description='Reroot a tree based on node containing a keyword.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--tree',
        required=True,
        metavar='NEWICK TREE',
        help='Input tree to root, in newick format.')
    parser.add_argument(
        '--reroot_tree',
        required=True,
        metavar='FILENAME',
        help='Filename for the output rerooted tree.')
    parser.add_argument(
        '--pattern',
        metavar='REGEX PATTERN',
        required=True,
        help="Pattern to search for the node to root on.")
    parser.add_argument(
        '--outgroup',
        required=False,
        default=0,
        metavar="0/1",
        help="Set as outgroup instead of tree root.")
    args = parser.parse_args()

    try:
        tree = Tree(args.tree)
    except:
        raise TreeFileParsingError('Could not read the input tree. Is this really newick format?')

    rerotted_tree = reroot_tree(tree, pattern=args.pattern, outgroup=args.outgroup)
    rerotted_tree.write(outfile=args.reroot_tree)


if __name__ == "__main__":
    main()

