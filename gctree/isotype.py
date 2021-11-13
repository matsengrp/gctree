#! /usr/bin/env python
# -*- coding: utf-8 -*-

import random
import ete3
import pickle
import argparse


class Isotype:
    # From https://doi.org/10.7554/eLife.16578.012, for humans:
    order = ["M", "G3", "A1", "G2", "G4", "E", "A2"]

    def __init__(self, isotype_name):
        if isotype_name == "?":
            self.isotype = None
        else:
            self.isotype = Isotype.order.index(isotype_name)

    def isbefore(self, t):
        return self.isotype <= t.isotype

    def __repr__(self):
        return f"Isotype({Isotype.order[self.isotype]})"

    def copy(self):
        return Isotype(Isotype.order[self.isotype])

    def __eq__(self, other):
        return self.isotype == other.isotype

    def __hash__(self):
        return hash(self.__repr__())


def resolutions(t: Isotype):
    """Returns list of all possible isotypes if passed an ambiguous isotype.
    Otherwise, returns a list containing a copy of the passed isotype.
    """
    if t.isotype is None:
        return [Isotype(name) for name in Isotype.order]
    else:
        return [t.copy()]


def isotype_distance(t1, t2):
    """This function is not symmetric on its arguments"""
    if not t1.isbefore(t2):
        return float("inf")
    else:
        return float(t1.isotype != t2.isotype)


def disambiguate_isotype(
    tree: ete3.Tree, random_state=None, dist_func=isotype_distance
) -> ete3.Tree:
    def is_ambiguous(isotype):
        return isotype.isotype is None

    if random_state is None:
        random.seed(tree.write(format=1))
    else:
        random.setstate(random_state)

    # First pass of Sankoff: compute cost vectors
    for node2 in tree.traverse(strategy="postorder"):
        node2.add_feature("costs", [[seq, 0] for seq in resolutions(node2.isotype)])
        if not node2.is_leaf():
            for seq_cost in node2.costs:
                for child in node2.children:
                    seq_cost[1] += min(
                        [
                            dist_func(seq_cost[0], child_seq) + child_cost
                            for child_seq, child_cost in child.costs
                        ]
                    )
    # Second pass: Choose base and adjust children's cost vectors.
    # Not necessary if we only want the minimum weight:
    for node2 in tree.traverse(strategy="preorder"):
        if not is_ambiguous(node2.isotype):
            continue
        min_cost = min([cost for _, cost in node2.costs])
        resolved_isotype = random.choice(
            [isotype for isotype, cost in node2.costs if cost == min_cost]
        )
        # Adjust child cost vectors
        if not node2.is_leaf():
            for child in node2.children:
                child.costs = [
                    [
                        isotype,
                        cost + dist_func(resolved_isotype, isotype),
                    ]
                    for isotype, cost in child.costs
                ]
        node2.isotype = resolved_isotype
    for node in tree.traverse():
        try:
            node.del_feature("costs")
        except (AttributeError, KeyError):
            pass
    return tree


def make_newidmap(idmap, isotype_map):
    newidmap = {}
    for id, cell_ids in idmap.items():
        newidmap[id] = {isotype_map[cell_id]: set() for cell_id in isotype_map}
        for cell_id in isotype_map:
            newidmap[id][isotype_map[cell_id]].add(cell_id)
    return newidmap


def add_observed_isotypes(tree: ete3.Tree, newidmap):
    # Drop observed nodes as leaves and explode by observed isotype:
    # Descend internal observed nodes as leaves:
    for node in list(tree.iter_descendants()):
        if node.abundance > 0 and not node.is_leaf():
            newchild = ete3.TreeNode(name=node.name)
            newchild.add_feature("sequence", node.sequence)
            newchild.add_feature("abundance", node.abundance)
            node.abundance = 0
            node.add_child(child=newchild)
    # Now duplicate nodes which represent multiple isotypes
    for node in list(tree.get_leaves()):
        if node.abundance == 0:
            node.add_feature("isotype", Isotype("?"))
        else:
            thisnode_isotypemap = newidmap[node.name]
            # node.name had better be in newidmap, since this is an observed node
            if len(thisnode_isotypemap) > 1:
                for isotype, cell_ids in thisnode_isotypemap.items():
                    # add new node below this leaf node. Must be below, and not child
                    # of parent, to preserve max parsimony in case that node.up has
                    # different sequence from node.
                    newlabel = node.name + "Ig" + isotype.order[isotype.isotype]
                    newchild = ete3.TreeNode(name=newlabel)
                    newchild.add_feature("abundance", len(cell_ids))
                    newchild.add_feature("sequence", node.sequence)
                    newchild.add_feature("isotype", isotype)
                    node.add_child(child=newchild)
                node.abundance = 0
            else:
                node.isotype = next(thisnode_isotypemap.keys())
    # Now add root and internal ambiguous isotypes
    tree.isotype = Isotype("M")
    for node in tree.iter_descendants():
        if not node.is_leaf():
            node.add_feature("isotype", Isotype("?"))


def collapse_tree_by_sequence_and_isotype(tree):
    """May be able to combine this function in utils with the collapsing infrastructure in
    CollapsedTree."""
    for node in tree.iter_descendants():
        node.dist = node.up.sequence != node.sequence or node.up.isotype != node.isotype
    for node in tree.iter_descendants():
        if node.dist == 0:
            node.up.abundance += node.abundance
            node.up.name = node.name
            node.delete()


def get_parser():
    # TODO Rework all this
    parser = argparse.ArgumentParser(
        description="Deduplicate sequences in a fasta file, write to stdout in phylip format, and creates a few other files (see arguments). Headers must be a unique ID of less than "
        "or equal to 10 ASCII characters."
        "An additional sequence representing the outgroup/root must be included (even if one or more observed sequences are identical to it)."
    )
    parser.add_argument(
        "infile",
        type=str,
        nargs="+",
        help="Fasta file with less than or equal to 10 characters unique header ID. "
        "Because dnapars will name internal nodes by integers a node name must include"
        "at least one non number character (unless using the id_abundances option).",
    )
    parser.add_argument(
        "--abundance_file",
        type=str,
        default=None,
        help="filename for the output file containing the counts.",
    )
    parser.add_argument(
        "--idmapfile",
        type=str,
        default=None,
        help="filename for the output file containing the map of new unique ids to original seq ids.",
    )
    parser.add_argument(
        "--id_abundances",
        action="store_true",
        help="flag to interpret integer ids in input as abundances",
    )
    parser.add_argument("--root", type=str, default="root", help="root sequence id")
    parser.add_argument(
        "--frame", type=int, default=None, help="codon frame", choices=(1, 2, 3)
    )
    parser.add_argument(
        "--colorfile",
        type=str,
        default=None,
        help="optional input csv filename for colors of each cell.",
    )
    parser.add_argument(
        "--colormap",
        type=str,
        default=None,
        help="optional output filename for colors map.",
    )
    return parser


def main(arg_list=None):
    with open(idmapfile, 'r') as fh:
        idmap = {}
        for line in fh:
            seqid, cell_ids = line.rstrip().split(",")
            cell_idset = {
                    cell_id for cell_id in cell_ids.split(":")
                    }
                idmap[seqid] = cell_idset
                }
    with open(isotypemap_file, 'r') as fh:
        isotypemap = {cell_id: isotype for line in fh for cell_id, isotype in line.split(",")}
    idmap = None
    isotypemap = None
    with open("parsimonyforest.p", "rb") as fh:
        forest = pickle.load(fh)
    # parse the idmap file and the isotypemap file
    for ctree in forest.forest:
        add_observed_isotypes(ctree.tree, idmap, isotypemap)
        collapse_tree_by_sequence_and_isotype(ctree.tree)
    # Do some isotype parsimony ranking and plotting.
