#! /usr/bin/env python
# -*- coding: utf-8 -*-

import re
import random
import ete3
from gctree.utils import hamming_distance
import pickle
import argparse
from pathlib import Path
import warnings


class IsotypeTemplate:
    def __init__(self, isotype_order):
        self.order = isotype_order

    def new(self, isotype_name):
        return Isotype(self.order, isotype_name)


class Isotype:
    # From https://doi.org/10.7554/eLife.16578.012, for humans:
    # order = ["M", "G3", "A1", "G2", "G4", "E", "A2"]

    def __init__(self, order, isotype_name):
        self.order = order
        if isotype_name == "?":
            self.isotype = None
        else:
            self.isotype = self.order.index(isotype_name)

    def isbefore(self, t):
        return self.isotype <= t.isotype

    def __repr__(self):
        return f"Isotype('{self.order[self.isotype]}')"

    def __str__(self):
        return f"Ig{self.order[self.isotype]}"

    def copy(self):
        return Isotype(self.order, self.order[self.isotype])

    def __eq__(self, other):
        return self.isotype == other.isotype

    def __hash__(self):
        return hash(self.__repr__())

    def resolutions(self):
        """Returns list of all possible isotypes if passed an ambiguous isotype.
        Otherwise, returns a list containing a copy of the passed isotype.
        """
        if self.isotype is None:
            return [Isotype(self.order, name) for name in self.order]
        else:
            return [self.copy()]


def isotype_distance(t1, t2):
    """This function is not symmetric on its arguments"""
    if not t1.isbefore(t2):
        return float("inf")
    else:
        return float(t1.isotype != t2.isotype)


def isotype_parsimony(tree: ete3.TreeNode):
    return sum(
        isotype_distance(node.up.isotype, node.isotype)
        for node in tree.iter_descendants()
    )


def disambiguate_isotype(
    tree: ete3.Tree, random_state=None, dist_func=isotype_distance
):
    def is_ambiguous(isotype):
        return isotype.isotype is None

    if random_state is None:
        random.seed(tree.write(format=1))
    else:
        random.setstate(random_state)

    # First pass of Sankoff: compute cost vectors
    for node2 in tree.traverse(strategy="postorder"):
        node2.add_feature("costs", [[seq, 0] for seq in node2.isotype.resolutions()])
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


def make_newidmap(idmap: dict, isotype_map: dict):
    newidmap = {}
    for id, cell_ids in idmap.items():
        isotypeset = set()
        for cell_id in cell_ids:
            try:
                isotypeset.add(isotype_map[cell_id])
            except KeyError as e:
                warnings.warn(
                    f"Sequence ID {id} has original sequence id {cell_id} "
                    "for which no observed isotype was provided. "
                    "Isotype will be assumed ambiguous if observed."
                )
                isotype_map[cell_id] = "?"
                isotypeset.add("?")
        newidmap[id] = {
            isotype: {
                cell_id for cell_id in cell_ids if isotype_map[cell_id] == isotype
            }
            for isotype in isotypeset
        }
    return newidmap


def add_observed_isotypes(tree: ete3.Tree, newidmap: dict, isotype_order: list):
    # Drop observed nodes as leaves and explode by observed isotype:
    # Descend internal observed nodes as leaves:
    newisotype = IsotypeTemplate(isotype_order).new
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
            node.add_feature("isotype", newisotype("?"))
        else:
            try:
                thisnode_isotypemap = newidmap[node.name]
            except KeyError as e:
                warnings.warn(
                    f"The sequence name {e} labels an observed node, but no mapping to an original sequence ID was found."
                    " Isotype will be assumed ambiguous."
                )
                thisnode_isotypemap = {
                    "?": {f"Unknown_id_{n+1}" for n in range(node.abundance)}
                }
            if "?" in thisnode_isotypemap:
                warnings.warn(
                    f"The sequence name {node.name} labels an observed node, and corresponds to sequence IDs for "
                    "which no observed isotype was provided. "
                    f" Isotype will be assumed ambiguous for: {', '.join(thisnode_isotypemap['?'])}"
                )
            # node.name had better be in newidmap, since this is an observed node
            if len(thisnode_isotypemap) > 1:
                for isotype, cell_ids in thisnode_isotypemap.items():
                    # add new node below this leaf node. Must be below, and not child
                    # of parent, to preserve max parsimony in case that node.up has
                    # different sequence from node.
                    newchild = ete3.TreeNode(name=node.name)
                    newchild.add_feature("abundance", len(cell_ids))
                    newchild.add_feature("sequence", node.sequence)
                    newchild.add_feature("isotype", newisotype(isotype))
                    node.add_child(child=newchild)
                node.abundance = 0
            else:
                node.isotype = newisotype(list(thisnode_isotypemap.keys())[0])
    # Now add root and internal ambiguous isotypes
    tree.add_feature("isotype", newisotype("M"))
    for node in tree.iter_descendants():
        if not node.is_leaf():
            node.add_feature("isotype", newisotype("?"))


def collapse_tree_by_sequence_and_isotype(tree):
    """May be able to combine this function in utils with the collapsing infrastructure in
    CollapsedTree."""
    for node in tree.iter_descendants():
        node.dist = node.up.sequence != node.sequence or node.up.isotype != node.isotype
    for node in tree.iter_descendants():
        if node.dist == 0:
            node.up.abundance += node.abundance
            node.up.name = node.name
            node.delete(prevent_nondicotomic=False)


def get_parser():
    parser = argparse.ArgumentParser(
        description=(
            "To be run on a directory containing gctree inference results."
            "Observed isotypes will be added to trees previously output by gctree,"
            "with unobserved isotypes inferred so as to minimize isotype switching"
            "while adhering to switching order."
        )
    )
    parser.add_argument(
        "--parsimony_forest",
        type=str,
        required=True,
        help="filename for parsimony_forest pickle file output by gctree inference",
    )
    parser.add_argument(
        "--inference_log",
        type=str,
        required=True,
        help="filename for gctree inference log file which contains branching process parameters",
    )
    parser.add_argument(
        "--isotype_mapfile",
        type=str,
        required=True,
        help="filename for a csv file mapping original sequence ids to observed isotypes"
        ". For example, each line should have the format 'somesequence_id, some_isotype'.",
    )
    parser.add_argument(
        "--idmapfile",
        type=str,
        required=True,
        help="filename for a csv file mapping sequence names to original sequence ids, like the one output by deduplicate.",
    )
    parser.add_argument(
        "--isotype_names",
        type=str,
        default=None,
        help="A list of isotype names used in isotype_mapfile, in order of most naive to most differentiated."
        """ Default is equivalent to passing 'M,G3,A1,G2,G4,E,A2'. """,
    )
    parser.add_argument(
        "--out_directory",
        type=str,
        default=None,
        help="Directory in which to place output. Default is working directory.",
    )
    return parser


def main(arg_list=None):
    isotype_palette = [
        "#a6cee3",
        "#1f78b4",
        "#b2df8a",
        "#33a02c",
        "#fb9a99",
        "#e31a1c",
        "#fdbf6f",
        "#ff7f00",
        "#cab2d6",
        "#6a3d9a",
        "#ffff99",
        "#b15928",
    ]
    args = get_parser().parse_args(arg_list)
    if args.out_directory:
        out_directory = args.out_directory + "/"
        p = Path(out_directory)
        if not p.exists():
            p.mkdir()
    else:
        out_directory = ""
    with open(args.isotype_mapfile, "r") as fh:
        isotypemap = dict(map(lambda x: x.strip(), line.split(",")) for line in fh)

    with open(args.idmapfile, "r") as fh:
        idmap = {}
        for line in fh:
            seqid, cell_ids = line.rstrip().split(",")
            cell_idset = {cell_id for cell_id in cell_ids.split(":") if cell_id}
            if len(cell_idset) > 0:
                idmap[seqid] = cell_idset

    with open(args.inference_log, "r") as fh:
        for line in fh:
            if re.match(r"params:", line):
                p = float(re.search(r"(?<=[\(])\S+(?=\,)", line).group())
                q = float(re.search(r"(?<=\,\s)\S+(?=[\)])", line).group())
                parameters = (p, q)
                break

    with open(args.parsimony_forest, "rb") as fh:
        forest = pickle.load(fh)
    # parse the idmap file and the isotypemap file
    forest.forest = tuple(
        sorted(
            forest.forest, key=lambda tree: -tree.ll(*parameters, build_cache=False)[0]
        )
    )
    tree_stats = [
        [
            ctree,
            ctree.ll(*parameters, build_cache=False)[0],
            idx,
            sum(1 for _ in ctree.tree.traverse()),
        ]
        for idx, ctree in enumerate(forest.forest)
    ]
    if not args.isotype_names:
        isotype_names = ["M", "G3", "A1", "G2", "G4", "E", "A2"]
    else:
        isotype_names = str(args.isotype_names).split(",")

    newidmap = make_newidmap(idmap, isotypemap)
    for ctree in forest.forest:
        add_observed_isotypes(ctree.tree, newidmap, isotype_names)
        disambiguate_isotype(ctree.tree)
        collapse_tree_by_sequence_and_isotype(ctree.tree)
        for node in ctree.tree.traverse():
            node.name = node.name + " " + str(node.isotype)
        for node in ctree.tree.iter_descendants():
            node.dist = hamming_distance(node.up.sequence, node.sequence)

    flattened_newidmap = {
        name + " " + str(isotype): cell_idset
        for name, cellid_map in newidmap.items()
        for isotype, cell_idset in cellid_map.items()
    }

    with open(out_directory + "isotyped.idmap", "w") as fh:
        for name, cellidset in flattened_newidmap.items():
            print(f"{name},{':'.join(cellidset)}", file=fh)

    for sublist in tree_stats:
        sublist.append(isotype_parsimony(sublist[0].tree))

    # Compute parsimony indices
    for index, sublist in enumerate(sorted(tree_stats, key=lambda slist: slist[2])):
        sublist.append(index)

    # Pickle forest (still sorted by likelihood)
    with open(out_directory + "parsimonyforest.isotyped.p", "wb") as fh:
        fh.write(pickle.dumps(forest))

    print(
        f"Parameters:\t{parameters}\n"
        "index\t ll\t\t\t original node count\t isotype parsimony\t new node count"
    )
    for (
        ctree,
        likelihood,
        likelihood_idx,
        original_numnodes,
        parsimony,
        parsimony_idx,
    ) in tree_stats:
        print(
            f"{likelihood_idx + 1}\t {likelihood}\t {original_numnodes}\t\t\t {parsimony}\t\t\t {sum(1 for _ in ctree.tree.traverse())}"
        )
        colormap = {
            node.name: isotype_palette[node.isotype.isotype]
            for node in ctree.tree.traverse()
        }
        ctree.render(
            outfile=out_directory
            + f"gctree.out.{likelihood_idx + 1}.isotype_parsimony.{int(parsimony)}.svg",
            colormap=colormap,
            idlabel=True,
        )
