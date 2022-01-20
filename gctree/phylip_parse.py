#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given an outputfile from one of the PHYLIP tools - `dnaml` or `dnapars` - produce an alignment (including
ancestral sequences), a newick tree (with matching internal node lables), and an svg rendering of said tree.
"""

import gctree.branching_processes as bp
import gctree.utils

import historydag as hdag
from ete3 import Tree
import re
import random
from collections import defaultdict

import pickle
import argparse
import warnings
from typing import Dict

code_vectors = {
    code: [
        0 if base in gctree.utils.ambiguous_dna_values[code] else float("inf")
        for base in gctree.utils.bases
    ]
    for code in gctree.utils.ambiguous_dna_values
}
cost_adjust = {
    base: [int(not i == j) for j in range(5)]
    for i, base in enumerate(gctree.utils.bases)
}


# iterate over recognized sections in the phylip output file.
def sections(fh):
    patterns = {
        "parents": r"\s+between\s+and\s+length",
        ("sequences", "dnaml"): r"\s*node\s+reconstructed\s+sequence",
        ("sequences", "dnapars"): r"from\s+to\s+any steps",
        "seqboot_dataset": r"Data\s*set",
    }
    patterns = {k: re.compile(v, re.IGNORECASE) for (k, v) in patterns.items()}
    for line in fh:
        for k, pat in patterns.items():
            if pat.match(line):
                yield k
                break


# iterate over entries in the distance section
def iter_edges(fh):
    #  152          >root2           0.01208     (     zero,     0.02525) **
    pat = re.compile(
        r"\s*(?P<parent>\w+)\s+(?P<child>[\w>_.-]+)\s+(?P<distance>\d+\.\d+)"
    )
    # drop the header underline
    fh.readline()
    matches = 0
    for line in fh:
        m = pat.match(line)
        if m:
            matches += 1
            yield (m.group("child"), m.group("parent"))
        # We only want to break at the very end of the block of matches; dnaml has an extra blank line between
        # header and rows that dnapars doesn't
        elif not matches:
            continue
        else:
            break


# iterate over entries in the sequences section
def parse_seqdict(fh, mode="dnaml"):
    #  152        sssssssssG AGGTGCAGCT GTTGGAGTCT GGGGGAGGCT TGGTACAGCC TGGGGGGTCC
    seqs = defaultdict(str)
    if mode == "dnaml":
        patterns = re.compile(r"^\s*(?P<id>[a-zA-Z0-9>_.-]*)\s+(?P<seq>[a-zA-Z \-]+)")
    elif mode == "dnapars":
        patterns = re.compile(
            r"^\s*\S+\s+(?P<id>[a-zA-Z0-9>_.-]*)\s+(yes\s+|no\s+|maybe\s+)?(?P<seq>[a-zA-Z \-\?]+)"
        )
    else:
        raise ValueError("invalid mode " + mode)
    fh.readline()
    for line in fh:
        m = patterns.match(line)
        if m and m.group("id") != "":
            last_blank = False
            seqs[m.group("id")] += m.group("seq").replace(" ", "").upper()
        elif line.rstrip() == "":
            if last_blank:
                break
            else:
                last_blank = True
                continue
        else:
            break
    return seqs


# parse the dnaml output file and return data structures containing a
# list biopython.SeqRecords and a dict containing adjacency
# relationships and distances between nodes.
def parse_outfile(outfile, abundance_file=None, root="root", disambiguate=False):
    """parse phylip outfile, and return dnapars trees, a dictionary mapping
    node names to sequences, and a dictionary mapping node names to observed
    abundances."""
    if abundance_file is not None:
        counts = {
            line.split(",")[0]: int(line.split(",")[1]) for line in open(abundance_file)
        }
    # No count, just make an empty count dictionary:
    else:
        counts = None
    trees = []
    bootstrap = False
    # Ugg... for compilation need to let python know that these will definely both be defined :-/
    sequences, parents = {}, {}
    with open(outfile, "rU") as fh:
        for sect in sections(fh):
            if sect == "parents":
                parents = {child: parent for child, parent in iter_edges(fh)}
            elif sect[0] == "sequences":
                sequences = parse_seqdict(fh, sect[1])
                # sanity check;  a valid tree should have exactly one node that is parentless
                if not len(parents) == len(sequences) - 1:
                    raise RuntimeError(
                        "invalid results attempting to parse {}: there are {} parentless sequences".format(
                            outfile, len(sequences) - len(parents)
                        )
                    )
                if bootstrap:
                    trees[-1].append(build_tree(sequences, parents, counts, root))
                else:
                    trees.append(build_tree(sequences, parents, counts, root))
            elif sect == "seqboot_dataset":
                bootstrap = True
                trees.append([])
            else:
                raise RuntimeError("unrecognized phylip section = {}".format(sect))
    if disambiguate:
        # Disambiguate sets node.dist for all nodes in disambiguated trees
        trees = [disambiguate(tree) for tree in trees]
    sequence_counts = {sequences[name]: count for name, count in counts.items()}
    return (trees, sequence_counts)


def make_dag(trees, sequence_counts={}, from_copy=True):
    """Build a history DAG from ambiguous or disambiguated trees, whose nodes
    have abundance, name, and sequence attributes."""
    # preprocess trees so they're acceptable inputs
    # Assume all trees have fixed root sequence and fixed leaf sequences
    leaf_seqs = {node.sequence for node in trees[0].get_leaves()}

    sequence_counts = sequence_counts.copy()
    if from_copy:
        trees = [tree.copy() for tree in trees]
    if all(len(tree.children) > 1 for tree in trees):
        pass  # We're all good!
    elif trees[0].sequence not in leaf_seqs:
        if trees[0].sequence not in sequence_counts:
            sequence_counts[trees[0].sequence] = 0
        for tree in trees:
            newleaf = tree.add_child(name="", dist=0)
            newleaf.add_feature("sequence", trees[0].sequence)
            if tree.sequence != newleaf.sequence:
                raise ValueError(
                    "At least some trees unifurcate at root, but root sequence is not fixed."
                )
    else:
        # This should never happen in parsimony setting, when internal edges
        # are collapsed by sequence
        raise RuntimeError(
            "Root sequence observed, but the corresponding leaf is not a child of the root node. "
            "Gctree inference may give nonsensical results. Are you sure these are parsimony trees?"
        )

    dag = hdag.history_dag_from_etes(
        trees,
        ["sequence"],
        attr_func=lambda n: {
            "name": n.name,
        },
    )
    # If there are too many ambiguities at too many nodes, disambiguation will
    # hang. Need to have an alternative (disambiguate each tree before putting in dag):
    if (
        dag.count_trees(expand_count_func=hdag.utils.sequence_resolutions_count)
        / dag.count_trees()
        > 500000000
    ):
        warnings.warn(
            "Parsimony trees have too many ambiguities for disambiguation in history DAG. "
            "Disambiguating trees individually. History DAG may find fewer new parsimony trees."
        )
        trees = [disambiguate(tree) for tree in trees]
        dag = hdag.history_dag_from_etes(
            trees,
            ["sequence"],
            attr_func=lambda n: {
                "name": n.name,
            },
        )
    dag.explode_nodes(expand_func=hdag.utils.sequence_resolutions)
    # Look for (even) more trees:
    dag.add_all_allowed_edges(adjacent_labels=True)
    dag.trim_optimal_weight()
    dag.convert_to_collapsed()
    # Add abundances to attrs:
    for node in dag.preorder(skip_root=True):
        if node.label.sequence in sequence_counts:
            node.attr["abundance"] = sequence_counts[node.label.sequence]
        else:
            node.attr["abundance"] = 0

    if len(dag.hamming_parsimony_count()) > 1:
        raise RuntimeError(
            f"History DAG parsimony search resulted in parsimony trees of unexpected weights:\n {dag.hamming_parsimony_count()}"
        )
    for node in dag.preorder(skip_root=True):
        if (
            node.attr["abundance"] != 0
            and not node.is_leaf()
            and frozenset({node.label}) not in node.clades
        ):
            raise RuntimeError(
                "An internal node not adjacent to a leaf with the same label was found with nonzero abundance."
            )

    # give disambiguated sequences unique names
    sequences = {
        node.attr["name"]: node.label.sequence for node in dag.preorder(skip_root=True)
    }
    n_max = max([int(name) for name in sequences.keys() if name.isdigit()])
    namedict = {sequence: name for name, sequence in sequences.items()}
    for node in dag.preorder(skip_root=True):
        if node.label.sequence not in namedict:
            n_max += 1
            namedict[node.label.sequence] = str(n_max)
            node.attr["name"] = str(n_max)
    return dag


def disambiguate(tree: Tree, random_state=None) -> Tree:
    """Randomly resolve ambiguous bases using a two-pass Sankoff Algorithm on
    subtrees of consecutive ambiguity codes."""
    if random_state is None:
        random.seed(tree.write(format=1))
    else:
        random.setstate(random_state)
    for node in tree.traverse():
        for site, base in enumerate(node.sequence):
            if base not in gctree.utils.bases:

                def is_leaf(node):
                    return (node.is_leaf()) or (
                        node.sequence[site] in gctree.utils.bases
                    )

                # First pass of Sankoff: compute cost vectors
                for node2 in node.traverse(strategy="postorder", is_leaf_fn=is_leaf):
                    base2 = node2.sequence[site]
                    node2.add_feature("cv", code_vectors[base2].copy())
                    if not is_leaf(node2):
                        for i in range(5):
                            for child in node2.children:
                                node2.cv[i] += min(
                                    [
                                        sum(v)
                                        for v in zip(
                                            child.cv, cost_adjust[gctree.utils.bases[i]]
                                        )
                                    ]
                                )
                # Second pass: Choose base and adjust children's cost vectors
                if not node.is_root():
                    node.cv = [
                        sum(v)
                        for v in zip(node.cv, cost_adjust[node.up.sequence[site]])
                    ]
                # traverse evaluates is_leaf(node) after yielding node.
                # Resolving base makes is_leaf true; must get order before
                # making changes.
                preorder = list(node.traverse(strategy="preorder", is_leaf_fn=is_leaf))
                for node2 in preorder:
                    if node2.sequence[site] in gctree.utils.bases:
                        continue
                    min_cost = min(node2.cv)
                    base_index = random.choice(
                        [i for i, val in enumerate(node2.cv) if val == min_cost]
                    )
                    new_base = gctree.utils.bases[base_index]
                    # Adjust child cost vectors
                    if not is_leaf(node2):
                        for child in node2.children:
                            child.cv = [
                                sum(v) for v in zip(child.cv, cost_adjust[new_base])
                            ]
                    node2.sequence = (
                        node2.sequence[:site] + new_base + node2.sequence[(site + 1) :]
                    )
    for node in tree.traverse():
        try:
            node.del_feature("cv")
        except (AttributeError, KeyError):
            pass
    tree.dist = 0
    for node in tree.iter_descendants():
        node.dist = gctree.utils.hamming_distance(node.up.sequence, node.sequence)
    return tree


# build a tree from a set of sequences and an adjacency dict.
def build_tree(
    sequences: Dict[str, str], parents: Dict[str, str], counts=None, root="root"
):
    """Build an ete tree from sequences and parents dictionaries.

    Args:
        sequences: a dictionary mapping names to sequences
        parents: a dictionary mapping parent node names to child node names
        counts: a dictionary mapping node names to observed abundances. This argument
            is no longer used in the main gctree inference pipeline (counts are assigned in DAG)
            but remains for compatibility.
        root: the name of the root node
    """
    # build an ete tree
    # first a dictionary of disconnected nodes
    nodes = {}
    for name in sequences:
        node = Tree()
        node.name = name
        node.add_feature("sequence", sequences[node.name])
        if counts is not None:
            if node.name in counts:
                node.add_feature("abundance", counts[node.name])
            else:
                node.add_feature("abundance", 0)
        nodes[name] = node
    for name in sequences:
        if name in parents:
            nodes[parents[name]].add_child(nodes[name])
        else:
            tree = nodes[name]
    # reroot on root
    if root is not None:
        root_id = [node for node in nodes if root in node][0]
        assert len(nodes[root_id].children) == 0
        root_parent = nodes[root_id].up
        root_parent.remove_child(nodes[root_id])
        nodes[root_id].add_child(root_parent)
        # remove possible unecessary unifurcation after rerooting
        if len(root_parent.children) == 1:
            root_parent.delete(prevent_nondicotomic=False)
            root_parent.children[0].dist = gctree.utils.hamming_distance(
                root_parent.children[0].sequence, nodes[root_id].sequence
            )
        tree = nodes[root_id]

    return tree


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "phylip_outfile",
        type=str,
        help="dnaml outfile (verbose output with inferred ancestral sequences, option 5).",
    )
    parser.add_argument("abundance_file", type=str, help="count file")
    parser.add_argument("--outputfile", default=None, help="output file.")
    parser.add_argument("--root", default="root", help="root sequence id")
    return parser


def main(arg_list=None):
    args = get_parser().parse_args(arg_list)

    if args.outputfile is None:
        args.outputfile = args.phylip_outfile + ".collapsed_forest.p"
    trees = parse_outfile(args.phylip_outfile, args.abundance_file, args.root)
    if isinstance(trees[0], list):
        print(trees[0][0])
        print(bp.CollapsedTree(tree=trees[0][0]))
        bootstraps = [
            [bp.CollapsedTree(tree=tree) for tree in bootstrap] for bootstrap in trees
        ]
        pickle.dump(
            [bp.CollapsedForest(forest=trees) for trees in bootstraps],
            open(args.outputfile, "w"),
        )
    else:
        trees = [bp.CollapsedTree(tree=tree) for tree in trees]
        pickle.dump(bp.CollapsedForest(forest=trees), open(args.outputfile, "w"))


if __name__ == "__main__":
    main()
