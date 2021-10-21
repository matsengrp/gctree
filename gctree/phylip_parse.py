#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Given an outputfile from one of the PHYLIP tools - `dnaml` or `dnapars` - produce an alignment (including
ancestral sequences), a newick tree (with matching internal node lables), and an svg rendering of said tree.
"""

import gctree.phylip_parse as pp
import gctree.branching_processes as bp
from gctree.utils import hamming_distance

from ete3 import Tree
import re
import random
from collections import defaultdict
from Bio.Data.IUPACData import ambiguous_dna_values

import pickle
import argparse

bases = "AGCT-"
ambiguous_dna_values.update({"?": "GATC-", "-": "-"})
code_vectors = {
    code: [0 if base in ambiguous_dna_values[code] else float("inf") for base in bases]
    for code in ambiguous_dna_values
}
cost_adjust = {
    base: [int(not i == j) for j in range(5)] for i, base in enumerate(bases)
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
def parse_outfile(outfile, abundance_file=None, root="root", **kwargs):
    """parse phylip outfile."""
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
                    trees[-1].append(
                        build_tree(sequences, parents, counts, root, **kwargs)
                    )
                else:
                    trees.append(build_tree(sequences, parents, counts, root, **kwargs))
            elif sect == "seqboot_dataset":
                bootstrap = True
                trees.append([])
            else:
                raise RuntimeError("unrecognized phylip section = {}".format(sect))
    return trees


def disambiguate(
    tree: Tree, random_state=None, dist_func=hamming_distance, distance_dependence=0
) -> Tree:
    """Randomly resolve ambiguous bases using a two-pass Sankoff Algorithm on
    subtrees of consecutive ambiguity codes."""

    def sequences(sequence, accum=""):
        """Iterates through possible disambiguations of sequence, recursively.

        Recursion-depth-limited by number of ambiguity codes in
        sequence, not sequence length.
        """
        if sequence:
            for index, base in enumerate(sequence):
                if base in bases:
                    accum += base
                else:
                    for newbase in ambiguous_dna_values[base]:
                        yield from sequences(
                            sequence[index + 1 :], accum=(accum + newbase)
                        )
                    return
        yield accum

    def iterate_independent_indices(sequencelist, dependence_window):
        """Sequencelist must contain sequences of the same length.

        dependence_window is max distance a motif extends past focused
        base, so a centered motif of length 5 has dependence_window 2.
        """
        # contains True at indices for columns which have an ambiguity, False otherwise.
        is_ambiguous_collapsed = [
            any((sequence[index] not in bases for sequence in sequencelist))
            for index in range(len(sequencelist[0]))
        ]
        ambiguous_indices = [
            index for index, val in enumerate(is_ambiguous_collapsed) if val is True
        ]

        def is_window_beginning(ambig_index):
            before_index = max([0, ambig_index - dependence_window])
            return not any(is_ambiguous_collapsed[before_index:ambig_index])

        def is_window_end(ambig_index):
            return not any(
                is_ambiguous_collapsed[
                    ambig_index + 1 : ambig_index + dependence_window + 1
                ]
            )

        window_beginnings = (
            index for index in ambiguous_indices if is_window_beginning(index)
        )
        window_ends = (index for index in ambiguous_indices if is_window_end(index))
        return (
            (max([0, a - dependence_window]), b + dependence_window + 1)
            for a, b in zip(window_beginnings, window_ends)
        )

    def is_ambiguous(sequence):
        return any((code not in bases for code in sequence))

    if random_state is None:
        random.seed(tree.write(format=1))
    else:
        random.setstate(random_state)

    seqlist = [node.sequence for node in tree.traverse()]
    windows_to_disambiguate = list(
        iterate_independent_indices(seqlist, distance_dependence)
    )
    for node in tree.traverse():
        for startidx, endidx in windows_to_disambiguate:
            if any((base not in bases for base in node.sequence[startidx:endidx])):

                def is_leaf(node):
                    return (node.is_leaf()) or not is_ambiguous(
                        node.sequence[startidx:endidx]
                    )

                # First pass of Sankoff: compute cost vectors
                for node2 in node.traverse(strategy="postorder", is_leaf_fn=is_leaf):
                    seq_fragment = node2.sequence[startidx:endidx]
                    node2.add_feature(
                        "costs", [[seq, 0] for seq in sequences(seq_fragment)]
                    )
                    if not is_leaf(node2):
                        for seq_cost in node2.costs:
                            for child in node2.children:
                                seq_cost[1] += min(
                                    [
                                        dist_func(seq_cost[0], child_seq) + child_cost
                                        for child_seq, child_cost in child.costs
                                    ]
                                )
                # Second pass: Choose base and adjust children's cost vectors
                # first adjust costs by above-node sequence
                if not node.is_root():
                    upseq_frag = node.up.sequence[startidx:endidx]
                    node.costs = [
                        [sequence, cost + dist_func(upseq_frag, sequence)]
                        for sequence, cost in node.costs
                    ]
                # traverse evaluates is_leaf(node) after yielding node.
                # Resolving base makes is_leaf true; must get order before
                # making changes.
                preorder = list(node.traverse(strategy="preorder", is_leaf_fn=is_leaf))
                for node2 in preorder:
                    if not is_ambiguous(node2.sequence[startidx:endidx]):
                        continue
                    min_cost = min([cost for _, cost in node2.costs])
                    resolved_fragment = random.choice(
                        [sequence for sequence, cost in node2.costs if cost == min_cost]
                    )
                    # Adjust child cost vectors
                    if not is_leaf(node2):
                        for child in node2.children:
                            child.costs = [
                                [
                                    sequence,
                                    cost + dist_func(resolved_fragment, sequence),
                                ]
                                for sequence, cost in child.costs
                            ]
                    node2.sequence = (
                        node2.sequence[:startidx]
                        + resolved_fragment
                        + node2.sequence[endidx:]
                    )
    for node in tree.traverse():
        try:
            node.del_feature("costs")
        except (AttributeError, KeyError):
            pass
    return tree


# build a tree from a set of sequences and an adjacency dict.
def build_tree(
    sequences, parents, counts=None, root="root", dist_func=hamming_distance, **kwargs
):
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
            root_parent.children[0].dist = dist_func(
                nodes[root_id].sequence, root_parent.children[0].sequence
            )
        tree = nodes[root_id]

    # make random choices for ambiguous bases
    tree = disambiguate(tree, dist_func=dist_func, **kwargs)

    # compute branch lengths
    tree.dist = 0  # no branch above root
    for node in tree.iter_descendants():
        node.dist = dist_func(node.up.sequence, node.sequence)

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
    trees = pp.parse_outfile(args.phylip_outfile, args.abundance_file, args.root)
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
