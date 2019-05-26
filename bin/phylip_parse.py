#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Given an outputfile from one of the PHYLIP tools - `dnaml` or `dnapars` - produce an alignment (including
ancestral sequences), a newick tree (with matching internal node lables), and an svg rendering of said tree.
"""
from __future__ import print_function
from ete3 import Tree
import re, random
from collections import defaultdict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data.IUPACData import ambiguous_dna_values
import gctree


# iterate over recognized sections in the phylip output file.
def sections(fh):
    patterns = {
        "parents": "\s+between\s+and\s+length",
        ("sequences", 'dnaml'): "\s*node\s+reconstructed\s+sequence",
        ('sequences', 'dnapars'): "from\s+to\s+any steps",
        "seqboot_dataset": "Data\s*set"}
    patterns = {k: re.compile(v, re.IGNORECASE) for (k,v) in patterns.items()}
    for line in fh:
        for k, pat in patterns.items():
            if pat.match(line):
                yield k
                break


# iterate over entries in the distance section
def iter_edges(fh):
    #  152          >naive2           0.01208     (     zero,     0.02525) **
    pat = re.compile("\s*(?P<parent>\w+)\s+(?P<child>[\w>_.-]+)\s+(?P<distance>\d+\.\d+)")
    # drop the header underline
    fh.next()
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
def parse_seqdict(fh, mode='dnaml'):
    #  152        sssssssssG AGGTGCAGCT GTTGGAGTCT GGGGGAGGCT TGGTACAGCC TGGGGGGTCC
    seqs = defaultdict(str)
    if mode == 'dnaml':
        patterns = re.compile("^\s*(?P<id>[a-zA-Z0-9>_.-]*)\s+(?P<seq>[a-zA-Z \-]+)")
    elif mode == 'dnapars':
        patterns = re.compile("^\s*\S+\s+(?P<id>[a-zA-Z0-9>_.-]*)\s+(yes\s+|no\s+|maybe\s+)?(?P<seq>[a-zA-Z \-\?]+)")
    else:
        raise ValueError('invalid mode '+mode)
    fh.next()
    for line in fh:
        m = patterns.match(line)
        if m and m.group("id") is not '':
            last_blank = False
            seqs[m.group("id")] += m.group("seq").replace(" ", "").upper()
        elif line.rstrip() == '':
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
def parse_outfile(outfile, countfile=None, naive='naive'):
    '''parse phylip outfile'''
    if countfile is not None:
        counts = {l.split(',')[0]:int(l.split(',')[1]) for l in open(countfile)}
    # No count, just make an empty count dictionary:
    else:
        counts = None
    trees = []
    bootstrap = False
    # Ugg... for compilation need to let python know that these will definely both be defined :-/
    sequences, parents = {}, {}
    with open(outfile, 'rU') as fh:
        for sect in sections(fh):
            if sect == 'parents':
                parents = { child:parent for child, parent in iter_edges(fh) }
            elif sect[0] == 'sequences':
                sequences = parse_seqdict(fh, sect[1])
                # sanity check;  a valid tree should have exactly one node that is parentless
                if not len(parents) == len(sequences) - 1:
                    raise RuntimeError('invalid results attempting to parse {}: there are {} parentless sequences'.format(outfile, len(sequences) - len(parents)))
                if bootstrap:
                    trees[-1].append(build_tree(sequences, parents, counts, naive))
                else:
                    trees.append(build_tree(sequences, parents, counts, naive))
            elif sect == 'seqboot_dataset':
                bootstrap = True
                trees.append([])
            else:
                raise RuntimeError("unrecognized phylip section = {}".format(sect))
    return trees


def disambiguate(tree):
    '''make random choices for ambiguous bases, respecting tree inheritance'''
    ambiguous_dna_values['?'] = 'GATC-'
    sequence_length = len(tree.sequence)
    for node in tree.traverse():
        for site in range(sequence_length):
            base = node.sequence[site]
            if base not in 'ACGT-':
                new_base = random.choice(ambiguous_dna_values[base])
                for node2 in node.traverse(is_leaf_fn=lambda n: False if base in [n2.sequence[site] for n2 in n.children] else True):
                    if node2.sequence[site] == base:
                        node2.sequence = node2.sequence[:site] + new_base + node2.sequence[(site+1):]
    return tree


# build a tree from a set of sequences and an adjacency dict.
def build_tree(sequences, parents, counts=None, naive='naive'):
    # build an ete tree
    # first a dictionary of disconnected nodes
    nodes = {}
    for name in sequences:
        node = Tree()
        node.name = name
        node.add_feature('sequence', sequences[node.name])
        if counts is not None:
            if node.name in counts:
                node.add_feature('frequency', counts[node.name])
            else:
                node.add_feature('frequency', 0)
        nodes[name] = node
    for name in sequences:
        if name in parents:
            nodes[parents[name]].add_child(nodes[name])
        else:
            tree = nodes[name]
    # reroot on naive
    if naive is not None:
        naive_id = [node for node in nodes if naive in node][0]
        assert len(nodes[naive_id].children) == 0
        naive_parent = nodes[naive_id].up
        naive_parent.remove_child(nodes[naive_id])
        nodes[naive_id].add_child(naive_parent)
        # remove possible unecessary unifurcation after rerooting
        if len(naive_parent.children) == 1:
            naive_parent.delete(prevent_nondicotomic=False)
            naive_parent.children[0].dist = gctree.hamming_distance(naive_parent.children[0].sequence, nodes[naive_id].sequence)
        tree = nodes[naive_id]

    # make random choices for ambiguous bases
    tree = disambiguate(tree)

    # compute branch lengths
    tree.dist = 0 # no branch above root
    for node in tree.iter_descendants():
        node.dist = gctree.hamming_distance(node.sequence, node.up.sequence)

    return tree

def main():

    import pickle, argparse, os

    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        return fname
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'phylip_outfile', type=existing_file, help='dnaml outfile (verbose output with inferred ancestral sequences, option 5).')
    parser.add_argument(
        'countfile', type=existing_file, help="count file")
    parser.add_argument(
        '--outputfile', default=None, help="output file.")
    parser.add_argument(
        '--naive', default='naive', help="naive sequence id")

    args = parser.parse_args()

    if args.outputfile is None:
        args.outputfile = args.phylip_outfile + '.collapsed_forest.p'
    trees = parse_outfile(args.phylip_outfile, args.countfile, args.naive)
    if isinstance(trees[0], list):
        print(trees[0][0])
        print(gctree.CollapsedTree(tree=trees[0][0]))
        bootstraps = [[gctree.CollapsedTree(tree=tree) for tree in bootstrap] for bootstrap in trees]
        pickle.dump([gctree.CollapsedForest(forest=trees) for trees in bootstraps], open(args.outputfile, 'w'))
    else:
        trees = [gctree.CollapsedTree(tree=tree) for tree in trees]
        pickle.dump(gctree.CollapsedForest(forest=trees), open(args.outputfile, 'w'))

if __name__ == "__main__":
    main()
