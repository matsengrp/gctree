#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
Parse the output ancestral state reconstruction from IgPhyML together with the topology
to create an ete3 tree with the ancestral sequences.
'''

from __future__ import division, print_function
import sys
import os
import re
import warnings
from Bio import SeqIO
from ete3 import Tree, TreeNode, NodeStyle, TreeStyle, TextFace, add_face_to_node, CircleFace, faces, AttrFace
sys.path.append(os.path.abspath('/'.join(os.path.realpath(__file__).split('/')[:-2]) + "/bin"))
from Bio import AlignIO

class FastaInputError(Exception):
    '''When the fasta file in not reflecting amino acid DNA coding for protein.'''


class TreeFileParsingError(Exception):
    '''When ete3 fails to read the input tree.'''


class TreeFileParsingError(Exception):
    '''When ete3 fails to read the input tree.'''


def ASR_parser(args):
    try:
        import cPickle as pickle
    except:
        import pickle
    from gctree import CollapsedForest, CollapsedTree, hamming_distance

    try:
        tree = Tree(args.tree)
    except:
        raise TreeFileParsingError('Could not read the input tree. Is this really newick format?')

    counts = {l.split(',')[0]:int(l.split(',')[1]) for l in open(args.counts)}
    tree.add_feature('frequency', 0)       # Placeholder will be deleted when rerooting
    tree.add_feature('sequence', 'DUMMY')  # Placeholder will be deleted when rerooting
    tree = map_asr_to_tree(args.asr_seq, tree, args.naive, counts)

    # Reroot to make the naive sequence the real root instead of just an outgroup:
    tree = reroot_tree(tree)

    # Recompute branch lengths as hamming distances:
    tree.dist = 0  # No branch above root
    for node in tree.iter_descendants():
        node.dist = hamming_distance(node.sequence, node.up.sequence)

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

    print('Done parsing IgPhyML tree')


def map_asr_to_tree(asr_seq, tree, naiveID, counts):
    '''Takes a IgPhyML fasta header and returns the matching ete3 tree node.'''
    for record in SeqIO.parse(asr_seq, "fasta"):
        descendants = record.id.split(';')[1].split(',')
        # Fix IgPhyML ASR turning fasta headers to upper case:
        descendants = [d.lower() for d in descendants]
        assert(descendants)
        if len(descendants) > 1:  # ASR infer node
            ancestor = tree.get_common_ancestor(descendants)
            frequency = 0
        else:  # Observed leaf or naive
            ancestor = tree.get_leaves_by_name(descendants[0])
            ancestor = ancestor[0]
            if descendants[0] in counts:
                frequency = counts[descendants[0]]
            else:
                frequency = 0

        # Add the features:
        assert(ancestor)
        ancestor.add_feature('frequency', frequency)
        ancestor.add_feature('sequence', str(record.seq))

    return tree




def make_igphyml_config(args):

    igphyml_path = which(args.igphyml_exe.rstrip('/'))
    assert(igphyml_path is not None)
    IGPHYML_DIR = re.sub(r'/src/\w+$', '', igphyml_path)
    MODEL = args.model
    assert(MODEL in ['gy94', 'hlp16'])
    # Find the length of the translated DNA sequence:
    len_nt = AlignIO.read(args.fasta_file, 'fasta').get_alignment_length()
    if len_nt % 3 != 0:
        raise FastaInputError('Problem with the input fasta file. Either is the not a multiple of three or it has indels.')
    LEN_AA = len_nt/3

    # Replace the keyword in the template file:
    with open(args.template) as fh:
        template = fh.read()

    cur_dir = os.getcwd()
    AMBIG = cur_dir + '/igphyml_files/dummy.ambig'

    template = template.replace('LEN_AA', str(LEN_AA))
    template = template.replace('IGPHYML_DIR', IGPHYML_DIR)
    template = template.replace('MODEL', MODEL)
    template = template.replace('AMBIG', AMBIG)


    # Write the new config file:
    with open(args.outfile, 'w') as fh_out:
        fh_out.write(template)

    print('Done making IgPhyML config file.')


def which(executable):
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, executable)):
            return os.path.realpath(os.path.join(path, executable))
    return None




def dedup_fasta(args):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet import generic_dna
    from Bio import AlignIO
    from Bio.Phylo.TreeConstruction import MultipleSeqAlignment

    aln = AlignIO.read(args.infile, 'fasta')

    seqs_unique_counts = {}
    for seq in aln:
        seq.id = seq.id.strip()
        # if id is just an integer, assume it represents count of that sequence
        if seq.id == args.naive:
            naive_seq = str(seq.seq)  # We want to keep the identity of the naive sequence
        elif seq.id.isdigit():
            seqs_unique_counts[str(seq.seq)] = int(seq.id)
        elif str(seq.seq) not in seqs_unique_counts:
            seqs_unique_counts[str(seq.seq)] = 1
        else:
            seqs_unique_counts[str(seq.seq)] += 1

    with open(args.outfile, 'w') as fh_out:
        print('>' + args.naive, file=fh_out)
        print(naive_seq, file=fh_out)
        for i, seq in enumerate(seqs_unique_counts):
            record = SeqRecord(Seq(seq, generic_dna), id=str(i+1)+'_'+str(seqs_unique_counts[seq]))
            print('>' + str(record.id), file=fh_out)
            print(str(record.seq), file=fh_out)

    print('Done deduplicating input fasta file.')


def reroot(args):
    from ete3 import Tree, NodeStyle, TreeStyle, TextFace, add_face_to_node

    try:
        tree = Tree(args.tree)
    except:
        raise TreeFileParsingError('Could not read the input tree. Is this really newick format?')

    rerotted_tree = reroot_tree(tree, pattern=args.pattern, outgroup=args.outgroup)
    rerotted_tree.write(outfile=args.reroot_tree)

    print('Done rerooting input tree.')


def find_node(tree, pattern):
    regex = re.compile(pattern).search
    nodes =  [node for node in tree.traverse() for m in [regex(node.name)] if m]
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
    if tree != node and outgroup == 1:
        s = node.get_sisters() # KBH - want the outgroup with a branch length of (near) zero
        s[0].dist = node.dist * 2
        node.dist = 0.0000001  # KBH - actual zero length branches cause problems
        tree.swap_children()   # KBH - want root to be at the last taxon in the newick file.
    elif tree != node:
        tree.remove_child(node)
        # Notice that an internal node between the outgroup (naive seq.)
        # and the the root is removed because ASR for IgPhyML is done using
        # the naive as an actual root.
        node.add_child(tree.get_children()[0].get_children()[0])
        node.add_child(tree.get_children()[0].get_children()[1])
        tree.dist = node.dist
        node.dist = 0
        tree = node
    return tree



def main():
    import argparse

    parser = argparse.ArgumentParser(description='Tools for running ASR with IgPhyML.')
    subparsers = parser.add_subparsers(help='Which program to run')

    # Parser for ASR_parser subprogram:
    parser_asr = subparsers.add_parser('ASR_parser',
                                       help='Reroot a tree based on node containing a keyword.',
                                       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_asr.add_argument('--tree', required=True, metavar='NEWICK TREE', help='Input tree used for topology.')
    parser_asr.add_argument('--counts', required=True, metavar='ALLELE_FREQUENCY', help='File containing allele frequencies (sequence counts) in the format: "SeqID,Nobs"')
    parser_asr.add_argument('--asr_seq', required=True, help='Input ancestral sequences.')
    parser_asr.add_argument('--outbase', required=True, metavar='FILENAME', help='Filename for the output ASR tree.')
    parser_asr.add_argument('--naive', type=str, default='naive', help='naive sequence id')
    parser_asr.set_defaults(func=ASR_parser)

    # Parser for make_igphyml_config subprogram:
    parser_conf = subparsers.add_parser('make_igphyml_config',
                                        help='Prepare config file for ASR under the HLP16 model.',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_conf.add_argument('--template', help='Config file template', type=str, metavar='FILE', required=True)
    parser_conf.add_argument('--outfile', help='Name of the config file to write.', type=str, metavar='FILENAME', required=True)
    parser_conf.add_argument('--igphyml_exe', help='IgPhyML executable. Will search for the default localt like UNIX `which`.', type=str, required=True)
    parser_conf.add_argument('--model', help='Which model to run? [gy94, hlp16]', type=str, required=True)
    parser_conf.add_argument('--fasta_file', help='To find the length of the amino acid sequence.', type=str, required=True)
    parser_conf.set_defaults(func=make_igphyml_config)

    # Parser for dedup_fasta subprogram:
    parser_dedup = subparsers.add_parser('dedup_fasta',
                                         help='Deduplicate a fasta file.',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_dedup.add_argument('--infile', type=str, help='fasta file with any integer ids indicating frequency', required=True)
    parser_dedup.add_argument('--outfile', type=str, help='Output filename.', required=True)
    parser_dedup.add_argument('--naive', type=str, default='naive', help='naive sequence id')
    parser_dedup.set_defaults(func=dedup_fasta)

    # Parser for reroot subprogram:
    parser_reroot = subparsers.add_parser('reroot',
                                          description='Reroot a tree based on node containing a keyword.',
                                          formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_reroot.add_argument('--tree', required=True, metavar='NEWICK TREE', help='Input tree to root, in newick format.')
    parser_reroot.add_argument('--reroot_tree', required=True, metavar='FILENAME', help='Filename for the output rerooted tree.')
    parser_reroot.add_argument('--pattern', metavar='REGEX PATTERN', required=True, help="Pattern to search for the node to root on.")
    parser_reroot.add_argument('--outgroup', required=False, type=int, default=0, metavar="0/1", help="Set as outgroup instead of tree root.")
    parser_reroot.set_defaults(func=reroot)


    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
