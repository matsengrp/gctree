#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
Parse the output ancestral state reconstruction from IgPhyML together with the topology
to create an ete3 tree with the ancestral sequences.
'''

from __future__ import division, print_function
import sys
import os
import warnings
sys.path.append(os.path.abspath("/fh/fast/matsen_e/kdavidse/gctree/bin"))


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
    from ete3 import Tree, TreeNode, NodeStyle, TreeStyle, TextFace, add_face_to_node, CircleFace, faces, AttrFace
    from Bio import SeqIO
    from gctree import CollapsedForest, CollapsedTree

    try:
        tree = Tree(args.tree)
    except:
        raise TreeFileParsingError('Could not read the input tree. Is this really newick format?')

    tree.add_feature('frequency', 0)
    tree = map_asr_to_tree(args.asr_seq, tree)

    # Reroot to make the naive sequence the real root instead of just an outgroup:
    tree = reroot_tree(tree)

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




def make_igphyml_config(args):
    import re

    igphyml_path = which(args.igphyml_exe.rstrip('/'))
    assert(igphyml_path is not None)
    IGPHYML_DIR = re.sub(r'/src/\w+$', '', igphyml_path)
    MODEL = args.model
    assert(MODEL in ['gy94', 'hlp16'])
    # Find the length of the translated DNA sequence:
    LEN_AA = None
    with open(args.fasta_file) as fh:
        for l in fh:
            if not l.startswith('>') and l != '':
                this_LEN_AA = len(l.strip()) / 3.0
                if int(this_LEN_AA) != this_LEN_AA or (LEN_AA is not None and this_LEN_AA != LEN_AA):
                    raise FastaInputError('Problem with the input fasta file. Either is the not a multiple of three or it has indels.')
                elif LEN_AA is None:
                    LEN_AA = int(this_LEN_AA)
    assert(LEN_AA is not None)

    # Replace the keyword in the template file:
    with open(args.template) as fh:
        template = fh.read()

    template = template.replace('LEN_AA', str(LEN_AA))
    template = template.replace('IGPHYML_DIR', IGPHYML_DIR)
    template = template.replace('MODEL', MODEL)

    # Write the new config file:
    with open(args.outfile, 'w') as fh_out:
        fh_out.write(template)


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




def reroot(args):
    from ete3 import Tree, NodeStyle, TreeStyle, TextFace, add_face_to_node
    import os
    import re

    try:
        tree = Tree(args.tree)
    except:
        raise TreeFileParsingError('Could not read the input tree. Is this really newick format?')

    rerotted_tree = reroot_tree(tree, pattern=args.pattern, outgroup=args.outgroup)
    rerotted_tree.write(outfile=args.reroot_tree)


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
    if tree != node and outgroup == 1:
        s = node.get_sisters() # KBH - want the outgroup with a branch length of (near) zero
        s[0].dist = node.dist * 2
        node.dist = 0.0000001  # KBH - actual zero length branches cause problems
        tree.swap_children()   # KBH - want root to be at the last taxon in the newick file.
    elif tree != node:
        tree.remove_child(node)
        node.add_child(tree)
        tree.dist = node.dist
        node.dist = 0
        tree = node        

    return tree



def main():
    import argparse

    parser = argparse.ArgumentParser(description='Tools for running ASR with IgPhyML.')
    subparsers = parser.add_subparsers(help='Which program to run')

    # Parser for ASR_parser subprogram:
    parser = subparsers.add_parser('ASR_parser',
                                   help='Reroot a tree based on node containing a keyword.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tree', required=True, metavar='NEWICK TREE', help='Input tree used for topology.')
    parser.add_argument('--asr_seq', required=True, help='Input ancestral sequences.')
    parser.add_argument('--outbase', required=True, metavar='FILENAME', help='Filename for the output ASR tree.')
    parser_test.set_defaults(func=ASR_parser)

    # Parser for make_igphyml_config subprogram:
    parser = subparsers.add_parser('make_igphyml_config',
                                   help='Prepare config file for ASR under the HLP16 model.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--template', help='Config file template', type=str, metavar='FILE', required=True)
    parser.add_argument('--outfile', help='Name of the config file to write.', type=str, metavar='FILENAME', required=True)
    parser.add_argument('--igphyml_exe', help='IgPhyML executable. Will search for the default localt like UNIX `which`.', type=str, required=True)
    parser.add_argument('--model', help='Which model to run? [gy94, hlp16]', type=str, required=True)
    parser.add_argument('--fasta_file', help='To find the length of the amino acid sequence.', type=str, required=True)
    parser_test.set_defaults(func=make_igphyml_config)

    # Parser for dedup_fasta subprogram:
    parser = subparsers.add_parser('dedup_fasta',
                                   help='Deduplicate a fasta file.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--infile', type=str, help='fasta file with any integer ids indicating frequency', required=True)
    parser.add_argument('--outfile', type=str, help='Output filename.', required=True)
    parser.add_argument('--naive', type=str, default='naive', help='naive sequence id')
    parser_test.set_defaults(func=dedup_fasta)

    # Parser for reroot subprogram:
    parser = subparsers.add_parser('reroot',
                                   description='Reroot a tree based on node containing a keyword.',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tree', required=True, metavar='NEWICK TREE', help='Input tree to root, in newick format.')
    parser.add_argument('--reroot_tree', required=True, metavar='FILENAME', help='Filename for the output rerooted tree.')
    parser.add_argument('--pattern', metavar='REGEX PATTERN', required=True, help="Pattern to search for the node to root on.")
    parser.add_argument('--outgroup', required=False, type=int, default=0, metavar="0/1", help="Set as outgroup instead of tree root.")
    parser_test.set_defaults(func=reroot)


    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
