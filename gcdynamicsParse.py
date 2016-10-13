#! usr/bin/env python

import sqlite3, argparse
from ete3 import Tree, NodeStyle, TreeStyle, TextFace, add_face_to_node, CircleFace, faces, AttrFace
from math import sqrt, pi
from scipy import mean
from recurse import hamming_distance
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO

parser = argparse.ArgumentParser(description='parse output db from gcdynamics simulation')
parser.add_argument('db', type=str, help='path to simulation output database')
args = parser.parse_args()

conn = sqlite3.connect(args.db)
c = conn.cursor()
# dictionary of sequence ids and sequence strings
sequences = {sequence_id:sequence for sequence_id, sequence in c.execute('SELECT * FROM sequences')}

# b cell data
bcells = {x[0]:x[1:] for x in c.execute('SELECT * FROM bcells')}

bcell_affinities = {x[0]:x[5] for x in c.execute('SELECT * FROM effectiveAffinity')}

print set(bcells.keys()) - set(bcell_affinities.keys())

# compute sequence affinities as average among cells with that sequence
sequence_affinities = {seq_id:mean([bcell_affinities[bcell] for bcell in bcells if bcells[bcell][2] == seq_id]) for seq_id in sequences}
    

#assert set([bcells[x][2] for x in bcells]) == set(sequences.keys())

nodes = {bcell:Tree(name=bcell) for bcell in bcells}

# first parse out individual lineages
# the roots don't have parents
trees = [nodes[node] for node in nodes if bcells[nodes[node].name][0] is None]

# attach the other nodes using parent data
for bcell in bcells:
    if bcells[bcell][0] is not None:
        nodes[bcells[bcell][0]].add_child(nodes[bcell])
        nodes[bcell].dist = hamming_distance(sequences[bcells[nodes[bcell].name][2]], sequences[bcells[bcells[bcell][0]][2]])

# collapse the trees
for tree_i, tree in enumerate(trees):
    collapsed_tree = []
    # the number of clonal leaf descendents is number of leaves we can get to on zero-length edges
    # root first
    clone_leaves = sum(1 for leaf in tree.iter_leaves() if tree.get_distance(leaf) == 0)

    # to get mutant offspring, first consider all nodes that are distance zero
    # the mutant offspring are children of these nodes with nonzero edge length
    mutant_offspring_nodes = [child for node in tree.traverse() if tree.get_distance(node) == 0 for child in node.children if child.dist != 0]
    mutant_offspring_edge_lengths = [node.dist for node in mutant_offspring_nodes]
    assert 0 not in mutant_offspring_edge_lengths
    collapsed_tree.append((clone_leaves, mutant_offspring_edge_lengths))
    # recurse into the mutant offspring
    done = False
    while not done:
        new_mutant_offspring_nodes = []
        for mutant in mutant_offspring_nodes:
            clone_leaves = sum(1 for leaf in mutant.iter_leaves() if mutant.get_distance(leaf) == 0)
            new_mutant_offspring_nodes_from_this_mutant = [child for node in mutant.traverse() if mutant.get_distance(node) == 0 for child in node.children if child.dist != 0]
            mutant_offspring_edge_lengths = [node.dist for node in new_mutant_offspring_nodes_from_this_mutant]
            assert 0 not in mutant_offspring_edge_lengths
            collapsed_tree.append((clone_leaves, mutant_offspring_edge_lengths))
            new_mutant_offspring_nodes.extend(new_mutant_offspring_nodes_from_this_mutant)
        mutant_offspring_nodes = new_mutant_offspring_nodes
        if len(new_mutant_offspring_nodes) == 0:
            done = True

    # make an ete version of collapsed tree for plotting
    nodes = [Tree(name=clone_leaves) for clone_leaves, mutant_offspring in collapsed_tree]
    nodes[0].dist = 0 # zero length edge for root node
    gen_size = 1
    gen_start_index = 0
    terminated = False
    while not terminated:
        k = 0
        for i in range(gen_start_index,gen_start_index + gen_size):
            for j in range(len(collapsed_tree[i][1])):
                nodes[i].add_child(nodes[gen_start_index+gen_size+k], dist=collapsed_tree[i][1][j])
                k += 1
        new_gen_start_index = gen_start_index + gen_size
        gen_size = sum(len(x[1]) for x in collapsed_tree[gen_start_index:(gen_start_index + gen_size)])
        gen_start_index = new_gen_start_index
        if gen_size == 0:
            terminated = True

    for node in nodes:
        nstyle = NodeStyle()
        if node.name == 0:
            nstyle['size'] = 5
            nstyle['fgcolor'] = 'grey'
        else:
            nstyle['size'] = 3*2*sqrt(pi*node.name)
            nstyle['fgcolor'] = 'black'
        if node.up is not None:
            nstyle['hz_line_color'] = 'black'
        node.set_style(nstyle)

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.rotation = 90
    def my_layout(node):
        if node.name > 1:
            N = TextFace(node.name, fsize=14, fgcolor='black')
            N.rotation = -90
            faces.add_face_to_node(N, node, 0, position='branch-top')
    ts.layout_fn = my_layout
    nodes[0].render(args.db+'.'+str(tree_i+1)+'.pdf', tree_style=ts)




# for spoof data we only want cells in the last round of simulation

#new_aln = MultipleSeqAlignment([germline])
#for i, seq in enumerate(seqs_unique_counts):
#    new_aln.append(SeqRecord(seq, id=str(i+1)+'_'+str(seqs_unique_counts[seq])))
#
#print new_aln.format('phylip')
