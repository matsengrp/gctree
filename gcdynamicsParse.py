#! usr/bin/env python

import sqlite3, argparse, gctree
from math import sqrt, pi
from scipy import mean
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from ete3 import Tree

parser = argparse.ArgumentParser(description='parse output db from gcdynamics simulation')
parser.add_argument('db', type=str, help='path to simulation output database')
parser.add_argument('--outbase', type=str, default=None, help='basename for output fasta and pdf')
args = parser.parse_args()

if args.outbase is None:
    args.outbase = args.db

conn = sqlite3.connect(args.db)
c = conn.cursor()
# dictionary of sequence ids and sequence strings
sequences = {sequence_id:sequence for sequence_id, sequence in c.execute('SELECT * FROM sequences')}

# b cell data
bcells = {x[0]:x[1:] for x in c.execute('SELECT * FROM bcells')}
#print set([bcells[bcell][1] for bcell in bcells])

bcell_affinities = {x[1]:x[5] for x in c.execute('SELECT * FROM affinity')}

print len(sequences), len(bcell_affinities)
raw_input(set(bcells.keys()) - set(bcell_affinities.keys()))


# compute sequence affinities as average among cells with that sequence
#sequence_affinities = {seq_id:mean([bcell_affinities[bcell] for bcell in bcells if bcells[bcell][2] == seq_id]) for seq_id in sequences}
    

#assert set([bcells[x][2] for x in bcells]) == set(sequences.keys())

nodes = {bcell:Tree(name=bcell) for bcell in bcells}
for bcell in nodes:
    nodes[bcell].add_feature('frequency', 1)

# first parse out individual lineages
# the roots don't have parents
trees = [nodes[node] for node in nodes if bcells[nodes[node].name][0] is None]

# attach the other nodes using parent data
for bcell in bcells:
    nodes[bcell].add_feature('sequence', sequences[bcells[bcell][2]])
    if bcells[bcell][0] is not None:
        nodes[bcells[bcell][0]].add_child(nodes[bcell])
        nodes[bcell].dist = gctree.hamming_distance(nodes[bcell].sequence, nodes[bcell].up.sequence)
        nodes[bcells[bcell][0]].frequency = 0 #add_feature('frequency', 0)
    else:
        nodes[bcell].dist = 0

# collapse the trees
for tree_i, tree in enumerate(trees):
    # write leaf data
    with open(args.outbase+'.'+str(tree_i)+'.leafdata.fa', 'w') as f:
        f.write('> GL\n')
        f.write(tree.sequence+'\n')
        i = 0
        for leaf in tree.iter_leaves():
            if leaf.frequency != 0:
                i += 1
                f.write('> seq%d\n' % i)
                f.write(leaf.sequence+'\n')
                leaf.name = 'seq%d' % i
    print i, 'simulated observed sequences'
    tree.render(args.outbase+'.tree'+str(tree_i)+'.png')

    # get collapsed tree
    collapsed_tree = gctree.CollapsedTree(tree=tree)
    collapsed_tree.render(args.outbase+'.'+str(tree_i)+'.collapsed_tree.png')



# for spoof data we only want cells in the last round of simulation

#new_aln = MultipleSeqAlignment([germline])
#for i, seq in enumerate(seqs_unique_counts):
#    new_aln.append(SeqRecord(seq, id=str(i+1)+'_'+str(seqs_unique_counts[seq])))
#
#print new_aln.format('phylip')
