#! usr/bin/env python

from __future__ import print_function
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import MultipleSeqAlignment

def main():
    parser = argparse.ArgumentParser(description='convert a Victora lab GC fasta file to phylip')
    parser.add_argument('infile', type=str, help='fasta file with any integer ids indicating frequency')
    args = parser.parse_args()

    aln = AlignIO.read(args.infile, 'fasta')

    seqs_unique_counts = {}
    for seq in aln:
        # if id is just an integer, assume it represents count of that sequence
        if seq.id.isdigit():
            seqs_unique_counts[str(seq.seq)] = int(seq.id)
        elif str(seq.seq) not in seqs_unique_counts:
            seqs_unique_counts[str(seq.seq)] = 1
        else:
            seqs_unique_counts[str(seq.seq)] += 1

    print(MultipleSeqAlignment([SeqRecord(Seq(seq, generic_dna), id=str(i+1)+'_'+str(seqs_unique_counts[seq])) for i, seq in enumerate(seqs_unique_counts)]).format('phylip'))

if __name__ == '__main__':
    main()
