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
    parser.add_argument('--germline', type=str, default='GL', help='germline sequence id')
    args = parser.parse_args()

    aln = AlignIO.read(args.infile, 'fasta')

    seqs_unique_counts = {}
    for seq in aln:
        # if id is just an integer, assume it represents count of that sequence
        if seq.id == args.germline:
            germline_seq = str(seq.seq)
        elif seq.id.isdigit():
            seqs_unique_counts[str(seq.seq)] = int(seq.id)
        elif str(seq.seq) not in seqs_unique_counts:
            seqs_unique_counts[str(seq.seq)] = 1
        else:
            seqs_unique_counts[str(seq.seq)] += 1

    new_aln = MultipleSeqAlignment([SeqRecord(Seq(germline_seq, generic_dna), id=args.germline)])
    for i, seq in enumerate(seqs_unique_counts):
        new_aln.append(SeqRecord(Seq(seq, generic_dna), id=str(i+1)+'_'+str(seqs_unique_counts[seq])))
    print(new_aln.format('phylip'))

if __name__ == '__main__':
    main()
