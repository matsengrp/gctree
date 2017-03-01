#! /usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import print_function
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import MultipleSeqAlignment

def Tas_parse(aln_file, naive, frame=None):
    aln = AlignIO.read(aln_file, 'fasta')
    sequence_length = aln.get_alignment_length()
    if frame is not None:
        start = frame-1
        end = start + 3*((sequence_length - start)//3)
    else:
        start = 0
        end = sequence_length

    seqs_unique_counts = {}
    for seq in aln:
        # if id is just an integer, assume it represents count of that sequence
        seqstr = str(seq.seq)[start:end]
        if seq.id == naive:
            naive_seq = seqstr
        elif seq.id.isdigit():
            seqs_unique_counts[seqstr] = int(seq.id)
        elif str(seq.seq) not in seqs_unique_counts:
            seqs_unique_counts[seqstr] = 1
        else:
            seqs_unique_counts[seqstr] += 1

    new_aln = MultipleSeqAlignment([SeqRecord(Seq(naive_seq, generic_dna), id=naive)])
    for i, seq in enumerate(seqs_unique_counts):
        new_aln.append(SeqRecord(Seq(seq, generic_dna), id=str(i+1)+'_'+str(seqs_unique_counts[seq])))

    return new_aln

def main():
    parser = argparse.ArgumentParser(description='convert a Victora lab GC fasta file to phylip')
    parser.add_argument('infile', type=str, help='fasta file with any integer ids indicating frequency')
    parser.add_argument('--naive', type=str, default='naive', help='naive sequence id')
    parser.add_argument('--frame', type=int, default=None, help='codon frame', choices=(1,2,3))
    args = parser.parse_args()

    new_aln = Tas_parse(args.infile, args.naive, frame=args.frame)

    print(new_aln.format('phylip'))

if __name__ == '__main__':
    main()
