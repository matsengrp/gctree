#! /usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import print_function
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import MultipleSeqAlignment

def Tas_parse(aln_file, count_fnam, naive, frame=None):
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
            if seqstr not in seqs_unique_counts:
                seqs_unique_counts[seqstr] = 0 # no observed naive unless we see it elsewhere
        elif seq.id.isdigit():
            seqs_unique_counts[seqstr] = int(seq.id)
        elif str(seq.seq) not in seqs_unique_counts:
            seqs_unique_counts[seqstr] = 1
        else:
            seqs_unique_counts[seqstr] += 1

    new_aln = MultipleSeqAlignment([SeqRecord(Seq(naive_seq, generic_dna), id=naive+'_'+str(seqs_unique_counts[naive_seq]))])
    for i, seq in enumerate(seqs_unique_counts):
        new_aln.append(SeqRecord(Seq(seq, generic_dna), id=str(i+1)+'_'+str(seqs_unique_counts[seq])))

    return new_aln


def check_header(header):
   try:
       header.decode('ascii')
   except UnicodeDecodeError as e:
       print('Sequence header must be an ascii-encoded string:', header)
       raise e
   if len(header) > 10:
       print('Sequence headers must be shorter than 10 characters:', header)
       raise Exception


def default_parse(aln_file, count_fnam, naive, frame=None):
    aln = AlignIO.read(aln_file, 'fasta')
    sequence_length = aln.get_alignment_length()
    if frame is not None:
        start = frame-1
        end = start + 3*((sequence_length - start)//3)
    else:
        start = 0
        end = sequence_length

    seqs_unique_counts = {}
    seq2id = dict()
    id_set = set()
    for seq in aln:
        check_header(seq.id)
        # if id is just an integer, assume it represents count of that sequence
        seqstr = str(seq.seq)[start:end]
        if seq.id == naive:
            naive_seq = seqstr
            if seqstr not in seqs_unique_counts:
                seqs_unique_counts[seqstr] = 0 # no observed naive unless we see it elsewhere
        elif str(seq.seq) not in seqs_unique_counts:
            seqs_unique_counts[seqstr] = 1
        else:
            seqs_unique_counts[seqstr] += 1

        if seq.id in id_set:
            print('Sequence ID found multiple times:', seq.id)
            raise Exception
        else:
            id_set.add(seq.id)

        if seqstr not in seq2id:
            seq2id[seqstr] = seq.id

    fh_out = open(basename+'.counts', 'w')
    new_aln = MultipleSeqAlignment([SeqRecord(Seq(naive_seq, generic_dna), id=naive)])
    for i, seq in enumerate(seqs_unique_counts):
        new_aln.append(SeqRecord(Seq(seq, generic_dna), id=seq2id[seq]))
        print('{},{}'.format(seq2id[seq], str(seqs_unique_counts[seq])), file=fh_out)
    fh_out.close()
    return new_aln


def main():
    parser = argparse.ArgumentParser(description='Convert a fasta file to philyp format. Headers must be a unique ID of less than '
                                                 'or equal to 10 ASCII characters. A special option for converting a Victora lab '
                                                 'GC fasta file to phylip is also included.')
    parser.add_argument('infile', type=str, help='Fasta file with less than or equal to 10 characters unique header ID. '
                                                 'For Vitora data any integer ids indicats frequency.')
    parser.add_argument('countfile', type=str, help='Filename for the output file containing the counts.')
    parser.add_argument('--converter', type=str, help='Use a special format convertion scheme e.g. for a Vitora lab GC fasta file. Options: [tas]')
    specified_coverters = ['tas']
    parser.add_argument('--naive', type=str, default='naive', help='naive sequence id')
    parser.add_argument('--frame', type=int, default=None, help='codon frame', choices=(1,2,3))
    args = parser.parse_args()

    if args.converter is not None and args.converter.lower() in specified_coverters:        
        new_aln = Tas_parse(args.infile, args.countfile, args.naive, frame=args.frame)
    elif args.converter is not None and args.converter.lower() not in specified_coverters:
        print('Cannot find the specified converter:', args.converter)
        print('Allowed converters:', specified_coverters.join(','))
        raise Exception
    else:
        new_aln = default_parse(args.infile, args.countfile, args.naive, frame=args.frame)
    print(new_aln.format('phylip'))

if __name__ == '__main__':
    main()

