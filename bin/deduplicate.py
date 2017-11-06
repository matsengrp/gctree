#! /usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import print_function
import argparse
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import MultipleSeqAlignment
from collections import defaultdict, Counter

def fasta_parse(aln_file, naive, frame=None, aln_file2=None, converter=None):
    # naive = naive.lower()
    if aln_file.endswith('fasta') or aln_file.endswith('fa'):
        aln_format = 'fasta'
    elif aln_file.endswith('phylip') or aln_file.endswith('phy'):
        aln_format = 'phylip'
    else:
        raise ValueError('unrecognized alignment file type: ' + aln_file)
    aln = AlignIO.read(aln_file, aln_format)
    if aln_file2 is not None:
        assert frame is None
        aln_combined = MultipleSeqAlignment([])
        if aln_file2.endswith('fasta') or aln_file2.endswith('fa'):
            aln_format2 = 'fasta'
        elif aln_file2.endswith('phylip') or aln_file2.endswith('phy'):
            aln_format2 = 'phylip'
        else:
            raise ValueError('unrecognized alignment file type: ' + aln_file2)
        aln2 = AlignIO.read(aln_file2, aln_format2)
        for seq in aln:
            cell = (seq.id[:-1] if seq.id != naive else naive)
            for seq2 in aln2:
                cell2 = (seq2.id[:-1] if seq2.id != naive else naive)
                if cell2 == cell:
                    aln_combined.append(SeqRecord(Seq(str(seq.seq) + str(seq2.seq), generic_dna), id=cell))
        aln = aln_combined

    sequence_length = aln.get_alignment_length()
    if frame is not None:
        start = frame-1
        end = start + 3*((sequence_length - start)//3)
    else:
        start = 0
        end = sequence_length

    seqs_unique_counts = defaultdict(list)
    id_set = set()
    naive_seq = None
    for seq in aln:
        # seq.id = seq.id.lower()
        # if id is just an integer, assume it represents count of that sequence
        if seq.id in id_set:
            raise ValueError('Sequence ID found multiple times:', seq.id)
        else:
            id_set.add(seq.id)
        seqstr = str(seq.seq)[start:end]
        if seq.id == naive:
            naive_seq = seqstr
            if seqstr not in seqs_unique_counts:
                seqs_unique_counts[seqstr] = [] # no observed naive unless we see it elsewhere
        elif seq.id.isdigit() and converter is not None:
            if converter.lower() == 'tas':
                seqs_unique_counts[seqstr] = [seq.id for _ in range(int(seq.id))]
            else:
                raise ValueError('invalid converter: '+converter)
        else:
            seqs_unique_counts[seqstr].append(seq.id)

    if naive_seq is None:
        raise ValueError('naive seq id {} not found'.format(naive))

    new_aln = MultipleSeqAlignment([SeqRecord(Seq(naive_seq, generic_dna), id=naive.lower())])
    counts = {naive.lower(): len(seqs_unique_counts[naive_seq])}  # Add the count for the naive sequence
    id_map = {naive.lower(): [x for x in seqs_unique_counts[naive_seq] if x != naive]}
    del seqs_unique_counts[naive_seq]  # Now delete the naive so it does not appear twice
    for i, seq in enumerate(seqs_unique_counts, 1):
        new_id = 'seq'+str(i)
        new_aln.append(SeqRecord(Seq(seq, generic_dna), id=new_id))
        counts[new_id] = len(seqs_unique_counts[seq])
        id_map[new_id] = seqs_unique_counts[seq]

    return new_aln, counts, id_map

def check_header(header):
    try:
        header.decode('ascii')
    except UnicodeDecodeError as e:
        print('Sequence header must be an ascii-encoded string:', header)
        raise e
    if len(header) > 10:
        print('Sequence headers must be shorter than 10 characters:', header)
        raise Exception
    try:
        int(header)
        raise Exception('Sequence headers must be distinguishable from an integer. Please add a non number character.')
    except:
        pass


def main():
    parser = argparse.ArgumentParser(description='Convert a fasta file to philyp format. Headers must be a unique ID of less than '
                                                 'or equal to 10 ASCII characters. A special option for converting a Victora lab '
                                                 'GC fasta file to phylip is also included. All headers are converted to lower case.')
    parser.add_argument('infile', type=str, nargs='+', help='Fasta file with less than or equal to 10 characters unique header ID. '
                                                 'For Victora data any integer ids indicates frequency.'
                                                 'Because dnapars will name internal nodes by intergers a node name must include'
                                                 'at least one non number character.')
    parser.add_argument('--countfile', type=str, default=None, help='filename for the output file containing the counts.')
    parser.add_argument('--idmapfile', type=str, default=None, help='filename for the output file containing the map of new unique ids to original seq ids.')
    parser.add_argument('--converter', type=str, default=None, help='use a special format convertion scheme e.g. for a Vitora lab GC fasta file. Options: [tas]')
    specified_coverters = ['tas']
    parser.add_argument('--naive', type=str, default='naive', help='naive sequence id')
    parser.add_argument('--frame', type=int, default=None, help='codon frame', choices=(1,2,3))
    parser.add_argument('--colorfile', type=str, default=None, help='optional input csv filename for colors of each cell.')
    parser.add_argument('--colormap', type=str, default=None, help='optional output filename for colors map.')
    args = parser.parse_args()

    if args.converter is not None and args.converter.lower() not in specified_coverters:
        print('Cannot find the specified converter:', args.converter)
        print('Allowed converters:', specified_coverters.join(','))
        raise Exception
    new_aln, counts, id_map = fasta_parse(args.infile[0],
                                          args.naive,
                                          frame=args.frame,
                                          aln_file2=args.infile[1] if len(args.infile) == 2 else None,
                                          converter=args.converter)
    print(new_aln.format('phylip'))
    if args.countfile is not None:
        with open(args.countfile, 'w') as f:
            for seqID, count in counts.items():
                print('{},{}'.format(seqID, count), file=f)
    if args.idmapfile is not None:
        with open(args.idmapfile, 'w') as f:
            for seq_id, cell_ids in id_map.items():
                print('{},{}'.format(seq_id, ':'.join(cell_ids)), file=f)
    if args.colorfile is not None:
        with open(args.colorfile, 'r') as f:
            cell_colors = dict(line.rstrip().split(',') for line in f)
        with open(args.colormap, 'w') as f:
            for seq_id, cell_ids in id_map.items():
                colors_counts = Counter(cell_colors[cell_id] if cell_id in cell_colors else None for cell_id in cell_ids).items()
                if len(colors_counts) == 1 and colors_counts[0][0] != None and colors_counts != []:
                    colors_str = colors_counts[0][0]
                    print('{}\t{}'.format(seq_id, colors_str), file=f)
                elif len(colors_counts) > 1:
                    colors_strs = ['{}:{}'.format(color, count) for color, count in colors_counts]
                    print('{}\t{}'.format(seq_id, ','.join(colors_strs)), file=f)


if __name__ == '__main__':
    main()
