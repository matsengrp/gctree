#! usr/bin/env python

#import pandas
import sys, itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import *
#from Bio.Emboss.Applications import FDNAParsCommandline

# fasta
infile = sys.argv[1]
germline = sys.argv[2]
## outfile basename
#outfile = sys.argv[2]

## parse to pandas dataframe
#df = pandas.read_excel(open(infile,'rb'), sheetname='Vibratome', skiprows=1, keep_default_na=False)

## get the IgH sequences of big blue cluster in fig 4b
#seqs = df[(df.GC==1) & (df.LN==2) & (df['Igh sequence']!='NA')]['Igh sequence'].tolist()

aln = AlignIO.read(infile, 'fasta')

seqs_unique_counts = {}
for seq in aln:
    if seq.id == germline:
        germline_seq = seq
    # if id is just an integer, assume it represents count of that sequence
    elif seq.id.isdigit():
        seqs_unique_counts[str(seq.seq)] = int(seq.id)
    elif str(seq.seq) not in seqs_unique_counts:
        seqs_unique_counts[str(seq.seq)] = 1
    else:
        seqs_unique_counts[str(seq.seq)] += 1

new_aln = MultipleSeqAlignment([germline_seq])
for i, seq in enumerate(seqs_unique_counts):
    new_aln.append(SeqRecord(Seq(seq, generic_dna), id=str(i+1)+'_'+str(seqs_unique_counts[seq])))

print new_aln.format('phylip')
sys.exit()
