#! usr/bin/env python

import pandas, sys, itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

# Tas data excel file
infile = sys.argv[1]
# .phylip outfile name
outfile = sys.argv[2]

# parse to pandas dataframe
df = pandas.read_excel(open(infile,'rb'), sheetname='Vibratome', skiprows=1, keep_default_na=False)

# get the IgH sequences of big blue cluster in fig 4b
seqs = df[(df.GC==1) & (df.LN==2) & (df['Igh sequence']!='NA')]['Igh sequence'].tolist()
# unique seqs and corresponding abundances
seqs = [(g[0], len(list(g[1]))) for g in itertools.groupby(sorted(seqs))]
# as SeqRecords, with abundance in description
seqs = (SeqRecord(Seq(x[0], generic_dna), id=str(i), description=str(x[1])) for (i, x) in enumerate(seqs))
# put seqs in fasta format, with abundance on comment line
SeqIO.write(seqs, open(outfile, 'w'), "fasta")
# align seqs
muscle_cline = MuscleCommandline('muscle', input=outfile, out=outfile)
#sys.exit('asdf')
muscle_cline()
# load fasta alignment
AlignIO.write(AlignIO.read(outfile, 'fasta'), outfile, 'phylip')
