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

# Tas data file, fasta
infile = sys.argv[1]
## outfile basename
#outfile = sys.argv[2]

## parse to pandas dataframe
#df = pandas.read_excel(open(infile,'rb'), sheetname='Vibratome', skiprows=1, keep_default_na=False)

## get the IgH sequences of big blue cluster in fig 4b
#seqs = df[(df.GC==1) & (df.LN==2) & (df['Igh sequence']!='NA')]['Igh sequence'].tolist()

aln = AlignIO.read(infile, 'fasta')

seqs_unique_counts = {}
for seq in aln:
    if seq.id == 'GL':
        germline = seq
    elif seq.id == '17':
        seqs_unique_counts[seq.seq] = 17
    elif seq.seq not in seqs_unique_counts:
        seqs_unique_counts[seq.seq] = 1
    else:
        seqs_unique_counts[seq.seq] += 1

new_aln = MultipleSeqAlignment([germline])
for i, seq in enumerate(seqs_unique_counts):
    new_aln.append(SeqRecord(seq, id=str(i+1)+'_'+str(seqs_unique_counts[seq])))

print new_aln.format('phylip')
sys.exit()

# the big expanded clone is marked with a single entry ">17", we have to duplicate that
# we also need to remove their "GL" sequence (unmutated Vh ancestor)
new_aln = MultipleSeqAlignment([])
for seq in aln:
    if seq.id == '17':
        for i in range(1,18):
            seq = copy.copy(seq)
            seq.id='17_'+str(i)
            new_aln.append(seq)
    elif seq.id == 'GL':
        unmutated_ancestor = MultipleSeqAlignment([seq])
    else:
        new_aln.append(seq)

# put unmutated ancestor first
unmutated_ancestor.extend(new_aln)
new_aln = unmutated_ancestor

print new_aln.format('phylip')

sys.exit()

seqs = [seq.seq for seq in seqs]
# unique seqs and corresponding abundances
#seqs = [(g[0], len(list(g[1]))) for g in itertools.groupby(sorted(seqs))]
seqs = [(seq, 1) for seq in seqs]
# as SeqRecords, with abundance in id
seqs = (SeqRecord(x[0], id=str(i)+'_'+str(x[1])) for (i, x) in enumerate(seqs))
# put seqs in fasta format, with abundance on comment line
SeqIO.write(seqs, open(outfile, 'w'), "fasta")
# align seqs
muscle_cline = MuscleCommandline('muscle', input=outfile, out=outfile)
#sys.exit('asdf')
muscle_cline()
# load fasta alignment and convert to phylip, and load phylip
AlignIO.write(AlignIO.read(outfile, 'fasta'), outfile, 'phylip')
aln = AlignIO.read(outfile, 'phylip')

# parsimony!
scorer = ParsimonyScorer()
searcher = NNITreeSearcher(scorer)
constructor = ParsimonyTreeConstructor(searcher)
pars_tree = constructor.build_tree(aln)
Phylo.draw_ascii(pars_tree)
