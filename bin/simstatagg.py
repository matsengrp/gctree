#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
aggregation plots of metrics from several simulation runs with same parameters
'''

from __future__ import division, print_function
import scipy, matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib import rc, ticker
import pandas as pd
import argparse
from deduplicate import fasta_parse
from utils import hamming_distance
import seaborn as sns
sns.set(style="white", color_codes=True)
from matplotlib.backends.backend_pdf import PdfPages


parser = argparse.ArgumentParser(description='aggregate validation of repeated runs with same parameters')
parser.add_argument('input', type=str, nargs='+', help='gctree.simulation.stats.tsv files')
parser.add_argument('--outbase', type=str, help='output file base name')
parser.add_argument('--experimental', type=str, help='experimental .fasta (like Tas et al. data)', default=None)
args = parser.parse_args()

for i, fname in enumerate(args.input):
    df = pd.read_csv(fname, sep='\t')
    df['simulation'] = i + 1
    if i == 0:
        aggdat = df
    else:
        aggdat = aggdat.append(df, ignore_index=True)

# aggdat.to_csv(args.outbase+'.tsv', sep='\t', index=False)

# fields = list(aggdat.columns)
# fields.remove('simulation')
# nfields = len(fields)

sims = set(aggdat['simulation'])
nsims = len(sims)

if args.experimental is not None:
    new_aln, counts = fasta_parse(args.experimental, naive='GL', converter='tas')[:2]
    exp_dict = {seq.id:str(seq.seq) for seq in new_aln}
    naive_id = [seq for seq in exp_dict if 'gl' in seq][0]
    frequency, distance_from_naive, degree = zip(*[(counts[seq],
                                                    hamming_distance(exp_dict[seq], exp_dict[naive_id]),
                                                    sum(hamming_distance(exp_dict[seq], exp_dict[seq2]) == 1 for seq2 in exp_dict if seq2 is not seq and counts[seq2] != 0))
                                                   for seq in exp_dict if counts[seq] != 0])
    exp_stats = pd.DataFrame({'genotype abundance':frequency,
                              'Hamming distance to root genotype':distance_from_naive,
                              'Hamming neighbor genotypes':degree})

# bw = .3
alpha = min([.9, 20/nsims])
bins = range(max(aggdat['Hamming distance to root genotype'].max(), exp_stats['Hamming distance to root genotype'].max() if args.experimental is not None else 0) + 2)

plt.figure(figsize=(6, 3))
plt.subplot(1, 2, 1)
for simulation, simulation_aggdat in aggdat.groupby('simulation'):
    sns.distplot(simulation_aggdat['Hamming distance to root genotype'], bins=bins, kde=False, hist_kws={'histtype':'step', 'cumulative':True, 'alpha':alpha, 'lw':1})
if args.experimental is not None:
    sns.distplot(exp_stats['Hamming distance to root genotype'], bins=bins, kde=False, hist_kws={'histtype':'step', 'cumulative':True, 'color':'k', 'lw':3, 'alpha':.8})
plt.xlabel('Hamming distance to root genotype')
plt.xlim([0, bins[-1]])
plt.ylabel('observed genotypes')
plt.tight_layout()

plt.subplot(1, 2, 2)
# g = sns.JointGrid('genotype abundance', 'Hamming neighbor genotypes', aggdat, space=0)
# levels = scipy.logspace(-2, 2, 10)
xbins = range(max(aggdat['genotype abundance'].max(), exp_stats['genotype abundance'].max() if args.experimental is not None else 0) + 2)
ybins = range(max(aggdat['Hamming neighbor genotypes'].max(), exp_stats['Hamming neighbor genotypes'].max() if args.experimental is not None else 0) + 2)
for simulation, simulation_aggdat in aggdat.groupby('simulation'):
    plt.plot(simulation_aggdat['genotype abundance'], simulation_aggdat['Hamming neighbor genotypes'], '+', mew=1, alpha=alpha/2)
if args.experimental is not None:
    plt.plot(exp_stats['genotype abundance'], exp_stats['Hamming neighbor genotypes'], 'o', mew=2, alpha=.8, markerfacecolor='none', color='k')
plt.xlabel('genotype abundance')
plt.ylabel('Hamming neighbor genotypes')
plt.xscale('symlog')
plt.yscale('symlog')
plt.xlim([.9, None])
plt.ylim([-.1, None])
plt.tight_layout()
plt.savefig(args.outbase+'.pdf')

# for sim in sims:
#     sns.jointplot(aggdat[aggdat['simulation']==sim][fields[0]],aggdat[aggdat['simulation']==sim][fields[1]], kind='kde', space=0, kwargs=dict(bw=.5, shaded=False))

# # g = sns.pairplot(aggdat, diag_kind='kde', hue='simulation', vars=fields)
# g = sns.PairGrid(aggdat, hue='simulation', vars=fields)
# g = g.map_diag(sns.kdeplot, bw=.4, alpha=20/nsims, lw=1)
# g = g.map_upper(sns.kdeplot, levels=scipy.linspace(0, 1, .1))
# g = g.map_lower(plt.scatter)
# # g = g.map_lower(sns.kdeplot)
# g.set(xscale='symlog', yscale='symlog', xlim=[0, None], ylim=[0, None])
# # plt.xscale('symlog')
# # ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())


# for i, field in enumerate(fields, 1):
#     ax = fig.add_subplot(1, nfields, i)
#     for simulation in sims:
#         sns.kdeplot(aggdat[aggdat['simulation']==simulation][field], legend=False, alpha=20/nsims, cumulative=True, lw=1)
#         plt.xlabel(field)
#         plt.xlim([min(aggdat[field]), max(aggdat[field])])
#         plt.ylim([0, 1])
#         plt.xscale('symlog')
#         ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

# plt.savefig(args.outbase+'.pdf')
