#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
aggregation plots of metrics from several simulation runs with same parameters
'''

from __future__ import division, print_function
import scipy, matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
from matplotlib import rc, ticker
import pandas as pd
import argparse
from fasta2phylip import Tas_parse
# from gctree import hamming_distance
import seaborn as sns
sns.set(style="ticks", color_codes=True)
from matplotlib.backends.backend_pdf import PdfPages

# something in numpy/seaborn breaks if I import this from gctree
def hamming_distance(seq1, seq2):
    '''Hamming distance between two sequences of equal length'''
    return sum(x != y for x, y in zip(seq1, seq2))

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

fields = list(aggdat.columns)
fields.remove('simulation')
nfields = len(fields)

sims = set(aggdat['simulation'])
nsims = len(sims)

if args.experimental is not None:
    new_aln, counts = Tas_parse(args.experimental, naive='GL')
    exp_dict = {seq.id:str(seq.seq) for seq in new_aln}
    naive_id = [seq for seq in exp_dict if 'gl' in seq][0]
    frequency, distance_from_naive, degree = zip(*[(counts[seq],
                                                    hamming_distance(exp_dict[seq], exp_dict[naive_id]),
                                                    sum(hamming_distance(exp_dict[seq], exp_dict[seq2]) == 1 for seq2 in exp_dict if seq2 is not seq and counts[seq2] != 0))
                                                   for seq in exp_dict if counts[seq] != 0])
    exp_stats = pd.DataFrame({'allele frequency':frequency,
                              'Hamming distance to naive sequence':distance_from_naive,
                              'degree':degree})

pp = PdfPages(args.outbase+'.pdf')
bw = .3
alpha = min([.9, 20/nsims])

fig = plt.figure()
for simulation, simulation_aggdat in aggdat.groupby('simulation'):
    sns.kdeplot(simulation_aggdat[fields[0]], cumulative=True, bw=bw, alpha=alpha, legend=False, lw=1)
if args.experimental is not None:
    sns.kdeplot(exp_stats[fields[0]], cumulative=True, bw=bw, color='k', legend=False, lw=3, alpha=.8)
plt.xlabel(fields[0])
plt.xlim([0, None])
pp.savefig()

fig = plt.figure()
g = sns.JointGrid(fields[1], fields[2], aggdat, space=0)
levels = scipy.logspace(-2, 2, 10)
for simulation, simulation_aggdat in aggdat.groupby('simulation'):
    sns.kdeplot(simulation_aggdat[fields[1]], bw=bw, ax=g.ax_marg_x, legend=False, alpha=alpha, lw=1)
    sns.kdeplot(simulation_aggdat[fields[2]], bw=bw, ax=g.ax_marg_y, vertical=True, legend=False, alpha=alpha, lw=1)
    g.ax_joint.plot(simulation_aggdat[fields[1]], simulation_aggdat[fields[2]], '+', mew=1, alpha=alpha/2)
    # sns.kdeplot(simulation_aggdat[fields[1]], simulation_aggdat[fields[2]], bw=3*bw, ax=g.ax_joint, levels=levels, alpha=alpha, shade=False)
if args.experimental is not None:
    sns.kdeplot(exp_stats[fields[1]], bw=bw, ax=g.ax_marg_x, color='k', legend=False, alpha=.8, lw=3)
    sns.kdeplot(exp_stats[fields[2]], bw=bw, ax=g.ax_marg_y, color='k', vertical=True, legend=False, alpha=.8, lw=3)
    g.ax_joint.plot(exp_stats[fields[1]], exp_stats[fields[2]], 'o', mew=2, alpha=.8, markerfacecolor='none', color='k')


g.ax_joint.set_xscale('symlog')
g.ax_joint.set_xlim([0, None])
g.ax_joint.set_yscale('symlog')
g.ax_joint.set_ylim([0, None])
g.ax_marg_x.set_xscale('symlog')
g.ax_marg_x.set_xlim([0, None])
g.ax_marg_y.set_yscale('symlog')
g.ax_marg_y.set_ylim([0, None])

pp.savefig()

pp.close()

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
