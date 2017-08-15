#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
aggregation plots of across parameters
'''

from __future__ import division, print_function
import scipy, matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import argparse
import seaborn as sns
from os import path, sep
import numpy as np

parser = argparse.ArgumentParser(description='outer aggregation')
parser.add_argument('input', type=str, nargs='+', help='validaggreg.tsv files')
parser.add_argument('--outbase', type=str, help='output file base name')
args = parser.parse_args()

aggdat = pd.DataFrame()

for i, fname in enumerate(args.input):
    df = pd.read_csv(fname, sep='\t')
    lambda_, lambda0 = map(float, path.dirname(fname).split(sep)[-2:])
    df['lambda'] = lambda_
    df['lambda0'] = lambda0
    aggdat = aggdat.append(df)

aggdat.to_csv(args.outbase+'.tsv', sep='\t', index=False)
columns = aggdat.columns.values.tolist()
# columns.remove('lambda0')

df = pd.melt(aggdat, id_vars=['method', 'N_taxa', 'lambda', 'lambda0'], value_vars=['RF','MRCA'], var_name='metric')
df.ix[df.ix[:,'metric'] == 'MRCA', 'value'] = np.log(1 + df.ix[df.ix[:,'metric'] == 'MRCA', 'value'])
df.ix[df.ix[:,'metric'] == 'MRCA', 'metric'] = 'logMRCA'


plot2var = [sns.swarmplot, sns.swarmplot, sns.lmplot]
variables = ['lambda', 'lambda0', 'N_taxa']
options = [dict(split=True, size=3), dict(split=True, size=3), dict(size=6)]
#options = [dict(size=3), dict(size=3), dict(size=3), dict(size=12)]


#options = [dict(split=True, size=3), dict(split=True, size=3), dict(split=True, size=3), dict()]
numb_variables = sum(len(set(df[var])) > 1 for var in variables) - 1  ### NOTICE -1 because of N_taxa cannot be plotted on same page because of "tight layout"

plot_func = {v:f for v, f in zip(variables, plot2var)}
plot_options = {v:o for v, o in zip(variables, options)}


gs = gridspec.GridSpec(numb_variables, 2)
fig = plt.figure(figsize=(8, 4*numb_variables))


i = 0
for var in variables:
    if not len(set(df[var])) > 1 or var == 'N_taxa':
        continue

    plf = plot_func[var]
    kwargs = plot_options[var]
    ax = plt.subplot(gs[i, 0])
    plot_data = df.ix[df.ix[:,'metric'] == 'RF', :]
    sns.boxplot(x=var, y="value", hue="method", data=plot_data, showfliers=False)
    plf(x=var, y="value", hue="method", data=plot_data, **kwargs)
    ax.set(ylabel='RF distance')

    ax = plt.subplot(gs[i, 1])
    plot_data = df.ix[df.ix[:,'metric'] == 'logMRCA', :]
    sns.boxplot(x=var, y="value", hue="method", data=plot_data, showfliers=False)
    plf(x=var, y="value", hue="method", data=plot_data, **kwargs)
    ax.set(ylabel='MRCA distance')

    i += 1


# fig.tight_layout()
plt.savefig(args.outbase+'.pdf')


var = 'N_taxa'
gs = gridspec.GridSpec(1, 1)
ax = plt.subplot(gs[0, 0])
fig = plt.figure(figsize=(8, 8))
plf = plot_func[var]
kwargs = plot_options[var]

plot_data = df.ix[df.ix[:,'metric'] == 'RF', :]
plf(x=var, y="value", hue="method", data=plot_data, legend_out=True, **kwargs)
#colors = {'gctree':'red', 'igphyml':'blue'}
#plot_data.plot(kind='scatter', x=var, y='value', c=plot_data['method'].apply(lambda x: colors[x]))
ax.set(ylabel='RF distance')
plt.title('Plot for tree size: {}. For RF.'.format(var))
#fig.tight_layout()
plt.savefig(args.outbase+'1.pdf')


plot_data = df.ix[df.ix[:,'metric'] == 'logMRCA', :]
plf(x=var, y="value", hue="method", data=plot_data, **kwargs)
ax.set(ylabel='MRCA distance')
plt.title('Plot for tree size: {}. For MRCA.'.format(var))
#fig.tight_layout()
plt.savefig(args.outbase+'2.pdf')






#if len(set(aggdat['lambda'])) == 1 and len(set(aggdat['r'])) == 1:
if False:

    gs = gridspec.GridSpec(5, 1, height_ratios=[1, 1, 2, 2, 2])
    gs2 = gridspec.GridSpec(5, 1, height_ratios=[1, 1, 2, 2, 2])
    gs.update(hspace=0.02)

    cut = 100
    fig = plt.figure(figsize=(10,10))
    ax = plt.subplot(gs[0])
    plot_data = aggdat.copy()
    plot_data.loc[plot_data['parsimony forest size'] <= cut, columns] = None
    sns.swarmplot(x='lambda0', y='parsimony forest size', data=plot_data, palette=sns.color_palette("GnBu_d"))
    ax.set_ylim(cut, None)
    ax.spines['bottom'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop='off')
    ax.get_xaxis().set_visible(False)
    ax.yaxis.set_label_coords(-.06,-.1)

    d = .005  # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False, linewidth=.5)
    ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
    ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal

    ax = plt.subplot(gs[1])
    plot_data = aggdat.copy()
    plot_data.loc[plot_data['parsimony forest size'] > cut, columns] = None
    sns.swarmplot(x='lambda0', y='parsimony forest size', data=plot_data, palette=sns.color_palette("GnBu_d"))
    ax.set_ylim(1, cut)
    ax.spines['top'].set_visible(False)
    ax.xaxis.tick_bottom()
    ax.get_xaxis().set_visible(False)
    ax.yaxis.label.set_visible(False)

    kwargs.update(transform=ax.transAxes)  # switch to the bottom axes
    ax.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    ax = plt.subplot(gs2[2])
    sns.swarmplot(x='lambda0', y='mean allele frequency', data=aggdat, palette=sns.color_palette("GnBu_d"))
    ax.get_xaxis().set_visible(False)

    ax = plt.subplot(gs2[3])
    sns.swarmplot(x='lambda0', y='mean branch length', data=aggdat, palette=sns.color_palette("GnBu_d"))
    ax.get_xaxis().set_visible(False)

    ax = plt.subplot(gs2[4])
    sns.swarmplot(x='lambda0', y='RF distance to true tree', data=aggdat, palette=sns.color_palette("GnBu_d"))
    # ax = fig.add_subplot(4,1,4)
    # sns.boxplot(x='lambda0', y='MRCA distance', data=aggdat)
    # ax.get_xaxis().set_visible(False)
    ax.set_xlabel('baseline mutation rate')
    plt.savefig(args.outbase+'.pdf')
