#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
aggregation plots of across parameters
'''

from __future__ import division, print_function
import scipy, matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import rc, ticker
import pandas as pd
import argparse
import seaborn as sns
from os import path, sep

parser = argparse.ArgumentParser(description='outer aggregation')
parser.add_argument('input', type=str, nargs='+', help='validaggreg.tsv files')
parser.add_argument('--outbase', type=str, help='output file base name')
args = parser.parse_args()

aggdat = pd.DataFrame()

for i, fname in enumerate(args.input):
    df = pd.read_csv(fname, sep='\t')
    lambda_, lambda0, r = map(float, path.dirname(fname).split(sep)[-3:])
    df['lambda'] = lambda_
    df['lambda0'] = lambda0
    df['r'] = r
    aggdat = aggdat.append(df)
    # aggdat.loc[i] = (forest_size, df['mean_frequency'][0], df['MRCA'][0], sum(x <= df['MRCA'][0] for x in df['MRCA']))

aggdat.to_csv(args.outbase+'.tsv', sep='\t', index=False)

if len(set(aggdat['lambda'])) == 1 and len(set(aggdat['r'])) == 1:
    plt.rc('text', usetex=True)
    fig = plt.figure()
    ax = fig.add_subplot(3,1,1)
    sns.boxplot(x='lambda0', y='parsimony forest size', data=aggdat)
    ax.set_yscale('log')
    ax = fig.add_subplot(3,1,2)
    sns.boxplot(x='lambda0', y='mean allele frequency', data=aggdat)
    ax = fig.add_subplot(3,1,3)
    sns.boxplot(x='lambda0', y='RF distance to true tree', data=aggdat)
    # ax = fig.add_subplot(4,1,4)
    # sns.boxplot(x='lambda0', y='MRCA distance', data=aggdat)

    ax.set_xlabel(r'$\lambda_0$')
    plt.savefig(args.outbase+'.pdf')
