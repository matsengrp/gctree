#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
aggregation plots of validation output from several simulation/validation runs
'''

from __future__ import division, print_function
import scipy, matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
from matplotlib import rc, ticker
import pandas as pd
import argparse
import seaborn as sns

parser = argparse.ArgumentParser(description='aggregate validation of repeated runs with same parameters')
parser.add_argument('input', type=str, nargs='+', help='validation.tsv files')
parser.add_argument('--outbase', type=str, help='output file base name')
args = parser.parse_args()

aggdat = pd.DataFrame()
for i, fname in enumerate(args.input):
    df = pd.read_csv(fname, sep='\t')
    aggdat = aggdat.append(df)

aggdat.to_csv(args.outbase+'.tsv', sep='\t', index=False)

plt.figure()
df1 = pd.melt(aggdat.ix[:, aggdat.columns != 'N_taxa'], id_vars=['method'], var_name='metric')
sns.factorplot(x="method", y="value", col="metric", col_wrap=2,
                   data=df1, kind="swarm", size=2.5, aspect=.8, sharey=False)
plt.savefig(args.outbase+'.pdf')
