#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
comparison of heavy and light trees
'''

from __future__ import division, print_function
from gctree import CollapsedTree, CollapsedForest
import pandas as pd
import scipy
import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
import seaborn as sns
sns.set(style='white', color_codes=True)
import pickle

def main():

    import argparse

    parser = argparse.ArgumentParser(description='compare heavy and light chain trees',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('heavy_forest', type=str, help='.p forest file')
    parser.add_argument('light_forest', type=str, help='.p forest file')
    # parser.add_argument('--outbase', type=str, required=True, help='output file base name')
    args = parser.parse_args()

    with open(args.heavy_forest, 'rb') as f:
        heavy_forest = pickle.load(f)
    with open(args.light_forest, 'rb') as f:
        light_forest = pickle.load(f)

    for node in heavy_forest.forest[0].tree.traverse():
        print(node.name)

    for node in light_forest.forest[0].tree.traverse():
        print(node.name)

if __name__ == '__main__':
    main()
