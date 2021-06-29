#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gctree.phylip_parse as pp
import gctree.branching_processes as bp

def main():

    import pickle, argparse, os

    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        return fname
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'phylip_outfile', type=existing_file, help='dnaml outfile (verbose output with inferred ancestral sequences, option 5).')
    parser.add_argument(
        'countfile', type=existing_file, help="count file")
    parser.add_argument(
        '--outputfile', default=None, help="output file.")
    parser.add_argument(
        '--naive', default='naive', help="naive sequence id")

    args = parser.parse_args()

    if args.outputfile is None:
        args.outputfile = args.phylip_outfile + '.collapsed_forest.p'
    trees = pp.parse_outfile(args.phylip_outfile, args.countfile, args.naive)
    if isinstance(trees[0], list):
        print(trees[0][0])
        print(bp.CollapsedTree(tree=trees[0][0]))
        bootstraps = [[bp.CollapsedTree(tree=tree) for tree in bootstrap] for bootstrap in trees]
        pickle.dump([bp.CollapsedForest(forest=trees) for trees in bootstraps], open(args.outputfile, 'w'))
    else:
        trees = [bp.CollapsedTree(tree=tree) for tree in trees]
        pickle.dump(bp.CollapsedForest(forest=trees), open(args.outputfile, 'w'))

if __name__ == "__main__":
    main()
