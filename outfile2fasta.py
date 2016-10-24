#! /bin/env python

import argparse, gctree

def main():
    parser = argparse.ArgumentParser(description='parse ancestral states from dnaml')
    parser.add_argument('phylipfile', type=str, help='mlpars outfile (verbose output with sequences at each site)')
    args = parser.parse_args()

    # parse phylip outfile
    outfiledat = [block.split('\n\n\n')[0].split('\n\n')[1:] for block in open(args.phylipfile, 'r').read().split('Probable sequences at interior nodes:\n')[1:]]

    for i, tree in enumerate(outfiledat):
        tree_sequence_dict = {}
        for j, block in enumerate(tree):
            if j == 0:
                for line in block.split('\n'):
                    fields = line.split()
                    if len(fields) == 0:
                        continue
                    name = fields[0]
                    seq = ''.join(fields[1:])
                    tree_sequence_dict[name] = seq
            else:
                for line in block.split('\n'):
                    fields = line.split()
                    if len(fields) == 0: continue
                    name = fields[0]
                    seq = ''.join(fields[1:])
                    tree_sequence_dict[name] += seq

        for name in tree_sequence_dict:
            print '> '+name
            print tree_sequence_dict[name]

if __name__ == "__main__":
    main()

