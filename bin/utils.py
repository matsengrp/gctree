#! /usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import jellyfish
    def hamming_distance(s1, s2):
        if s1 == s2:
            return 0
        else:
            return jellyfish.hamming_distance(unicode(s1), unicode(s2))
except:
    def hamming_distance(seq1, seq2):
        '''Hamming distance between two sequences of equal length'''
        return sum(x != y for x, y in zip(seq1, seq2))
    print('Couldn\'t find the python module "jellyfish" which is used for fast string comparison. Falling back to pure python function.')
