#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""Utility functions."""


def hamming_distance(seq1, seq2):
    """Hamming distance between two sequences of equal length."""

    if len(seq1) != len(seq2):
        raise ValueError(
            f"sequences must have equal length, got {len(seq1)} and {len(seq2)}"
        )

    return sum(x != y for x, y in zip(seq1, seq2))
