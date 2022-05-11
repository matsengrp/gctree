r"""Utility functions."""
from functools import wraps, reduce
import Bio.Data.IUPACData
import operator
from typing import Sequence, Any

Multiplicable = Any

bases = "AGCT-"
ambiguous_dna_values = Bio.Data.IUPACData.ambiguous_dna_values.copy()
ambiguous_dna_values.update({"?": "GATC-", "-": "-"})
ambiguous_dna_keys = {frozenset(val): key for key, val in ambiguous_dna_values.items()}


def _check_distance_arguments(distance):
    @wraps(distance)
    def new_distance(seq1: str, seq2: str, *args, **kwargs):
        if len(seq1) != len(seq2):
            raise ValueError(
                f"sequences must have equal length, got {len(seq1)} and {len(seq2)}"
            )
        return distance(seq1, seq2, *args, **kwargs)

    return new_distance


@_check_distance_arguments
def hamming_distance(seq1: str, seq2: str) -> int:
    r"""Hamming distance between two sequences of equal length.

    Args:
        seq1: sequence 1
        seq2: sequence 2
    """
    return sum(x != y for x, y in zip(seq1, seq2))


def product(
    factors: Sequence[Multiplicable], f=operator.mul, identity=1
) -> Multiplicable:
    return reduce(f, factors, identity)
