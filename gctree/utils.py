r"""Utility functions."""


def hamming_distance(seq1: str, seq2: str) -> int:
    r"""Hamming distance between two sequences of equal length.

    Args:
        seq1: sequence 1
        seq2: sequence 2
    """

    if len(seq1) != len(seq2):
        raise ValueError(
            f"sequences must have equal length, got {len(seq1)} and {len(seq2)}"
        )

    return sum(x != y for x, y in zip(seq1, seq2))
