r"""Utility functions."""


def check_distance_arguments(distance):

    def new_distance(seq1: str, seq2: str, *args, **kwargs):
        if len(seq1) != len(seq2):
            raise ValueError(
                f"sequences must have equal length, got {len(seq1)} and {len(seq2)}"
            )
        return distance(seq1, seq2, *args, **kwargs)

    return new_distance

@check_distance_arguments
def hamming_distance(seq1: str, seq2: str) -> int:
    r"""Hamming distance between two sequences of equal length.

    Args:
        seq1: sequence 1
        seq2: sequence 2
    """
    return sum(x != y for x, y in zip(seq1, seq2))

@check_distance_arguments
def mutability_distance(seq1: str, seq2: str, mutability_model) -> float:
    """Assume that sequences being compared are already padded, so this will be a sum of
    negative log biases over bases within the padding margin...but that will be a problem
    when there's an ambiguity at the beginning of the sequence"""
    mutabilities = mutability_model.mutabilities(seq1)
    normalizing_constant = sum(mut for mut, _ in mutabilities)
    normalized_mutabilities = [(mut / normalizing_constant, biases) for mut, biases in mutabilities]
    transition_costs = [-log(mut[0] * mut[1][seq2[index]]) if seq1[index] == seq2[index]
                        else -log(1-mut[0])
                        for index, mutability in enumerate(mutabilities)]
    return sum(transition_costs)
