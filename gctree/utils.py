from math import log
from Bio.Data.IUPACData import ambiguous_dna_values

r"""Utility functions."""

bases = "AGCT-"
ambiguous_dna_values.update({"?": "GATC-", "-": "-"})

def disambiguations(sequence, accum=""):
    """Iterates through possible disambiguations of sequence, recursively.
    Recursion-depth-limited by number of ambiguity codes in
    sequence, not sequence length.
    """
    if sequence:
        for index, base in enumerate(sequence):
            if base in bases:
                accum += base
            else:
                for newbase in ambiguous_dna_values[base]:
                    yield from disambiguations(
                        sequence[index + 1 :], accum=(accum + newbase)
                    )
                return
    yield accum

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


def mutability_distance(mutation_model):
    """Returns a fast distance function based on mutability_model.
    First, caches computed mutabilities for k-mers with k // 2 N's on either
    end. This is pretty fast for k=5, but the distance function should be created
    once and reused."""
    # Caching could be moved to the MutationModel class instead.
    context_model = mutation_model.context_model.copy()
    k = mutation_model.k
    h = k // 2
    # Build all sequences with (when k=5) one or two Ns on either end
    templates = [("N"* left, "N" * (k - left - right), "N" * right)
                 for left in range(h + 1)
                 for right in range(h + 1)
                 if left != 0 or right != 0]

    kmers_to_compute = [
        leftns + stub + rightns
        for leftns, ambig_stub, rightns in templates
        for stub in disambiguations(ambig_stub)
    ]
    # Cache all these mutabilities in context_model also
    context_model.update({kmer: mutation_model.mutability(kmer) for kmer in kmers_to_compute})

    def mutabilities(seq):
        newseq = "N" * h + seq + "N" * h
        return [
            context_model[newseq[i - h : i + h + 1]] for i in range(h, len(seq) + h)
        ]

    @check_distance_arguments
    def distance(seq1: str, seq2: str) -> float:
        """Assume that sequences being compared are already padded, so this will be
        a sum of negative log biases over bases within the padding margin...but
        that will be a problem when there's an ambiguity at the beginning of the
        sequence."""
        muts = mutabilities(seq1)
        nc = sum(mut for mut, _ in muts)
        return sum(
            -log((mut[0] / nc) * mut[1][seq2[index]])
            if seq1[index] != seq2[index]
            else -log(1 - (mut[0] / nc))
            for index, mut in enumerate(muts)
        )

    return distance
