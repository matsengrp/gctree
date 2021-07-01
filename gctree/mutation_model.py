#! /usr/bin/env python

"""Mutation models."""

from gctree.utils import hamming_distance

from ete3 import TreeNode
from scipy.stats import poisson
import random
import scipy
from Bio.Seq import Seq


class MutationModel:
    """a class for a mutation model, and functions to mutate sequences."""

    def __init__(
        self,
        mutability_file=None,
        substitution_file=None,
        mutation_order=True,
        with_replacement=True,
    ):
        """initialized with input files of the S5F format.

        @param mutation_order: whether or not to mutate sequences using a context sensitive manner
                                where mutation order matters
        @param with_replacement: allow the same position to mutate multiple times on a single branch
        """
        self.mutation_order = mutation_order
        self.with_replacement = with_replacement
        if mutability_file is not None and substitution_file is not None:
            self.context_model = {}
            with open(mutability_file, "r") as f:
                # eat header
                f.readline()
                for line in f:
                    motif, score = line.replace('"', "").split()[:2]
                    self.context_model[motif] = float(score)

            # kmer k
            self.k = None
            with open(substitution_file, "r") as f:
                # eat header
                f.readline()
                for line in f:
                    fields = line.replace('"', "").split()
                    motif = fields[0]
                    if self.k is None:
                        self.k = len(motif)
                        assert self.k % 2 == 1
                    else:
                        assert len(motif) == self.k
                    self.context_model[motif] = (
                        self.context_model[motif],
                        {b: float(x) for b, x in zip("ACGT", fields[1:5])},
                    )
        else:
            self.context_model = None

    @staticmethod
    def disambiguate(sequence):
        """generator of all possible nt sequences implied by a sequence
        containing Ns."""
        # find the first N nucleotide
        N_index = sequence.find("N")
        # if there is no N nucleotide, yield the input sequence
        if N_index == -1:
            yield sequence
        else:
            for n_replace in "ACGT":
                # ooooh, recursion
                # NOTE: in python3 we could simply use "yield from..." instead of this loop
                for sequence_recurse in MutationModel.disambiguate(
                    sequence[:N_index] + n_replace + sequence[N_index + 1 :]
                ):
                    yield sequence_recurse

    def mutability(self, kmer):
        """returns the mutability of a central base of kmer, along with
        nucleotide bias averages over N nucleotide identities."""
        if self.context_model is None:
            raise ValueError("kmer mutability only defined for context models")
        if len(kmer) != self.k:
            raise ValueError(
                "kmer of length {} inconsistent with context model kmer length {}".format(
                    len(kmer), self.k
                )
            )
        if not all(n in "ACGTN" for n in kmer):
            raise ValueError(
                "sequence {} must contain only characters A, C, G, T, or N".format(kmer)
            )

        mutabilities_to_average, substitutions_to_average = zip(
            *[self.context_model[x] for x in MutationModel.disambiguate(kmer)]
        )

        average_mutability = scipy.mean(mutabilities_to_average)
        average_substitution = {
            b: sum(
                substitution_dict[b] for substitution_dict in substitutions_to_average
            )
            / len(substitutions_to_average)
            for b in "ACGT"
        }

        return average_mutability, average_substitution

    def mutabilities(self, sequence):
        """returns the mutability of a sequence at each site, along with
        nucleotide biases."""
        if self.context_model is None:
            return [
                (1, dict((n2, 1 / 3) if n2 is not n else (n2, 0.0) for n2 in "ACGT"))
                for n in sequence
            ]
        else:
            # pad with Ns to allow averaged edge effects
            sequence = "N" * (self.k // 2) + sequence + "N" * (self.k // 2)
            # mutabilities of each nucleotide
            return [
                self.mutability(sequence[(i - self.k // 2) : (i + self.k // 2 + 1)])
                for i in range(self.k // 2, len(sequence) - self.k // 2)
            ]

    def mutate(self, sequence, lambda0=1, frame=None):
        """Mutate a sequence, with lamdba0 the baseline mutability Cannot
        mutate the same position multiple times.

        @param sequence: the original sequence to mutate
        @param lambda0: a "baseline" mutation rate
        @param frame: the reading frame index
        """
        sequence_length = len(sequence)
        if frame is not None:
            codon_start = frame - 1
            codon_end = codon_start + 3 * ((sequence_length - codon_start) // 3)
            if "*" in Seq(sequence[codon_start:codon_end]).translate():
                raise RuntimeError("sequence contains stop codon!")

        mutabilities = self.mutabilities(sequence)
        sequence_mutability = (
            sum(mutability[0] for mutability in mutabilities) / sequence_length
        )
        # poisson rate for this sequence (given its relative mutability)
        lambda_sequence = sequence_mutability * lambda0
        # number of mutations m
        trials = 20
        for trial in range(1, trials + 1):
            m = scipy.random.poisson(lambda_sequence)
            if m <= sequence_length or self.with_replacement:
                break
            if trial == trials:
                raise RuntimeError("mutations saturating, consider reducing lambda0")

        # mutate the sites with mutations
        # if frame is not None and sequence contains stop codon, try again, up to 10 times
        unmutated_positions = range(sequence_length)
        for i in range(m):
            sequence_list = list(sequence)  # make string a list so we can modify it
            # Determine the position to mutate from the mutability matrix
            mutability_p = scipy.array(
                [mutabilities[pos][0] for pos in unmutated_positions]
            )
            for trial in range(1, trials + 1):
                mut_pos = scipy.random.choice(
                    unmutated_positions, p=mutability_p / mutability_p.sum()
                )
                # Now draw the target nucleotide using the substitution matrix
                substitution_p = [mutabilities[mut_pos][1][n] for n in "ACGT"]
                assert 0 <= abs(sum(substitution_p) - 1.0) < 1e-10
                chosen_target = scipy.random.choice(4, p=substitution_p)
                original_base = sequence_list[mut_pos]
                sequence_list[mut_pos] = "ACGT"[chosen_target]
                sequence = "".join(sequence_list)  # reconstruct our sequence
                if (
                    frame is None
                    or "*" not in Seq(sequence[codon_start:codon_end]).translate()
                ):
                    if self.mutation_order:
                        # if mutation order matters, the mutabilities of the sequence need to be updated
                        mutabilities = self.mutabilities(sequence)
                    if not self.with_replacement:
                        # Remove this position so we don't mutate it again
                        unmutated_positions.remove(mut_pos)
                    break
                if trial == trials:
                    raise RuntimeError(
                        "stop codon in simulated sequence on "
                        + str(trials)
                        + " consecutive attempts"
                    )
                sequence_list[
                    mut_pos
                ] = original_base  # <-- we only get here if we are retrying

        return sequence

    def one_mutant(self, sequence, Nmuts, frame=1, lambda0=0.1):
        """Make a single mutant with a distance, in amino acid sequence, of
        Nmuts away from the starting point."""
        trial = 100  # Allow 100 trials before quitting
        while trial > 0:
            mut_seq = sequence[:]
            aa = str(
                Seq(
                    sequence[
                        (frame - 1) : (
                            frame - 1 + (3 * (((len(sequence) - (frame - 1)) // 3)))
                        )
                    ]
                ).translate()
            )
            aa_mut = Seq(
                mut_seq[
                    (frame - 1) : (
                        frame - 1 + (3 * (((len(mut_seq) - (frame - 1)) // 3)))
                    )
                ]
            ).translate()
            dist = hamming_distance(aa, aa_mut)
            while dist < Nmuts:
                mut_seq = self.mutate(mut_seq, lambda0=lambda0, frame=frame)
                aa_mut = str(
                    Seq(
                        mut_seq[
                            (frame - 1) : (
                                frame - 1 + (3 * (((len(mut_seq) - (frame - 1)) // 3)))
                            )
                        ]
                    ).translate()
                )
                dist = hamming_distance(aa, aa_mut)
            if dist == Nmuts:
                return aa_mut
            else:
                trial -= 1
        raise RuntimeError(
            "100 consecutive attempts for creating a target sequence failed."
        )

    def simulate(
        self,
        sequence,
        seq_bounds=None,
        progeny=poisson(0.9),
        lambda0=[1],
        frame=None,
        N=None,
        T=None,
        n=None,
        verbose=False,
    ):
        """simulate neutral binary branching process with mutation model
        progeny must be like a scipy.stats distribution, with rvs() and mean()
        methods."""
        # Checking the validity of the input parameters:
        if N is not None and T is not None:
            raise ValueError("Only one of N and T can be used. One must be None.")
        elif N is None and T is None:
            raise ValueError("Either N or T must be specified.")
        if N is not None and n is not None and n > N:
            raise ValueError("n ({}) must not larger than N ({})".format(n, N))

        # Planting the tree:
        tree = TreeNode()
        tree.dist = 0
        tree.add_feature("sequence", sequence)
        tree.add_feature("terminated", False)
        tree.add_feature("frequency", 0)
        tree.add_feature("time", 0)

        t = 0  # <-- time
        leaves_unterminated = 1
        while (
            leaves_unterminated > 0
            and (leaves_unterminated < N if N is not None else True)
            and (t < max(T) if T is not None else True)
        ):
            if verbose:
                print("At time:", t)
            t += 1
            list_of_leaves = list(tree.iter_leaves())
            random.shuffle(list_of_leaves)
            for leaf in list_of_leaves:
                if not leaf.terminated:
                    n_children = progeny.rvs()
                    leaves_unterminated += (
                        n_children - 1
                    )  # <-- this kills the parent if we drew a zero
                    if not n_children:
                        leaf.terminated = True
                    for child_count in range(n_children):
                        # If sequence pair mutate them separately with their own mutation rate:
                        if seq_bounds is not None:
                            mutated_sequence1 = self.mutate(
                                leaf.sequence[seq_bounds[0][0] : seq_bounds[0][1]],
                                lambda0=lambda0[0],
                                frame=frame,
                            )
                            mutated_sequence2 = self.mutate(
                                leaf.sequence[seq_bounds[1][0] : seq_bounds[1][1]],
                                lambda0=lambda0[1],
                                frame=frame,
                            )
                            mutated_sequence = mutated_sequence1 + mutated_sequence2
                        else:
                            mutated_sequence = self.mutate(
                                leaf.sequence, lambda0=lambda0[0], frame=frame
                            )
                        child = TreeNode()
                        child.dist = sum(
                            x != y for x, y in zip(mutated_sequence, leaf.sequence)
                        )
                        child.add_feature("sequence", mutated_sequence)
                        child.add_feature("frequency", 0)
                        child.add_feature("terminated", False)
                        child.add_feature("time", t)
                        leaf.add_child(child)

        if leaves_unterminated < N:
            raise RuntimeError(
                "tree terminated with {} leaves, {} desired".format(
                    leaves_unterminated, N
                )
            )

        # each leaf in final generation gets an observation frequency of 1, unless downsampled
        if T is not None and len(T) > 1:
            # Iterate the intermediate time steps:
            for Ti in sorted(T)[:-1]:
                # Only sample those that have been 'sampled' at intermediate sampling times:
                final_leaves = [
                    leaf
                    for leaf in tree.iter_descendants()
                    if leaf.time == Ti and leaf.sampled
                ]
                if len(final_leaves) < n:
                    raise RuntimeError(
                        "tree terminated with {} leaves, less than what desired after downsampling {}".format(
                            leaves_unterminated, n
                        )
                    )
                for (
                    leaf
                ) in (
                    final_leaves
                ):  # No need to down-sample, this was already done in the simulation loop
                    leaf.frequency = 1
        # Do the normal sampling of the last time step:
        final_leaves = [leaf for leaf in tree.iter_leaves() if leaf.time == t]
        # by default, downsample to the target simulation size
        if n is not None and len(final_leaves) >= n:
            for leaf in random.sample(final_leaves, n):
                leaf.frequency = 1
        elif n is None and N is not None:
            for leaf in random.sample(final_leaves, N):
                leaf.frequency = 1
        elif N is None and T is not None:
            for leaf in final_leaves:
                leaf.frequency = 1
        elif n is not None and len(final_leaves) < n:
            raise RuntimeError(
                "tree terminated with {} leaves, less than what desired after downsampling {}".format(
                    leaves_unterminated, n
                )
            )
        else:
            raise RuntimeError("Unknown option.")

        # prune away lineages that are unobserved
        for node in tree.iter_descendants():
            if sum(node2.frequency for node2 in node.traverse()) == 0:
                node.detach()

        # # remove unobserved unifurcations
        # for node in tree.iter_descendants():
        #     parent = node.up
        #     if node.frequency == 0 and len(node.children) == 1:
        #         node.delete(prevent_nondicotomic=False)
        #         node.children[0].dist = hamming_distance(node.children[0].sequence, parent.sequence)

        # assign unique names to each node
        for i, node in enumerate(tree.traverse(), 1):
            node.name = "simcell_{}".format(i)

        # return the uncollapsed tree
        return tree
