#!/bin/env python

import argparse, scipy, gctree
from ete3 import TreeNode
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna

class MutationModel():
    """a class for a mutation model, and functions to mutate sequences"""
    def __init__(self, mutability_file, substitution_file):
        """initialized with input files of the S5F format"""
        self._mutation_model = {}
        with open(mutability_file, 'r') as f:
            # eat header
            f.readline()
            for line in f:
                motif, score = line.replace('"', '').split()[:2]
                self._mutation_model[motif] = float(score)

        # kmer k
        self.k = None
        with open(substitution_file, 'r') as f:
            # eat header
            f.readline()
            for line in f:
                fields = line.replace('"', '').split()
                motif = fields[0]
                if self.k is None:
                    self.k = len(motif)
                    assert self.k % 2 == 1
                else:
                    assert len(motif) == self.k
                self._mutation_model[motif] = (self._mutation_model[motif], {b:float(x) for b, x in zip('ACGT', fields[1:5])})

    def mutability(self, kmer):
        """"returns the mutability of a kmer, along with nucleotide biases"""
        assert len(kmer) == self.k
        return self._mutation_model[kmer]

    def mutate(self, sequence, lambda0=1):
        """mutate a sequence, with q the baseline mutability"""
        assert all(n in 'ACGT' for n in sequence)
        # mutabilities of each nucleotide
        mutabilities = []
        sequence_length = len(sequence)
        assert sequence_length >= 5
        # ambiguous left end motifs
        for i in range(self.k//2 + 1):
            kmer_suffix = sequence[:(i+self.k//2+1)]
            matches = [value for key, value in self._mutation_model.iteritems() if key.endswith(kmer_suffix)]
            len_matches = len(matches)
            assert len_matches == 4**(self.k - len(kmer_suffix))
            # use mean over matches
            mutability = sum(match[0] for match in matches)/float(len_matches)
            substitution = {n:sum(d[1][n] for d in matches)/float(len_matches) for n in 'ACGT'}
            mutabilities.append((mutability, substitution))
        # unambiguous internal kmers
        for i in range(self.k//2, sequence_length - self.k//2):
            mutabilities.append(self.mutability(sequence[(i-self.k//2):(i+self.k//2+1)]))
        # ambiguous right end motifs
        for i in range(sequence_length - self.k//2 + 1, sequence_length):
            kmer_prefix = sequence[(i-self.k//2):]
            matches = [value for key, value in self._mutation_model.iteritems() if key.startswith(kmer_prefix)]
            len_matches = len(matches)
            assert len_matches == 4**(self.k - len(kmer_prefix))
            # use mean over matches
            mutability = sum(match[0] for match in matches)/float(len_matches)
            substitution = {n:sum(d[1][n] for d in matches)/float(len_matches) for n in 'ACGT'}
            mutabilities.append((mutability, substitution))

        assert len(mutabilities) == sequence_length

        # mean mutability
        sequence_mutability = sum(mutability[0] for mutability in mutabilities)/float(sequence_length)
        # baseline Piosson 
        #lambda_0 = -scipy.log(1-q)
        # poisson rate for this sequence (given its relative mutability)
        lambda_sequence = sequence_mutability*lambda0
        # number of mutations 
        m = scipy.random.poisson(lambda_sequence)

        if m > 0:
            # now we choose random positions for the m mutations, weighted by mutabilities
            # invoking a long sequence limit, we don't allow back mutations
            # draw a multinomial rv for the number of mutations in each site 
            p = [mutability[0]/(sequence_length*sequence_mutability) for mutability in mutabilities]
            assert 0 <= abs(sum(p) - 1.) < 1e-10
            mutated_sites = scipy.random.multinomial(m, p)
            trial = 0
            while max(mutated_sites) > 1:
                print 'repeated mutations, trying again'
                trial += 1
                if trial > 5:
                    raise RuntimeError('mutations saturating')
                mutated_sites = scipy.random.multinomial(m, p)
            sequence = list(sequence) # mutable
            for i in range(sequence_length):
                if mutated_sites[i]:
                    p = [mutabilities[i][1][n] for n in 'ACGT']
                    assert 0 <= abs(sum(p) - 1.) < 1e-10
                    sequence[i] = 'ACGT'[scipy.nonzero(scipy.random.multinomial(1, p))[0][0]]
            sequence = "".join(sequence)

        return sequence


    def simulate(self, sequence, outbase, p=.4, lambda0=1, r=1.):
        """"simulate neutral binary branching process with mutation model"""
        if p >= .5:
            raw_input('WARNING: p = %f is not subcritical, tree termination not garanteed! [ENTER] to proceed')
        self.tree = TreeNode(dist=0)
        self.tree.add_feature('sequence', sequence)
        self.tree.add_feature('terminated', False)
        self.tree.add_feature('frequency', 0)
        nodes_unterminated = 1
        while nodes_unterminated > 0:
            for leaf in self.tree.iter_leaves():
                if not leaf.terminated:
                    if scipy.random.random() < p:
                        for child_count in range(2):
                            mutated_sequence = self.mutate(leaf.sequence, lambda0=lambda0)
                            child = TreeNode(dist=sum(x!=y for x,y in zip(mutated_sequence, leaf.sequence)))
                            child.add_feature('sequence', mutated_sequence)
                            child.add_feature('frequency', 0)
                            leaf.add_child(child)
                            child.add_feature('terminated' ,False)
                        nodes_unterminated += 1
                    else:
                        leaf.terminated = True
                        nodes_unterminated -= 1

        # each leaf gets an observation frequency of 1
        for node in self.tree.iter_leaves():
            if scipy.random.random() < r:
                node.frequency = 1

        with open(outbase+'.leafdata.fa', 'w') as f:
            f.write('> GL\n')
            f.write(sequence+'\n')
            i = 0
            for leaf in self.tree.iter_leaves():
                if leaf.frequency != 0:# and '*' not in Seq(leaf.sequence, generic_dna).translate():
                    i += 1
                    f.write('> seq%d\n' % i)
                    f.write(leaf.sequence+'\n')
                    leaf.name = 'seq%d' % i
        print i, 'simulated observed sequences'
        #self.tree.link_to_alignment(alignment=outbase+'.leafdata.fa', alg_format='fasta')
        self.tree.render(outbase+'.tree.png')


        # get collapsed tree
        self.collapsed_tree = gctree.CollapsedTree(tree=self.tree)
        self.collapsed_tree.render(outbase+'.collapsed_tree.png')

        return self


def main():
    parser = argparse.ArgumentParser(description='neutral model gctree simulation')
    parser.add_argument('sequence', type=str, help='seed germline nucleotide sequence')
    parser.add_argument('mutability', type=str, help='mutability model')
    parser.add_argument('substitution', type=str, help='substitution model')
    parser.add_argument('outbase', type=str, help='base name for output')
    parser.add_argument('--p', type=float, default=.4, help='branching probability')
    parser.add_argument('--lambda0', type=float, default=None, help='baseline mutation rate')
    parser.add_argument('--r', type=float, default=1., help='sampling probability')
    args = parser.parse_args()

    if args.lambda0 is None:
        args.lambda0 = max([1, int(.01*len(sequence))])

    args.sequence = args.sequence.upper()

    mutation_model = MutationModel(args.mutability, args.substitution)
    mutation_model.simulate(args.sequence, args.outbase, p=args.p, lambda0=args.lambda0, r=args.r)

if __name__ == "__main__":
    main()
