#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
This module contains classes for simulation and inference for a binary branching process with mutation
in which the tree is collapsed to nodes that count the number of clonal leaves of each type
'''

from __future__ import division, print_function
import scipy, warnings, random
try:
    import cPickle as pickle
except:
    import pickle
from scipy.misc import logsumexp
from scipy.optimize import minimize, check_grad, fsolve
from itertools import cycle
import random
import numpy as np
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib import rc, ticker
import pandas as pd
#import seaborn as sns; sns.set(style="white", color_codes=True)
from scipy.stats import poisson
from ete3 import TreeNode, NodeStyle, TreeStyle, TextFace, add_face_to_node, CircleFace, PieChartFace, faces, AttrFace, SVG_COLORS
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
# from Bio.Data.IUPACData import ambiguous_dna_values
from utils import SVG2HEX
import phylip_parse
try:
    # Last time I checked this module was the fastest kid on the block,
    # 5-10x faster than pure python, however 2x slower than a simple Cython function
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


scipy.seterr(all='raise')


class LeavesAndClades():
    '''
    This is a base class for simulating, and computing likelihood for, an infinite type branching
    process with branching probability p, mutation probability q, and we collapse mutant clades off the
    root type and consider just the number of clone leaves, c, and mutant clades, m.

      /\
     /\ ^          (3)
      /\     ==>   / \\
       /\\
        ^
    '''
    def __init__(self, params=None, c=None, m=None):
        '''initialize with branching probability p and mutation probability q, both in the unit interval'''
        if params is not None:
            p, q = params
            if not (0 <= p <= 1 and 0 <= q <= 1):
                raise ValueError('p and q must be in the unit interval')
        self._nparams = 2#len(params)
        self.params = params
        if c is not None or m is not None:
            if not (c >= 0) and (m >= 0) and (c+m > 0):
                raise ValueError('c and m must be nonnegative integers summing greater than zero')
            self.c = c
            self.m = m

    def simulate(self):
        '''simulate the number of clone leaves and mutant clades off a root node'''
        if self.params[0]>=.5:
            warnings.warn('p >= .5 is not subcritical, tree simulations not garanteed to terminate')
        if self.params is None:
            raise ValueError('params must be defined for simulation\n')

        # let's track the tree in breadth first order, listing number clone and mutant descendants of each node
        # mutant clades terminate in this view
        cumsum_clones = 0
        len_tree = 0
        self.c = 0
        self.m = 0
        # while termination condition not met
        while cumsum_clones > len_tree - 1:
            if random.random() < self.params[0]:
                mutants = sum(random.random() < self.params[1] for child in range(2))
                clones = 2 - mutants
                self.m += mutants
            else:
                mutants = 0
                clones = 0
                self.c += 1
            cumsum_clones += clones
            len_tree += 1
        assert cumsum_clones == len_tree - 1

    f_hash = {} # <--- class variable for hashing calls to the following function
    def f(self, params):
        '''
        Probability of getting c leaves that are clones of the root and m mutant clades off
        the root lineage, given branching probability p and mutation probability q
        Also returns gradient wrt (p, q)
        Computed by dynamic programming
        '''
        p, q = params
        c, m = self.c, self.m
        if (p, q, c, m) not in LeavesAndClades.f_hash:
            if c==m==0 or (c==0 and m==1):
                f_result = 0
                dfdp_result = 0
                dfdq_result = 0
            elif c==1 and m==0:
                f_result = 1-p
                dfdp_result = -1
                dfdq_result = 0
            elif c==0 and m==2:
                f_result = p*q**2
                dfdp_result = q**2
                dfdq_result = 2*p*q
            else:
                if m >= 1:
                    neighbor = LeavesAndClades(params=params, c=c, m=m-1)
                    neighbor_f, (neighbor_dfdp, neighbor_dfdq) = neighbor.f(params)
                    f_result = 2*p*q*(1-q)*neighbor_f
                    dfdp_result =   2*q*(1-q) * neighbor_f + \
                                  2*p*q*(1-q) * neighbor_dfdp
                    dfdq_result = (2*p - 4*p*q) * neighbor_f + \
                                   2*p*q*(1-q)  * neighbor_dfdq
                else:
                    f_result = 0.
                    dfdp_result = 0.
                    dfdq_result = 0.
                for cx in range(c+1):
                    for mx in range(m+1):
                        if (not (cx==0 and mx==0)) and (not (cx==c and mx==m)):
                            neighbor1 = LeavesAndClades(params=params, c=cx, m=mx)
                            neighbor2 = LeavesAndClades(params=params, c=c-cx, m=m-mx)
                            neighbor1_f, (neighbor1_dfdp, neighbor1_dfdq) = neighbor1.f(params)
                            neighbor2_f, (neighbor2_dfdp, neighbor2_dfdq) = neighbor2.f(params)
                            f_result += p*(1-q)**2*neighbor1_f*neighbor2_f
                            dfdp_result +=   (1-q)**2 * neighbor1_f    * neighbor2_f + \
                                           p*(1-q)**2 * neighbor1_dfdp * neighbor2_f + \
                                           p*(1-q)**2 * neighbor1_f    * neighbor2_dfdp
                            dfdq_result += -2*p*(1-q) * neighbor1_f    * neighbor2_f + \
                                           p*(1-q)**2 * neighbor1_dfdq * neighbor2_f + \
                                           p*(1-q)**2 * neighbor1_f    * neighbor2_dfdq
            LeavesAndClades.f_hash[(p, q, c, m)] = (f_result, scipy.array([dfdp_result, dfdq_result]))
        return LeavesAndClades.f_hash[(p, q, c, m)]


class CollapsedTree(LeavesAndClades):
    '''
    Here's a derived class for a collapsed tree, where we recurse into the mutant clades
          (4)
         / | \\
       (3)(1)(2)
           |   \\
          (2)  (1)
    '''
    def __init__(self, params=None, tree=None, frame=None, collapse_syn=False, allow_repeats=False):
        '''
        For intialization, either params or tree (or both) must be provided
        params: offspring distribution parameters
        tree: ete tree with frequency node feature. If uncollapsed, it will be collapsed
        frame: tranlation frame, with default None, no tranlation attempted
        '''
        LeavesAndClades.__init__(self, params=params)
        if frame is not None and frame not in (1, 2, 3):
            raise RuntimeError('frame must be 1, 2, 3, or None')
        self.frame = frame

        if collapse_syn is True:
            tree.dist = 0  # no branch above root
            for node in tree.iter_descendants():
                aa = Seq(node.sequence[(frame-1):(frame-1+(3*(((len(node.sequence)-(frame-1))//3))))],
                         generic_dna).translate()
                aa_parent = Seq(node.up.sequence[(frame-1):(frame-1+(3*(((len(node.sequence)-(frame-1))//3))))],
                                generic_dna).translate()
                node.dist = hamming_distance(aa, aa_parent)


        if tree is not None:
            self.tree = tree.copy()
            # iterate over the tree below root and collapse edges of zero length
            for node in self.tree.get_descendants():
                if node.dist == 0:
                    node.up.frequency += node.frequency
                    # if node.name and not node.up.name:
                    if node.up != self.tree:
                        node.up.name = node.name
                    node.delete(prevent_nondicotomic=False)

            assert sum(node.frequency for node in tree.traverse()) == sum(node.frequency for node in self.tree.traverse())
            if 'sequence' in tree.features:
                rep_seq = sum(node.frequency > 0 for node in self.tree.traverse()) - len(set([node.sequence for node in self.tree.traverse() if node.frequency > 0]))
                if not allow_repeats and rep_seq:
                    raise RuntimeError('Repeated observed sequences in collapsed tree. {} sequences were found repeated.'.format(rep_seq))
                elif allow_repeats and rep_seq:
                    rep_seq = sum(node.frequency > 0 for node in self.tree.traverse()) - len(set([node.sequence for node in self.tree.traverse() if node.frequency > 0]))
                    print('Repeated observed sequences in collapsed tree. {} sequences were found repeated.'.format(rep_seq))
            # now we do a custom ladderize accounting for abundance and sequence to break ties in abundance
            for node in self.tree.traverse(strategy='postorder'):
                # add a partition feature and compute it recursively up the tree
                node.add_feature('partition', node.frequency + sum(node2.partition for node2 in node.children))
                # sort children of this node based on partion and sequence
                node.children.sort(key=lambda node: (node.partition, node.sequence if 'sequence' in tree.features else None))
        else:
            self.tree = tree

    def l(self, params, sign=1):
        '''
        log likelihood of params, conditioned on collapsed tree, and its gradient wrt params
        optional parameter sign must be 1 or -1, with the latter useful for MLE by minimization
        '''
        if self.tree is None:
            raise ValueError('tree data must be defined to compute likelihood')
        if sign not in (-1, 1):
            raise ValueError('sign must be 1 or -1')
        leaves_and_clades_list = [LeavesAndClades(c=node.frequency, m=len(node.children)) for node in self.tree.traverse()]
        if leaves_and_clades_list[0].c == 0 and leaves_and_clades_list[0].m == 1 and leaves_and_clades_list[0].f(params)[0] == 0:
            # if unifurcation not possible under current model, add a psuedocount for the naive
            leaves_and_clades_list[0].c = 1
        # extract vector of function values and gradient components
        f_data = [leaves_and_clades.f(params) for leaves_and_clades in leaves_and_clades_list]
        fs = scipy.array([[x[0]] for x in f_data])
        logf = scipy.log(fs).sum()
        grad_fs = scipy.array([x[1] for x in f_data])
        grad_logf = (grad_fs/fs).sum(axis=0)
        return sign*logf, sign*grad_logf

    def mle(self, **kwargs):
        '''
        Maximum likelihood estimate for params given tree
        updates params if not None
        returns optimization result
        '''
        # random initalization
        x_0 = (random.random(), random.random())
        #x_0 = (.5, .5)
        bounds = ((.01, .99), (.001, .999))
        kwargs['sign'] = -1
        grad_check = check_grad(lambda x: self.l(x, **kwargs)[0], lambda x: self.l(x, **kwargs)[1], (.4, .5))
        if grad_check > 1e-3:
            warnings.warn('gradient mismatches finite difference approximation by {}'.format(grad_check), RuntimeWarning)
        result = minimize(lambda x: self.l(x, **kwargs), x0=x_0, jac=True, method='L-BFGS-B', options={'ftol':1e-10}, bounds=bounds)
        # update params if None and optimization successful
        if not result.success:
            warnings.warn('optimization not sucessful, '+result.message, RuntimeWarning)
        elif self.params is None:
            self.params = result.x.tolist()
        return result

    def simulate(self):
        '''
        simulate a collapsed tree given params
        replaces existing tree data member with simulation result, and returns self
        '''
        if self.params is None:
            raise ValueError('params must be defined for simulation')

        # initiate by running a LeavesAndClades simulation to get the number of clones and mutants
        # in the root node of the collapsed tree
        LeavesAndClades.simulate(self)
        self.tree = TreeNode()
        self.tree.add_feature('frequency', self.c)
        if self.m == 0:
            return self
        for _ in range(self.m):
            # ooooh, recursion
            child = CollapsedTree(params=self.params, frame=self.frame).simulate().tree
            child.dist = 1
            self.tree.add_child(child)

        return self


    def __str__(self):
        '''return a string representation for printing'''
        return 'params = ' + str(self.params)+ '\ntree:\n' + str(self.tree)

    def render(self, outfile, colormap=None, leftright_split=None):
        '''render to image file, filetype inferred from suffix, svg for color images'''
        def my_layout(node):
            #if node.frequency > 0:
                circle_color = 'lightgray' if colormap is None or node.name not in colormap else colormap[node.name]
                text_color = 'black'
                if isinstance(circle_color, str):
                    C = CircleFace(radius=max(3, 10*scipy.sqrt(node.frequency)), color=circle_color, label={'text':str(node.frequency), 'color':text_color} if node.frequency > 0 else None)
                    C.rotation = -90
                    C.hz_align = 1
                    faces.add_face_to_node(C, node, 0)
                else:
                    #for color in circle_color:
                    #    C = CircleFace(radius=max(.1, 10*scipy.sqrt(circle_color[color])), color=color, label={'text':str(circle_color[color]), 'color':text_color})
                    #    C.rotation = -90
                    #    faces.add_face_to_node(C, node, 0, position='float')
                    P = PieChartFace([100*x/node.frequency for x in circle_color.values()], 2*10*scipy.sqrt(node.frequency), 2*10*scipy.sqrt(node.frequency), colors=list(circle_color.keys()), line_color=None)
                    T = TextFace(' '.join([str(x) for x in list(circle_color.values())]))
                    T.hz_align = 1
                    T.vt_align = 1
                    T.rotation = -90
                    faces.add_face_to_node(P, node, 0, position='float')#0)
                    faces.add_face_to_node(T, node, 0, position='float')#0)
        for node in self.tree.traverse():
            nstyle = NodeStyle()
            #if node.frequency == 0:
            #    nstyle['size'] = 5
            #    nstyle['fgcolor'] = 'grey'
            #else:
            nstyle['size'] = 0
            #     nstyle['size'] = 3*2*scipy.sqrt(scipy.pi*node.frequency)
            #     nstyle['fgcolor'] = 'black'
            if node.up is not None:
                if set(node.sequence.upper()) == set('ACGT'):
                    if leftright_split is not None:
                        assert self.frame is None
                        if node.frequency > 0:
                            print(node.sequence[:leftright_split], node.sequence[leftright_split:])
                        leftseq_mutated = hamming_distance(node.sequence[:leftright_split], node.up.sequence[:leftright_split]) > 0
                        rightseq_mutated = hamming_distance(node.sequence[leftright_split:], node.up.sequence[leftright_split:]) > 0
                        if leftseq_mutated and rightseq_mutated:
                            nstyle['hz_line_color'] = 'purple'
                            nstyle['hz_line_width'] = 3
                        elif leftseq_mutated:
                            nstyle['hz_line_color'] = 'red'
                            nstyle['hz_line_width'] = 2
                        elif rightseq_mutated:
                            nstyle['hz_line_color'] = 'blue'
                            nstyle['hz_line_width'] = 2
                    if self.frame is not None:
                        aa = Seq(node.sequence[(self.frame-1):(self.frame-1+(3*(((len(node.sequence)-(self.frame-1))//3))))],
                                 generic_dna).translate()
                        aa_parent = Seq(node.up.sequence[(self.frame-1):(self.frame-1+(3*(((len(node.sequence)-(self.frame-1))//3))))],
                                        generic_dna).translate()
                        nonsyn = hamming_distance(aa, aa_parent)
                        if '*' in aa:
                            nstyle['bgcolor'] = 'red'
                        if nonsyn > 0:
                            nstyle['hz_line_color'] = 'black'
                            nstyle['hz_line_width'] = nonsyn
                        else:
                            nstyle['hz_line_type'] = 1
            node.set_style(nstyle)

        ts = TreeStyle()
        ts.show_leaf_name = False
        ts.rotation = 90
        ts.draw_aligned_faces_as_table = False
        ts.allow_face_overlap = True
        ts.layout_fn = my_layout
        ts.show_scale = False
        #self.tree.ladderize()
        self.tree.render(outfile, tree_style=ts)

    def write(self, file_name):
        '''serialize tree to file'''
        with open(file_name, 'wb') as f:
            pickle.dump(self, f)

    def compare(self, tree2, method='identity'):
        '''compare this tree to the other tree'''
        if method == 'identity':
            # we compare lists of seq, parent, abundance
            # return true if these lists are identical, else false
            list1 = sorted((node.sequence, node.frequency, node.up.sequence if node.up is not None else None) for node in self.tree.traverse())
            list2 = sorted((node.sequence, node.frequency, node.up.sequence if node.up is not None else None) for node in tree2.tree.traverse())
            return list1 == list2
        elif method == 'MRCA':
            # here's Erick's idea of matrix of hamming distance of common ancestors of taxa
            # takes a true and inferred tree as CollapsedTree objects
            taxa = [node.sequence for node in self.tree.traverse() if node.frequency]
            n_taxa = len(taxa)
            d = scipy.zeros(shape=(n_taxa, n_taxa))
            for i in range(n_taxa):
                nodei_true = self.tree.iter_search_nodes(sequence=taxa[i]).next()
                nodei      =      tree2.tree.iter_search_nodes(sequence=taxa[i]).next()
                for j in range(i + 1, n_taxa):
                    nodej_true = self.tree.iter_search_nodes(sequence=taxa[j]).next()
                    nodej      =      tree2.tree.iter_search_nodes(sequence=taxa[j]).next()
                    MRCA_true = self.tree.get_common_ancestor((nodei_true, nodej_true)).sequence
                    MRCA =           tree2.tree.get_common_ancestor((nodei, nodej)).sequence
                    d[i, j] = hamming_distance(MRCA_true, MRCA)
            return d.sum()
        elif method == 'RF':
            tree1_copy = self.tree.copy(method='deepcopy')
            tree2_copy = tree2.tree.copy(method='deepcopy')
            for treex in (tree1_copy, tree2_copy):
                for node in list(treex.traverse()):
                    # for _ in range(node.frequency):
                    if node.frequency > 0:
                        child = TreeNode()
                        child.add_feature('sequence', node.sequence)
                        node.add_child(child)
            return tree1_copy.robinson_foulds(tree2_copy, attr_t1='sequence', attr_t2='sequence', unrooted_trees=True)[0]
        else:
            raise ValueError('invalid distance method: '+method)


class CollapsedForest(CollapsedTree):
    '''
    simply a set of CollapsedTrees, with the same p and q parameters
          (4)          (3)
         / | \\         / \\
       (3)(1)(2)     (1) (2)
           |   \\  ,          , ...
          (2)  (1)
    '''
    def __init__(self, params=None, n_trees=None, forest=None):
        '''
        in addition to p and q, we need number of trees
        can also intialize with forest, a list of trees, each an instance of CollapsedTree
        '''
        CollapsedTree.__init__(self, params=params)
        if forest is None and params is None:
            raise ValueError('either params or forest (or both) must be provided')
        if forest is not None:
            if len(forest) == 0:
                raise ValueError('passed empty tree list')
            if n_trees is not None and len(forest) != n_trees:
                raise ValueError('n_trees not consistent with forest')
            self.forest = forest
        if n_trees is not None:
            if type(n_trees) is not int or n_trees < 1:
                raise ValueError('number of trees must be at least one')
            self.n_trees = n_trees
        if n_trees is None and forest is not None:
            self.n_trees = len(forest)

    def simulate(self):
        '''
        simulate a forest of collapsed trees given params and number of trees
        replaces existing forest data member with simulation result, and returns self
        '''
        if self.params is None or self.n_trees is None:
            raise ValueError('params and n_trees parameters must be defined for simulation')
        self.forest = [CollapsedTree(self.params).simulate() for x in range(self.n_trees)]
        return self

    def l(self, params, sign=1, Vlad_sum=False):
        '''
        likelihood of params, given forest, and it's gradient wrt params
        optional parameter sign must be 1 or -1, with the latter useful for MLE by minimization
        if optional parameter Vlad_sum is true, we're doing the Vlad sum for estimating params for
        as set of parsimony trees
        '''
        if self.forest is None:
            raise ValueError('forest data must be defined to compute likelihood')
        if sign not in (-1, 1):
            raise ValueError('sign must be 1 or -1')
        # since the l method on the CollapsedTree class returns l and grad_l...
        terms = [tree.l(params) for tree in self.forest]
        ls = scipy.array([term[0] for term in terms])
        grad_ls = scipy.array([term[1] for term in terms])
        if Vlad_sum:
            # we need to find the smallest derivative component for each
            # coordinate, then subtract off to get positive things to logsumexp
            grad_l = []
            for j in range(len(params)):
                i_prime = grad_ls[:,j].argmin()
                grad_l.append(grad_ls[i_prime,j] +
                              scipy.exp(logsumexp(ls - ls[i_prime],
                                                  b=grad_ls[:,j]-grad_ls[i_prime,j]) -
                                        logsumexp(ls - ls[i_prime])))
            return sign*(-scipy.log(len(ls)) + logsumexp(ls)), sign*scipy.array(grad_l)
        else:
            return sign*ls.sum(), sign*grad_ls.sum(axis=0)


    # NOTE: we get mle() method for free by inheritance/polymorphism magic


    def __str__(self):
        '''return a string representation for printing'''
        return 'params = {}, n_trees = {}\n'.format(self.params, self.n_trees) + \
                '\n'.join([str(tree) for tree in self.forest])


def disambiguate(tree):
    '''make random choices for ambiguous bases, respecting tree inheritance'''
    sequence_length = len(tree.sequence)
    for node in tree.traverse():
        for site in range(sequence_length):
            base = node.sequence[site]
            if base not in 'ACGT':
                new_base = random.choice(ambiguous_dna_values[base])
                for node2 in node.traverse(is_leaf_fn=lambda n: False if base in [n2.sequence[site] for n2 in n.children] else True):
                    if node2.sequence[site] == base:
                        node2.sequence = node2.sequence[:site] + new_base + node2.sequence[(site+1):]
    return tree


class MutationModel():
    '''a class for a mutation model, and functions to mutate sequences'''
    def __init__(self, mutability_file=None, substitution_file=None, mutation_order=True, with_replacement=True):
        """
        initialized with input files of the S5F format
        @param mutation_order: whether or not to mutate sequences using a context sensitive manner
                                where mutation order matters
        @param with_replacement: allow the same position to mutate multiple times on a single branch
        """
        self.mutation_order = mutation_order
        self.with_replacement = with_replacement
        if mutability_file is not None and substitution_file is not None:
            self.context_model = {}
            with open(mutability_file, 'r') as f:
                # eat header
                f.readline()
                for line in f:
                    motif, score = line.replace('"', '').split()[:2]
                    self.context_model[motif] = float(score)

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
                    self.context_model[motif] = (self.context_model[motif], {b:float(x) for b, x in zip('ACGT', fields[1:5])})
        else:
            self.context_model = None


    def mutabilities(self, sequence):
        '''returns the mutability of a sequence at each site, along with nucleotide biases'''
        assert all(n in 'ACGT' for n in sequence)
        sequence_length = len(sequence)
        if self.context_model is not None:
            # mutabilities of each nucleotide
            mutabilities = []
            assert sequence_length >= 5
            # ambiguous left end motifs
            for i in range(self.k//2):
                kmer_suffix = sequence[:(i+self.k//2+1)]
                matches = [value for key, value in self.context_model.iteritems() if key.endswith(kmer_suffix)]
                len_matches = len(matches)
                assert len_matches == 4**(self.k - len(kmer_suffix))
                # use mean over matches
                mutability = sum(match[0] for match in matches)/len_matches
                substitution = {n:sum(d[1][n] for d in matches)/len_matches for n in 'ACGT'}
                mutabilities.append((mutability, substitution))
            # unambiguous internal kmers
            for i in range(self.k//2, sequence_length - self.k//2):
                mutabilities.append(self.context_model[sequence[(i-self.k//2):(i+self.k//2+1)]])
            # ambiguous right end motifs
            for i in range(sequence_length - self.k//2, sequence_length):
                kmer_prefix = sequence[(i-self.k//2):]
                matches = [value for key, value in self.context_model.iteritems() if key.startswith(kmer_prefix)]
                len_matches = len(matches)
                assert len_matches == 4**(self.k - len(kmer_prefix))
                # use mean over matches
                mutability = sum(match[0] for match in matches)/len_matches
                substitution = {n:sum(d[1][n] for d in matches)/len_matches for n in 'ACGT'}
                mutabilities.append((mutability, substitution))
            return mutabilities
        else:
            return [(1, dict((n2, 1/3) if n2 is not n else (n2, 0.) for n2 in 'ACGT')) for n in sequence]

    def mutate(self, sequence, lambda0=1, frame=None):
        """
        Mutate a sequence, with lamdba0 the baseline mutability
        Cannot mutate the same position multiple times
        @param sequence: the original sequence to mutate
        @param lambda0: a "baseline" mutation rate
        @param frame: the reading frame index
        """
        sequence_length = len(sequence)
        if frame is not None:
            codon_start = frame-1
            codon_end = codon_start + 3*((sequence_length - codon_start)//3)
            if '*' in Seq(sequence[codon_start:codon_end], generic_dna).translate():
                raise RuntimeError('sequence contains stop codon!')

        mutabilities = self.mutabilities(sequence)
        sequence_mutability = sum(mutability[0] for mutability in mutabilities)/sequence_length
        # baseline Poisson
        # poisson rate for this sequence (given its relative mutability)
        lambda_sequence = sequence_mutability*lambda0
        # number of mutations
        trials = 20
        for trial in range(1, trials+1):
            m = scipy.random.poisson(lambda_sequence)
            if m <= sequence_length or self.with_replacement:
                break
            if trial == trials:
                raise RuntimeError('mutations saturating, consider reducing lambda0')

        # mutate the sites with mutations
        # if contains stop codon, try again, up to 10 times
        unmutated_positions = range(sequence_length)
        for i in range(m):
            sequence_list = list(sequence) # make string a list so we can modify it
            # Determine the position to mutate from the mutability matrix
            mutability_p = scipy.array([mutabilities[pos][0] for pos in unmutated_positions])
            for trial in range(1, trials+1):
                mut_pos = scipy.random.choice(unmutated_positions, p=mutability_p/mutability_p.sum())
                # Now draw the target nucleotide using the substitution matrix
                substitution_p = [mutabilities[mut_pos][1][n] for n in 'ACGT']
                assert 0 <= abs(sum(substitution_p) - 1.) < 1e-10
                chosen_target = scipy.random.choice(4, p=substitution_p)
                original_base = sequence_list[mut_pos]
                sequence_list[mut_pos] = 'ACGT'[chosen_target]
                sequence = ''.join(sequence_list) # reconstruct our sequence
                if frame is None or '*' not in Seq(sequence[codon_start:codon_end], generic_dna).translate():
                    if self.mutation_order:
                        # if mutation order matters, the mutabilities of the sequence need to be updated
                        mutabilities = self.mutabilities(sequence)
                    if not self.with_replacement:
                        # Remove this position so we don't mutate it again
                        unmutated_positions.remove(mut_pos)
                    break
                if trial == trials:
                    raise RuntimeError('stop codon in simulated sequence on '+str(trials)+' consecutive attempts')
                sequence_list[mut_pos] = original_base # <-- we only get here if we are retrying

        return sequence


    def one_mutant(self, sequence, Nmuts, frame=1, lambda0=0.1):
        '''
        Make a single mutant with a distance, in amino acid sequence, of Nmuts away from the starting point.
        '''
        trial = 100  # Allow 100 trials before quitting
        while trial > 0:
            mut_seq = sequence[:]
            aa = str(Seq(sequence[(frame-1):(frame-1+(3*(((len(sequence)-(frame-1))//3))))], generic_dna).translate())
            aa_mut = Seq(mut_seq[(frame-1):(frame-1+(3*(((len(mut_seq)-(frame-1))//3))))], generic_dna).translate()
            dist = hamming_distance(aa, aa_mut)
            while dist < Nmuts:
                mut_seq = self.mutate(mut_seq, lambda0=lambda0, frame=frame)
                aa_mut = str(Seq(mut_seq[(frame-1):(frame-1+(3*(((len(mut_seq)-(frame-1))//3))))], generic_dna).translate())
                dist = hamming_distance(aa, aa_mut)
            if dist == Nmuts:
                return aa_mut
            else:
                trial -= 1
        raise RuntimeError('100 consecutive attempts for creating a target sequence failed.')


    def simulate(self, sequence, seq_bounds=None, progeny=poisson(.9), lambda0=[1], frame=None,
                 N=None, T=None, n=None, verbose=False, selection_params=None):
        '''
        simulate neutral binary branching process with mutation model
        progeny must be like a scipy.stats distribution, with rvs() and mean() methods
        '''
        # Checking the validity of the input parameters:
        if N is not None and T is not None:
            raise ValueError('Only one of N and T can be used. One must be None.')
        if selection_params is not None and T is None:
            raise ValueError('Simulation with selection was chosen. A time, T, must be specified.')
        elif N is None and T is None:
            raise ValueError('Either N or T must be specified.')
            # expected_progeny = progeny.mean()
            # if expected_progeny >= 1:
            #     raise ValueError('E[progeny] = {} is not subcritical, tree termination not gauranteed!'.format(expected_progeny))
        if N is not None and n > N:
            raise ValueError('n ({}) must not larger than N ({})'.format(n, N))
        if selection_params is not None and frame is None:
            raise ValueError('Simulation with selection was chosen. A frame must must be specified.')
        if T is not None and len(T) > 1 and n is None:
            raise ValueError('When sampling intermediate time points it must be a subsample of the full population specified by "n".')

        # Planting the tree:
        tree = TreeNode()
        tree.dist = 0
        tree.add_feature('sequence', sequence)
        tree.add_feature('terminated', False)
        tree.add_feature('sampled', False)
        tree.add_feature('frequency', 0)
        tree.add_feature('time', 0)

        # <---- Selection:
        def calc_Kd(seqAA, targetAAseqs, hd2affy):
            '''Find the closest target sequence to and apply the "hamming distance to affinity" transformation function.'''
            hd = min([hamming_distance(seqAA, t) for t in targetAAseqs])
            return(hd2affy(hd))

        if selection_params is not None:
            hd_generation = list()  # Collect an array of the counts of each hamming distance at each time step
            mature_affy, naive_affy, target_dist, skip_update, targetAAseqs, A_total, B_total, Lp, k, outbase = selection_params
            # Assert that the target sequences are comparable to the naive sequence:
            aa = Seq(tree.sequence[(frame-1):(frame-1+(3*(((len(tree.sequence)-(frame-1))//3))))], generic_dna).translate()
            assert(sum([1 for t in targetAAseqs if len(t) != len(aa)]) == 0)  # All targets are same length
            assert(sum([1 for t in targetAAseqs if hamming_distance(aa, t) == target_dist]))  # All target are "target_dist" away from the naive sequence
            # Affinity is an exponential function of hamming distance:
            if target_dist > 0:
                def hd2affy(hd): return(mature_affy + hd**k * (naive_affy - mature_affy) / target_dist**k)
            else:
                def hd2affy(hd): return(mature_affy)
            # We store both the amino acid sequence and the affinity as tree features:
            tree.add_feature('AAseq', str(aa))
            tree.add_feature('Kd', calc_Kd(tree.AAseq, targetAAseqs, hd2affy))
            tree.add_feature('target_dist', min([hamming_distance(tree.AAseq, taa) for taa in targetAAseqs]))

        def lambda_selection(node, tree, targetAAseqs, hd2affy, A_total, B_total, Lp):
            '''
            Given a node and its tree and a "hamming distance to affinity" transformation function
            reutrn the poisson lambda parameter for the progeny distribution.
            '''
            def calc_BnA(Kd_n, A, B_total):
                '''
                This calculated the fraction B:A (B bound to A), at equilibrium also referred to as "binding time",
                of all the different Bs in the population given the number of free As in solution.
                '''
                BnA = B_total/(1+Kd_n/A)
                return(BnA)

            def return_objective_A(Kd_n, A_total, B_total):
                '''
                The objective function that solves the set of differential equations setup to find the number of free As,
                at equilibrium, given a number of Bs with some affinity listed in Kd_n.
                '''
                def obj(A): return((A_total - (A + np.sum(B_total/(1+Kd_n/A))))**2)
                return(obj)

            def calc_binding_time(Kd_n, A_total, B_total):
                '''
                Solves the objective function to find the number of free As and then uses this,
                to calculate the fraction B:A (B bound to A) for all the different Bs.
                '''
                obj = return_objective_A(Kd_n, A_total, B_total)
                # Different minimizers have been tested and 'L-BFGS-B' was significant faster than anything else:
                obj_min = minimize(obj, A_total, bounds=[[1e-10, A_total]], method='L-BFGS-B', tol=1e-20)
                BnA = calc_BnA(Kd_n, obj_min.x[0], B_total)
                # Terminate if the precision is not good enough:
                assert(BnA.sum()+obj_min.x[0]-A_total < A_total/100)
                return(BnA)

            def trans_BA(BA, Lp):
                '''Transform the fraction B:A (B bound to A) to a poisson lambda between 0 and 2.'''
                # We keep alpha to enable the possibility that there is a minimum lambda_:
                alpha, beta, Q = Lp
                lambda_ = alpha + (2 - alpha) / (1 + Q*np.exp(-beta*BA))
                return(lambda_)

            # Update the list of affinities for all the live nodes:
            Kd_n = np.array([n.Kd for n in tree.iter_leaves() if not n.terminated])
            BnA = calc_binding_time(Kd_n, A_total, B_total)
            lambdas = trans_BA(BnA, Lp)
            i = 0
            for n in tree.iter_leaves():
                if n.terminated:
                    continue
                n.add_feature('lambda_', lambdas[i])
                i += 1
            return(tree)
        # ----/> Selection

        t = 0  # <-- time
        leaves_unterminated = 1
        # Small lambdas are causing problems so make a minimum:
        lambda_min = 10e-10
        while leaves_unterminated > 0 and (leaves_unterminated < N if N is not None else True) and (t < max(T) if T is not None else True):
            t += 1
            if verbose:
                print('At time:', t)
            skip_lambda_n = 0  # At every new round reset all the lambdas
            unterminated_leaves = [l for l in tree.iter_leaves() if not l.terminated]
            random.shuffle(unterminated_leaves)
            # Sample intermediate time point:
            if T is not None and len(T) > 1 and (t-1) in T:
                if len(unterminated_leaves) < n:
                    raise RuntimeError('tree terminated with {} leaves, less than what desired after downsampling {}'.format(leaves_unterminated, n))
                # Make the sample and kill the cells sampled:
                for leaf in random.sample(unterminated_leaves, n):
                    leaves_unterminated -= 1
                    leaf.sampled = True
                    leaf.terminated = True
                if verbose:
                    print('Made an intermediate sample at time:', t-1)
                # Update the list of unterminated leafs:
                unterminated_leaves = [l for l in tree.iter_leaves() if not l.terminated]

            for leaf in unterminated_leaves:
                # <---- Selection:
                if selection_params is not None:
                    if skip_lambda_n == 0:
                        skip_lambda_n = skip_update + 1  # Add one so skip_update=0 is no skip
                        tree = lambda_selection(leaf, tree, targetAAseqs, hd2affy, A_total, B_total, Lp)
                    if leaf.lambda_ > lambda_min:
                        progeny = poisson(leaf.lambda_)
                    else:
                        progeny = poisson(lambda_min)
                    skip_lambda_n -= 1
                # ----/> Selection
                n_children = progeny.rvs()
                leaves_unterminated += n_children - 1 # <-- this kills the parent if we drew a zero
                if not n_children:
                    leaf.terminated = True
                for child_count in range(n_children):
                    # If sequence pair mutate them separately with their own mutation rate:
                    if seq_bounds is not None:
                        mutated_sequence1 = self.mutate(leaf.sequence[seq_bounds[0][0]:seq_bounds[0][1]], lambda0=lambda0[0], frame=frame)
                        mutated_sequence2 = self.mutate(leaf.sequence[seq_bounds[1][0]:seq_bounds[1][1]], lambda0=lambda0[1], frame=frame)
                        mutated_sequence = mutated_sequence1 + mutated_sequence2
                    else:
                        mutated_sequence = self.mutate(leaf.sequence, lambda0=lambda0[0], frame=frame)
                    child = TreeNode()
                    child.dist = sum(x!=y for x,y in zip(mutated_sequence, leaf.sequence))
                    child.add_feature('sequence', mutated_sequence)
                    # <---- Selection:
                    if selection_params is not None:
                        aa = Seq(child.sequence[(frame-1):(frame-1+(3*(((len(child.sequence)-(frame-1))//3))))], generic_dna).translate()
                        child.add_feature('AAseq', str(aa))
                        child.add_feature('Kd', calc_Kd(child.AAseq, targetAAseqs, hd2affy))
                        child.add_feature('target_dist', min([hamming_distance(child.AAseq, taa) for taa in targetAAseqs]))
                    # ----/> Selection
                    child.add_feature('frequency', 0)
                    child.add_feature('terminated', False)
                    child.add_feature('sampled', False)
                    child.add_feature('time', t)
                    leaf.add_child(child)
            # <---- Selection:
            if selection_params is not None:
                hd_distrib = [min([hamming_distance(tn.AAseq, ta) for ta in targetAAseqs]) for tn in tree.iter_leaves() if not tn.terminated]
                if target_dist > 0:
                    hist = np.histogram(hd_distrib, bins=list(range(target_dist*10)))
                else:  # Just make a minimum of 10 bins
                    hist = np.histogram(hd_distrib, bins=list(range(10)))
                hd_generation.append(hist)
                if verbose and hd_distrib:
                    print('Total cell population:', sum(hist[0]))
                    print('Majority hamming distance:', np.argmax(hist[0]))
                    print('Affinity of latest sampled leaf:', leaf.Kd)
                    print('Progeny distribution lambda for the latest sampled leaf:', leaf.lambda_)
            # ----/> Selection

        # <---- Selection:
        if selection_params is not None:
            # Keep a histogram of the hamming distances at each generation:
            with open(outbase + 'selection_sim.runstats.p', 'wb') as f:
                pickle.dump(hd_generation, f)
        # ----/> Selection

        if leaves_unterminated < N:
            raise RuntimeError('tree terminated with {} leaves, {} desired'.format(leaves_unterminated, N))

        # each leaf in final generation gets an observation frequency of 1, unless downsampled
        if T is not None and len(T) > 1:
            # Iterate the intermediate time steps:
            for Ti in sorted(T)[:-1]:
                # Only sample those that have been 'sampled' at intermediate sampling times:
                final_leaves = [leaf for leaf in tree.iter_descendants() if leaf.time == Ti and leaf.sampled]
                if len(final_leaves) < n:
                    raise RuntimeError('tree terminated with {} leaves, less than what desired after downsampling {}'.format(leaves_unterminated, n))
                for leaf in final_leaves:  # No need to down-sample, this was already done in the simulation loop
                    leaf.frequency = 1
        if selection_params and max(T) != t:
            raise RuntimeError('tree terminated with before the requested sample time.')

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
            raise RuntimeError('tree terminated with {} leaves, less than what desired after downsampling {}'.format(leaves_unterminated, n))
        else:
            raise RuntimeError('Unknown option.')

        # prune away lineages that are unobserved
        for node in tree.iter_descendants():
            if sum(node2.frequency for node2 in node.traverse()) == 0:
                node.detach()

        # remove unobserved unifurcations
        for node in tree.iter_descendants():
            parent = node.up
            if node.frequency == 0 and len(node.children) == 1:
                node.delete(prevent_nondicotomic=False)
                node.children[0].dist = hamming_distance(node.children[0].sequence, parent.sequence)

        # assign unique names to each node
        for i, node in enumerate(tree.traverse(), 1):
            node.name = 'simcell_gctreeinternal_{}'.format(i)

        # return the fine (uncollapsed) tree
        return tree


def test(args):
    '''
    test subprogram
    checks likelihood against a by-hand calculation for a simple tree, simulates a forest, computes MLE parameters, and plots some sanity check figures to plot_file
    command line arguments are p, q, number of trees to simulate, and plot file name
    '''

    import seaborn as sns
    sns.set(style='white', color_codes=True)
    sns.set_style('ticks')
    plt.rc('text', usetex=True)

    # compare likelihood to empirical likelihood (might need large n)
    n = 10000
    df = pd.DataFrame(columns=('p', 'q', 'parameters', 'f', 'L'))
    ct = 0
    ps = (.1, .2, .3, .4)
    qs = (.2, .4, .6, .8)
    scipy.seterr(all='ignore')
    for p in ps:
        for q in qs:
            forest = CollapsedForest((p, q), n)
            print('parameters: p = {}, q = {}'.format(p, q))
            forest.simulate()
            tree_dict = {}
            for tree in forest.forest:
                tree_hash = tuple((node.frequency, len(node.children)) for node in tree.tree.traverse())
                if tree_hash not in tree_dict:
                    tree_dict[tree_hash] = [tree, tree.l((p, q))[0], 1]
                else:
                    tree_dict[tree_hash][-1] += 1
            L_empirical, L_theoretical = zip(*[(tree_dict[tree_hash][2], scipy.exp(tree_dict[tree_hash][1])) for tree_hash in tree_dict if tree_dict[tree_hash][2]])
            for tree_hash in tree_dict:
                df.loc[ct] = (p, q, 'p={}, q={}'.format(p, q), tree_dict[tree_hash][2], scipy.exp(tree_dict[tree_hash][1]))
                ct += 1
    print()

    plt.figure()
    limx = (1/n, 1.1)
    limy = (1, 1.1*n)
    g = sns.lmplot(x='L', y='f', col='p', row='q', hue='parameters', data=df,
                   fit_reg=False, scatter_kws={'alpha':.3}, size=1.2, legend=False,
                   row_order=reversed(qs))
    g.set(xscale='log', yscale='log', xlim=limx, ylim=limy)
    g.fig.subplots_adjust(hspace=0.4, wspace=0.4)
    # g.ax.legend(title='parameters', loc='lower right', frameon=True, bbox_to_anchor=(1.1, 0))
    for i in range(len(ps)):
        for j in range(len(qs)):
            g.axes[i, j].plot(limx, limy, ls='--', c='black', lw=.5, zorder=0, markeredgewidth=.1)
            g.axes[i, j].set_title('$p={}$\n$q={}$'.format(ps[j], list(reversed(qs))[i]), x=.05, y=.9, size='x-small', ha='left', va='top')
    g.set_axis_labels('', '')
    # g.axes[-1, len(ps)//2].set_xlabel('GCtree likelihood')
    # g.axes[len(qs)//2, 0].set_ylabel('frequency among {} simulations'.format(n))
    g.fig.text(0.45, .02, s='GCtree likelihood', multialignment='center')
    g.fig.text(.03, 0.7, s='frequency among {} simulations'.format(n), rotation=90, multialignment='center')
    plt.savefig(args.outbase+'.pdf')

    # MLE check
    n = 10
    n2 = 1000
    df = pd.DataFrame(columns=('true parameters', '$\hat{p}$', '$\hat{q}$'))
    ct = 0
    for p in ps:
        for q in qs:
            for _ in range(n):
                forest = CollapsedForest((p, q), n2)
                print('parameters: p = {}, q = {}'.format(p, q))
                forest.simulate()
                result = forest.mle()
                df.loc[ct] = ('p={}, q={}'.format(p, q), result.x[0], result.x[1])
                ct += 1

    plt.figure()
    g = sns.lmplot(x='$\hat{p}$', y='$\hat{q}$', hue='true parameters', data=df,
                   fit_reg=False, scatter_kws={'alpha':.5}, size=4.5, legend=False)
    g.set(xlim=(0.05, .45), xticks=scipy.arange(0., .6, .1), ylim=(.1, .9), yticks=scipy.arange(0., 1.2, .2))
    for i in range(len(ps)):
        for j in range(len(qs)):
            plt.scatter([ps[i]], [qs[j]], c='black', marker='+', zorder=0)
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # plt.tight_layout()
    plt.savefig(args.outbase+'.2.pdf')

    return



    p = args.p
    q = args.q
    n = args.n

    plot_file = args.outbase

    if plot_file[-4:] != '.pdf':
        plot_file += '.pdf'

    print('Let''s check our likelihood against a by-hand calculation for the following simple tree')
    # ete tree
    parent = TreeNode(format=1)
    parent.add_feature('frequency', 2)
    parent.name = parent.frequency
    child = TreeNode()
    child.add_feature('frequency', 1)
    child.dist = 1
    child.name = child.frequency
    parent.add_child(child)
    tree = CollapsedTree(tree=parent)
    f = 6*p**2*(1-p)**3*q*(1-q)**3
    dfdp = 6*(1 - p)**2*p*(-2 + 5*p)*(-1 + q)**3*q #6*q*(1-q)**3*(2*p*(1-p)**3-3*p**2*(1-p)**2)
    dfdq = 6*(-1 + p)**3*p**2*(1 - q)**2*(-1 + 4*q) #6*p**2*(1-p)**3*((1-q)**3-3*q*(1-q)**2)
    print( '    T =', tree.tree.get_ascii(show_internal=True))
    print( '    Summing the probabilities of the two possible fine structures, we have')
    print( '    logP =', scipy.log(f))
    print(u'    \u2207logP = ', (dfdp/f, dfdq/f))
    print( '    Now, our dynamic programming algorithm gives')
    print(u'    logP , \u2207logP =', tree.l((p, q)))
    print('')

    # total leaf counts
    total_data = sorted([sum(node.frequency for node in tree.tree.traverse()) for tree in forest.forest])
    max_total = max(total_data)
    len_total = len(total_data)

    totals = []
    freq = []
    log_prob = []
    for x in range(1, max_total+1):
        totals.append(x)
        freq.append(total_data.count(x))
        tmp_tree = TreeNode(format=1)
        tmp_tree.add_feature('frequency', x)
        log_prob.append(CollapsedTree(tree=tmp_tree).l((p, 0))[0])
    theoretical_cdf = scipy.cumsum(scipy.exp(log_prob))
    empirical_cdf = scipy.cumsum(freq)/len_total

    #sns.reset_orig() # <-- don't use seaborn
    fig = plt.figure()
    fig.set_tight_layout(True)
    plt.rc('text', usetex=True)

    # plot the empirical and theoretical distribution of total leaf counts

    ax = fig.add_subplot(2,2,1)
    ax.plot(totals, scipy.exp(log_prob), 'ko', markerfacecolor='none', alpha=.5, label='theoretical PMF')
    ax.plot(totals, scipy.array(freq)/len_total, 'k.', label='empirical PMF')
    ax.legend(numpoints=1, loc=1, fontsize='small')
    ax.set_xlabel('total leaves')
    ax.set_ylabel('$\Pr($total leaves$)$')
    ax.set_ylim([0, 1.1])
    #ax.set_xscale('log')
    #ax.set_yscale('symlog')

# uncomment this if you want the CDF
#    ax = fig.add_subplot(2,2,2)
#    ax.plot(totals, theoretical_cdf, 'ko', markerfacecolor='none', alpha=.5, label='theoretical CDF')
#    ax.plot(totals, empirical_cdf, 'k.', label='empirical CDF')
#    ax.legend(numpoints=1, loc=4, fontsize='small')
#    ax.set_xlabel('number of leaves')
#    ax.set_ylim([0, 1.1])


    empirical_quantiles = []
    theoretical_quantiles = []
    for x in total_data:
        empirical_quantiles.append(sum(y <= x for y in total_data)/len_total)
        to_add = 0.
        for y in range(1, x+1):
            tmp_tree = TreeNode(format=1)
            tmp_tree.add_feature('frequency', y)
            to_add += scipy.exp(CollapsedTree(tree=tmp_tree).l((p, 0))[0])
        theoretical_quantiles.append(to_add)

    ax = fig.add_subplot(2,2,2)
    ax.plot(theoretical_quantiles, empirical_quantiles, 'ko', alpha=.1)
    ax.plot([0, 1], [0, 1], 'k')
    ax.set_title('total leaves')
    ax.set_xlabel('theoretical quantiles')
    ax.set_ylabel('empirical quantiles')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_aspect('equal')

    mle = forest.mle()

    #for tree in forest.forest:
    #    print(tree)
    print('    MLE parameters:  p = {}, q = {}'.format(*mle.x.tolist()))

    # plot the 2-norm of the difference between the gradient and its finite difference approximation
    print('computing plot data...')
    X, Y = scipy.mgrid[slice(.05, 1, .05),
                       slice(.05, 1, .05)]
    Z = scipy.zeros((X.shape[0], X.shape[1]))
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            Z[i, j] = check_grad(lambda x: forest.l(x)[0], lambda x: forest.l(x)[1], (X[i, j], Y[i, j]))

    print('done')
    ax = fig.add_subplot(2,2,3)
    ax.set_title(r'$||\nabla \ell(p, q) - \Delta \ell(p, q)||_2$')
    im = ax.contourf(X, Y, Z, locator=ticker.LogLocator(), cmap='Greys')
    ax.set_xlabel(r'$p$')
    ax.set_ylabel(r'$q$')
    ax.set_aspect('equal')
    fig.colorbar(im, ax=ax)


    # plot likelihood surface, with true and MLE parameters shown
    X, Y = scipy.mgrid[slice(.02, 1, .02),
                       slice(.02, 1, .02)]
    Z = scipy.zeros((X.shape[0], X.shape[1]))
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            l, grad_l = forest.l((X[i, j], Y[i, j]))
            z = l
            Z[i, j] = z
    ax = fig.add_subplot(2,2,4)
    ax.set_title(r'$\ell(p, q)$')
    contour = ax.contour(X, Y, Z, 10, colors='k', label='likelihood contours')
    for c in contour.collections:
        c.set_linestyle('solid')

    ax.clabel(contour, fontsize=8, inline=1)
    ax.plot([p], [q], 'k+', label='true parameters')
    ax.plot(mle.x[0], mle.x[1], 'ko', markerfacecolor='none', label='MLE parameters')
    ax.set_xlabel(r'$p$')
    ax.set_ylabel(r'$q$')
    ax.set_aspect('equal')
    ax.legend(numpoints = 1, fontsize='small')

    plt.savefig(plot_file)
    print('plot saved to', plot_file)







def infer(args):
    '''inference subprogram'''
    phylip_collapsed = [CollapsedTree(tree=tree, frame=args.frame) for tree in phylip_parse.parse_outfile(args.phylipfile, args.countfile, args.naive)]
    phylip_collapsed_unique = []
    for tree in phylip_collapsed:
        if sum(tree.compare(tree2, method='identity') for tree2 in phylip_collapsed_unique) == 0:
            phylip_collapsed_unique.append(tree)

    parsimony_forest = CollapsedForest(forest=phylip_collapsed_unique)

    if parsimony_forest.n_trees == 1:
        warnings.warn('only one parsimony tree reported from dnapars')

    print('number of trees with integer branch lengths:', parsimony_forest.n_trees)

    # check for unifurcations at root
    unifurcations = sum(tree.tree.frequency == 0 and len(tree.tree.children) == 1 for tree in parsimony_forest.forest)
    if unifurcations:
        print('WARNING: {} trees exhibit unobserved unifurcation from root, which is not possible under current model. Adding psuedocounts to these nodes'.format(unifurcations))

    # fit p and q using all trees
    # if we get floating point errors, try a few more times (starting params are random)
    max_tries = 10
    for tries in range(max_tries):
        try:
            parsimony_forest.mle(Vlad_sum=True)
            break
        except FloatingPointError as e:
            if tries + 1 < max_tries:
                print('floating point error in MLE: {}. Attempt {} of {}. Rerunning with new random start.'.format(e, tries+1, max_tries))
            else:
                raise
        else:
            raise

    print('params = {}'.format(parsimony_forest.params))

    # get likelihoods and sort by them
    ls = [tree.l(parsimony_forest.params)[0] for tree in parsimony_forest.forest]
    ls, parsimony_forest.forest = zip(*sorted(zip(ls, parsimony_forest.forest), reverse=True))

    with open(args.outbase+'.inference.parsimony_forest.p', 'wb') as f:
        pickle.dump(parsimony_forest, f)

    if args.colormapfile is not None:
        with open(args.colormapfile, 'r') as f:
            colormap = {}
            for line in f:
                seqid, color = line.rstrip().split('\t')
                if ',' in color:
                    colors = {x.split(':')[0]:int(x.split(':')[1]) for x in color.split(',')}
                    colormap[seqid] = colors
                else:
                    colormap[seqid] = color
    else:
        colormap = None

    print('tree\talleles\tlogLikelihood')
    for i, (l, collapsed_tree) in enumerate(zip(ls, parsimony_forest.forest), 1):
        alleles = sum(1 for _ in collapsed_tree.tree.traverse())
        print('{}\t{}\t{}'.format(i, alleles, l))
        collapsed_tree.render(args.outbase+'.inference.{}.svg'.format(i), colormap, args.leftright_split)

    # rank plot of likelihoods
    plt.figure(figsize=(6.5,2))
    try:
        plt.plot(scipy.exp(ls), 'ko', clip_on=False, markersize=4)
        plt.ylabel('GCtree likelihood')
        plt.yscale('log')
        plt.ylim([None, 1.1*max(scipy.exp(ls))])
    except FloatingPointError:
        plt.plot(ls, 'ko', clip_on=False, markersize=4)
        plt.ylabel('GCtree log-likelihood')
        plt.ylim([None, 1.1*max(ls)])
    plt.xlabel('parsimony tree')
    plt.xlim([-1, len(ls)])
    plt.tick_params(axis='y', direction='out', which='both')
    plt.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off') # labels along the bottom edge are off
    plt.savefig(args.outbase + '.inference.likelihood_rank.pdf')

    # rank plot of observed allele frequencies
    y = sorted((node.frequency for node in parsimony_forest.forest[0].tree.traverse() if node.frequency != 0), reverse=True)
    plt.figure()
    plt.bar(range(1, len(y) + 1), y, color='black')
    plt.xlabel('genotype')
    plt.ylabel('abundance')
    plt.savefig(args.outbase + '.inference.abundance_rank.pdf')


def find_A_total(carry_cap, B_total, f_full, mature_affy, U):
    def A_total_fun(A, B_total, Kd_n): return(A + np.sum(B_total/(1+Kd_n/A)))

    def C_A(A, A_total, f_full, U): return(U * (A_total - A) / f_full)

    def A_obj(carry_cap, B_total, f_full, Kd_n, U):
        def obj(A): return((carry_cap - C_A(A, A_total_fun(A, B_total, Kd_n), f_full, U))**2)
        return(obj)

    Kd_n = np.array([mature_affy] * carry_cap)
    obj = A_obj(carry_cap, B_total, f_full, Kd_n, U)
    # Some funny "zero encountered in true_divide" errors are not affecting results so ignore them:
    old_settings = np.seterr(all='ignore')  # Keep old settings
    np.seterr(divide='ignore')
    obj_min = minimize(obj, 1e-20, bounds=[[1e-20, carry_cap]], method='L-BFGS-B', tol=1e-20)
    np.seterr(**old_settings)  # Reset to default
    A = obj_min.x[0]
    A_total = A_total_fun(A, B_total, Kd_n)
    assert(C_A(A, A_total, f_full, U) > carry_cap * 99/100)
    return(A_total)


def find_Lp(f_full, U):
    assert(U > 1)
    def T_BA(BA, p):
        # We keep alpha to enable the possibility
        # that there is a minimum lambda_
        alpha, beta, Q = p
        lambda_ = alpha + (2 - alpha) / (1 + Q*np.exp(-beta*BA))
        return(lambda_)

    def solve_T_BA(p, f_full, U):
        epsilon = 1/1000
        C1 = (T_BA(0, p) - 0)**2
        C2 = (T_BA(f_full/U, p) - 1)**2
        C3 = (T_BA(1*f_full, p) - (2 - 2*epsilon))**2
        return(C1, C2, C3)

    def solve_T_BA_low_epsilon(p, f_full, U):
        epsilon = 1/1000
        C1 = (T_BA(0, p) - 0)**2
        C2 = (T_BA(f_full/U, p) - 1)**2
        C3 = (T_BA(1*f_full, p) - (2 - 2*epsilon))**2 * ((2 - T_BA(1*f_full, p)) < 2*epsilon)
        return(C1, C2, C3)

    # Some funny "FloatingPointError" errors are not affecting results so ignore them:
    old_settings = np.seterr(all='ignore')  # Keep old settings
    np.seterr(over='ignore')
    try:
        def obj_T_A(p): return(solve_T_BA(p, f_full, U))
        p = fsolve(obj_T_A, (0, 10e-5, 1), xtol=1e-20, maxfev=1000)
        assert(sum(solve_T_BA(p, f_full, U)) < f_full * 1/1000)
    except:
        print('The U parameter is large and therefore the epsilon parameter has to be adjusted to find a valid solution.')
        def obj_T_A(p): return(solve_T_BA_low_epsilon(p, f_full, U))
        p = fsolve(obj_T_A, (0, 10e-5, 1), xtol=1e-20, maxfev=1000)
        assert(sum(solve_T_BA(p, f_full, U)) < f_full * 1/1000)
    np.seterr(**old_settings)  # Reset to default
    return(p)


def simulate(args):
    '''
    Simulation subprogram. Can simulate in two modes.
    a) Neutral mode. A Galton–Watson process, with mutation probabilities according to a user defined motif model e.g. S5F
    b) Selection mode. Using the same mutation process as in a), but in selection mode the poisson progeny distribution's lambda
    is variable accordring to the hamming distance to a list of target sequences. The closer a sequence gets to one of the targets
    the higher fitness and the closer lambda will approach 2, vice versa when the sequence is far away lambda approaches 0.
    '''
    mutation_model = MutationModel(args.mutability, args.substitution)
    if args.lambda0 is None:
        args.lambda0 = [max([1, int(.01*len(args.sequence))])]
    args.sequence = args.sequence.upper()
    if args.sequence2 is not None:
        if len(args.lambda0) == 1:  # Use the same mutation rate on both sequences
            args.lambda0 = [args.lambda0[0], args.lambda0[0]]
        elif len(args.lambda0) != 2:
            raise Exception('Only one or two lambda0 can be defined for a two sequence simulation.')
        # Require both sequences to be in frame 1:
        if args.frame is not None and args.frame != 1:
            print('Warning: When simulating with two sequences they are truncated to be beginning at frame 1.')
            args.sequence = args.sequence[(args.frame-1):(args.frame-1+(3*(((len(args.sequence)-(args.frame-1))//3))))]
            args.sequence2 = args.sequence2[(args.frame-1):(args.frame-1+(3*(((len(args.sequence2)-(args.frame-1))//3))))]
        # Extract the bounds between sequence 1 and 2:
        seq_bounds = ((0, len(args.sequence)), (len(args.sequence), len(args.sequence)+len(args.sequence2)))
        # Merge the two seqeunces to simplify future dealing with the pair:
        args.sequence += args.sequence2.upper()
    else:
        seq_bounds = None
    # <---- Selection:
    if args.selection:
        if args.frame is None:
            raise Exception('Frame must be defined when simulating with selection.')
        assert(args.B_total >= args.f_full)  # the fully activating fraction on BA must be possible to reach within B_total
        # Make a list of target sequences:
        targetAAseqs = [mutation_model.one_mutant(args.sequence, args.target_dist, frame=args.frame) for i in range(args.target_count)]
        # Find the total amount of A necessary for sustaining the inputted carrying capacity:
        print((args.carry_cap, args.B_total, args.f_full, args.mature_affy))
        A_total = find_A_total(args.carry_cap, args.B_total, args.f_full, args.mature_affy, args.U)
        # Calculate the parameters for the logistic function:
        Lp = find_Lp(args.f_full, args.U)
        selection_params = [args.mature_affy, args.naive_affy, args.target_dist, args.skip_update, targetAAseqs, A_total, args.B_total, Lp, args.k, args.outbase]
    else:
        selection_params = None
    # ----/> Selection

    trials = 1000
    # this loop makes us resimulate if size too small, or backmutation
    for trial in range(trials):
        try:
            tree = mutation_model.simulate(args.sequence,
                                           seq_bounds=seq_bounds,
                                           progeny=poisson(args.lambda_),
                                           lambda0=args.lambda0,
                                           n=args.n,
                                           N=args.N,
                                           T=args.T,
                                           frame=args.frame,
                                           verbose=args.verbose,
                                           selection_params=selection_params)
            if args.selection:
                collapsed_tree = CollapsedTree(tree=tree, frame=args.frame, collapse_syn=True, allow_repeats=True)
            else:
                collapsed_tree = CollapsedTree(tree=tree, frame=args.frame) # <-- this will fail if backmutations
            tree.ladderize()
            uniques = sum(node.frequency > 0 for node in collapsed_tree.tree.traverse())
            if uniques < 2:
                raise RuntimeError('collapsed tree contains {} sampled sequences'.format(uniques))
            break
        except RuntimeError as e:
            print('{}, trying again'.format(e))
        else:
            raise
    if trial == trials - 1:
        raise RuntimeError('{} attempts exceeded'.format(trials))

    # In the case of a sequence pair print them to separate files:
    if args.sequence2 is not None:
        fh1 = open(args.outbase+'.simulation_seq1.fasta', 'w')
        fh2 = open(args.outbase+'.simulation_seq2.fasta', 'w')
        fh1.write('>naive\n')
        fh1.write(args.sequence[seq_bounds[0][0]:seq_bounds[0][1]]+'\n')
        fh2.write('>naive\n')
        fh2.write(args.sequence[seq_bounds[1][0]:seq_bounds[1][1]]+'\n')
        i = 0
        for leaf in tree.iter_leaves():
            if leaf.frequency != 0:# and '*' not in Seq(leaf.sequence, generic_dna).translate():
                i += 1
                fh1.write('>simcell{}\n'.format(i))
                fh1.write(leaf.sequence[seq_bounds[0][0]:seq_bounds[0][1]]+'\n')
                fh2.write('>simcell{}\n'.format(i))
                fh2.write(leaf.sequence[seq_bounds[1][0]:seq_bounds[1][1]]+'\n')
                leaf.name = 'simcell{}'.format(i)
    else:
        with open(args.outbase+'.simulation.fasta', 'w') as f:
            f.write('>naive\n')
            f.write(args.sequence+'\n')
            i = 0
            for leaf in tree.iter_leaves():
                if leaf.frequency != 0:# and '*' not in Seq(leaf.sequence, generic_dna).translate():
                    i += 1
                    f.write('>simcell{}\n'.format(i))
                    f.write(leaf.sequence+'\n')
                    leaf.name = 'simcell{}'.format(i)

    # some observable simulation stats to write
    frequency, distance_from_naive, degree = zip(*[(node.frequency,
                                                    hamming_distance(node.sequence, args.sequence),
                                                    sum(hamming_distance(node.sequence, node2.sequence) == 1 for node2 in collapsed_tree.tree.traverse() if node2.frequency and node2 is not node))
                                                   for node in collapsed_tree.tree.traverse() if node.frequency])
    stats = pd.DataFrame({'genotype abundance':frequency,
                          'Hamming distance to root genotype':distance_from_naive,
                          'Hamming neighbor genotypes':degree})
    stats.to_csv(args.outbase+'.simulation.stats.tsv', sep='\t', index=False)

    print('{} simulated observed sequences'.format(sum(leaf.frequency for leaf in collapsed_tree.tree.traverse())))

    # render the full lineage tree
    ts = TreeStyle()
    ts.rotation = 90
    ts.show_leaf_name = False
    ts.show_scale = False

    colors = {}
    palette = SVG_COLORS
    palette -= set(['black', 'white', 'gray'])
    palette = cycle(list(palette))  # <-- circular iterator

    # Either plot by DNA sequence or amino acid sequence:
    if args.plotAA and args.selection:
        colors[tree.AAseq] = 'gray'
    else:
        colors[tree.sequence] = 'gray'

    for n in tree.traverse():
        nstyle = NodeStyle()
        nstyle["size"] = 10
        if args.plotAA:
            if n.AAseq not in colors:
                colors[n.AAseq] = next(palette)
            nstyle['fgcolor'] = colors[n.AAseq]
        else:
            if n.sequence not in colors:
                colors[n.sequence] = next(palette)
            nstyle['fgcolor'] = colors[n.sequence]
        n.set_style(nstyle)

    tree.render(args.outbase+'.simulation.lineage_tree.svg', tree_style=ts)

    # render collapsed tree
    # create an id-wise colormap
    if args.plotAA and args.selection:
        colormap = {node.name:colors[node.AAseq] for node in collapsed_tree.tree.traverse()}
    else:
        colormap = {node.name:colors[node.sequence] for node in collapsed_tree.tree.traverse()}
    collapsed_tree.write( args.outbase+'.simulation.collapsed_tree.p')
    collapsed_tree.render(args.outbase+'.simulation.collapsed_tree.svg', colormap=colormap)
    # print colormap to file
    with open(args.outbase+'.simulation.collapsed_tree.colormap.tsv', 'w') as f:
        for name, color in colormap.items():
            f.write(name + '\t' + color + '\n')


    if args.selection:
        # Define a list a suitable colors that are easy to distinguish:
        palette = ['crimson', 'purple', 'hotpink', 'limegreen', 'darkorange', 'darkkhaki', 'brown', 'lightsalmon', 'darkgreen', 'darkseagreen', 'darkslateblue', 'teal', 'olive', 'wheat', 'magenta', 'lightsteelblue', 'plum', 'gold']
        palette = cycle(list(palette)) # <-- circular iterator
        colors = {i: next(palette) for i in range(int(len(args.sequence) // 3))}
        # The minimum distance to the target is colored:
        colormap = {node.name:colors[node.target_dist] for node in collapsed_tree.tree.traverse()}
        collapsed_tree.write( args.outbase+'.simulation.collapsed_runstat_color_tree.p')
        collapsed_tree.render(args.outbase+'.simulation.collapsed_runstat_color_tree.svg', colormap=colormap)
        # Write a file with the selection run stats. These are also plotted:
        with open(args.outbase + 'selection_sim.runstats.p', 'rb') as fh:
            runstats = pickle.load(fh)
            plot_runstats(runstats, args.outbase, colors)


def plot_runstats(runstats, outbase, colors):
    def make_bounds(runstats):
        all_counts = runstats[0][0].copy()
        for l in runstats:
            all_counts += l[0]
        i = None
        ii = None
        for j, c in enumerate(all_counts):
            if i is None and c > 0:
                i = j
            elif c > 0:
                ii = j
        return(i, ii)
    # Total population size:
    pop_size = np.array([sum(r[0]) for r in runstats])
    # min:max of the hamming distances to plot:
    bounds = make_bounds(runstats)

    fig = plt.figure()
    ax = plt.subplot(111)
    t = np.array(list(range(len(pop_size))))  # The x-axis are generations
    ax.plot(t, pop_size, lw=2, label='All cells')  # Total population size is plotted
    # Then plot the counts for each hamming distance as a function on generation:
    for k in list(range(*bounds)):
        color = SVG2HEX[colors[k]]
        ax.plot(t, np.array([r[0][k] for r in runstats]), lw=2, color=color, label='Dist {}'.format(k))

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)

    # Shrink current axis by 20% to make the legend fit:
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    plt.ylabel('Count')
    plt.xlabel('GC generation')
    plt.title('Cell count as function of GC generation')
    fig.savefig(outbase + '.selection_sim.runstats.pdf')


def main():
    import argparse
    parser = argparse.ArgumentParser(description='germinal center tree inference and simulation')
    subparsers = parser.add_subparsers(help='which program to run')

    # parser for test subprogram
    parser_test = subparsers.add_parser('test',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        help='run tests on library functions')
    parser_test.set_defaults(func=test)

    # parser for inference subprogram
    parser_infer = subparsers.add_parser('infer',
                                         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                         help='likelihood ranking of parsimony trees')
    parser_infer.add_argument('--naive', type=str, default=None, help='name of naive sequence (outgroup root)')
    parser_infer.add_argument('phylipfile', type=str, help='dnapars outfile (verbose output with sequences at each site)')
    parser_infer.add_argument('countfile', type=str, help='File containing allele frequencies (sequence counts) in the format: "SeqID,Nobs"')
    parser_infer.add_argument('--colormapfile', type=str, default=None, help='File containing color map in the format: "SeqID\tcolor"')
    parser_infer.add_argument('--leftright_split', type=int, default=None, help='split between heavy and light for combined seqs')
    parser_infer.set_defaults(func=infer)

    # parser for simulation subprogram
    parser_sim = subparsers.add_parser('simulate',
                                       formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                       help='Neutral, and target selective, model gctree simulation')
    parser_sim.add_argument('sequence', type=str, help='seed naive nucleotide sequence')
    parser_sim.add_argument('mutability', type=str, help='path to mutability model file')
    parser_sim.add_argument('substitution', type=str, help='path to substitution model file')
    parser_sim.add_argument('--sequence2', type=str, default=None, help='Second seed naive nucleotide sequence. For simulating heavy/light chain co-evolution.')
    parser_sim.add_argument('--lambda', dest='lambda_', type=float, default=.9, help='poisson branching parameter')
    parser_sim.add_argument('--lambda0', type=float, default=None, nargs='*', help='List of one or two elements with the baseline mutation rates. Space separated input values.'
                            'First element belonging to seed sequence one and optionally the next to sequence 2. If only one rate is provided for two sequences,'
                            'this rate will be used on both.')
    parser_sim.add_argument('--n', type=int, default=None, help='cells downsampled')
    parser_sim.add_argument('--N', type=int, default=None, help='target simulation size')
    parser_sim.add_argument('--T', type=int, default=None, help='observation time, if None we run until termination and take all leaves')
    parser_sim.add_argument('--selection', type=bool, default=False, help='Simulation with selection? true/false. When doing simulation with selection an observation time cut must be set.')
    parser_sim.add_argument('--carry_cap', type=int, default=1000, help='The carrying capacity of the simulation with selection. This number affects the fixation time of a new mutation.'
                            'Fixation time is approx. log2(carry_cap), e.g. log2(1000) ~= 10.')
    parser_sim.add_argument('--target_count', type=int, default=10, help='The number of targets to generate.')
    parser_sim.add_argument('--target_dist', type=int, default=10, help='The number of non-synonymous mutations the target should be away from the naive.')
    parser_sim.add_argument('--naive_affy', type=float, default=100, help='Affinity of the naive sequence in nano molar.')
    parser_sim.add_argument('--mature_affy', type=float, default=1, help='Affinity of the mature sequences in nano molar.')
    parser_sim.add_argument('--skip_update', type=int, default=100, help='When iterating through the leafs the B:A fraction is recalculated every time.'
                            'It is possible though to update less often and get the same approximate results. This parameter sets the number of iterations to skip,'
                            'before updating the B:A results. skip_update < carry_cap/10 recommended.')
    parser_sim.add_argument('--B_total', type=float, default=1, help='Total number of BCRs per B cell normalized to 10e4. So 1 equals 10e4, 100 equals 10e6 etc.'
                            'It is recommended to keep this as the default.')
    parser_sim.add_argument('--U', type=float, default=5, help='Controls the fraction of BCRs binding antigen necessary to only sustain the life of the B cell'
                            'It is recommended to keep this as the default.')
    parser_sim.add_argument('--f_full', type=float, default=1, help='The fraction of antigen bound BCRs on a B cell that is needed to elicit close to maximum reponse.'
                            'Cannot be smaller than B_total. It is recommended to keep this as the default.')
    parser_sim.add_argument('--k', type=float, default=2, help='The exponent in the function to map hamming distance to affinity.'
                            'It is recommended to keep this as the default.')
    parser_sim.add_argument('--plotAA', type=bool, default=False, help='Plot trees with collapsing and coloring on amino acid level.')
    parser_sim.add_argument('--verbose', type=bool, default=False, help='Print progress during simulation. Mostly useful for simulation with selection since this can take a while.')
    parser_sim.set_defaults(func=simulate)

    # a common outbase parameter
    for subparser in [parser_test, parser_infer, parser_sim]:
        subparser.add_argument('--outbase', type=str, default='gctree.out', help='output file base name')

    # a common parameter for the inference and simulation subprograms
    for subparser in [parser_infer, parser_sim]:
        subparser.add_argument('--frame', type=int, default=None, choices=(1, 2, 3), help='codon frame')

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
