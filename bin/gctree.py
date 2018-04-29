#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
This module contains classes for simulation and inference for a binary branching process with mutation
in which the tree is collapsed to nodes that count the number of clonal leaves of each type
'''

from __future__ import division, print_function

import phylip_parse, scipy, warnings, random, collections, pandas as pd, os
from scipy.misc import logsumexp
from scipy.optimize import minimize, check_grad, fsolve
from itertools import cycle
from scipy.stats import poisson
from ete3 import TreeNode, NodeStyle, TreeStyle, TextFace, add_face_to_node, CircleFace, PieChartFace, faces, SVG_COLORS
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import MultipleSeqAlignment
import matplotlib; matplotlib.use('agg')
from matplotlib import pyplot as plt, ticker
try:
    import cPickle as pickle
except:
    import pickle

from utils import hamming_distance
import selection_utils
from mutation_model import MutationModel

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

        # let's track the tree in breadth first order, listing number of clonal and mutant descendants of each node
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
            # remove unobserved internal unifurcations
            for node in self.tree.iter_descendants():
                parent = node.up
                if node.frequency == 0 and len(node.children) == 1:
                    node.delete(prevent_nondicotomic=False)
                    node.children[0].dist = hamming_distance(node.children[0].sequence, parent.sequence)

            # iterate over the tree below root and collapse edges of zero length
            # if the node is a leaf and it's parent has nonzero frequency we combine taxa names to a set
            # this acommodates bootstrap samples that result in repeated genotypes
            observed_genotypes = set((leaf.name for leaf in self.tree))
            observed_genotypes.add(self.tree.name)
            for node in self.tree.get_descendants(strategy='postorder'):
                if node.dist == 0:
                    node.up.frequency += node.frequency
                    node_set = set([node.name]) if isinstance(node.name, str) else set(node.name)
                    node_up_set = set([node.up.name]) if isinstance(node.up.name, str) else set(node.up.name)
                    if node_up_set < observed_genotypes:
                        if node_set < observed_genotypes:
                            node.up.name = tuple(node_set | node_up_set)
                            if len(node.up.name) == 1:
                                node.up.name = node.up.name[0]
                    elif node_set < observed_genotypes:
                        node.up.name = tuple(node_set)
                        if len(node.up.name) == 1:
                            node.up.name = node.up.name[0]
                    node.delete(prevent_nondicotomic=False)

            final_observed_genotypes = set([name for node in self.tree.traverse() if node.frequency > 0 or node == self.tree for name in ((node.name,) if isinstance(node.name, str) else node.name)])
            if final_observed_genotypes != observed_genotypes:
                raise RuntimeError('observed genotypes don\'t match after collapse\n\tbefore: {}\n\tafter: {}\n\tsymmetric diff: {}'.format(observed_genotypes, final_observed_genotypes, observed_genotypes ^ final_observed_genotypes))
            assert sum(node.frequency for node in tree.traverse()) == sum(node.frequency for node in self.tree.traverse())

            rep_seq = sum(node.frequency > 0 for node in self.tree.traverse()) - len(set([node.sequence for node in self.tree.traverse() if node.frequency > 0]))
            if not allow_repeats and rep_seq:
                raise RuntimeError('Repeated observed sequences in collapsed tree. {} sequences were found repeated.'.format(rep_seq))
            elif allow_repeats and rep_seq:
                rep_seq = sum(node.frequency > 0 for node in self.tree.traverse()) - len(set([node.sequence for node in self.tree.traverse() if node.frequency > 0]))
                print('Repeated observed sequences in collapsed tree. {} sequences were found repeated.'.format(rep_seq))
            # a custom ladderize accounting for abundance and sequence to break ties in abundance
            for node in self.tree.traverse(strategy='postorder'):
                # add a partition feature and compute it recursively up the tree
                node.add_feature('partition', node.frequency + sum(node2.partition for node2 in node.children))
                # sort children of this node based on partion and sequence
                node.children.sort(key=lambda node: (node.partition, node.sequence))
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

    def render(self, outfile, idlabel=False, colormap=None, show_support=False, chain_split=None):
        '''render to image file, filetype inferred from suffix, svg for color images'''
        def my_layout(node):
            circle_color = 'lightgray' if colormap is None or node.name not in colormap else colormap[node.name]
            text_color = 'black'
            if isinstance(circle_color, str):
                C = CircleFace(radius=max(3, 10*scipy.sqrt(node.frequency)), color=circle_color, label={'text':str(node.frequency), 'color':text_color} if node.frequency > 0 else None)
                C.rotation = -90
                C.hz_align = 1
                faces.add_face_to_node(C, node, 0)
            else:
                P = PieChartFace([100*x/node.frequency for x in circle_color.values()], 2*10*scipy.sqrt(node.frequency), 2*10*scipy.sqrt(node.frequency), colors=[(color if color != 'None' else 'lightgray') for color in list(circle_color.keys())], line_color=None)
                T = TextFace(' '.join([str(x) for x in list(circle_color.values())]), tight_text=True)
                T.hz_align = 1
                T.rotation = -90
                faces.add_face_to_node(P, node, 0, position='branch-right')
                faces.add_face_to_node(T, node, 1, position='branch-right')
            if idlabel:
                T = TextFace(node.name, tight_text=True, fsize=6)
                T.rotation = -90
                T.hz_align = 1
                faces.add_face_to_node(T, node, 1 if isinstance(circle_color, str) else 2, position='branch-right')
        for node in self.tree.traverse():
            nstyle = NodeStyle()
            nstyle['size'] = 0
            if node.up is not None:
                if set(node.sequence.upper()) == set('ACGT'):
                    if chain_split is not None:
                        if self.frame is not None:
                            raise NotImplementedError('frame not implemented with chain_split')
                        leftseq_mutated = hamming_distance(node.sequence[:chain_split], node.up.sequence[:chain_split]) > 0
                        rightseq_mutated = hamming_distance(node.sequence[chain_split:], node.up.sequence[chain_split:]) > 0
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
        ts.show_branch_support = show_support
        self.tree.render(outfile, tree_style=ts)
        # if we labelled seqs, let's also write the alignment out so we have the sequences (including of internal nodes)
        if idlabel:
            aln = MultipleSeqAlignment([])
            for node in self.tree.traverse():
                aln.append(SeqRecord(Seq(str(node.sequence), generic_dna), id=str(node.name), description='abundance={}'.format(node.frequency)))
            AlignIO.write(aln, open(os.path.splitext(outfile)[0] + '.fasta', 'w'), 'fasta')

    def write(self, file_name):
        '''serialize tree to file'''
        with open(file_name, 'wb') as f:
            pickle.dump(self, f)

    def newick(self, file_name):
        '''write to newick file'''
        self.tree.write(format=1, outfile=file_name)

    def compare(self, tree2, method='identity'):
        '''compare this tree to the other tree'''
        if method == 'identity':
            # we compare lists of seq, parent, abundance
            # return true if these lists are identical, else false
            list1 = sorted((node.sequence, node.frequency, node.up.sequence if node.up is not None else None) for node in self.tree.traverse())
            list2 = sorted((node.sequence, node.frequency, node.up.sequence if node.up is not None else None) for node in tree2.tree.traverse())
            return list1 == list2
        elif method == 'MRCA':
            # matrix of hamming distance of common ancestors of taxa
            # takes a true and inferred tree as CollapsedTree objects
            taxa = [node.sequence for node in self.tree.traverse() if node.frequency]
            n_taxa = len(taxa)
            d = scipy.zeros(shape=(n_taxa, n_taxa))
            sum_sites = scipy.zeros(shape=(n_taxa, n_taxa))
            for i in range(n_taxa):
                nodei_true = self.tree.iter_search_nodes(sequence=taxa[i]).next()
                nodei      =      tree2.tree.iter_search_nodes(sequence=taxa[i]).next()
                for j in range(i + 1, n_taxa):
                    nodej_true = self.tree.iter_search_nodes(sequence=taxa[j]).next()
                    nodej      =      tree2.tree.iter_search_nodes(sequence=taxa[j]).next()
                    MRCA_true = self.tree.get_common_ancestor((nodei_true, nodej_true)).sequence
                    MRCA =           tree2.tree.get_common_ancestor((nodei, nodej)).sequence
                    d[i, j] = hamming_distance(MRCA_true, MRCA)
                    sum_sites[i, j] = len(MRCA_true)
            return d.sum() / sum_sites.sum()
        elif method == 'RF':
            tree1_copy = self.tree.copy(method='deepcopy')
            tree2_copy = tree2.tree.copy(method='deepcopy')
            for treex in (tree1_copy, tree2_copy):
                for node in list(treex.traverse()):
                    if node.frequency > 0:
                        child = TreeNode()
                        child.add_feature('sequence', node.sequence)
                        node.add_child(child)
            try:
                return tree1_copy.robinson_foulds(tree2_copy, attr_t1='sequence', attr_t2='sequence', unrooted_trees=True)[0]
            except:
                return tree1_copy.robinson_foulds(tree2_copy, attr_t1='sequence', attr_t2='sequence', unrooted_trees=True, allow_dup=True)[0]
        else:
            raise ValueError('invalid distance method: ' + method)

    def get_split(self, node):
        '''return the bipartition resulting from clipping this node's edge above'''
        if node.get_tree_root() != self.tree:
            raise ValueError('node not found')
        if node == self.tree:
            raise ValueError('this node is the root (no split above)')
        parent = node.up
        taxa1 = []
        for node2 in node.traverse():
            if node2.frequency > 0 or node2 == self.tree:
                if isinstance(node2.name, str):
                    taxa1.append(node2.name)
                else:
                    taxa1.extend(node2.name)
        taxa1 = set(taxa1)
        node.detach()
        taxa2 = []
        for node2 in self.tree.traverse():
            if node2.frequency > 0 or node2 == self.tree:
                if isinstance(node2.name, str):
                    taxa2.append(node2.name)
                else:
                    taxa2.extend(node2.name)
        taxa2 = set(taxa2)
        parent.add_child(node)
        assert taxa1.isdisjoint(taxa2)
        assert taxa1.union(taxa2) == set((name for node in self.tree.traverse() if node.frequency > 0 or node == self.tree for name in ((node.name,) if isinstance(node.name, str) else node.name)))
        return tuple(sorted([taxa1, taxa2]))

    @staticmethod
    def split_compatibility(split1, split2):
        diff = split1[0].union(split1[1]) ^ split2[0].union(split2[1])
        if diff:
            raise ValueError('splits do not cover the same taxa\n\ttaxa not in both: {}'.format(diff))
        for partition1 in split1:
            for partition2 in split2:
                if partition1.isdisjoint(partition2):
                    return True
        return False

    def support(self, bootstrap_trees_list, weights=None, compatibility=False):
        '''
        compute support from a list of bootstrap GCtrees
        weights (optional) is needed for weighting parsimony degenerate trees
        compatibility mode counts trees that don't disconfirm the split
        '''
        for node in self.tree.get_descendants():
            split = self.get_split(node)
            support = 0
            compatibility_ = 0
            for i, tree in enumerate(bootstrap_trees_list):
                compatible = True
                supported = False
                for boot_node in tree.tree.get_descendants():
                    boot_split = tree.get_split(boot_node)
                    if compatibility and compatible and not self.split_compatibility(split, boot_split):
                        compatible = False
                    if not compatibility and not supported and boot_split == split:
                        supported = True
                if supported:
                    support += weights[i] if weights is not None else 1
                if compatible:
                    compatibility_ += weights[i] if weights is not None else 1
            node.support = compatibility_ if compatibility else support

        return self



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

    def l(self, params, sign=1, empirical_bayes_sum=False):
        '''
        likelihood of params, given forest, and it's gradient wrt params
        optional parameter sign must be 1 or -1, with the latter useful for MLE by minimization
        if optional parameter empirical_bayes_sum is true, we're doing the Vlad sum for estimating params for
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
        if empirical_bayes_sum:
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

def test(args):
    '''
    test subprogram
    checks likelihood against a by-hand calculation for a simple tree, simulates a forest, computes MLE parameters, and plots some sanity check figures to plot_file
    command line arguments are p, q, number of trees to simulate, and plot file name
    '''

    import seaborn as sns
    sns.set(style='white', color_codes=True)
    sns.set_style('ticks')

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
                   fit_reg=False, scatter_kws={'alpha':.3}, size=1.5, legend=False,
                   row_order=reversed(qs))
    g.set(xscale='log', yscale='log', xlim=limx, ylim=limy)
    for i in range(len(ps)):
        for j in range(len(qs)):
            g.axes[i, j].plot(limx, limy, ls='--', c='black', lw=.5, zorder=0, markeredgewidth=.1)
            g.axes[i, j].set_title('p={}\nq={}'.format(ps[j], list(reversed(qs))[i]), x=.05, y=.7, size='x-small', ha='left')
    g.set_axis_labels('', '')
    g.fig.text(0.45, .02, s='GCtree likelihood', multialignment='center')
    g.fig.text(.02, 0.7, s='frequency among {} simulations'.format(n), rotation=90, multialignment='center')
    plt.savefig(args.outbase+'.pdf')

    # MLE check
    n = 20
    n2 = 500
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
                   fit_reg=False, scatter_kws={'alpha':.2}, size=6, legend=False)
    g.set(xlim=(0.05, .45), xticks=scipy.arange(0., .6, .1), ylim=(.1, .9), yticks=scipy.arange(0., 1.2, .2))
    for i in range(len(ps)):
        for j in range(len(qs)):
            plt.scatter([ps[i]], [qs[j]], c='black', marker='+')
    plt.savefig(args.outbase+'.2.pdf')

    return


def infer(args):
    '''inference subprogram'''
    outfiles = [phylip_parse.parse_outfile(args.phylipfile, args.countfile, args.naive)]
    if args.bootstrap_phylipfile is not None:
        outfiles.extend(phylip_parse.parse_outfile(args.bootstrap_phylipfile, args.countfile, args.naive))
    bootstrap = len(outfiles) > 1
    if bootstrap:
        # store bootstrap parameter data
        df = pd.DataFrame(columns=('$\hat{p}$', '$\hat{q}$'))
        gctrees = [] # we'll store the mle gctrees here for computing support later

    for i, content in enumerate(outfiles):
        if i > 0:
            print('bootstrap sample {}'.format(i))
            print('----------------------')
            outbase = args.outbase + '.bootstrap_{}'.format(i)
        else:
            outbase = args.outbase
        phylip_collapsed = [CollapsedTree(tree=tree, frame=args.frame, allow_repeats=(i>0)) for tree in content]
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
            print('{} trees exhibit unobserved unifurcation from root. Adding psuedocounts to these roots'.format(unifurcations))

        # fit p and q using all trees
        # if we get floating point errors, try a few more times (starting params are random)
        max_tries = 10
        for tries in range(max_tries):
            try:
                parsimony_forest.mle(empirical_bayes_sum=True)
                break
            except FloatingPointError as e:
                if tries + 1 < max_tries:
                    print('floating point error in MLE: {}. Attempt {} of {}. Rerunning with new random start.'.format(e, tries+1, max_tries))
                else:
                    raise
            else:
                raise

        print('params = {}'.format(parsimony_forest.params))

        if i > 0:
            df.loc[i-1] = parsimony_forest.params

        # get likelihoods and sort by them
        ls = [tree.l(parsimony_forest.params)[0] for tree in parsimony_forest.forest]
        ls, parsimony_forest.forest = zip(*sorted(zip(ls, parsimony_forest.forest), reverse=True))

        if bootstrap:
            gctrees.append(parsimony_forest.forest[0])

        with open(outbase+'.inference.parsimony_forest.p', 'wb') as f:
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
        for j, (l, collapsed_tree) in enumerate(zip(ls, parsimony_forest.forest), 1):
            alleles = sum(1 for _ in collapsed_tree.tree.traverse())
            print('{}\t{}\t{}'.format(j, alleles, l))
            collapsed_tree.render(outbase+'.inference.{}.svg'.format(j),
                                  idlabel=args.idlabel,
                                  colormap=colormap,
                                  chain_split=args.chain_split)
            collapsed_tree.newick(outbase+'.inference.{}.nk'.format(j))

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
        plt.savefig(outbase + '.inference.likelihood_rank.pdf')

        # rank plot of observed allele frequencies
        y = sorted((node.frequency for node in parsimony_forest.forest[0].tree.traverse() if node.frequency != 0), reverse=True)
        plt.figure()
        plt.bar(range(1, len(y) + 1), y, color='black')
        plt.xlabel('genotype')
        plt.ylabel('abundance')
        plt.savefig(outbase + '.inference.abundance_rank.pdf')

    if bootstrap:
        import seaborn as sns; sns.set(style="white", color_codes=True)
        sns.set_style("ticks")
        plt.figure()
        scipy.seterr(all='ignore')
        g = sns.jointplot('$\hat{p}$', '$\hat{q}$', data=df, joint_kws={'alpha':.1, 's':.1},
                          space=0, color='k', stat_func=None, xlim=(0, 1), ylim=(0, 1), size=3)
        plt.savefig(args.outbase+'.inference.bootstrap_theta.pdf')
        gctrees[0].support(gctrees[1:])
        gctrees[0].render(args.outbase+'.inference.bootstrap_support.svg',
                          colormap=colormap,
                          idlabel=args.idlabel,
                          chain_split=args.chain_split,
                          show_support=True)
        gctrees[0].support(gctrees[1:], compatibility=True)
        gctrees[0].render(args.outbase+'.inference.bootstrap_compatibility.svg',
                          colormap=colormap,
                          idlabel=args.idlabel,
                          chain_split=args.chain_split,
                          show_support=True)


def simulate(args):
    '''
    Simulation subprogram. Can simulate in two modes.
    a) Neutral mode. A Galton–Watson process, with mutation probabilities according to a user defined motif model e.g. S5F
    b) Selection mode. Using the same mutation process as in a), but in selection mode the poisson progeny distribution's lambda
    is variable accordring to the hamming distance to a list of target sequences. The closer a sequence gets to one of the targets
    the higher fitness and the closer lambda will approach 2, vice versa when the sequence is far away lambda approaches 0.
    '''
    random.seed(a=args.seed)
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
    if args.selection:
        if args.frame is None:
            raise Exception('Frame must be defined when simulating with selection.')
        assert(args.B_total >= args.f_full)  # the fully activating fraction on BA must be possible to reach within B_total
        # Make a list of target sequences:
        targetAAseqs = [mutation_model.one_mutant(args.sequence, args.target_dist, frame=args.frame) for i in range(args.target_count)]
        # Find the total amount of A necessary for sustaining the inputted carrying capacity:
        print((args.carry_cap, args.B_total, args.f_full, args.mature_affy))
        A_total = selection_utils.find_A_total(args.carry_cap, args.B_total, args.f_full, args.mature_affy, args.U)
        # Calculate the parameters for the logistic function:
        Lp = selection_utils.find_Lp(args.f_full, args.U)
        selection_params = [args.stop_dist, args.mature_affy, args.naive_affy, args.target_dist, args.skip_update, targetAAseqs, A_total, args.B_total, Lp, args.k, args.outbase]
    else:
        selection_params = None

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
                collapsed_tree = CollapsedTree(tree=tree, frame=args.frame, collapse_syn=False, allow_repeats=True)
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
        for leaf in tree.iter_leaves():
            if leaf.frequency != 0:# and '*' not in Seq(leaf.sequence, generic_dna).translate():
                fh1.write('>' + leaf.name + '\n')
                fh1.write(leaf.sequence[seq_bounds[0][0]:seq_bounds[0][1]]+'\n')
                fh2.write('>' + leaf.name + '\n')
                fh2.write(leaf.sequence[seq_bounds[1][0]:seq_bounds[1][1]]+'\n')
    else:
        with open(args.outbase+'.simulation.fasta', 'w') as f:
            f.write('>naive\n')
            f.write(args.sequence+'\n')
            for leaf in tree.iter_leaves():
                if leaf.frequency != 0:# and '*' not in Seq(leaf.sequence, generic_dna).translate():
                    f.write('>' + leaf.name + '\n')
                    f.write(leaf.sequence + '\n')

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

    # this makes the rendered branch lenths correspond to time
    for node in tree.iter_descendants():
        node.dist = node.time - node.up.time
    tree.render(args.outbase+'.simulation.lineage_tree.svg', tree_style=ts)

    # render collapsed tree
    # create an id-wise colormap
    # NOTE: node.name can be a set
    if args.plotAA and args.selection:
        colormap = {node.name:colors[node.AAseq] for node in collapsed_tree.tree.traverse()}
    else:
        colormap = {node.name:colors[node.sequence] for node in collapsed_tree.tree.traverse()}
    collapsed_tree.write( args.outbase+'.simulation.collapsed_tree.p')
    collapsed_tree.render(args.outbase+'.simulation.collapsed_tree.svg',
                          idlabel=args.idlabel,
                          colormap=colormap)
    # print colormap to file
    with open(args.outbase+'.simulation.collapsed_tree.colormap.tsv', 'w') as f:
        for name, color in colormap.items():
            f.write((name if isinstance(name, str) else ','.join(name)) + '\t' + color + '\n')


    if args.selection:
        # Define a list a suitable colors that are easy to distinguish:
        palette = ['crimson', 'purple', 'hotpink', 'limegreen', 'darkorange', 'darkkhaki', 'brown', 'lightsalmon', 'darkgreen', 'darkseagreen', 'darkslateblue', 'teal', 'olive', 'wheat', 'magenta', 'lightsteelblue', 'plum', 'gold']
        palette = cycle(list(palette)) # <-- circular iterator
        colors = {i: next(palette) for i in range(int(len(args.sequence) // 3))}
        # The minimum distance to the target is colored:
        colormap = {node.name:colors[node.target_dist] for node in collapsed_tree.tree.traverse()}
        collapsed_tree.write( args.outbase+'.simulation.collapsed_runstat_color_tree.p')
        collapsed_tree.render(args.outbase+'.simulation.collapsed_runstat_color_tree.svg',
                              idlabel=args.idlabel,
                              colormap=colormap)
        # Write a file with the selection run stats. These are also plotted:
        with open(args.outbase + 'selection_sim.runstats.p', 'rb') as fh:
            runstats = pickle.load(fh)
            selection_utils.plot_runstats(runstats, args.outbase, colors)


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
    parser_infer.add_argument('--bootstrap_phylipfile', type=str, help='dnapars outfile from seqboot (multiple data sets)')
    parser_infer.add_argument('--colormapfile', type=str, default=None, help='File containing color map in the format: "SeqID\tcolor"')
    parser_infer.add_argument('--chain_split', type=int, default=None, help='split between heavy and light for combined seqs')
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
    parser_sim.add_argument('--T', type=int, nargs='+', default=None, help='observation time, if None we run until termination and take all leaves')
    parser_sim.add_argument('--seed', type=int, default=None, help='integer random seed')
    parser_sim.add_argument('--selection', type=bool, default=False, help='Simulation with selection? true/false. When doing simulation with selection an observation time cut must be set.')
    parser_sim.add_argument('--stop_dist', type=int, default=None, help='Stop when this distance has been reached in the selection model.')
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

    # common parameters for the inference and simulation subprograms
    for subparser in [parser_infer, parser_sim]:
        subparser.add_argument('--frame', type=int, default=None, choices=(1, 2, 3), help='codon frame')
        subparser.add_argument('--idlabel', action='store_true', help='flag for labeling the sequence ids of the nodes in the output tree images, also write associated fasta alignment if True')

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
