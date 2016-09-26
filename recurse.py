#! /bin/env python

import scipy, warnings, random
from scipy.optimize import minimize, check_grad

"""
This module contains classes for simulation and inference for a binary branching process with mutation
in which the tree is collapsed to nodes that count the number of clonal leaves of each type
"""

class LeavesAndClades():
    """
    This is a base class for simulating, and computing likelihood for, a binary infinite type branching
    process with branching probability p, mutation probability q, and we collapse mutant clades off the
    root type and consider just the number of clone leaves, c, and mutant clades, m.

      /\            
     /\ ^          (3)
      /\     ==>   / \\
       /\\
        ^
    """
    def __init__(self, p=None, q=None, c=None, m=None):
        """initialize with branching probability p and mutation probability q, both in the unit interval"""
        if p is not None or q is not None:
            if not (0 <= p <= 1 and 0 <= q <= 1):
                raise ValueError('p and q must be in the unit interval')
        self._p = p
        self._q = q
        if c is not None or m is not None:
            if not (c >= 0) and (m >= 0) and (c+m > 0):
                raise ValueError('c and m must be nonnegative integers summing greater than zero')
            self._c = c
            self._m = m
        self._extinction_time = None # <-- number of generations clone line survives for

    def simulate(self):
        """simulate the number of clone leaves and mutant clades off a root node"""
        if self._p>=.5:
            warnings.warn('p >= .5 is not subcritical, tree simulations not garanteed to terminate')
        if self._p is None or self._q is None:
            raise ValueError('p and q parameters must be defined for simulation\n')

        # let's track the tree in breadth first order, listing number clone and mutant descendants of each node
        # mutant clades terminate in this view
        cumsum_clones = 0
        len_tree = 0
        self._c = 0
        self._m = 0
        # while termination condition not met
        while cumsum_clones > len_tree - 1:
            if random.random() < self._p:
                mutants = sum(random.random() < self._q for child in range(2))
                clones = 2 - mutants 
                self._m += mutants
            else:
                mutants = 0
                clones = 0
                self._c += 1
            cumsum_clones += clones
            len_tree += 1
        assert cumsum_clones == len_tree - 1

    f_hash = {} # <--- class variable for hashing calls to the following function
    def f(self, p, q, sign=1):
        """
        Probability of getting c leaves that are clones of the root and m mutant clades off
        the root line, given branching probability p and mutation probability q 
        Also returns gradient wrt (p, q)
        Computed by dynamic programming
        """
        c, m = self._c, self._m
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
                    neighbor = LeavesAndClades(p=p, q=q, c=c, m=m-1)
                    neighbor_f, (neighbor_dfdp, neighbor_dfdq) = neighbor.f(p, q)
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
                            neighbor1 = LeavesAndClades(p=p, q=q, c=cx, m=mx)
                            neighbor2 = LeavesAndClades(p=p, q=q, c=c-cx, m=m-mx)
                            neighbor1_f, (neighbor1_dfdp, neighbor1_dfdq) = neighbor1.f(p, q)
                            neighbor2_f, (neighbor2_dfdp, neighbor2_dfdq) = neighbor2.f(p, q)
                            f_result += p*(1-q)**2*neighbor1_f*neighbor2_f
                            dfdp_result +=   (1-q)**2 * neighbor1_f    * neighbor2_f + \
                                           p*(1-q)**2 * neighbor1_dfdp * neighbor2_f + \
                                           p*(1-q)**2 * neighbor1_f    * neighbor2_dfdp
                            dfdq_result += -2*p*(1-q) * neighbor1_f    * neighbor2_f + \
                                           p*(1-q)**2 * neighbor1_dfdq * neighbor2_f + \
                                           p*(1-q)**2 * neighbor1_f    * neighbor2_dfdq
            LeavesAndClades.f_hash[(p, q, c, m)] = (f_result, scipy.array([dfdp_result, dfdq_result]))
        return LeavesAndClades.f_hash[(p, q, c, m)]

    def get(self, param_name=None):
        """
        return a dictionary of member variables, or a single parameter indicated by param_name
        param_name may equal 'p', 'q', or 'tree', or None.
        """
        if param_name is None:
            return {'p':self._p, 'q':self._q, 'c':self._c, 'm':self._m}
        elif param_name is 'p':
            return self._p
        elif param_name is 'q':
            return self._q
        elif param_name is 'c':
            return self._c
        elif param_name is 'm':
            return self._m
        else:
            raise ValueError("param_name may equal 'p', 'q', 'c', 'm', or None.")


class CollapsedTree(LeavesAndClades):
    """
    Here's a derived class for a collapsed tree, where we recurse into the mutant clades
          (4)
         / | \\
       (3)(1)(2)
           |   \\
          (2)  (1)
    """
    def __init__(self, p=None, q=None, tree=None):
        """
        For intialization, either p and q or tree (or all three) must be provided
        p: branching probability
        q: mutation probability
        tree: Clonal leaf count and count of mutant clades are provided as tuples in
        breadth first order.
        """
        if p is None and q is None and tree is None:
            raise ValueError('either p and q or tree (or all three) must be provided')
        LeavesAndClades.__init__(self, p=p, q=q)
        # check that tree is valid
        if tree is not None:
            k = len(tree)
            if k==0 or \
               set(map(len, tree))!=set([2]) or \
               not all(scipy.greater_equal([x for y in tree for x in y], 0)):
                raise ValueError('"tree" must be a nonempty list of 2-element tuples of nonnegative integers')
            cs, ms = zip(*tree)
            if not all(scipy.greater(scipy.cumsum([x[1] for x in tree])[:-1], scipy.arange(1, k)-1)) \
            or sum(x[1] for x in tree) != k - 1:
                raise ValueError('inconsistent breadth first tree data')
        self._tree = tree
        self._extinction_time = None

    def phi(self, x, n):
        """
        The nth composition of the generating function of the offspring distribution
        This is the generating function of the total number of (uncollapsed) nodes in the nth generation
        Note: since collapsed tree simulations don't currently capture fine structure, this is of limited use
        """
        if n == 1:
            return (1-self._p) + self._p*x**2
        elif n > 1:
            return phi(x, n-1)
        else:
            raise ValueError('n must be a natural number')

    def sf(self, n):
        """
        The survival function of the extinction time, n (integer number of generations), of the uncollapsed tree
        This is computed analytically in terms of the generating funtion, phi, of the offsprint distribution
        Note: since collapsed tree simulations don't currently capture fine structure, this is of limited use
        """
        return 1 - phi(self, 0, n)

    def l(self, (p, q), sign=1):
        """
        log likelihood of p and q, conditioned on collapsed tree, and its gradient wrt (p, q)
        optional parameter sign must be 1 or -1, with the latter useful for MLE by minimization
        """
        if self._tree is None:
            raise ValueError('tree data must be defined to compute likelihood')
        if sign not in (-1, 1):
            raise ValueError('sign must be 1 or -1')
        f_data = [LeavesAndClades(c=c, m=m).f(p, q) for (c, m) in self._tree]
        # extract vector of function values and gradient components
        fs = scipy.array([x[0] for x in f_data])
        dfdps = scipy.array([x[1][0] for x in f_data])
        dfdqs = scipy.array([x[1][1] for x in f_data])
        return sign*scipy.log(fs).sum(), sign*scipy.array([(dfdps/fs).sum(), (dfdqs/fs).sum()])

    def mle(self):
        """
        Maximum likelihood estimate for p and q given tree
        updates p and q if not None
        returns optimization result
        """
        # random initalization
        x_0 = (random.random(), random.random())
        bounds = ((.001, .999), (.001, .999))
        result = minimize(self.l, x0=x_0, args=(-1,), jac=True, method='L-BFGS-B', bounds=bounds)
        # update p and q if None and optimization successful
        if not result.success:
            warnings.warn('optimization not sucessful, '+result.message, RuntimeWarning)
        elif self._p is None and self._q is None:
            self._p, self._q = result.x
        return result

    def simulate(self):
        """
        simulate a collapsed tree given p and q
        replaces existing tree data member with simulation result, and returns self
        """
        if self._p is None or self._q is None:
            raise ValueError('p and q parameters must be defined for simulation')

#        This code simulates the fine structure of the tree, possibly useful if we want extinction time
#        # initiate with an unmutated root node (type 0)
#        # this tree is the full tree
#        tree = [(2 if random.random() < self._p else 0, 0)]
#        # while termination condition not met
#        cumsum_offspring = tree[0][0]
#        len_tree = 0
#        while cumsum_offspring > len_tree - 1:
#            offspring = 2 if random.random() < self._p else 0
#            a_type = random.random() < self._q
#            tree.append((offspring, is_mutant))
#            cumsum_offspring += offspring
#            len_tree += 1
#        assert cumsum_offspring == len_tree - 1
#
#        # use breadth first structure to figure out generation boundaries
#        self._extinction_time = 1
#        nextgen_size = tree[0][0]
#        i = 1
#        while i + nextgen_size < len_tree:
#            next_nextgen_size = sum(offspring for offspring, is_mutant in tree[i:(i+nextgen_size)])
#            i += nextgen_size
#            nextgen_size = next_nextgen_size
#            self._extinction_time += 1
#
#        # now we would need to collapse the tree given the simulated fine structure...
        

        # initiate by running a LeavesAndClades simulation to get the number of clones and mutants
        # in the root node of the collapsed tree
        LeavesAndClades.simulate(self)
        tree = [(self._c, self._m)] # <-- accessing member variable in base class (updated by base class method)
        # now for each mutant off the root, we do a LeavesAndClades simulation, recursing
        more_mutants = tree[0][1] # aka self._m
        while more_mutants > 0:
            new_nodes = []
            for m in range(more_mutants):
                LeavesAndClades.simulate(self)
                new_nodes.append((self._c, self._m))
            more_mutants = sum(x[1] for x in new_nodes) # mutant clades from this generation
            tree += new_nodes
        self._tree = tree # replace tree data member
        return self
                
    def get(self, param_name=None):
        """
        return a dictionary of member variables, or a single parameter indicated by param_name
        param_name may equal 'p', 'q', or 'tree', or None.
        """
        if param_name is None:
            return {'p':self._p, 'q':self._q, 'tree':self._tree}
        elif param_name is 'p':
            return self._p
        elif param_name is 'q':
            return self._q
        elif param_name is 'tree':
            return self._tree
        else:
            raise ValueError("param_name may equal 'p', 'q', or 'tree', or None.")

    def __str__(self):
        """return a string representation for printing"""
        return 'p = %f, q = %f\ntree: ' % (self._p, self._q) + str(self._tree)
        
class CollapsedForest(CollapsedTree):
    """
    simply a set of CollapsedTrees, with the same p and q parameters
          (4)          (3)
         / | \\         / \\
       (3)(1)(2)     (1) (2)
           |   \\  ,          , ...
          (2)  (1)
    """
    def __init__(self, p=None, q=None, n_trees=None, forest=None):
        """
        in addition to p and q, we need number of trees
        can also intialize with forest, a list of trees, each same format as tree member of CollapsedTree
        """
        if p is not None and q is not None:
            CollapsedTree.__init__(self, p=p, q=q)
        elif forest is None:
            raise ValueError('either p and q or forest (or all three) must be provided')
        if forest is not None:
            if len(forest) == 0:
                raise ValueError('passed empty tree list')
            if n_trees is not None and len(forest) != n_trees:
                raise ValueError('n_trees not consistent with forest')
            self._forest = forest
        if n_trees is not None and n_trees < 1:
            raise ValueError('number of trees must be at least one')
        if n_trees is None and forest is not None:
            self._n_trees = len(forest)
        self._n_trees = n_trees
        
    def simulate(self):
        """
        simulate a forest of collapsed trees given p and q and number of trees
        replaces existing forest data member with simulation result, and returns self
        """
        if self._p is None or self._q is None or self._n_trees is None:
            raise ValueError('p, q, and n_trees parameters must be defined for simulation')
        tree = CollapsedTree(self._p, self._q)
        self._forest = [tree.simulate().get('tree') for x in range(self._n_trees)]
        return self

    def l(self, (p, q), sign=1):
        """
        likelihood of (p, q), given forest, and it's gradient wrt (p, q)
        optional parameter sign must be 1 or -1, with the latter useful for MLE by minimization
        """
        if self._forest is None:
            raise ValueError('forest data must be defined to compute likelihood')
        if sign not in (-1, 1):
            raise ValueError('sign must be 1 or -1')
        # since the l method on the CollapsedTree class returns l and grad_l...
        terms = [CollapsedTree(tree=tree).l((p, q), sign=sign) for tree in self._forest]
        return sum(x[0] for x in terms), scipy.array([sum(x[1][0] for x in terms), sum(x[1][1] for x in terms)])

    # I would have coded a method for Maximum likelihood method for p and q given forest of independent trees
    # but we get this for free from inheritance and polymorphism magic.

    def get(self, param_name=None):
        """
        return a dictionary of member variables (None argument), or a single parameter indicated by param_name
        param_name may equal 'p', 'q', 'n_trees', or 'forest'.
        """
        if param_name is None:
            return {'p':self._p, 'q':self._q, 'n_trees':self._n_trees, 'forest':self._forest}
        elif param_name is 'p':
            return self._p
        elif param_name is 'q':
            return self._q
        elif param_name is 'n_trees':
            return self._n_trees
        elif param_name is 'forest':
            return self._forest
        else:
            raise ValueError("param_name may equal 'p', 'q', or 'tree', or None.")

    def __str__(self):
        """return a string representation for printing"""
        return ('p = %f, q = %f, n_trees = %d\n'+
                '\n'.join([str(tree) for tree in self._forest])) % (self._p, self._q, self._n_trees)
        
def test(p, q, n, plot_file):
    """
    checks likelihood against a by-hand calculation for a simple tree, simulates a forest, computes MLE parameters, and plots some sanity check figures to plot_file
    command line arguments are p, q, number of trees to simulate, and plot file name
    """

    if plot_file[-4:] != '.pdf':
        plot_file += '.pdf'

    import matplotlib
    matplotlib.use('PDF')
    from matplotlib import pyplot as plt
    from matplotlib import rc, ticker
    from scipy.stats import probplot

    print 'Let''s check our likelihood against a by-hand calculation for the following simple tree'
    tree = CollapsedTree(tree=[(2,1), (1,0)])
    print '    T =', str(tree.get('tree'))
    print '    Summing the probabilities of the two possible fine structures, we have'
    print '    Pr(T) = 6 p^2 (1-p)^3 q (1-q)^3 =', 6*p**2*(1-p)**3*q*(1-q)**3
    print '    Now, our dynamic programming algorithm gives'
    print '    Pr(T) =', scipy.exp(tree.l((p, q))[0])
    print ''

    print 'Simulating a forest of %d trees' % n
    forest = CollapsedForest(p, q, n)
    print '    true parameters: p = %f, q = %f' % (p, q)
    forest.simulate()

    # total leaf counts
    total_data = sorted([sum(x[0] for x in tree) for tree in forest.get('forest')])
    max_total = max(total_data)
    len_total = len(total_data)

    totals, freq, log_prob = zip(*[(x, total_data.count(x), CollapsedTree(tree=[(x, 0)]).l((p, 0))[0]) for x in range(1, max_total+1)])
    theoretical_cdf = scipy.cumsum(scipy.exp(log_prob))
    empirical_cdf = scipy.cumsum(freq)/float(len_total)

    fig = plt.figure()
    fig.set_tight_layout(True)
    plt.rc('text', usetex=True)

    # plot the empirical and theoretical distribution of total leaf counts

    ax = fig.add_subplot(2,2,1)
    ax.plot(totals, scipy.exp(log_prob), 'ko', markerfacecolor='none', alpha=.5, label='theoretical PMF')
    ax.plot(totals, scipy.array(freq)/float(len_total), 'k.', label='empirical PMF')
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
        empirical_quantiles.append(sum(y <= x for y in total_data)/float(len_total))
        theoretical_quantiles.append(scipy.sum(scipy.exp([CollapsedTree(tree=[(y, 0)]).l((p, 0))[0] for y in range(1, x+1)])))

    ax = fig.add_subplot(2,2,2)
    ax.plot(theoretical_quantiles, empirical_quantiles, 'ko', alpha=.1)
    ax.plot([0, 1], [0, 1], 'k')
    ax.set_title('total leaves')
    ax.set_xlabel('theoretical quantiles')
    ax.set_ylabel('empirical quantiles')
    ax.set_xlim([0, 1])
    ax.set_ylim([0, 1])
    ax.set_aspect('equal')

    mle = forest.mle().x
    print '    MLE parameters:  p = %f, q = %f' % tuple(mle.tolist())

    # plot the 2-norm of the difference between the gradient and its finite difference approximation
    print 'computing plot data...'
    X, Y = scipy.mgrid[slice(.05, 1, .05),
                       slice(.05, 1, .05)]
    Z = scipy.zeros((X.shape[0], X.shape[1]))
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            Z[i, j] = check_grad(lambda x: forest.l(x)[0], lambda x: forest.l(x)[1], (X[i, j], Y[i, j]))

    print 'done'
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
    contour = ax.contour(X, Y, Z, colors='k', label='likelihood contours')
    for c in contour.collections:
        c.set_linestyle('solid')

    ax.clabel(contour, fontsize=9, inline=1)
    ax.plot([p], [q], 'k+', label='true parameters')
    ax.plot(mle[0], mle[1], 'ko', markerfacecolor='none', label='MLE parameters')
    ax.set_xlabel(r'$p$')
    ax.set_ylabel(r'$q$')
    ax.set_aspect('equal')
    ax.legend(numpoints = 1, fontsize='small')

    plt.savefig(plot_file)
    print 'plot saved to', plot_file

def main():
    """if "--test" option is passed, run the test suite, else load newick file and do MLEs for each tree"""
    import sys, argparse
    parser = argparse.ArgumentParser(description='multitype tree modeling')
    parser.add_argument('--test', action='store_true', default=False, help='run tests on library functions')
    parser.add_argument('--p', type=float, default=.4, help='branching probability for test mode')
    parser.add_argument('--q', type=float, default=.5, help='mutation probability for test mode')
    parser.add_argument('--n', type=int, default=100, help='forest size for test mode')
    parser.add_argument('--plot_file', type=str, default='foo.pdf', help='output file for plots from test mode')
    parser.add_argument('--germline', type=str, default=None, help='name of germline sequence (outgroup root)')

    #parser.add_argument('outfile', type=str, help='dnapars outfile (verbose output)')
    parser.add_argument('outtree', type=str, help='newick file of trees (dnapars outtree)')
    
    args = parser.parse_args()

    if args.test:
        test(args.p, args.q, args.n, args.plot_file)
        return
    
    #with open(args.outfile, 'r') as f:
    #    node_sequences = f.read().split('\n\n\n\n')
    #for i in range(len(node_sequences)):
    #    blocks = node_sequences[i].split('From    To     Any Steps?    State at upper node\n                            \n')[1].split('\n\n')
    #    for j in range(len(blocks)):
    #        if j == 0:
    #            # get edge ids
    #            edge_ids = []
    #            for x in blocks[j].rstrip().lstrip().split('\n')[1:]:
    #                edge_ids.append(x.split()[:2])
    #        blocks[j] = [x.split()[2] for x in blocks[j].split('\n')[1:]]
    #    print blocks
    #    assert all(len(block) == len(edge_ids) for block in blocks)
    #    node_sequences[i] = {tuple(edge_ids[j]):''.join([block[j] for block in blocks]) for j in range(len(edge_ids))}
        
        #print node_sequences[i]
    #return



    import dendropy, copy
    trees = dendropy.TreeList.get(path=args.outtree, schema='newick')

    # let's infer the seq length from the per site branch lengths
    lengths = [node.edge.length for node in trees[0] if node.edge.length is not None]
    min_delta = 100
    for length in range(1,1000):
        delta = sum(abs(round(x*length) - x*length) for x in lengths)
        if delta < min_delta:
            min_delta = delta
            best_length = length

    print 'inferred sequence length:', best_length

    # adjust branch lens and remove trees with fractional branches
    i = 0
    while i < len(trees):
        # if germline seq is defined, we reroot on it
        # the reroot_at_node function doesn't seem to work right because it expects an internal node
        if args.germline is not None:
            GL = trees[i].find_node_with_taxon_label(args.germline)
            assert len(GL.child_nodes()) == 0
            for child in trees[i].seed_node.child_node_iter():
                child.edge_length += GL.edge_length
            trees[i].seed_node.remove_child(GL)#, suppress_unifurcations=True)
            trees[i].seed_node.label = GL
        # we need to adjust the branch lengths, since they're per site
        for edge in trees[i].postorder_edge_iter():
            if edge.length is not None:
                edge.length = round(edge.length*best_length,1)
                # if branch length is fractional, forget about this tree
                if edge.length % 1 != 0:
                    trees.remove(trees[i])
                    bad_tree = True
                    break
                else:
                    bad_tree = False
        if not bad_tree:
            i += 1

    n_trees = len(trees)

    print 'number of trees with integer branch lengths:', n_trees
    # now we need to get collapsed trees, is there a less ugly way to do this in dendropy?
    best_likelihood_sofar = None
    for i, tree in enumerate(trees):
        collapsed_tree = []
        # the number of clonal leaf descendents is number of leaves we can get to on zero-length edges
        # root first
        clone_leaves = sum(node.distance_from_root() == 0 for node in tree.leaf_nodes())
        # to get mutant offspring, first consider all nodes that are distance zero
        # the mutant offspring are children of these nodes with nonzero edge length
        mutant_offspring_nodes = [child_node for node in tree if node.distance_from_root() == 0 for child_node in node.child_node_iter() if child_node.edge.length != 0]
        mutant_offspring_edge_lengths = [node.edge.length for node in mutant_offspring_nodes]
        #mutant_offspring = len(mutant_offspring_nodes)
        collapsed_tree.append((clone_leaves, mutant_offspring_edge_lengths))#mutant_offspring))
        # recurse into the mutant offspring
        done = False
        while not done:
            new_mutant_offspring_nodes = []
            for mutant in mutant_offspring_nodes:
                new_mutant_offspring_nodes_from_this_mutant = []
                if mutant.num_child_nodes() == 0:
                    clone_leaves = 1
                else:
                    clone_leaves = 0 # <-- initialize
                    # loop over nodes of the subtree rooted at this mutant
                    for clonal_descendent in mutant.preorder_iter():
                        # skip the root of the subtree
                        if clonal_descendent == mutant:
                            continue
                        # let's get the sequence of branch lengths connecting clonal_descendent with this mutant
                        # Looking at that sequence we can see if a node is a clone of the mutant, or the beginning 
                        # of a new mutant clade
                        distances = [clonal_descendent.edge_length]
                        clonal_descendent_parent = clonal_descendent.parent_node
                        while clonal_descendent_parent is not mutant:
                            distances.append(clonal_descendent_parent.edge_length)
                            clonal_descendent_parent = clonal_descendent_parent.parent_node
                        if sum(distances) == 0 and clonal_descendent.is_leaf():
                            clone_leaves += 1
                        elif distances[0] != 0 and all(distances[i] == 0 for i in range(1, len(distances))):
                            new_mutant_offspring_nodes_from_this_mutant.append(clonal_descendent)
                #mutant_offspring = len(new_mutant_offspring_nodes_from_this_mutant)
                mutant_offspring_edge_lengths = [node.edge.length for node in new_mutant_offspring_nodes_from_this_mutant]
                collapsed_tree.append((clone_leaves, mutant_offspring_edge_lengths))#mutant_offspring))
                new_mutant_offspring_nodes.extend(new_mutant_offspring_nodes_from_this_mutant)
            mutant_offspring_nodes = new_mutant_offspring_nodes
            if len(new_mutant_offspring_nodes) == 0:
                done = True

        # make sure number of leaves is consistent
        assert sum(clone_leaves for clone_leaves, mutant_offspring_edge_lengths in collapsed_tree) == len(tree.leaf_nodes())

        # ok, so now we have the parsimony tree in dendropy format and the corresponding collapsed tree in 
        # CollapsedTree format (but with explicit edge lengths). 

        result = CollapsedTree(tree=[(clone_leaves, len(mutant_offspring_edge_lengths)) for clone_leaves, mutant_offspring_edge_lengths in collapsed_tree]).mle()
        assert result.success

        if best_likelihood_sofar is None or -result.fun > best_likelihood_sofar:
            best_likelihood_sofar = -result.fun
            best_likelihood_sofar_params = result.x
            best_i = i
            best_tree = collapsed_tree
        print 'tree %d: l = %f, p = %f, q = %f' % (i+1, -result.fun, result.x[0], result.x[1])
        sys.stdout.flush()

        #break

    print 'best tree: tree %d, l = %f, p = %f, q = %f' % (best_i + 1,
                                                          best_likelihood_sofar,
                                                          best_likelihood_sofar_params[0],
                                                          best_likelihood_sofar_params[1])

    from ete3 import Tree, NodeStyle, TreeStyle, TextFace, add_face_to_node, CircleFace, faces, AttrFace

    # make an ete version of collapsed tree for plotting
    nodes = [Tree(name=clone_leaves) for clone_leaves, mutant_offspring in collapsed_tree]
    nodes[0].dist = 0 # zero length edge for root node
    #nodes = [dendropy.Node(label=str(clone_leaves)) for clone_leaves, mutant_offspring in collapsed_tree]
    #tree = dendropy.Tree(seed_node=nodes[0]) # seed_node=nodes[0])        
    gen_size = 1
    gen_start_index = 0
    terminated = False
    while not terminated:
        k = 0
        for i in range(gen_start_index,gen_start_index + gen_size):
            for j in range(len(collapsed_tree[i][1])):
                nodes[i].add_child(nodes[gen_start_index+gen_size+k], dist=collapsed_tree[i][1][j])
                k += 1
        new_gen_start_index = gen_start_index + gen_size
        gen_size = sum(len(x[1]) for x in collapsed_tree[gen_start_index:(gen_start_index + gen_size)])
        gen_start_index = new_gen_start_index
        if gen_size == 0:
            terminated = True

    for node in nodes:
        nstyle = NodeStyle()
        nstyle['size'] = 10*scipy.sqrt(node.name)
        node.set_style(nstyle)

    ts = TreeStyle()
    ts.show_leaf_name = False
    #ts.show_branch_length = True
    ts.rotation = 90
    def my_layout(node):
        N = AttrFace('name', fsize=14, fgcolor='black')
        N.rotation = -90
        faces.add_face_to_node(N, node, 0, position='branch-top')
        #C = CircleFace(radius=5*scipy.sqrt(node.name), color='RoyalBlue', style='sphere')
        #C.opacity = 1#0.5
        # And place as a float face over the tree
        #faces.add_face_to_node(C, node, 0, position="float")
    ts.layout_fn = my_layout
    nodes[0].render(args.plot_file, tree_style=ts)
    print 'tree plot saved to', args.plot_file
    

if __name__ == "__main__":
    main()





