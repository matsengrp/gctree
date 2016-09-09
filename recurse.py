#! /bin/env python

import sys, scipy, warnings, random
from scipy.optimize import minimize, check_grad
import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
from matplotlib import rc

"""
This module contains classes for simulation and inference for a binary branching process with mutation
in which the tree is collapsed to nodes that count the number of clonal leaves of each type
"""

class LeavesAndClades():
    """
    This is a base class for simulating, and computing likelihood for, a binary infinite type branching process with branching
    probability p, mutation probability q, and we collapse mutant clades off the root type
    and consider just the number of clone leaves, c, and mutant clades, m.

      /\            
     /\ ^          (3)
      /\     ==>   / \\
       /\\
        ^
    """
    def __init__(self, p=None, q=None, c=None, m=None):
        """initialize with branching probability p and mutation probability q, both in the unit interval"""
        if p is not None or q is not None:
            if not (0 < p < 1 and 0 < q < 1):
                raise ValueError('p and q must be in the unit interval')
            self._p = p
            self._q = q
        if c is not None or m is not None:
            if not (c >= 0) and (m >= 0) and (c+m > 0):
                raise ValueError('c and m must be nonnegative integers summing greater than zero')
            self._c = c
            self._m = m
    def simulate(self, max_gen=20):
        """simulate the number of clone leaves and mutant clades off a root node"""
        if self._p>=.5:
            warnings.warn('p >= .5 is not subcritical, tree simulations not garanteed to terminate')
        if self._p is None or self._q is None:
            raise ValueError('p and q parameters must be defined for simulation\n')
        if max_gen == 0:
            raise RuntimeError('tree not terminated before max_gen generations\n')
        # if there are no children, there is just one leaf and no mutant clades
        if random.random() > self._p:
            self._c = 1
            self._m = 0
            return self
        clone_leaves = 0
        mutant_clades = 0
        # if not a leaf, there are two children
        for child in range(2):
            # if mutant, it just counts as a clade
            if random.random() < self._q:
                mutant_clades += 1
            else:
                LeavesAndClades.simulate(self, max_gen=max_gen-1)
                new_clone_leaves = self._c
                new_mutant_clades = self._m
                if new_clone_leaves is not None and new_mutant_clades is not None:
                    clone_leaves += new_clone_leaves
                    mutant_clades += new_mutant_clades
                else:
                    clone_leaves = mutant_clades = None
                    break
        self._c = clone_leaves
        self._m = mutant_clades
        return self

    f_hash = {} # <--- class variable for hashing these
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
                return (0., (0., 0.))
            if c==1 and m==0:
                return (1-p, (-1, 0.))
            if c==1 and m==1:
                return (2*p*q*(1-q), (2*q*(1-q), 2*p - 4*p*q))
            if c==0 and m==2:
                return (p*q**2, (q**2, 2*p*q))
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
            CollapsedTree.f_hash[(p, q, c, m)] = (f_result, scipy.array([dfdp_result, dfdq_result]))
        return CollapsedTree.f_hash[(p, q, c, m)]

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
        if p is not None and q is not None:
            LeavesAndClades.__init__(self, p=p, q=q)
        else:
            if tree is None:
                raise ValueError('either p and q or tree (or all three) must be provided')
        # check that tree is valid
        if tree is not None:
            k = len(tree)
            if k==0 or \
               set(map(len, tree))!=set([2]) or \
               not all(scipy.greater_equal([x for y in tree for x in y], 0)):
                raise ValueError('"tree" must be a nonempty list of 2-element tuples of nonnegative integers')
            cs, ms = zip(*tree)
            if not all(scipy.greater(scipy.cumsum([x[1] for x in tree])[:-1], scipy.arange(1, k)-1)):
                raise ValueError('inconsistent breadth first tree data')
        self._tree = tree


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
        

def main():
    """some tests"""
    p = random.random()/2
    q = random.random()
    n_trees = 100
    print 'initializing a forest of %d trees with random (subcritical) parameters' % n_trees
    forest = CollapsedForest(p, q, n_trees)
    print 'true parameters: p = %f, q = %f' % (p, q)

    print 'simulating the forest'
    forest.simulate()

    print 'computing MLE'
    mle = forest.mle().x
    print 'MLE parameters:  p = %f, q = %f' % tuple(mle.tolist())

    # plot the 2-norm of the difference between the gradient and a finite difference approximation
    # ideally small everywhere
    print 'computing plot data (may take a sec)...'
    X, Y = scipy.mgrid[slice(.02, 1, .02),
                       slice(.02, 1, .02)]
    Z = scipy.zeros((X.shape[0], X.shape[1]))
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            Z[i, j] = check_grad(lambda x: forest.l(x)[0], lambda x: forest.l(x)[1], (X[i, j], Y[i, j]))
    fig = plt.figure()
    plt.rc('text', usetex=True)
    plt.title(r'$||\nabla \ell(p, q) - \Delta \ell(p, q)||_2$')
    plt.pcolor(X, Y, Z)
    plt.xlabel(r'$p$')
    plt.ylabel(r'$q$')
    plt.axes().set_aspect('equal')
    plt.colorbar()
    plt.savefig('foo.pdf')


    # plot likelihood surface
    X, Y = scipy.mgrid[slice(.02, 1, .02),
                       slice(.02, 1, .02)]
    Z = scipy.zeros((X.shape[0], X.shape[1]))
    #U = scipy.zeros((X.shape[0], X.shape[1]))
    #V = scipy.zeros((X.shape[0], X.shape[1]))
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            l, grad_l = forest.l((X[i, j], Y[i, j]))
            z = l #scipy.exp(l)
            Z[i, j] = z
    #        U[i, j] = z*grad_l[0]
    #        V[i, j] = z*grad_l[1]
    fig = plt.figure()
    plt.rc('text', usetex=True)
    plt.title(r'$\ell(p, q)$')
    #plt.title(r'$P(T \mid p, q), \nabla_{p q} P(T \mid p, q)$')
    contour = plt.contour(X, Y, Z, colors='k', label='likelihood contours')#cmap='Blues', vmin=vmin, vmax=1.1*vmax)
    for c in contour.collections:
        c.set_linestyle('solid')

    plt.clabel(contour, fontsize=9, inline=1)
    #plt.quiver(X, Y, U, V, units='x', color='w',edgecolors=('grey'))
    plt.plot([p], [q], 'k+', label='true parameters')
    plt.plot(mle[0], mle[1], 'ko', markerfacecolor='none', label='MLE parameters')
    plt.xlabel(r'$p$')
    plt.ylabel(r'$q$')
    plt.axes().set_aspect('equal')
    plt.legend(numpoints = 1)
    plt.savefig('bar.pdf')


if __name__ == "__main__":
    main()
