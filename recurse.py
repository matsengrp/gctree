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
    This is a base class for simulating a binary infinite type branching process with branching
    probability p, mutation probability q, and we collapse mutant clades off the root type
    and consider just the number of clone leaves and mutant clades.

      /\            
     /\ ^          (3)
      /\     ==>   / \\
       /\\
        ^
    """
    def __init__(self, p=None, q=None):
        """initialize with branching probability p and mutation probability q, both in the unit interval"""
        if not (0 < p < 1 and 0 < q < 1):
            raise ValueError('p and q must be in the unit interval')
        if p>=.5:
            warnings.warn('p >= .5 is not subcritical, tree simulations not garanteed to terminate')
        self._p = p
        self._q = q
    def simulate(self, max_gen=20):
        """simulate the number of clone leaves and mutant clades off a root node"""
        if self._p is None or self._q is None:
            raise ValueError('p and q parameters must be defined for simulation\n')
        if max_gen == 0:
            raise RuntimeError('tree not terminated before max_gen generations\n')
        # if there are no children, there is just one leaf and no mutant clades
        if random.random() > self._p:
            return (1, 0)
        clone_leaves = 0
        mutant_clades = 0
        # if not a leaf, there are two children
        for child in range(2):
            # if mutant, it just counts as a clade
            if random.random() < self._q:
                mutant_clades += 1
            else:
                new_clone_leaves, new_mutant_clades = LeavesAndClades.simulate(self, max_gen=max_gen-1)
                if new_clone_leaves is not None and new_mutant_clades is not  None:
                    clone_leaves += new_clone_leaves
                    mutant_clades += new_mutant_clades
                else:
                    clone_leaves = mutant_clades = None
                    break
        return (clone_leaves, mutant_clades)

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

    f_hash = {} # <--- class variable for hashing these
    @staticmethod
    def __f(p, q, c, m):
        """
        probability of getting c leaves that are clones of the root and m mutant clades off
        the root line, given branching probability p and mutation probability q 
        """
        if (p, q, c, m) not in CollapsedTree.f_hash:
            if c==m==0 or (c==0 and m==1):
                return 0.
            if c==1 and m==0:
                return 1-p
            if c==1 and m==1:
                return 2*p*q*(1-q)
            if c==0 and m==2:
                return p*q**2
            if m >= 1:
                result = 2*p*q*(1-q)*CollapsedTree.__f(p, q, c, m-1)
            else:
                result = 0.
            summation = 0.
            for cx in range(c+1):
                for mx in range(m+1):
                    if (not (cx==0 and mx==0)) and (not (cx==c and mx==m)):
                        summation += CollapsedTree.__f(p, q, cx, mx)*CollapsedTree.__f(p, q, c-cx, m-mx)
            CollapsedTree.f_hash[(p, q, c, m)] = result + p*(1-q)**2*summation
        return CollapsedTree.f_hash[(p, q, c, m)]

    dfdp_hash = {}
    @staticmethod
    def __dfdp(p, q, c, m):
        """derivative of f wrt p"""
        if (p, q, c, m) not in CollapsedTree.dfdp_hash:
            if c==m==0 or (c==0 and m==1):
                return 0.
            if c==1 and m==0:
                return -1
            if c==1 and m==1:
                return 2*q*(1-q)
            if c==0 and m==2:
                return q**2
            if m >= 1:
                result = 2*p*q*(1-q)*CollapsedTree.__dfdp(p, q, c, m-1) + \
                         2*q*(1-q)*CollapsedTree.__f(p, q, c, m-1)
            else:
                result = 0.
            summation = 0.
            for cx in range(c+1):
                for mx in range(m+1):
                    if (not (cx==0 and mx==0)) and (not (cx==c and mx==m)):
                        summation += CollapsedTree.__dfdp(p, q, cx, mx)* \
                                      CollapsedTree.__f(p, q, c-cx, m-mx)*p*(1-q)**2 \
                                   + CollapsedTree.__f(p, q, cx, mx)* \
                                      CollapsedTree.__dfdp(p, q, c-cx, m-mx)*p*(1-q)**2 \
                                   + CollapsedTree.__f(p, q, cx, mx)* \
                                      CollapsedTree.__f(p, q, c-cx, m-mx)*(1-q)**2
            CollapsedTree.dfdp_hash[(p, q, c, m)] = result + summation
        return CollapsedTree.dfdp_hash[(p, q, c, m)]

    dfdq_hash = {}
    @staticmethod
    def __dfdq(p, q, c, m):
        """derivative of f wrt q"""
        if (p, q, c, m) not in CollapsedTree.dfdq_hash:
            if c==m==0 or (c==0 and m==1):
                return 0.
            if c==1 and m==0:
                return 0
            if c==1 and m==1:
                return 2*p - 4*p*q
            if c==0 and m==2:
                return 2*p*q
            if m >= 1:
                result = (2*p - 4*p*q)*CollapsedTree.__f(p, q, c, m-1) + \
                         2*p*q*(1-q)*CollapsedTree.__dfdq(p, q, c, m-1)
            else:
                result = 0.
            summation = 0.
            for cx in range(c+1):
                for mx in range(m+1):
                    if (not (cx==0 and mx==0)) and (not (cx==c and mx==m)):
                        summation += CollapsedTree.__dfdq(p, q, cx, mx)* \
                                      CollapsedTree.__f(p, q, c-cx, m-mx)*p*(1-q)**2 \
                                   + CollapsedTree.__f(p, q, cx, mx)* \
                                      CollapsedTree.__dfdq(p, q, c-cx, m-mx)*p*(1-q)**2 \
                                   + CollapsedTree.__f(p, q, cx, mx)* \
                                      CollapsedTree.__f(p, q, c-cx, m-mx)*p*2*(q-1)
            CollapsedTree.dfdq_hash[(p, q, c, m)] = result + summation
        return CollapsedTree.dfdq_hash[(p, q, c, m)]

    def l(self, (p, q), sign=1):
        """
        log likelihood of p and q, conditioned on collapsed tree, and its gradient wrt (p, q)
        optional parameter sign must be 1 or -1, with the latter useful for MLE by minimization
        """
        if self._tree is None:
            raise ValueError('tree data must be defined to compute likelihood')
        if sign not in (-1, 1):
            raise ValueError('sign must be 1 or -1')
        fs = scipy.array([CollapsedTree.__f(p, q, c, m) for (c, m) in self._tree])
        dfdps = scipy.array([CollapsedTree.__dfdp(p, q, c, m) for (c, m) in self._tree])
        dfdqs = scipy.array([CollapsedTree.__dfdq(p, q, c, m) for (c, m) in self._tree])
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
        tree = [LeavesAndClades.simulate(self)] # initialize with root
        more_mutants = sum(x[1] for x in tree) # mutant clades off root
        while more_mutants > 0:
            new_nodes = [LeavesAndClades.simulate(self) for m in range(more_mutants)]
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
    n_trees = 50
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
    X, Y = scipy.mgrid[slice(.05, 1, .05),
                       slice(.05, 1, .05)]
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
            z = scipy.exp(l)
            Z[i, j] = z
    #        U[i, j] = z*grad_l[0]
    #        V[i, j] = z*grad_l[1]
    fig = plt.figure()
    plt.rc('text', usetex=True)
    plt.title(r'$\ell(p, q)$')
    #plt.title(r'$P(T \mid p, q), \nabla_{p q} P(T \mid p, q)$')
    plt.contour(X, Y, Z, 10, colors='0.5')#cmap='Blues', vmin=vmin, vmax=1.1*vmax)
    #plt.quiver(X, Y, U, V, units='x', color='w',edgecolors=('grey'))
    plt.plot([p], [q], 'ro')
    plt.plot(mle[0], mle[1], 'bo')
    plt.xlabel(r'$p$')
    plt.ylabel(r'$q$')
    plt.axes().set_aspect('equal')
    plt.show()
    plt.savefig('bar.pdf')


if __name__ == "__main__":
    main()
