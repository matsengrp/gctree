#! /bin/env python

import sys, scipy, time, random
from scipy.optimize import minimize, check_grad
import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
from matplotlib import rc

# here's a base class for simulating a binary ininite type branching process with branching
# probability p and mutation probability q, and we collapse mutant clades off the root type
# and consider just the number of clone leaves and mutant clades.
#
#   /\            
#  /\ ^          (3)
#   /\     ==>   / \
#    /\
#     ^
#
class leavesAndClades():
    def __init__(self, p=None, q=None):
        assert 0 < p < 1 and 0 < q < 1
        self._p = p
        self._q = q
    # simulate the number of clone leaves and mutant clades off a root node
    def rootSim(self, maxGen=20):
        if self._p is None or self._q is None:
            sys.stderr.write('ERROR: p and q parameters must be defined for simulation\n')
            return
        if maxGen == 0:
            sys.stderr.write('tree not terminated before maxGen generations\n')
            return None, None
        # if there are no children, there is just one leaf and no mutant clades
        if random.random() > self._p:
            return (1, 0)
        cloneLeaves = 0
        mutantClades = 0
        # if not a leaf, there are two children
        for child in range(2):
            # if mutant, it just counts as a clade
            if random.random() < self._q:
                mutantClades += 1
            else:
                newCloneLeaves, newMutantClades = self.rootSim(maxGen=maxGen-1)
                if newCloneLeaves is not None and newMutantClades is not  None:
                    cloneLeaves += newCloneLeaves
                    mutantClades += newMutantClades
                else:
                    cloneLeaves = mutantClades = None
                    break
        return (cloneLeaves, mutantClades)


# Here's a derived class for a collapsed tree, where we recurse into the mutant clades
#
#      (4)
#     / | \
#   (3)(1)(2)
#       |   \
#      (2)  (1)
#
class collapsedTree(leavesAndClades):
    #Clonal leaf count and count of mutant clades are provided in
    # breadth first order, assumes inputs are valid.
    def __init__(self, p=None, q=None, treeData=None):
        leavesAndClades.__init__(self, p=p, q=q)
        # check that treeData valid
        if treeData is not None:
            cs, ms = zip(*treeData)
            k = len(cs)
            assert k >= 1 and len(ms) == k
            assert all(scipy.greater_equal(cs, 0)) and all(scipy.greater_equal(ms, 0))
            assert all(scipy.greater(scipy.cumsum(ms)[:-1], scipy.arange(1, k)-1))
        self._treeData = treeData

    # probability of getting c leaves that are clones of the root and m mutant clades off
    # the root line, dynamic programming
    fHash = {} # <--- class variable for hashing these
    @staticmethod
    def __f(p, q, c, m):
        if (p, q, c, m) not in collapsedTree.fHash:
            if c==m==0 or (c==0 and m==1):
                return 0.
            if c==1 and m==0:
                return 1-p
            if c==1 and m==1:
                return 2*p*q*(1-q)
            if c==0 and m==2:
                return p*q**2
            if m >= 1:
                result = 2*p*q*(1-q)*collapsedTree.__f(p, q, c, m-1)
            else:
                result = 0.
            summation = 0.
            for cx in range(c+1):
                for mx in range(m+1):
                    if (not (cx==0 and mx==0)) and (not (cx==c and mx==m)):
                        summation += collapsedTree.__f(p, q, cx, mx)*collapsedTree.__f(p, q, c-cx, m-mx)
            collapsedTree.fHash[(p, q, c, m)] = result + p*(1-q)**2*summation
        return collapsedTree.fHash[(p, q, c, m)]

    # derivatives
    dfdpHash = {}
    @staticmethod
    def __dfdp(p, q, c, m):
        if (p, q, c, m) not in collapsedTree.dfdpHash:
            if c==m==0 or (c==0 and m==1):
                return 0.
            if c==1 and m==0:
                return -1
            if c==1 and m==1:
                return 2*q*(1-q)
            if c==0 and m==2:
                return q**2
            if m >= 1:
                result = 2*p*q*(1-q)*collapsedTree.__dfdp(p, q, c, m-1) + \
                         2*q*(1-q)*collapsedTree.__f(p, q, c, m-1)
            else:
                result = 0.
            summation = 0.
            for cx in range(c+1):
                for mx in range(m+1):
                    if (not (cx==0 and mx==0)) and (not (cx==c and mx==m)):
                        summation += collapsedTree.__dfdp(p, q, cx, mx)* \
                                      collapsedTree.__f(p, q, c-cx, m-mx)*p*(1-q)**2 \
                                   + collapsedTree.__f(p, q, cx, mx)* \
                                      collapsedTree.__dfdp(p, q, c-cx, m-mx)*p*(1-q)**2 \
                                   + collapsedTree.__f(p, q, cx, mx)* \
                                      collapsedTree.__f(p, q, c-cx, m-mx)*(1-q)**2
            collapsedTree.dfdpHash[(p, q, c, m)] = result + summation
        return collapsedTree.dfdpHash[(p, q, c, m)]

    dfdqHash = {}
    @staticmethod
    def __dfdq(p, q, c, m):
        if (p, q, c, m) not in collapsedTree.dfdqHash:
            if c==m==0 or (c==0 and m==1):
                return 0.
            if c==1 and m==0:
                return 0
            if c==1 and m==1:
                return 2*p - 4*p*q
            if c==0 and m==2:
                return 2*p*q
            if m >= 1:
                result = (2*p - 4*p*q)*collapsedTree.__f(p, q, c, m-1) + \
                         2*p*q*(1-q)*collapsedTree.__dfdq(p, q, c, m-1)
            else:
                result = 0.
            summation = 0.
            for cx in range(c+1):
                for mx in range(m+1):
                    if (not (cx==0 and mx==0)) and (not (cx==c and mx==m)):
                        summation += collapsedTree.__dfdq(p, q, cx, mx)* \
                                      collapsedTree.__f(p, q, c-cx, m-mx)*p*(1-q)**2 \
                                   + collapsedTree.__f(p, q, cx, mx)* \
                                      collapsedTree.__dfdq(p, q, c-cx, m-mx)*p*(1-q)**2 \
                                   + collapsedTree.__f(p, q, cx, mx)* \
                                      collapsedTree.__f(p, q, c-cx, m-mx)*p*2*(q-1)
            collapsedTree.dfdqHash[(p, q, c, m)] = result + summation
        return collapsedTree.dfdqHash[(p, q, c, m)]

    # likelihood of collapsed tree. We also return the gradient wrt p, q
    def logL(self, (p, q), sign=1):
        # sign = -1 is useful for maximizing with scipy's minimize
        # let's avoid redundancy
        fs = scipy.array([collapsedTree.__f(p, q, c, m) for (c, m) in self._treeData])
        dfdps = scipy.array([collapsedTree.__dfdp(p, q, c, m) for (c, m) in self._treeData])
        dfdqs = scipy.array([collapsedTree.__dfdq(p, q, c, m) for (c, m) in self._treeData])
        return sign*scipy.log(fs).sum(), sign*scipy.array([(dfdps/fs).sum(), (dfdqs/fs).sum()])

    # MLE for p and q given tree
    def MLE(self):
        if self._treeData is None:
            sys.stderr.write('ERROR: tree must be defined for simulation\n')
            return
        # initalize optimization
        x0 = (random.random(), random.random())
        bounds = ((.001, .999), (.001, .999))
        result = minimize(self.logL, x0=x0, args=(-1,), jac=True, method='L-BFGS-B', bounds=bounds)
        # update p and q if None
        if not result.success:
            sys.stderr.write('ERROR: optiomization failed\n')
        elif self._p is None and self._q is None:
            self._p, self._q = result.x
        return result

    # simulate a tree
    def simulate(self):
        if self._p is None or self._q is None:
            sys.stderr.write('ERROR: p and q parameters must be defined for simulation\n')
            return
        # initialize with root
        tree = [self.rootSim()]
        moreMutants = sum(x[1] for x in tree) # mutant clades off root
        while moreMutants > 0:
            newNodes = [self.rootSim() for m in range(moreMutants)]
            moreMutants = sum(x[1] for x in newNodes) # mutant clades from this generation
            tree += newNodes
        #if self._treeData is None:
        # replace tree data
        self._treeData = tree
        return self
                
    def getparams(self):
        return {'p':self._p, 'q':self._q, 'tree':self._treeData}

    def __str__(self):
        return 'p = %f, q = %f\ntree: ' % (self._p, self._q) + str(self._treeData)
        
        
def main():
    # initialize
    p, q = [float(x) for x in sys.argv[1:]]
    tree = collapsedTree(p, q)
    MLEhash = {} # hash MLE results to avoid redundant comp
    MLEs = []
    ct = 0
    while ct < 10000:
        tree.simulate()
        if str(tree) not in MLEhash:
            MLEhash[str(tree)] = tree.MLE().x
        MLEs.append(MLEhash[str(tree)])
        ct += 1
        sys.stderr.write('simulate 10000 trees and compute MLE for each: %d  \r' % (ct))
    ps, qs = zip(*MLEs)
    sys.stderr.write('\n')

    # plot dis
    fig = plt.figure()
    plt.hist2d(ps, qs, bins=50, range=[[0,1], [0,1]], cmap='Greys')
    plt.xlim([0,1])
    plt.ylim([0,1])
    plt.xlabel(r'$p$')
    plt.ylabel(r'$q$')
    plt.plot([p], [q], 'ro')
    plt.axes().set_aspect('equal')
    plt.savefig('foo.pdf')

    ##this checks that the derivative is accurate, should print a small number
    #print 'gradient check:', check_grad(lambda x: tree.logL(x)[0], lambda x: tree.logL(x)[1], (p, q))

    print tree

    # plot MLE for last simulation
    X, Y = scipy.mgrid[slice(.02, 1, .02),
                       slice(.02, 1, .02)]
    Z = scipy.zeros((X.shape[0], X.shape[1]))
    U = scipy.zeros((X.shape[0], X.shape[1]))
    V = scipy.zeros((X.shape[0], X.shape[1]))
    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            L, gradL = tree.logL((X[i, j], Y[i, j]))
            z = scipy.exp(L)
            Z[i, j] = z
            U[i, j] = z*gradL[0]
            V[i, j] = z*gradL[1]
    fig = plt.figure()
    plt.rc('text', usetex=True)
    plt.title(r'$P(T \mid p, q), \nabla_{p q} P(T \mid p, q)$')
    plt.contour(X, Y, Z, 10, colors='0.5')#cmap='Blues', vmin=vmin, vmax=1.1*vmax)
    plt.quiver(X, Y, U, V, units='x', color='w',edgecolors=('grey'))
    plt.plot([p], [q], 'ro')
    plt.plot([ps[-1]], [qs[-1]], 'bo')
    plt.xlabel(r'$p$')
    plt.ylabel(r'$q$')
    plt.axes().set_aspect('equal')
    plt.show()
    plt.savefig('bar.pdf')


if __name__ == "__main__":
    main()
