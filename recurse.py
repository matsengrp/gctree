#! /bin/env python

import sys, scipy, time
from scipy.optimize import minimize, check_grad
import matplotlib
matplotlib.use('PDF')
from matplotlib import pyplot as plt
from matplotlib import rc

# probability of getting c leaves that are clones of the root and m mutant clades off
# the root type line
# NOTE: this could be used recursively on the mutant clades to compute the probability
#       of c, m1, m2, ...., i.e. all the mutant leaf abundances
def f(c, m, p, q):
    if c==m==0 or (c==0 and m==1):
        return 0.
    if c==1 and m==0:
        return 1-p
    if c==1 and m==1:
        return 2*p*q*(1-q)
    if c==0 and m==2:
        return p*q**2
    if m >= 1:
        result = 2*p*q*(1-q)*f(c, m-1, p, q)
    else:
        result = 0.
    Sum = 0.
    for cx in range(c+1):
        for mx in range(m+1):
            if (not (cx==0 and mx==0)) and (not (cx==c and mx==m)):
                Sum += f(cx, mx, p, q)*f(c-cx, m-mx, p, q)
    return result + p*(1-q)**2*Sum

# derivatives
def dfdp(c, m, p, q):
    if c==m==0 or (c==0 and m==1):
        return 0.
    if c==1 and m==0:
        return -1
    if c==1 and m==1:
        return 2*q*(1-q)
    if c==0 and m==2:
        return q**2
    if m >= 1:
        result = 2*p*q*(1-q)*dfdp(c, m-1, p, q) + 2*q*(1-q)*f(c, m-1, p, q)
    else:
        result = 0.
    Sum = 0.
    for cx in range(c+1):
        for mx in range(m+1):
            if (not (cx==0 and mx==0)) and (not (cx==c and mx==m)):
                Sum += dfdp(cx, mx, p, q)*   f(c-cx, m-mx, p, q)*p*(1-q)**2 \
                    +     f(cx, mx, p, q)*dfdp(c-cx, m-mx, p, q)*p*(1-q)**2 \
                    +     f(cx, mx, p, q)*   f(c-cx, m-mx, p, q)*(1-q)**2
    return result + Sum
    
def dfdq(c, m, p, q):
    if c==m==0 or (c==0 and m==1):
        return 0.
    if c==1 and m==0:
        return 0
    if c==1 and m==1:
        return 2*p - 4*p*q
    if c==0 and m==2:
        return 2*p*q
    if m >= 1:
        result = (2*p - 4*p*q)*f(c, m-1, p, q) + 2*p*q*(1-q)*dfdq(c, m-1, p, q)
    else:
        result = 0.
    Sum = 0.
    for cx in range(c+1):
        for mx in range(m+1):
            if (not (cx==0 and mx==0)) and (not (cx==c and mx==m)):
                Sum += dfdq(cx, mx, p, q)*   f(c-cx, m-mx, p, q)*p*(1-q)**2 \
                    +     f(cx, mx, p, q)*dfdq(c-cx, m-mx, p, q)*p*(1-q)**2 \
                    +     f(cx, mx, p, q)*   f(c-cx, m-mx, p, q)*p*2*(q-1)
    return result + Sum


# Ok now we recurse this to compute the probability of a collapsed tree, where each node has a count of
# clonal leaves, and a count of mutant clades, and are provided in breadth first order
# assumes inputs are valid
# we also return the gradient wrt p, q
def logL((p, q), cs, ms, sign=1):
    # sign = -1 is useful for maximizing with scipy's minimize
    # let's avoid redundancy
    fs = scipy.array([f(c, m, p, q) for (c, m) in zip(cs, ms)])
    dfdps = scipy.array([dfdp(c, m, p, q) for (c, m) in zip(cs, ms)])
    dfdqs = scipy.array([dfdq(c, m, p, q) for (c, m) in zip(cs, ms)])
    return sign*scipy.log(fs).sum(), sign*scipy.array([(dfdps/fs).sum(), (dfdqs/fs).sum()])

# get comma-separated pairs of c and m from command line args
cms = [tuple(int(x) for x in x.split(',')) for x in sys.argv[1:]]
cs, ms = zip(*cms)

# check that cs and ms are valid
k = len(cs)
assert k >= 1 and len(ms) == k
assert all(scipy.greater_equal(cs, 0)) and all(scipy.greater_equal(ms, 0))
assert all(scipy.greater(scipy.cumsum(ms)[:-1], scipy.arange(1, k)-1))

# initalize
x0 = (.5, .5)
print 'initial parameters: p0=%f, q0=%f' % x0

#this checks that the derivative is accurate, should print a small number
print 'gradient check:', check_grad(lambda x: logL(x, cs, ms)[0], lambda x: logL(x, cs, ms)[1], x0)

print '\nminimization results:'
t = time.time()
result = minimize(logL, x0=x0, args=(cs, ms, -1), jac=True, method='L-BFGS-B', bounds=((.01, .99), (.01, .99)))
print result
print '    elapsed time:', time.time()-t, 'sec'

# plot MLE
X, Y = scipy.mgrid[slice(.02, 1, .02),
                   slice(.02, 1, .02)]
Z = scipy.zeros((X.shape[0], X.shape[1]))
U = scipy.zeros((X.shape[0], X.shape[1]))
V = scipy.zeros((X.shape[0], X.shape[1]))
for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        L, gradL = logL((X[i, j], Y[i, j]), cs, ms)
        z = scipy.exp(L)
        Z[i, j] = scipy.exp(z)
        U[i, j] = z*gradL[0]
        V[i, j] = z*gradL[1]
#vmin, vmax = Z.min(), Z.max()
fig = plt.figure()
plt.rc('text', usetex=True)
plt.title(r'$P(T \mid p, q), \nabla_{p q} P(T \mid p, q)$')
plt.contour(X, Y, Z, 10, colors='0.5')#cmap='Blues', vmin=vmin, vmax=1.1*vmax)
#plt.plot((x0[0], result.x[0]), (x0[1], result.x[1]), 'ko')
plt.annotate('initial guess', x0, xytext=[x+.2 for x in x0],
            arrowprops=dict(facecolor='k', width=.01, headwidth=5),
            )
plt.quiver(X, Y, U, V, units='x', color='w',edgecolors=('grey'))
#cbar = plt.colorbar()
#cbar.ax.set_ylabel('likelihood', rotation=270)
plt.xlabel(r'$p$')
plt.ylabel(r'$q$')
plt.annotate('MLE', result.x, xytext=[x+.3 for x in result.x],
            arrowprops=dict(facecolor='k', width=.01, headwidth=5),
            )
plt.show()
plt.savefig('foo.pdf')

