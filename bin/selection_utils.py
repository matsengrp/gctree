#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
utility functions for selection simulation
'''

from __future__ import division, print_function

import scipy
from scipy.optimize import minimize, fsolve
import matplotlib; matplotlib.use('agg')
from matplotlib import pyplot as plt
from utils import hamming_distance

def calc_Kd(seqAA, targetAAseqs, hd2affy):
    '''Find the closest target sequence to and apply the "hamming distance to affinity" transformation function.'''
    hd = min([hamming_distance(seqAA, t) for t in targetAAseqs])
    return(hd2affy(hd))

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
        return lambda A: (A_total - (A + scipy.sum(B_total/(1+Kd_n/A))))**2

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
        lambda_ = alpha + (2 - alpha) / (1 + Q*scipy.exp(-beta*BA))
        return(lambda_)

    # Update the list of affinities for all the live nodes:
    Kd_n = scipy.array([n.Kd for n in tree.iter_leaves() if not n.terminated])
    BnA = calc_binding_time(Kd_n, A_total, B_total)
    lambdas = trans_BA(BnA, Lp)
    i = 0
    for n in tree.iter_leaves():
        if n.terminated:
            continue
        n.add_feature('lambda_', lambdas[i])
        i += 1
    return(tree)

def find_A_total(carry_cap, B_total, f_full, mature_affy, U):
    def A_total_fun(A, B_total, Kd_n): return(A + scipy.sum(B_total/(1+Kd_n/A)))

    def C_A(A, A_total, f_full, U): return(U * (A_total - A) / f_full)

    def A_obj(carry_cap, B_total, f_full, Kd_n, U):
        def obj(A): return((carry_cap - C_A(A, A_total_fun(A, B_total, Kd_n), f_full, U))**2)
        return(obj)

    Kd_n = scipy.array([mature_affy] * carry_cap)
    obj = A_obj(carry_cap, B_total, f_full, Kd_n, U)
    # Some funny "zero encountered in true_divide" errors are not affecting results so ignore them:
    old_settings = scipy.seterr(all='ignore')  # Keep old settings
    scipy.seterr(divide='ignore')
    obj_min = minimize(obj, 1e-20, bounds=[[1e-20, carry_cap]], method='L-BFGS-B', tol=1e-20)
    scipy.seterr(**old_settings)  # Reset to default
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
        lambda_ = alpha + (2 - alpha) / (1 + Q*scipy.exp(-beta*BA))
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

    # FloatingPointError errors are not affecting results so ignore them:
    old_settings = scipy.seterr(all='ignore')  # Keep old settings
    scipy.seterr(over='ignore')
    try:
        def obj_T_A(p): return(solve_T_BA(p, f_full, U))
        p = fsolve(obj_T_A, (0, 10e-5, 1), xtol=1e-20, maxfev=1000)
        assert(sum(solve_T_BA(p, f_full, U)) < f_full * 1/1000)
    except:
        print('The U parameter is large and therefore the epsilon parameter has to be adjusted to find a valid solution.')
        def obj_T_A(p): return(solve_T_BA_low_epsilon(p, f_full, U))
        p = fsolve(obj_T_A, (0, 10e-5, 1), xtol=1e-20, maxfev=1000)
        assert(sum(solve_T_BA(p, f_full, U)) < f_full * 1/1000)
    scipy.seterr(**old_settings)  # Reset to default
    return(p)


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
    pop_size = scipy.array([sum(r[0]) for r in runstats])
    # min:max of the hamming distances to plot:
    bounds = make_bounds(runstats)

    fig = plt.figure()
    ax = plt.subplot(111)
    t = scipy.array(list(range(len(pop_size))))  # The x-axis are generations
    ax.plot(t, pop_size, lw=2, label='All cells')  # Total population size is plotted
    # Then plot the counts for each hamming distance as a function on generation:
    for k in list(range(*bounds)):
        color = colors[k]
        ax.plot(t, scipy.array([r[0][k] for r in runstats]), lw=2, color=color, label='Dist {}'.format(k))

    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)

    # Shrink current axis by 20% to make the legend fit:
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    plt.ylabel('Count')
    plt.xlabel('GC generation')
    plt.title('Cell count as function of GC generation')
    fig.savefig(outbase + '.selection_sim.runstats.pdf')
