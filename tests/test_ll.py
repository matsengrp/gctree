import gctree.branching_processes as bp
import gctree.phylip_parse as pp
import numpy as np
from scipy.special import logsumexp, softmax
from functools import lru_cache


trees = pp.parse_outfile('tests/small_outfile', root='GL', abundance_file='tests/abundances.csv')
forest = bp.CollapsedForest([bp.CollapsedTree(tree) for tree in trees])

def test_ll_genotype_cache():
    p, q = 0.4, 0.6
    c_max, m_max = 10, 10
    bp.CollapsedTree._build_ll_genotype_cache(c_max, m_max, p, q)

def test_forest_likelihoods():
    # make sure ll doesn't have underflow
    p, q = 0.4, 0.6
    forest.ll(p, q, marginal=True)

def validate_ll_genotype():
    c_max, m_max = 10, 10
    for params in [(0.4, 0.6), (0.3, 0.5)]:
        for c in range(c_max):
            for m in range(m_max):
                if c > 0 or m > 1:
                    true_res = ll_validate(c, m, *params)
                    res = bp.CollapsedTree._ll_genotype(c, m, *params)
                    assert np.isclose(true_res[0], res[0])
                    assert np.isclose(true_res[1][0], res[1][0])
                    assert np.isclose(true_res[1][1], res[1][1])


@lru_cache(maxsize=None)
def ll_validate(c, m, p, q):
    r"""Log-probability of getting :math:`c` leaves that are clones of the root and
    :math:`m` mutant clades off the root lineage, given branching probability :math:`p` and
    mutation probability :math:`q`.
    AKA the spaceship distribution. Also returns gradient wrt p and q
    (p, q). Computed by dynamic programming.
    Args:
        c: clonal leaves
        m: mutant clades
        p: branching probability
        q: mutation probability
    Returns:
        log-likelihood and gradient wrt :math:`p` and :math:`q`.
    """
    if c == m == 0 or (c == 0 and m == 1):
        logf_result = -np.inf
        dlogfdp_result = 0
        dlogfdq_result = 0
    elif c == 1 and m == 0:
        logf_result = np.log(1 - p)
        dlogfdp_result = -1 / (1 - p)
        dlogfdq_result = 0
    elif c == 0 and m == 2:
        logf_result = np.log(p) + 2 * np.log(q)
        dlogfdp_result = 1 / p
        dlogfdq_result = 2 / q
    else:
        if m >= 1:
            (
                neighbor_ll_genotype,
                (
                    neighbor_dlogfdp,
                    neighbor_dlogfdq,
                ),
            ) = ll_validate(c, m - 1, p, q)
            logf_result = (
                np.log(2) + np.log(p) + np.log(q) + np.log(1 - q) + neighbor_ll_genotype
            )
            dlogfdp_result = 1 / p + neighbor_dlogfdp
            dlogfdq_result = 1 / q - 1 / (1 - q) + neighbor_dlogfdq
        else:
            logf_result = -np.inf
            dlogfdp_result = 0.0
            dlogfdq_result = 0.0
        logg_array = [logf_result]
        dloggdp_array = [dlogfdp_result]
        dloggdq_array = [dlogfdq_result]
        for cx in range(c + 1):
            for mx in range(m + 1):
                if (not (cx == mx == 0)) and (not (cx == c and mx == m)):
                    (
                        neighbor1_ll_genotype,
                        (neighbor1_dlogfdp, neighbor1_dlogfdq),
                    ) = ll_validate(cx, mx, p, q)
                    (
                        neighbor2_ll_genotype,
                        (neighbor2_dlogfdp, neighbor2_dlogfdq),
                    ) = ll_validate(c - cx, m - mx, p, q)
                    logg = (
                        np.log(p)
                        + 2 * np.log(1 - q)
                        + neighbor1_ll_genotype
                        + neighbor2_ll_genotype
                    )
                    dloggdp = 1 / p + neighbor1_dlogfdp + neighbor2_dlogfdp
                    dloggdq = -2 / (1 - q) + neighbor1_dlogfdq + neighbor2_dlogfdq
                    logg_array.append(logg)
                    dloggdp_array.append(dloggdp)
                    dloggdq_array.append(dloggdq)
        logf_result = logsumexp(logg_array)
        softmax_logg_array = softmax(logg_array)
        dlogfdp_result = np.multiply(softmax_logg_array, dloggdp_array).sum()
        dlogfdq_result = np.multiply(softmax_logg_array, dloggdq_array).sum()

    return (logf_result, np.array([dlogfdp_result, dlogfdq_result]))
