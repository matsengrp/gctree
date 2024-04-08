import gctree.branching_processes as bp
import gctree.phylip_parse as pp
import gctree.utils as utils
import gctree.mutation_model as mm
from math import log

import numpy as np
from multiset import FrozenMultiset
from oldcode import OldCollapsedTree, OldCollapsedForest


def make_oldctree(tree):
    """Make an old CollapsedTree from an hDAG clade tree"""
    etetree = tree.to_ete(
        name_func=lambda n: n.attr["name"],
        features=["sequence"],
        feature_funcs={"abundance": lambda n: n.label.abundance},
    )
    for node in etetree.traverse():
        if not node.is_leaf():
            node.abundance = 0
    etetree.dist = 0
    for node in etetree.iter_descendants():
        node.dist = utils.hamming_distance(node.up.sequence, node.sequence)
    return OldCollapsedTree(etetree)


def make_oldcforest(newforest):
    """Make an old CollapsedForest"""
    return OldCollapsedForest(
        [make_oldctree(tree) for tree in newforest._forest.get_trees()]
    )


trees1 = pp.parse_outfile(
    "tests/example_output/original/small_outfile",
    abundance_file="tests/example_output/original/abundances.csv",
    root="GL",
)
# Sample trees disambiguated
trees1dis = [pp.disambiguate(tree.copy()) for tree in trees1]
trees2 = pp.parse_outfile(
    "tests/example_output/observed_root/small_outfile",
    abundance_file="tests/example_output/observed_root/abundances.csv",
    root="GL",
)
# Sample trees with observed root disambiguated
trees2dis = [pp.disambiguate(tree.copy()) for tree in trees2]

# Ambiguous and disambiguated, lists of sample trees and their associated counts
testtrees = [trees1, trees2]
testtreesdis = [trees1dis, trees2dis]

# The three kinds of CollapsedForests we're comparing:
# new ones with hDAG
newforests = [bp.CollapsedForest(trees) for trees in testtrees]

# new ones with a list of ctrees
newforests_ctrees = []
for forest in newforests:
    f = bp.CollapsedForest()
    f._ctrees = list(forest)
    f.n_trees = len(f._ctrees)
    f._cm_countlist = None
    newforests_ctrees.append(f)

# old ones
oldforests = [make_oldcforest(newforest) for newforest in newforests]

# all zipped
allforests = zip(newforests, newforests_ctrees, oldforests)

cmcount_dagfuncs = bp._cmcounter_dagfuncs()


def first(it):
    """Return the first element of an arbitrary iterable"""
    return list(zip(it, range(1)))[0][0]


def ll_isclose(ll1, ll2):
    """Check if likelihoods output by ll methods are the same"""
    return all(
        [
            np.isclose(ll1[0], ll2[0]),
            np.isclose(ll1[1][0], ll2[1][0]),
            np.isclose(ll1[1][1], ll2[1][1]),
        ]
    )


def test_newcounters():
    """Make sure the cmcounts found by new CollapsedTree init agree with old"""

    def pseudocount(mset):
        if (0, 1) in mset:
            return mset - {(0, 1)} + {(1, 1)}
        else:
            return mset

    for newforest, newforest_ctrees, oldforest in allforests:
        # To build cached _cm_counts
        newforest.ll(0.5, 0.5)
        newforest_ctrees.ll(0.5, 0.5)

        # Test just a single tree (the first in each forest)
        oldtreecounts = pseudocount(FrozenMultiset(oldforest.forest[0]._cm_list))
        newtreecounts = pseudocount(
            first(first(newforest._forest.get_trees()).weight_count(**cmcount_dagfuncs))
        )
        newtreectreecounts = FrozenMultiset(dict(first(newforest_ctrees)._cm_counts))

        assert oldtreecounts == newtreecounts
        assert oldtreecounts == newtreectreecounts

        # Test for the whole forests all together
        oldcmcounts = FrozenMultiset(
            [pseudocount(FrozenMultiset(ctree._cm_list)) for ctree in oldforest.forest]
        )

        newcmcounts_ctrees = FrozenMultiset()
        for treecounts, a in newforest_ctrees._cm_countlist:
            newcmcounts_ctrees += FrozenMultiset({FrozenMultiset(dict(treecounts)): a})

        newcmcounts = FrozenMultiset()
        for treecounts, a in newforest._cm_countlist:
            newcmcounts += FrozenMultiset({FrozenMultiset(dict(treecounts)): a})

        assert newcmcounts_ctrees == newcmcounts
        assert oldcmcounts == newcmcounts_ctrees


def test_newlikelihoods():
    """Make sure the likelihoods found by new CollapsedTree.ll agree with old"""
    p, q = 0.4, 0.6
    ll_dagfuncs = bp._ll_genotype_dagfuncs(p, q)

    # new and old trees agree
    for newforest, newforest_ctrees, oldforest in allforests:
        oldtreell = oldforest.forest[0].ll(p, q)
        newtreell = first(newforest).ll(p, q)
        newtreellctree = first(newforest_ctrees).ll(p, q)

        assert ll_isclose(oldtreell, newtreell)
        assert ll_isclose(oldtreell, newtreellctree)

        for marginal in [True, False]:
            oldfll = oldforest.ll(p, q, marginal=marginal)
            newfll = newforest.ll(p, q, marginal=marginal)
            newfctreell = newforest_ctrees.ll(p, q, marginal=marginal)

            assert ll_isclose(oldfll, newfll)
            assert ll_isclose(oldfll, newfctreell)
        for dagll, treell in zip(
            sorted(newforest._forest.weight_count(**ll_dagfuncs).elements()),
            sorted(ctree.ll(p, q) for ctree in newforest),
        ):
            assert np.isclose(dagll, treell)


def test_fit():
    """mle is the same from old code and computed directly using the DAG.
    This test extends ``gctree test`` (which verifies that inferred parameters
    on simulated forests consisting of a list of ctrees match simulation parameters)
    to also validate parameter inference in forests using the DAG to store trees.
    ``newforest_ctrees`` is the forest constructed from all the ctrees extracted from
    the DAG in ``newforest``"""
    for newforest, newforest_ctrees, oldforest in allforests:
        newforest._cm_countlist = None
        newforest_ctrees._cm_countlist = None
        assert newforest.n_trees == newforest_ctrees.n_trees
        pnew, qnew = newforest.mle(marginal=True)
        pnewc, qnewc = newforest_ctrees.mle(marginal=True)
        pold, qold = oldforest.mle(marginal=True)
        assert np.isclose(pnew, pold)
        assert np.isclose(qnew, qold)
        assert np.isclose(pnew, pnewc)
        assert np.isclose(qnew, qnewc)


def test_validate_ll_genotype():
    """compare old ll_genotype to new ll_genotype directly"""
    c_max, m_max = 10, 10
    for params in [(0.4, 0.6), (0.3, 0.5)]:
        for c in range(c_max):
            for m in range(m_max):
                if c > 0 or m > 1:
                    with np.errstate(all="raise"):
                        true_res = OldCollapsedTree._ll_genotype(c, m, *params)
                        res = bp.CollapsedTree._ll_genotype(c, m, *params)
                    assert np.isclose(true_res[0], res[0])
                    assert np.isclose(true_res[1][0], res[1][0])
                    assert np.isclose(true_res[1][1], res[1][1])


def test_recursion_depth():
    """Be sure ahead-of-time caching is implemented correctly to avoid
    recursion depth issues"""
    bp.CollapsedTree._ll_genotype.cache_clear()
    bp.CollapsedTree._max_ll_cache = {}
    with np.errstate(all="raise"):
        bp.CollapsedTree._ll_genotype(2, 500, 0.4, 0.6)


def test_context_likelihood():
    # These files will be present if pytest is run through `make test`.
    mutation_model = mm.MutationModel(
        mutability_file="HS5F_Mutability.csv", substitution_file="HS5F_Substitution.csv"
    )
    log_likelihood = mm._context_poisson_likelihood(mutation_model, splits=[])

    parent_seq = "AAGAAA"
    child_seq = "AATCAA"

    term1 = sum(
        log(
            mutation_model.mutability(fivemer)[0]
            * mutation_model.mutability(fivemer)[1][target_base]
        )
        for fivemer, target_base in [("AAGAA", "T"), ("AGAAA", "C")]
    )
    sum_mutabilities = sum(
        mutation_model.mutability(fivemer)[0]
        for fivemer in ["NNAAG", "NAAGA", "AAGAA", "AGAAA", "GAAAN", "AAANN"]
    )
    true_val = term1 + 2 * log(2 / sum_mutabilities) - 2
    assert true_val == log_likelihood(parent_seq, child_seq)

    # Now test chain split:
    parent_seq = parent_seq + parent_seq
    child_seq = child_seq + child_seq
    # At index 6, the second concatenated sequence starts.
    log_likelihood = mm._context_poisson_likelihood(mutation_model, splits=[6])

    true_val = 2 * term1 + 4 * log(4 / (2 * sum_mutabilities)) - 4
    assert true_val == log_likelihood(parent_seq, child_seq)
