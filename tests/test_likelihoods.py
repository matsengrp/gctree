import gctree.branching_processes as bp
import gctree.phylip_parse as pp
import gctree.utils as utils

import warnings
from collections import Counter
import numpy as np
from scipy.special import logsumexp, softmax
from scipy.optimize import minimize, check_grad
import ete3
from multiset import FrozenMultiset
from typing import Tuple, List, Callable
from functools import lru_cache


def make_ctree(cladetree):
    etetree = cladetree.to_ete(
        name_func=lambda n: n.attr["name"],
        features=["sequence"],
        feature_funcs={"abundance": lambda n: n.attr["abundance"]},
    )
    return bp.CollapsedTree(etetree)


def make_oldctree(tree):
    etetree = tree.to_ete(
        name_func=lambda n: n.attr["name"],
        features=["sequence"],
        feature_funcs={"abundance": lambda n: n.attr["abundance"]},
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
    return OldCollapsedForest([make_oldctree(tree) for tree in newforest._forest.get_trees()])


trees_seqcounts1 = pp.parse_outfile(
    "tests/example_output/original/small_outfile",
    abundance_file="tests/example_output/original/abundances.csv",
    root="GL",
)
trees1 = trees_seqcounts1[0]
trees1dis = [pp.disambiguate(tree.copy()) for tree in trees1]
trees_seqcounts2 = pp.parse_outfile(
    "tests/example_output/observed_root/small_outfile",
    abundance_file="tests/example_output/observed_root/abundances.csv",
    root="GL",
)
trees2 = trees_seqcounts2[0]
trees2dis = [pp.disambiguate(tree.copy()) for tree in trees2]
singletrees = [pp.disambiguate(trees[0].copy()) for trees in [trees1, trees2]]

forests = [
    bp.CollapsedForest(*trees_seqcounts)
    for trees_seqcounts in [trees_seqcounts1, trees_seqcounts2]
]

# ctrees and new cforests containing only one tree:
singletree_ctrees = [bp.CollapsedTree(tree) for tree in singletrees]
singletree_cforests = [bp.CollapsedForest([pp.disambiguate(tree.copy())], sequence_counts=trees_seqcount) for tree, trees_seqcount in [[trees1[0], trees_seqcounts1[1]], [trees2[0], trees_seqcounts2[1]]]]
singletree_forests = [
    bp.CollapsedForest([tree], seqcount)
    for tree, seqcount in zip(singletrees, [trees_seqcounts1[1], trees_seqcounts2[1]])
]
# These are for testing if added leaf is dealt with correctly.

cmcount_dagfuncs = bp._cmcounter_dagfuncs()


def tree_likelihood(dag, p, q, marginal=True):
    cmcounters = dag.weight_count(**cmcount_dagfuncs)
    cmcountlist = [(tuple(mset.items()), mult) for mset, mult in cmcounters.items()]
    return bp._lltree(cmcountlist, p, q, marginal=marginal)


def cmset_to_tuple(mset):
    # _lltree checks first item in list for unobserved root unifurcation,
    # but here it could end up not being first.
    if (0, 1) in mset:
        assert mset[(0, 1)] == 1
        mset = mset - {(0, 1)} + {(1, 1)}
    return tuple(mset.items())


def dag_likelihood(dag, p, q, marginal=True):
    cmcounters = dag.weight_count(**cmcount_dagfuncs)
    cmcountlist = [(cmset_to_tuple(mset), mult) for mset, mult in cmcounters.items()]
    return bp._llforest(cmcountlist, p, q, marginal=marginal)


def test_newcounters():
    # Make sure the cmcounts found by new CollapsedTree init agree with old
    # CollapsedTree init:
    oldctrees = [OldCollapsedTree(tree.copy()) for tree in singletrees]
    oldcmcounts = [FrozenMultiset(ctree._cm_list) for ctree in oldctrees]
    for forest in singletree_cforests:
        forest.ll(.5, .5)
    newcmcounts = [
        FrozenMultiset(dict(ctree._cm_counts))
        for ctree in singletree_ctrees
    ]
    assert oldcmcounts == newcmcounts


def test_newlikelihoods():
    # Make sure the likelihoods found by new CollapsedTree.ll agree with old
    p, q = 0.4, 0.6
    # new and old trees agree
    for tree in trees1dis + trees2dis:
        oldctree = OldCollapsedTree(tree.copy())
        oldtreell = oldctree.ll(p, q)
        newctree = bp.CollapsedTree(tree.copy())
        newtreell = newctree.ll(p, q)
        assert np.isclose(oldtreell[0], newtreell[0])
        assert np.isclose(oldtreell[1][0], newtreell[1][0])
        assert np.isclose(oldtreell[1][1], newtreell[1][1])

    # new forest and new forest made with ctrees agree
    newforests = [
        bp.CollapsedForest(trees, sequence_counts=seqcounts)
        for trees, seqcounts in [[trees1dis, trees_seqcounts1[1]], [trees2dis, trees_seqcounts2[1]]]
    ]
    oldforests = [make_oldcforest(newcforest) for newcforest in newforests]
    for marginal in [True, False]:
        for newforest, oldforest in zip(newforests, oldforests):
            newfll = newforest.ll(p, q, marginal=marginal)
            oldfll = oldforest.ll(p, q, marginal=marginal)
            assert np.isclose(oldfll[0], newfll[0])
            assert np.isclose(oldfll[1][0], newfll[1][0])
            assert np.isclose(oldfll[1][1], newfll[1][1])


# def test_newcounters_after_dag():
#     # Make sure the cmcounts found by new CollapsedTree init,
#     # after going through DAG, agree with old
#     # CollapsedTree init:
#     oldctrees = [OldCollapsedTree(tree.copy()) for tree in trees1dis + trees2dis]
#     dags = [pp.make_dag([tree], trees_seqcounts1[1]) for tree in trees1dis]
#     dags.extend([pp.make_dag([tree], trees_seqcounts2[1]) for tree in trees2dis])
#     newctrees = [make_ctree(dag) for dag in dags]
#     oldcmlists = [ctree._cm_list for ctree in oldctrees]
#     oldcmcounts = [FrozenMultiset(oldcmlist) for oldcmlist in oldcmlists]
#     newcmcounts = [FrozenMultiset(dict(ctree._cm_counts)) for ctree in newctrees]
#     assert oldcmcounts == newcmcounts


# def test_newlikelihoods_singletree_after_dag():
#     # old likelihoods are the same as new likelihoods (computed by
#     # CollapsedTree) for single trees added to
#     # DAG then removed and made new CollapsedTree
#     p, q = 0.4, 0.6
#     for trees, seqcount in [
#         (trees1dis, trees_seqcounts1[1]),
#         (trees2dis, trees_seqcounts2[1]),
#     ]:
#         for tree in trees:
#             oldctree = OldCollapsedTree(tree.copy())
#             oldtreell = oldctree.ll(p, q)
#             newctree = bp.clade_tree_to_ctree(pp.make_dag([tree], seqcount))
#             newtreell = newctree.ll(p, q)
#             assert np.isclose(oldtreell[0], newtreell[0])
#             assert np.isclose(oldtreell[1][0], newtreell[1][0])
#             assert np.isclose(oldtreell[1][1], newtreell[1][1])


# def test_newcounters_singletree_in_dag():
#     # old cmcounters are the same as new cmcounters computed in DAG
#     # for single trees...
#     for trees, seqcount in [
#         (trees1dis, trees_seqcounts1[1]),
#         (trees2dis, trees_seqcounts2[1]),
#     ]:
#         for tree in trees:
#             tree.dist = 0
#             for node in tree.iter_descendants():
#                 node.dist = utils.hamming_distance(node.up.sequence, node.sequence)
#             oldctree = OldCollapsedTree(tree.copy())
#             oldcmlist = oldctree._cm_list
#             oldtreecounters = FrozenMultiset(oldcmlist)
#             newctree = bp.CollapsedTree(tree.copy())
#             newtreecounters = FrozenMultiset(dict(newctree._cm_counts))
#             dag = pp.make_dag([tree], seqcount)
#             dag.convert_to_collapsed()
#             cmcounters = dag.weight_count(**cmcount_dagfuncs)
#             assert len(cmcounters) == 1  # there's just one tree in the dag
#             cmcounts = list(cmcounters.keys())[0]
#             assert cmcounts == oldtreecounters
#             assert cmcounts == newtreecounters


# def test_problem():
#     tree = trees2dis[0]
#     print(tree.name)
#     print(len(tree.children))
#     print(tree.children[0].dist)
#     print(tree.sequence == tree.children[0].sequence)
#     print(tree.sequence in {node.sequence for node in tree.iter_leaves()})
#     print(tree.abundance)
#     seqcounts = trees_seqcounts2[1]
#     for node in tree:
#         if node.sequence in seqcounts:
#             node.abundance = seqcounts[node.sequence]
#         else:
#             node.abundance = 0
#     print(seqcounts[node.sequence])
#     oldctree = OldCollapsedTree(tree.copy())
#     oldtreecounters = FrozenMultiset(oldctree._cm_list)

#     dag = pp.make_dag([tree], seqcounts)
#     cmcounters = dag.weight_count(**cmcount_dagfuncs)
#     assert list(cmcounters.keys())[0] == oldtreecounters


# def test_newlikelihoods_singletree_in_dag():
#     # old likelihoods are the same as new likelihoods computed in DAG
#     # for single trees...
#     p, q = 0.4, 0.6
#     for trees, seqcount in [
#         (trees1dis, trees_seqcounts1[1]),
#         (trees2dis, trees_seqcounts2[1]),
#     ]:
#         for tree in trees:
#             oldctree = OldCollapsedTree(tree.copy())
#             oldtreell = oldctree.ll(p, q)
#             dag = pp.make_dag([tree], seqcount)
#             cmcounters = dag.weight_count(**cmcount_dagfuncs)
#             assert len(cmcounters) == 1  # there's just one tree in the dag
#             cmcounts = cmset_to_tuple(list(cmcounters.keys())[0])
#             newtreell = bp._lltree(cmcounts, p, q)
#             assert np.isclose(oldtreell[0], newtreell[0])
#             assert np.isclose(oldtreell[1][0], newtreell[1][0])
#             assert np.isclose(oldtreell[1][1], newtreell[1][1])


# def test_newcounters_in_dag():
#     # old cm counters are the same as new likelihoods computed in DAG
#     oldcforests = [make_oldcforest(dag) for dag in dags]
#     for dag, forest, oldforest in zip(dags, cforests, oldcforests):
#         cmcounters = dag.weight_count(**cmcount_dagfuncs)
#         forestcounters = Counter(
#             [FrozenMultiset(dict(tree._cm_counts)) for tree in forest.forest]
#         )

#         oldforestcounters = Counter(
#             [FrozenMultiset(tree._cm_list) for tree in oldforest.forest]
#         )
#         assert cmcounters == forestcounters
#         assert cmcounters == oldforestcounters


# def test_newll_in_dag():
#     # ll computed in DAG is the same as ll computed in old forest, or in new forest.
#     p, q = 0.4, 0.6

#     oldcforests = [make_oldcforest(dag) for dag in dags]
#     for dag, forest, oldforest in zip(dags, cforests, oldcforests):
#         for marginal in [True, False]:
#             dagll = dag_likelihood(dag, p, q, marginal=marginal)
#             forestll = forest.ll(p, q, marginal=marginal)
#             oldforestll = oldforest.ll(p, q, marginal=marginal)
#             for other in [forestll, oldforestll]:
#                 assert np.isclose(dagll[0], other[0])
#                 assert np.isclose(dagll[1][0], other[1][0])
#                 assert np.isclose(dagll[1][1], other[1][1])


# def test_fit():
#     # mle is the same from old code and computed directly using the DAG.
#     for dag in dags:
#         p, q = bp.fit_branching_process(dag, verbose=True, marginal=True)
#         oldforest = make_oldcforest(dag)
#         pold, qold = oldforest.mle(marginal=True)
#         assert np.isclose(p, pold)
#         assert np.isclose(q, qold)


# def test_validate_ll_genotype():
#     # compare old ll_genotype to new ll_genotype directly
#     c_max, m_max = 10, 10
#     for params in [(0.4, 0.6), (0.3, 0.5)]:
#         for c in range(c_max):
#             for m in range(m_max):
#                 if c > 0 or m > 1:
#                     true_res = OldCollapsedTree._ll_genotype(c, m, *params)
#                     res = bp.CollapsedTree._ll_genotype(c, m, *params)
#                     assert np.isclose(true_res[0], res[0])
#                     assert np.isclose(true_res[1][0], res[1][0])
#                     assert np.isclose(true_res[1][1], res[1][1])


# def test_recursion_depth():
#     # Be sure ahead-of-time caching is implemented correctly to avoid
#     # recursion depth issues
#     bp.CollapsedTree._ll_genotype.cache_clear()
#     bp.CollapsedTree._max_ll_cache = {}
#     bp.CollapsedTree._ll_genotype(2, 500, 0.4, 0.6)


class OldCollapsedTree:
    def __init__(self, tree: ete3.TreeNode = None, allow_repeats: bool = False):
        if tree is not None:
            self.tree = tree.copy()
            # remove unobserved internal unifurcations
            for node in self.tree.iter_descendants():
                parent = node.up
                if node.abundance == 0 and len(node.children) == 1:
                    node.delete(prevent_nondicotomic=False)
                    node.children[0].dist = utils.hamming_distance(
                        node.children[0].sequence, parent.sequence
                    )

            # iterate over the tree below root and collapse edges of zero
            # length if the node is a leaf and it's parent has nonzero
            # abundance we combine taxa names to a set to acommodate
            # bootstrap samples that result in repeated genotypes
            observed_genotypes = set((leaf.name for leaf in self.tree))
            observed_genotypes.add(self.tree.name)
            for node in self.tree.get_descendants(strategy="postorder"):
                if node.dist == 0:
                    node.up.abundance += node.abundance
                    if isinstance(node.name, str):
                        node_set = set([node.name])
                    else:
                        node_set = set(node.name)
                    if isinstance(node.up.name, str):
                        node_up_set = set([node.up.name])
                    else:
                        node_up_set = set(node.up.name)
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

            final_observed_genotypes = set()
            for node in self.tree.traverse():
                if node.abundance > 0 or node == self.tree:
                    for name in (
                        (node.name,) if isinstance(node.name, str) else node.name
                    ):
                        final_observed_genotypes.add(name)
            if final_observed_genotypes != observed_genotypes:
                raise RuntimeError(
                    "observed genotypes don't match after "
                    f"collapse\n\tbefore: {observed_genotypes}"
                    f"\n\tafter: {final_observed_genotypes}\n\t"
                    "symmetric diff: "
                    f"{observed_genotypes ^ final_observed_genotypes}"
                )
            assert sum(node.abundance for node in tree.traverse()) == sum(
                node.abundance for node in self.tree.traverse()
            )

            rep_seq = sum(node.abundance > 0 for node in self.tree.traverse()) - len(
                set(
                    [
                        node.sequence
                        for node in self.tree.traverse()
                        if node.abundance > 0
                    ]
                )
            )
            if not allow_repeats and rep_seq:
                raise RuntimeError(
                    "Repeated observed sequences in collapsed "
                    f"tree. {rep_seq} sequences were found repeated."
                )
            elif allow_repeats and rep_seq:
                rep_seq = sum(
                    node.abundance > 0 for node in self.tree.traverse()
                ) - len(
                    set(
                        [
                            node.sequence
                            for node in self.tree.traverse()
                            if node.abundance > 0
                        ]
                    )
                )
                print(
                    "Repeated observed sequences in collapsed tree. "
                    f"{rep_seq} sequences were found repeated."
                )
            # a custom ladderize accounting for abundance and sequence to break
            # ties in abundance
            for node in self.tree.traverse(strategy="postorder"):
                # add a partition feature and compute it recursively up tree
                node.add_feature(
                    "partition",
                    node.abundance + sum(node2.partition for node2 in node.children),
                )
                # sort children of this node based on partion and sequence
                node.children.sort(key=lambda node: (node.partition, node.sequence))

            # create list of (c, m) for each node
            self._cm_list = [
                (node.abundance, len(node.children)) for node in self.tree.traverse()
            ]
            # store max c and m
            self._c_max = max(node.abundance for node in self.tree.traverse())
            self._m_max = max(len(node.children) for node in self.tree.traverse())
        else:
            self.tree = tree

    def ll(
        self, p: np.float64, q: np.float64, build_cache: bool = True
    ) -> Tuple[np.float64, np.ndarray]:
        if self.tree is None:
            raise ValueError("tree data must be defined to compute likelihood")
        if self._cm_list[0][0] == 0 and self._cm_list[0][1] == 1:
            # if unifurcation not possible under current model, add a
            # psuedocount for the root
            self._cm_list[0] = (1, self._cm_list[0][1])
        # extract vector of function values and gradient components
        logf_data = [
            OldCollapsedTree._ll_genotype(c, m, p, q) for c, m in self._cm_list
        ]
        logf = np.array([[x[0]] for x in logf_data]).sum()
        grad_ll_genotype = np.array([x[1] for x in logf_data]).sum(axis=0)
        return logf, grad_ll_genotype

    @staticmethod
    @lru_cache(maxsize=None)
    def _ll_genotype(
        c: int, m: int, p: np.float64, q: np.float64
    ) -> Tuple[np.float64, np.ndarray]:
        if c == m == 0 or (c == 0 and m == 1):
            raise ValueError("Zero likelihood event")
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
                ) = OldCollapsedTree._ll_genotype(c, m - 1, p, q)
                logg_array = [
                    (
                        np.log(2)
                        + np.log(p)
                        + np.log(q)
                        + np.log(1 - q)
                        + neighbor_ll_genotype
                    )
                ]
                dloggdp_array = [1 / p + neighbor_dlogfdp]
                dloggdq_array = [1 / q - 1 / (1 - q) + neighbor_dlogfdq]
            else:
                logg_array = []
                dloggdp_array = []
                dloggdq_array = []
            for cx in range(c + 1):
                for mx in range(m + 1):
                    if (cx > 0 or mx > 1) and (c - cx > 0 or m - mx > 1):
                        (
                            neighbor1_ll_genotype,
                            (neighbor1_dlogfdp, neighbor1_dlogfdq),
                        ) = OldCollapsedTree._ll_genotype(cx, mx, p, q)
                        (
                            neighbor2_ll_genotype,
                            (neighbor2_dlogfdp, neighbor2_dlogfdq),
                        ) = OldCollapsedTree._ll_genotype(c - cx, m - mx, p, q)
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
            if not logg_array:
                raise ValueError("Zero likelihood event")
            else:
                logf_result = logsumexp(logg_array)
                softmax_logg_array = softmax(logg_array)
                dlogfdp_result = np.multiply(softmax_logg_array, dloggdp_array).sum()
                dlogfdq_result = np.multiply(softmax_logg_array, dloggdq_array).sum()

        return (logf_result, np.array([dlogfdp_result, dlogfdq_result]))


class OldCollapsedForest:
    r"""A collection of :class:`CollapsedTree`
    We can intialize with a list of trees, each an instance of :class:`CollapsedTree`, or we can simulate the forest later.
    Attributes:
        forest: list of trees
        n_trees: number of trees in forest
    Args:
        forest: list of :class:`CollapsedTree`
    """

    def __init__(self, forest: List[OldCollapsedTree] = None):
        if forest is not None:
            if len(forest) == 0:
                raise ValueError("passed empty tree list")
            self.forest = forest
            self.n_trees = len(forest)
            self._c_max = max(
                node.abundance for tree in self.forest for node in tree.tree.traverse()
            )
            self._m_max = max(
                len(node.children)
                for tree in self.forest
                for node in tree.tree.traverse()
            )
        else:
            self.forest = None
            self.n_trees = None

    def ll(
        self,
        p: np.float64,
        q: np.float64,
        marginal: bool = False,
    ) -> Tuple[np.float64, np.ndarray]:
        r"""Log likelihood of branching process parameters :math:`(p, q)` given tree topologies :math:`T_1, \dots, T_n` and corresponding genotype abundances vectors :math:`A_1, \dots, A_n` for each of :math:`n` trees in the forest.
        If ``marginal=False`` (the default), compute the joint log likelihood
        .. math::
            \ell(p, q; T, A) = \sum_{i=1}^n\log\mathbb{P}(T_i, A_i \mid p, q),
        otherwise compute the marginal log likelihood
        .. math::
            \ell(p, q; T, A) = \log\left(\sum_{i=1}^n\mathbb{P}(T_i, A_i \mid p, q)\right).
        Args:
            p: branching probability
            q: mutation probability
            marginal: compute the marginal likelihood over trees, otherwise compute the joint likelihood of trees
        Returns:
            Log likelihood :math:`\ell(p, q; T, A)` and its gradient :math:`\nabla\ell(p, q; T, A)`
        """
        if self.forest is None:
            raise ValueError("forest data must be defined to compute likelihood")
        # we don't want to build the cache again in each tree
        terms = [tree.ll(p, q, build_cache=False) for tree in self.forest]
        ls = np.array([term[0] for term in terms])
        grad_ls = np.array([term[1] for term in terms])
        if marginal:
            # we need to find the smallest derivative component for each
            # coordinate, then subtract off to get positive things to logsumexp
            grad_l = []
            for j in range(len((p, q))):
                i_prime = grad_ls[:, j].argmin()
                b = grad_ls[:, j] - grad_ls[i_prime, j]
                # believe it or not, logsumexp can't handle 0 in b
                # when np.seterr(underflow='raise') on newer processors:
                grad_l.append(
                    grad_ls[i_prime, j]
                    + np.exp(
                        logsumexp((ls - ls[i_prime])[b != 0], b=b[b != 0])
                        - logsumexp(ls - ls[i_prime])
                    )
                )
            return (-np.log(len(ls)) + logsumexp(ls)), np.array(grad_l)
        else:
            return ls.sum(), grad_ls.sum(axis=0)

    def mle(self, **kwargs) -> Tuple[np.float64, np.float64]:
        r"""Maximum likelihood estimate of :math:`(p, q)`.
        .. math::
            (p, q) = \arg\max_{p,q\in [0,1]}\ell(p, q)
        Args:
            kwargs: keyword arguments passed along to the log likelihood :meth:`CollapsedTree.ll`
        Returns:
            Tuple :math:`(p, q)` with estimated branching probability and estimated mutation probability
        """
        return _mle_helper(self.ll, **kwargs)


def _mle_helper(
    ll: Callable[[np.float64, np.float64], Tuple[np.float64, np.ndarray]], **kwargs
) -> Tuple[np.float64, np.float64]:
    # initialization
    x_0 = (0.5, 0.5)
    bounds = ((1e-6, 1 - 1e-6), (1e-6, 1 - 1e-6))

    def f(x):
        """negative log likelihood."""
        return tuple(-x for x in ll(*x, **kwargs))

    grad_check = check_grad(lambda x: f(x)[0], lambda x: f(x)[1], x_0)
    if grad_check > 1e-3:
        warnings.warn(
            "gradient mismatches finite difference " f"approximation by {grad_check}",
            RuntimeWarning,
        )
    result = minimize(
        f,
        jac=True,
        x0=x_0,
        method="L-BFGS-B",
        options={"ftol": 1e-10},
        bounds=bounds,
    )
    # update params if None and optimization successful
    if not result.success:
        warnings.warn("optimization not sucessful, " + result.message, RuntimeWarning)
    return result.x[0], result.x[1]
