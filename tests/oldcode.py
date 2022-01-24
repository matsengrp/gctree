import gctree.utils as utils

import warnings
import numpy as np
from scipy.special import logsumexp, softmax
from scipy.optimize import minimize, check_grad
import ete3
from typing import Tuple, List, Callable
from functools import lru_cache


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
