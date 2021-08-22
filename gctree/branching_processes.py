r"""
This module contains classes for simulation and inference for a binary
branching process with mutation in which the tree is collapsed to nodes that
count the number of clonal leaves of each type.
"""

from __future__ import annotations

from gctree.utils import hamming_distance

import numpy as np
import warnings
import random
import os
from scipy.special import logsumexp
from scipy.special import softmax
from scipy.optimize import minimize, check_grad
import ete3
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import MultipleSeqAlignment
import pickle
from functools import lru_cache
from typing import Tuple, Dict, List, Union, Set, Callable

np.seterr(all="raise")


class CollapsedTree:
    r"""A collapsed tree, modeled as an infinite type Galton-Watson process run
    to extinction.

    Attributes:
        tree: :class:`ete3.TreeNode` object with ``abundance`` node features

    Args:
        tree: ete3 tree with ``abundance`` node features. If uncollapsed, it will be collapsed along branches with no mutations. Can be ommitted on initializaion, and later simulated.
        allow_repeats: tolerate the existence of nodes with the same genotype after collapse, e.g. in sister clades.
    """

    def __init__(self, tree: ete3.TreeNode = None, allow_repeats: bool = False):
        if tree is not None:
            self.tree = tree.copy()
            # remove unobserved internal unifurcations
            for node in self.tree.iter_descendants():
                parent = node.up
                if node.abundance == 0 and len(node.children) == 1:
                    node.delete(prevent_nondicotomic=False)
                    node.children[0].dist = hamming_distance(
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

    @staticmethod
    def _simulate_genotype(p: np.float64, q: np.float64) -> Tuple[int, int]:
        r"""Simulate the number of clonal leaves :math:`c` and mutant clades
        :math:`m` as a Galton-Watson process with branching probability
        :math:`p` and mutation probability :math:`q`.

        Args:
            p: branching probability
            q: mutation probability

        Returns:
            Tuple :math:`(c, m)` of the number of clonal leaves and mutant clades
        """
        if not (0 <= p <= 1 and 0 <= q <= 1):
            raise ValueError(
                "p and q must be in the unit interval: " f"p = {p}, q = {q}"
            )
        if p >= 0.5:
            warnings.warn(
                f"p = {p} is not subcritical, tree simulations not garanteed to terminate!"
            )
        # let's track the tree in breadth first order, listing number of clonal
        # and mutant descendants of each node mutant clades terminate in this
        # view
        cumsum_clones = 0
        len_tree = 0
        c = 0
        m = 0
        # while termination condition not met
        while cumsum_clones > len_tree - 1:
            if random.random() < p:
                mutants = sum(random.random() < q for child in range(2))
                clones = 2 - mutants
                m += mutants
            else:
                mutants = 0
                clones = 0
                c += 1
            cumsum_clones += clones
            len_tree += 1
        assert cumsum_clones == len_tree - 1
        return c, m

    @staticmethod
    @lru_cache(maxsize=None)
    def _ll_genotype(
        c: int, m: int, p: np.float64, q: np.float64
    ) -> Tuple[np.float64, np.ndarray]:
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
                neighbor_ll_genotype, (
                    neighbor_dlogfdp,
                    neighbor_dlogfdq,
                ) = CollapsedTree._ll_genotype(c, m - 1, p, q)
                logf_result = (
                    np.log(2)
                    + np.log(p)
                    + np.log(q)
                    + np.log(1 - q)
                    + neighbor_ll_genotype
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
                        ) = CollapsedTree._ll_genotype(cx, mx, p, q)
                        (
                            neighbor2_ll_genotype,
                            (neighbor2_dlogfdp, neighbor2_dlogfdq),
                        ) = CollapsedTree._ll_genotype(c - cx, m - mx, p, q)
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

    @staticmethod
    def _build_ll_genotype_cache(c_max: int, m_max: int, p: np.float64, q: np.float64):
        r"""build up the lru_cache from the bottom to avoid recursion depth
        issues."""
        CollapsedTree._ll_genotype.cache_clear()
        for c in range(c_max + 1):
            for m in range(m_max + 1):
                CollapsedTree._ll_genotype(c, m, p, q)

    def ll(
        self, p: np.float64, q: np.float64, build_cache: bool = True
    ) -> Tuple[np.float64, np.ndarray]:
        r"""Log likelihood of branching process parameters :math:`(p, q)` given tree topology :math:`T` and genotype abundances :math:`A`.

        .. math::
            \ell(p, q; T, A) = \log\mathbb{P}(T, A \mid p, q)

        Args:
            p: branching probability
            q: mutation probability
            build_cache: build cache from the bottom up. Normally this should be left to its default ``True``.

        Returns:
            Log likelihood :math:`\ell(p, q; T, A)` and its gradient :math:`\nabla\ell(p, q; T, A)`
        """
        if self.tree is None:
            raise ValueError("tree data must be defined to compute likelihood")
        if build_cache:
            self._build_ll_genotype_cache(self._c_max, self._m_max, p, q)
        if (
            self._cm_list[0][0] == 0
            and self._cm_list[0][1] == 1
            and CollapsedTree._ll_genotype(
                self._cm_list[0][0], self._cm_list[0][1], p, q
            )[0]
            == -np.inf
        ):
            # if unifurcation not possible under current model, add a
            # psuedocount for the root
            self._cm_list[0] = (1, self._cm_list[0][1])
        # extract vector of function values and gradient components
        logf_data = [CollapsedTree._ll_genotype(c, m, p, q) for c, m in self._cm_list]
        logf = np.array([[x[0]] for x in logf_data]).sum()
        grad_ll_genotype = np.array([x[1] for x in logf_data]).sum(axis=0)
        return logf, grad_ll_genotype

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

    def simulate(self, p: np.float64, q: np.float64, root: bool = True):
        r"""Simulate a collapsed tree as an infinite type Galton-Watson process
        run to extintion, with branching probability :math:`p` and mutation
        probability :math:`q`. Overwrites existing tree attribute.

        Args:
            p: branching probability
            q: mutation probability
            root: flag indicating simulation is being run from the root of the tree, so we should update tree attributes (should usually be ``True``)
        """
        c, m = self._simulate_genotype(p, q)
        self.tree = ete3.TreeNode()
        self.tree.add_feature("abundance", c)
        for _ in range(m):
            # ooooh, recursion
            child = CollapsedTree()
            child.simulate(p, q, root=False)
            child = child.tree
            child.dist = 1
            self.tree.add_child(child)

        if root:
            # create list of (c, m) for each node
            self._cm_list = [
                (node.abundance, len(node.children)) for node in self.tree.traverse()
            ]
            # store max c and m
            self._c_max = max(node.abundance for node in self.tree.traverse())
            self._m_max = max(len(node.children) for node in self.tree.traverse())

    def __repr__(self):
        r"""return a string representation for printing."""
        return str(self.tree)

    def render(
        self,
        outfile: str,
        idlabel: bool = False,
        colormap: Dict = None,
        frame: int = None,
        position_map: List = None,
        chain_split: int = None,
        frame2: int = None,
        position_map2: List = None,
        show_support: bool = False,
    ):
        r"""Render to tree image file.

        Args:
            outfile: file name to render to, filetype inferred from suffix, .svg for color
            idlabel: label nodes with seq ids, and write sequences of all nodes to a fasta file with same base name as ``outfile``
            frame: coding frame for annotating amino acid substitutions
            position_map: mapping of position names for sequence indices, to be used with substitution annotations and the ``frame`` argument
            chain_split: if sequences are a concatenation two gene sequences, this is the index at which the 2nd one starts (requires ``frame`` and ``frame2`` arguments)
            frame2: coding frame for 2nd sequence when using ``chain_split``
            position_map2: like ``position_map``, but for 2nd sequence when using ``chain_split``
            show_support: annotate bootstrap support if available
        """
        if frame is not None and frame not in (1, 2, 3):
            raise RuntimeError("frame must be 1, 2, or 3")

        def my_layout(node):
            if colormap is None or node.name not in colormap:
                circle_color = "lightgray"
            else:
                circle_color = colormap[node.name]
            text_color = "black"
            if isinstance(circle_color, str):
                C = ete3.CircleFace(
                    radius=max(3, 10 * np.sqrt(node.abundance)),
                    color=circle_color,
                    label={"text": str(node.abundance), "color": text_color}
                    if node.abundance > 0
                    else None,
                )
                C.rotation = -90
                C.hz_align = 1
                ete3.faces.add_face_to_node(C, node, 0)
            else:
                P = ete3.PieChartFace(
                    [100 * x / node.abundance for x in circle_color.values()],
                    2 * 10 * np.sqrt(node.abundance),
                    2 * 10 * np.sqrt(node.abundance),
                    colors=[
                        (color if color != "None" else "lightgray")
                        for color in list(circle_color.keys())
                    ],
                    line_color=None,
                )
                T = ete3.TextFace(
                    " ".join([str(x) for x in list(circle_color.values())]),
                    tight_text=True,
                )
                T.hz_align = 1
                T.rotation = -90
                ete3.faces.add_face_to_node(P, node, 0, position="branch-right")
                ete3.faces.add_face_to_node(T, node, 1, position="branch-right")
            if idlabel:
                T = ete3.TextFace(node.name, tight_text=True, fsize=6)
                T.rotation = -90
                T.hz_align = 1
                ete3.faces.add_face_to_node(
                    T,
                    node,
                    1 if isinstance(circle_color, str) else 2,
                    position="branch-right",
                )

        # we render on a copy, so faces are not permanent
        tree_copy = self.tree.copy(method="deepcopy")
        for node in tree_copy.traverse():
            nstyle = ete3.NodeStyle()
            nstyle["size"] = 0
            if node.up is not None:
                if "sequence" in tree_copy.features and set(
                    node.sequence.upper()
                ) == set("ACGT"):
                    if frame is not None:
                        if chain_split is not None and frame2 is None:
                            raise ValueError(
                                "must define frame2 when using chain_split"
                            )
                        if frame2 is not None and chain_split is None:
                            raise ValueError(
                                "must define chain_split when using frame2"
                            )
                        # loop over split heavy/light chain subsequences
                        for start, end, framex, position_mapx in (
                            (0, chain_split, frame, position_map),
                            (chain_split, None, frame2, position_map2),
                        ):
                            if start is not None:
                                seq = node.sequence[start:end]
                                aa = Seq(
                                    seq[
                                        (framex - 1) : (
                                            framex
                                            - 1
                                            + (3 * ((len(seq) - (framex - 1)) // 3))
                                        )
                                    ]
                                ).translate()
                                seq = node.up.sequence[start:end]
                                aa_parent = Seq(
                                    seq[
                                        (framex - 1) : (
                                            framex
                                            - 1
                                            + (3 * ((len(seq) - (framex - 1)) // 3))
                                        )
                                    ]
                                ).translate()
                                mutations = [
                                    f"{aa1}{pos if position_mapx is None else position_mapx[pos]}{aa2}"
                                    for pos, (aa1, aa2) in enumerate(zip(aa_parent, aa))
                                    if aa1 != aa2
                                ]
                                if mutations:
                                    T = ete3.TextFace(
                                        "\n".join(mutations),
                                        fsize=6,
                                        tight_text=False,
                                        ftype="Courier",
                                    )
                                    if start == 0:
                                        T.margin_top = 6
                                    else:
                                        T.margin_bottom = 6
                                    T.rotation = -90
                                    node.add_face(
                                        T,
                                        0,
                                        position="branch-bottom"
                                        if start == 0
                                        else "branch-top",
                                    )
                                if "*" in aa:
                                    nstyle["hz_line_color"] = "red"
            node.set_style(nstyle)

        ts = ete3.TreeStyle()
        ts.show_leaf_name = False
        ts.rotation = 90
        ts.draw_aligned_faces_as_table = False
        ts.allow_face_overlap = True
        ts.layout_fn = my_layout
        ts.show_scale = False
        ts.show_branch_support = show_support
        # if we labelled seqs, let's also write the alignment out so we have
        # the sequences (including of internal nodes)
        if idlabel:
            aln = MultipleSeqAlignment([])
            for node in tree_copy.traverse():
                aln.append(
                    SeqRecord(
                        Seq(str(node.sequence)),
                        id=str(node.name),
                        description=f"abundance={node.abundance}",
                    )
                )
            AlignIO.write(
                aln, open(os.path.splitext(outfile)[0] + ".fasta", "w"), "fasta"
            )
        return tree_copy.render(outfile, tree_style=ts)

    def write(self, file_name: str):
        r"""Serialize to pickle file.

        Args:
            file_name: file name (.p suffix recommended)
        """
        with open(file_name, "wb") as f:
            pickle.dump(self, f)

    def newick(self, file_name: str):
        r"""Write to newick file.

        Args:
            file_name: file name (.nk suffix recommended)
        """
        self.tree.write(format=1, outfile=file_name)

    def compare(
        self, tree2: CollapsedTree, method: str = "identity"
    ) -> Union[bool, np.float64]:
        r"""Compare this tree to the other tree.

        Args:
            tree2: another object of this type
            method: comparison type (``identity``, ``MRCA``, or ``RF``)

        Returns:
            tree difference
        """
        if method == "identity":
            # we compare lists of seq, parent, abundance
            # return true if these lists are identical, else false
            list1 = sorted(
                (
                    node.sequence,
                    node.abundance,
                    node.up.sequence if node.up is not None else None,
                )
                for node in self.tree.traverse()
            )
            list2 = sorted(
                (
                    node.sequence,
                    node.abundance,
                    node.up.sequence if node.up is not None else None,
                )
                for node in tree2.tree.traverse()
            )
            return list1 == list2
        elif method == "MRCA":
            # matrix of hamming distance of common ancestors of taxa
            # takes a true and inferred tree as CollapsedTree objects
            taxa = [node.sequence for node in self.tree.traverse() if node.abundance]
            n_taxa = len(taxa)
            d = np.zeros(shape=(n_taxa, n_taxa))
            sum_sites = np.zeros(shape=(n_taxa, n_taxa))
            for i in range(n_taxa):
                nodei_true = self.tree.iter_search_nodes(sequence=taxa[i]).__next__()
                nodei = tree2.tree.iter_search_nodes(sequence=taxa[i]).__next__()
                for j in range(i + 1, n_taxa):
                    nodej_true = self.tree.iter_search_nodes(
                        sequence=taxa[j]
                    ).__next__()
                    nodej = tree2.tree.iter_search_nodes(sequence=taxa[j]).__next__()
                    MRCA_true = self.tree.get_common_ancestor(
                        (nodei_true, nodej_true)
                    ).sequence
                    MRCA = tree2.tree.get_common_ancestor((nodei, nodej)).sequence
                    d[i, j] = hamming_distance(MRCA_true, MRCA)
                    sum_sites[i, j] = len(MRCA_true)
            return d.sum() / sum_sites.sum()
        elif method == "RF":
            tree1_copy = self.tree.copy(method="deepcopy")
            tree2_copy = tree2.tree.copy(method="deepcopy")
            for treex in (tree1_copy, tree2_copy):
                for node in list(treex.traverse()):
                    if node.abundance > 0:
                        child = ete3.TreeNode()
                        child.add_feature("sequence", node.sequence)
                        node.add_child(child)
            return tree1_copy.robinson_foulds(
                tree2_copy,
                attr_t1="sequence",
                attr_t2="sequence",
                unrooted_trees=True,
            )[0]
        else:
            raise ValueError("invalid distance method: " + method)

    def _get_split(self, node: ete3.TreeNode) -> Tuple[Set, Set]:
        r"""Return the bipartition resulting from clipping this node's edge
        above.

        Args:
            node: tree node

        Returns:
            A tuple of two sets
        """
        if node.get_tree_root() != self.tree:
            raise ValueError("node not found")
        if node == self.tree:
            raise ValueError("this node is the root (no split above)")
        parent = node.up
        taxa1 = []
        for node2 in node.traverse():
            if node2.abundance > 0 or node2 == self.tree:
                if isinstance(node2.name, str):
                    taxa1.append(node2.name)
                else:
                    taxa1.extend(node2.name)
        taxa1 = set(taxa1)
        node.detach()
        taxa2 = []
        for node2 in self.tree.traverse():
            if node2.abundance > 0 or node2 == self.tree:
                if isinstance(node2.name, str):
                    taxa2.append(node2.name)
                else:
                    taxa2.extend(node2.name)
        taxa2 = set(taxa2)
        parent.add_child(node)
        assert taxa1.isdisjoint(taxa2)
        assert taxa1.union(taxa2) == set(
            (
                name
                for node in self.tree.traverse()
                if node.abundance > 0 or node == self.tree
                for name in ((node.name,) if isinstance(node.name, str) else node.name)
            )
        )
        return tuple(sorted([taxa1, taxa2]))

    @staticmethod
    def _split_compatibility(split1, split2):
        diff = split1[0].union(split1[1]) ^ split2[0].union(split2[1])
        if diff:
            raise ValueError(
                "splits do not cover the same taxa\n" f"\ttaxa not in both: {diff}"
            )
        for partition1 in split1:
            for partition2 in split2:
                if partition1.isdisjoint(partition2):
                    return True
        return False

    def support(
        self,
        bootstrap_trees_list: List[CollapsedTree],
        weights: List[np.float64] = None,
        compatibility: bool = False,
    ):
        r"""Compute support from a list of bootstrap :class:`CollapsedTree` objects, and add to tree attibute.

        Args:
            bootstrap_trees_list: List of trees
            weights: weights for each tree, perhaps for weighting parsimony degenerate trees
            compatibility: counts trees that don't disconfirm the split.
        """
        for node in self.tree.get_descendants():
            split = self._get_split(node)
            support = 0
            compatibility_ = 0
            for i, tree in enumerate(bootstrap_trees_list):
                compatible = True
                supported = False
                for boot_node in tree.tree.get_descendants():
                    boot_split = tree._get_split(boot_node)
                    if (
                        compatibility
                        and compatible
                        and not self._split_compatibility(split, boot_split)
                    ):
                        compatible = False
                    if not compatibility and not supported and boot_split == split:
                        supported = True
                if supported:
                    support += weights[i] if weights is not None else 1
                if compatible:
                    compatibility_ += weights[i] if weights is not None else 1
            node.support = compatibility_ if compatibility else support


class CollapsedForest:
    r"""A collection of :class:`CollapsedTree`

    We can intialize with a list of trees, each an instance of :class:`CollapsedTree`, or we can simulate the forest later.

    Attributes:
        forest: list of trees
        n_trees: number of trees in forest

    Args:
        forest: list of :class:`CollapsedTree`
    """

    def __init__(self, forest: List[CollapsedTree] = None):
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

    def simulate(self, p: np.float64, q: np.float64, n_trees: int):
        r"""Simulate a forest of collapsed trees. Overwrites existing forest attribute.

        Args:
            p: branching probability
            q: mutation probability
            n_trees: number of trees
        """
        self.forest = [CollapsedTree() for _ in range(n_trees)]
        self.n_trees = n_trees
        for tree in self.forest:
            tree.simulate(p, q)

        self._c_max = max(
            node.abundance for tree in self.forest for node in tree.tree.traverse()
        )
        self._m_max = max(
            len(node.children) for tree in self.forest for node in tree.tree.traverse()
        )

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
        CollapsedTree._build_ll_genotype_cache(self._c_max, self._m_max, p, q)
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
                grad_l.append(
                    grad_ls[i_prime, j]
                    + np.exp(
                        logsumexp(
                            ls - ls[i_prime], b=grad_ls[:, j] - grad_ls[i_prime, j]
                        )
                        - logsumexp(ls - ls[i_prime])
                    )
                )
            return (-np.log(len(ls)) + logsumexp(ls)), np.array(grad_l)
        else:
            return ls.sum(), grad_ls.sum(axis=0)

    def mle(self, **kwargs) -> Tuple[np.float64, np.float64]:
        return _mle_helper(self.ll, **kwargs)

    mle.__doc__ = CollapsedTree.mle.__doc__

    def __repr__(self):
        r"""return a string representation for printing."""
        return f"n_trees = {self.n_trees}\n" "\n".join(
            [str(tree) for tree in self.forest]
        )


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
