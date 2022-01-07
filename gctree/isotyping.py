r"""A module for incorporating isotype information in gctree inference."""

import random
import ete3
from gctree.utils import hamming_distance
import warnings
from typing import Dict, Callable, Optional, Set, Sequence
from functools import wraps


def _assert_switching_order_match(
    fn: Callable[["Isotype", "Isotype"], bool]
) -> Callable:
    @wraps(fn)
    def newfunc(t1: "Isotype", t2: "Isotype") -> bool:
        if t1.order != t2.order:
            raise TypeError(
                "Comparison attempted between isotypes"
                " with different assumed class switching order"
                f" first argument has order {t1.order}"
                f" second argument has order {t2.order}"
            )
        if t1.weight_matrix != t2.weight_matrix:
            raise TypeError(
                "Comparison attempted between isotypes with different weight matrices"
            )
        else:
            return fn(t1, t2)

    return newfunc


class IsotypeTemplate:
    """A template constructor class for :meth:`Isotype`

    Args:
        isotype_order: A list or other Sequence of isotype names, in the allowed class switching order
        weight_matrix: An optional matrix containing transition weights between isotypes. Default weight is 1 for all allowed isotype transitions.
    """

    def __init__(
        self,
        isotype_order: Sequence[str],
        weight_matrix: Optional[Sequence[Sequence[float]]] = None,
    ):
        n = len(isotype_order)
        if weight_matrix is None:
            weight_matrix = [
                [
                    0.0
                    if target == origin
                    else 1.0
                    if target > origin
                    else float("inf")
                    for target in range(n)
                ]
                for origin in range(n)
            ]
        else:
            if len(weight_matrix) != n or len(weight_matrix[0]) != n:
                raise TypeError(
                    "Weight matrix should be an n by n matrix, with n equal to the length of isotype_order."
                )
            # Ensure that illegal transitions are encoded as 'inf' in weight_matrix
            for i in range(n):
                for j in range(n):
                    if i > j:
                        weight_matrix[i][j] = float("inf")
                    elif i == j:
                        weight_matrix[i][j] = 0.0
        self.order = isotype_order
        self.weight_matrix = weight_matrix

    def new(self, isotype_name):
        """Return an :meth:`Isotype` object with the isotype order and weight
        matrix provided to this class's constructor."""
        return Isotype(self.order, self.weight_matrix, isotype_name)


class Isotype:
    """An isotype, and associated information about class switching order and
    transition weights.

    Attributes:
        order: A list or other sequence of isotype names, in their allowed class switching order
        weight_matrix: Transition weights between isotypes, with first index the original isotype index, and second index the new isotype index.
        isotype: Index of isotype name in ``self.order``

    Objects of this class shall be instantiated using :meth:`IsotypeTemplate.new`.
    """

    # From https://doi.org/10.7554/eLife.16578.012, for humans:
    # order = ["M", "G3", "A1", "G2", "G4", "E", "A2"]

    def __init__(
        self,
        order: Sequence[str],
        weight_matrix: Sequence[Sequence[float]],
        isotype_name: str,
    ):
        self.order = order
        self.weight_matrix = weight_matrix
        if isotype_name == "?":
            self.isotype = None
        else:
            try:
                self.isotype = self.order.index(isotype_name)
            except ValueError as e:
                raise ValueError(
                    "Unrecognized isotype name. Check that default"
                    " or provided isotype name list contains all"
                    " observed isotypes\n" + str(e)
                )

    @_assert_switching_order_match
    def isbefore(self, t: "Isotype") -> bool:
        return self.isotype <= t.isotype

    def __repr__(self) -> str:
        return f"Isotype('{self.order[self.isotype]}')"

    def __str__(self) -> str:
        return f"{self.order[self.isotype]}"

    def copy(self) -> "Isotype":
        return Isotype(self.order, self.weight_matrix, self.order[self.isotype])

    @_assert_switching_order_match
    def __eq__(self, other: "Isotype") -> bool:
        return self.isotype == other.isotype

    def __hash__(self) -> int:
        return hash(self.__repr__())

    def resolutions(self) -> Sequence["Isotype"]:
        """Returns list of all possible isotypes if passed an ambiguous
        isotype.

        Otherwise, returns a list containing a copy of the passed
        isotype.
        """
        if self.isotype is None:
            return [
                Isotype(self.order, self.weight_matrix, name) for name in self.order
            ]
        else:
            return [self.copy()]


def isotype_tree(
    tree: ete3.TreeNode,
    newidmap: Dict[str, Dict[str, str]],
    isotype_names: Sequence[str],
    weight_matrix: Optional[Sequence[Sequence[float]]] = None,
) -> ete3.TreeNode:
    """Method adds isotypes to ``tree``, minimizing isotype switching and
    obeying switching order.

    * Adds observed isotypes to each observed node in the collapsed
      trees output by gctree inference. If cells with the same sequence
      but different isotypes are observed, then collapsed tree nodes
      must be ‘exploded’ into new nodes with the appropriate isotypes
      and abundances. Each unique sequence ID generated by gctree is
      prepended to its observed isotype, and a new `isotyped.idmap`
      mapping these new sequence IDs to original sequence IDs is
      written in the output directory.
    * Resolves isotypes of unobserved ancestral genotypes in a way
      that minimizes isotype switching and obeys isotype switching
      order. If observed isotypes of an observed internal node and its
      children violate switching order, then the observed internal node
      is replaced with an unobserved node with the same sequence, and
      the observed internal node is placed as a child leaf. This
      procedure always allows switching order conflicts to be resolved,
      and should usually increase isotype transitions required in the
      resulting tree.

    Args:
        tree: ete3 Tree
        newidmap: mapping of sequence IDs to isotypes, such as that output by :meth:`explode_idmap`.
        isotype_names: list or other sequence of isotype names observed, in correct switching order.

    Returns:
        A new ete3 Tree whose nodes have isotype annotations in the attribute ``isotype``.
        Node names in this tree also contain isotype names.
    """
    tree = tree.copy()
    _add_observed_isotypes(tree, newidmap, isotype_names, weight_matrix=weight_matrix)
    _disambiguate_isotype(tree)
    _collapse_tree_by_sequence_and_isotype(tree)
    for node in tree.traverse():
        node.name = node.name + " " + str(node.isotype)
    for node in tree.iter_descendants():
        node.dist = hamming_distance(node.up.sequence, node.sequence)
    return tree


@_assert_switching_order_match
def isotype_distance(t1: Isotype, t2: Isotype) -> float:
    """Returns the weight of the transition from isotype ``t1`` to isotype
    ``t2``.

    This function is not symmetric on its arguments, so isn't a true
    distance.
    """
    return t1.weight_matrix[t1.isotype][t2.isotype]


def isotype_parsimony(tree: ete3.TreeNode) -> float:
    """Computes the sum of :meth:`isotype_distance` along each edge in an
    isotyped tree.

    If no weight matrix was provided during isotyping of the tree, then
    the return value of this function is the number of isotype
    transitions along edges in the tree.
    """
    return sum(
        isotype_distance(node.up.isotype, node.isotype)
        for node in tree.iter_descendants()
    )


def _disambiguate_isotype(
    tree: ete3.Tree, random_state=None, dist_func=isotype_distance
):
    def is_ambiguous(isotype):
        return isotype.isotype is None

    if random_state is None:
        random.seed(tree.write(format=1))
    else:
        random.setstate(random_state)

    # First pass of Sankoff: compute cost vectors
    for node2 in tree.traverse(strategy="postorder"):
        node2.add_feature("costs", [[seq, 0] for seq in node2.isotype.resolutions()])
        if not node2.is_leaf():
            for seq_cost in node2.costs:
                for child in node2.children:
                    seq_cost[1] += min(
                        [
                            dist_func(seq_cost[0], child_seq) + child_cost
                            for child_seq, child_cost in child.costs
                        ]
                    )
    # Second pass: Choose base and adjust children's cost vectors.
    # Not necessary if we only want the minimum weight:
    for node2 in tree.traverse(strategy="preorder"):
        if not is_ambiguous(node2.isotype):
            continue
        min_cost = min([cost for _, cost in node2.costs])
        resolved_isotype = random.choice(
            [isotype for isotype, cost in node2.costs if cost == min_cost]
        )
        # Adjust child cost vectors
        if not node2.is_leaf():
            for child in node2.children:
                child.costs = [
                    [
                        isotype,
                        cost + dist_func(resolved_isotype, isotype),
                    ]
                    for isotype, cost in child.costs
                ]
        node2.isotype = resolved_isotype
    for node in tree.traverse():
        try:
            node.del_feature("costs")
        except (AttributeError, KeyError):
            pass


def explode_idmap(
    idmap: Dict[str, Set[str]], isotype_map: Dict[str, str]
) -> Dict[str, Dict[str, Set[str]]]:
    """Uses provided mapping of original sequence IDs to observed isotypes to
    'explode' the provided mapping of unique sequence IDs to original sequence
    IDs.

    Args:
        idmap: A dictionary mapping unique sequence IDs to sets of original IDs of observed sequences
        isotype_map: A dictionary mapping original IDs to observed isotype names

    Returns:
        A dictionary mapping unique sequence IDs to dictionaries, mapping isotype names to sets of original sequence IDs.
    """
    newidmap = {}
    for id, cell_ids in idmap.items():
        isotypeset = set()
        for cell_id in cell_ids:
            try:
                isotypeset.add(isotype_map[cell_id])
            except KeyError:
                warnings.warn(
                    f"Sequence ID {id} has original sequence id {cell_id} "
                    "for which no observed isotype was provided. "
                    "Isotype will be assumed ambiguous if observed."
                )
                isotype_map[cell_id] = "?"
                isotypeset.add("?")
        newidmap[id] = {
            isotype: {
                cell_id for cell_id in cell_ids if isotype_map[cell_id] == isotype
            }
            for isotype in isotypeset
        }
    return newidmap


def _add_observed_isotypes(
    tree: ete3.Tree,
    newidmap: Dict[str, str],
    isotype_order: Sequence[str],
    weight_matrix: Optional[Sequence[Sequence[float]]] = None,
):
    # Drop observed nodes as leaves and explode by observed isotype:
    # Descend internal observed nodes as leaves:
    newisotype = IsotypeTemplate(isotype_order, weight_matrix=weight_matrix).new
    for node in list(tree.iter_descendants()):
        if node.abundance > 0 and not node.is_leaf():
            newchild = ete3.TreeNode(name=node.name)
            newchild.add_feature("sequence", node.sequence)
            newchild.add_feature("abundance", node.abundance)
            node.abundance = 0
            node.add_child(child=newchild)
    # Now duplicate nodes which represent multiple isotypes
    for node in list(tree.get_leaves()):
        if node.abundance == 0:
            node.add_feature("isotype", newisotype("?"))
        else:
            try:
                thisnode_isotypemap = newidmap[node.name]
            except KeyError as e:
                warnings.warn(
                    f"The sequence name {e} labels an observed node, but no mapping to an original sequence ID was found."
                    " Isotype will be assumed ambiguous."
                )
                thisnode_isotypemap = {
                    "?": {f"Unknown_id_{n+1}" for n in range(node.abundance)}
                }
            if "?" in thisnode_isotypemap:
                warnings.warn(
                    f"The sequence name {node.name} labels an observed node, and corresponds to sequence IDs for "
                    "which no observed isotype was provided. "
                    f" Isotype will be assumed ambiguous for: {', '.join(thisnode_isotypemap['?'])}"
                )
            # node.name had better be in newidmap, since this is an observed node
            if len(thisnode_isotypemap) > 1:
                for isotype, cell_ids in thisnode_isotypemap.items():
                    # add new node below this leaf node. Must be below, and not child
                    # of parent, to preserve max parsimony in case that node.up has
                    # different sequence from node.
                    newchild = ete3.TreeNode(name=node.name)
                    newchild.add_feature("abundance", len(cell_ids))
                    newchild.add_feature("sequence", node.sequence)
                    newchild.add_feature("isotype", newisotype(isotype))
                    node.add_child(child=newchild)
                node.abundance = 0
            else:
                node.isotype = newisotype(list(thisnode_isotypemap.keys())[0])
    # Now add ancestral ambiguous isotypes
    for node in tree.traverse():
        if not node.is_leaf():
            node.add_feature("isotype", newisotype("?"))


def _collapse_tree_by_sequence_and_isotype(tree: ete3.TreeNode):
    for node in tree.iter_descendants():
        node.dist = node.up.sequence != node.sequence or node.up.isotype != node.isotype
    for node in tree.iter_descendants():
        if node.dist == 0:
            node.up.abundance += node.abundance
            node.up.name = node.name
            node.delete(prevent_nondicotomic=False)
