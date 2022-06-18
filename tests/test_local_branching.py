import ete3
import gctree
import numpy as np
import pytest

tree = ete3.TreeNode(name="naive", dist=0)
tree.abundance = 0
tree.sequence = "A"
tree.isotype = {}

child1 = ete3.TreeNode(name="seq1", dist=1)
child1.abundance = 1
child1.sequence = "C"
child1.isotype = {}

child2 = ete3.TreeNode(name="seq2", dist=2)
child2.abundance = 2
child2.sequence = "G"
child2.isotype = {}

grandchild = ete3.TreeNode(name="seq3", dist=3)
grandchild.abundance = 3
grandchild.sequence = "T"
grandchild.isotype = {}

tree.add_child(child1)
tree.add_child(child2)
child1.add_child(grandchild)

ctree = gctree.CollapsedTree()
ctree.tree = tree

τ = τ0 = 1

ctree.local_branching(tau=τ, tau0=τ0)

# dummy integrated branch length for exploded clonal abundances
clone_contribution = τ * (1 - np.exp(-τ0 / τ))

LB_down = {
    tree: {
        tree: 0,
        child1: τ * (1 - np.exp(-(child1.dist + grandchild.dist) / τ))
        + np.exp(-child1.dist / τ) * child1.abundance * clone_contribution
        + np.exp(-(child1.dist + grandchild.dist) / τ)
        * grandchild.abundance
        * clone_contribution,
        child2: τ * (1 - np.exp(-child2.dist / τ))
        + np.exp(-child2.dist / τ) * child2.abundance * clone_contribution,
    },
    child1: {
        child1: child1.abundance * clone_contribution,
        grandchild: τ * (1 - np.exp(-grandchild.dist / τ))
        + np.exp(-grandchild.dist / τ) * grandchild.abundance * clone_contribution,
    },
    child2: {child2: child2.abundance * clone_contribution},
    grandchild: {grandchild: grandchild.abundance * clone_contribution},
}

LB_up = {tree: τ}
LB_up[child1] = τ * (1 - np.exp(-child1.dist / τ)) + np.exp(-child1.dist / τ) * (
    LB_up[tree]
    + sum(LB_down[tree][message] for message in LB_down[tree] if message != child1)
)
LB_up[child2] = τ * (1 - np.exp(-child2.dist / τ)) + np.exp(-child2.dist / τ) * (
    LB_up[tree]
    + sum(LB_down[tree][message] for message in LB_down[tree] if message != child2)
)
LB_up[grandchild] = τ * (1 - np.exp(-grandchild.dist / τ)) + np.exp(
    -grandchild.dist / τ
) * (
    LB_up[child1]
    + sum(
        LB_down[child1][message] for message in LB_down[child1] if message != grandchild
    )
)

LBI = {node: LB_up[node] + sum(LB_down[node].values()) for node in tree.traverse()}
LBR = {node: sum(LB_down[node].values()) / LB_up[node] for node in tree.traverse()}

print("name\tLBI\tLBR")
for node in ctree.tree.traverse():
    print(f"{node.name}\t{node.LBI:.3f}\t{node.LBR:.3f}")
    assert LB_up[node] == pytest.approx(node.LB_up)
    assert LB_down[node] == pytest.approx(node.LB_down)
    assert LBI[node] == pytest.approx(node.LBI)
    assert LBR[node] == pytest.approx(node.LBR) or (
        np.isnan(LBR[node]) and np.isnan(node.LBR)
    )
