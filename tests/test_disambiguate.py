import random
import ete3
import gctree.phylip_parse as pps

newick_tree1 = (
    "((((12[&&NHX:name=12:sequence=T])4[&&NHX:name=4:sequence=C],"
    "(6[&&NHX:name=6:sequence=C],"
    "7[&&NHX:name=7:sequence=A])5[&&NHX:name=5:sequence=M])"
    "3[&&NHX:name=3:sequence=M],8[&&NHX:name=8:sequence=A],"
    "(11[&&NHX:name=11:sequence=A],10[&&NHX:name=10:sequence=G])"
    "9[&&NHX:name=9:sequence=R])2[&&NHX:name=2:sequence=R])"
    "1[&&NHX:name=1:sequence=G];"
)

newick_tree2 = (
    "((((12[&&NHX:name=12:sequence=T])4[&&NHX:name=4:sequence=C],"
    "(6[&&NHX:name=6:sequence=C],"
    "7[&&NHX:name=7:sequence=A])5[&&NHX:name=5:sequence=?])"
    "3[&&NHX:name=3:sequence=?],8[&&NHX:name=8:sequence=A],"
    "(11[&&NHX:name=11:sequence=A],10[&&NHX:name=10:sequence=G])"
    "9[&&NHX:name=9:sequence=?])2[&&NHX:name=2:sequence=?])"
    "1[&&NHX:name=1:sequence=G];"
)


def treeprint(tree: ete3.TreeNode):
    tree = tree.copy()
    for node in tree.traverse():
        node.name = node.sequence
    return tree.write(format=8)


def sample(treestring, n=20):
    ogt = ete3.TreeNode(newick=newick_tree2, format=1)
    print(ogt.sequence)
    newickset = set()
    for i in range(n):
        t = ogt.copy()
        random.seed(i)
        t = pps.disambiguate(t, random_state=random.getstate())
        newickset.add(treeprint(t))
    return newickset


def test_full_ambiguity():
    newickset = sample(newick_tree2)
    correctset = {
        "((((T)C,(C,A)C)C,A,(A,G)G)G);",
        "((((T)C,(C,A)A)A,A,(A,G)A)A);",
        "((((T)C,(C,A)C)C,A,(A,G)A)A);",
    }
    if not newickset == correctset:
        missing = correctset - newickset
        wrong = newickset - correctset
        print(
            f"\nDisambiguate function is missing {missing}\n"
            f"and came up with these incorrect trees: {wrong}"
        )
        raise ValueError("Invalid Disambiguation")


def test_restricted_ambiguity():
    newickset = sample(newick_tree1)
    correctset = {
        "((((T)C,(C,A)C)C,A,(A,G)G)G);",
        "((((T)C,(C,A)A)A,A,(A,G)A)A);",
        "((((T)C,(C,A)C)C,A,(A,G)A)A);",
    }
    if not newickset == correctset:
        missing = correctset - newickset
        wrong = newickset - correctset
        print(
            f"\nDisambiguate function is missing {missing}\n"
            f"and came up with these incorrect trees: {wrong}"
        )
        raise ValueError("Invalid Disambiguation")
