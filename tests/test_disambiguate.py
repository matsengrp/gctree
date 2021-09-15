import random
import ete3
from gctree.utils import hamming_distance
import gctree.phylip_parse as pps

newick_tree1 = ("((((12[&&NHX:name=12:sequence=T])4[&&NHX:name=4:sequence=C],"
                "(6[&&NHX:name=6:sequence=C],"
                "7[&&NHX:name=7:sequence=A])5[&&NHX:name=5:sequence=M])"
                "3[&&NHX:name=3:sequence=M],8[&&NHX:name=8:sequence=A],"
                "(11[&&NHX:name=11:sequence=A],10[&&NHX:name=10:sequence=G])"
                "9[&&NHX:name=9:sequence=R])2[&&NHX:name=2:sequence=R])"
                "1[&&NHX:name=1:sequence=G];")

def treeprint(tree: ete3.TreeNode):
    for node in tree.traverse():
        node.name = node.sequence
    return(tree.write(format=8))

def total_weight(tree: ete3.TreeNode) -> float:
    return(sum(node.dist for node in tree.traverse()))

def test_optimal():
    ogt = ete3.TreeNode(newick=newick_tree1, format=1)
    for i in range(1000):
        t = ogt.copy()
        random.seed(i)
        print(random.randint(1, 5))
        t = pps.disambiguate(t, random_state=random.getstate())
        for node in t.iter_descendants():
            node.dist = hamming_distance(node.up.sequence, node.sequence)
        print(treeprint(t))
        if not total_weight(t) == 5.0:
            raise ValueError(f"A tree was found with weight {total_weight(t)}")
        
        # elif t.write(format=8) == '((((T)C,(C,A)C)C,A,(A,G)G)G);':
        #     raise ValueError("Illegal inheritance of disambiguation")
        #     pass
    raise ValueError("Hi!!")

def test_exhaustive():
    ogt = ete3.TreeNode(newick=newick_tree1, format=1)
    newickset = set()
    for i in range(100):
        t = ogt.copy()
        random.seed(i)
        t = pps.disambiguate(t, random_state=random.getstate())
        for node in t.traverse():
            try:
                cv = list(node.cv)
            except AttributeError:
                cv = 'None'
            print('\n', node.name, cv)
        newickset.add(treeprint(t))
        print(treeprint(t))
    correctset = {"((((T)C,(C,A)C)C,A,(A,G)G)G);",
                  "((((T)C,(C,A)A)A,A,(A,G)G)G);",
                  "((((T)C,(C,A)C)C,A,(A,G)A)A);"}
    if not newickset == correctset:
        missing = correctset - newickset
        wrong = newickset - correctset
        print(f"\nDisambiguate function is missing {missing}\n"
              f"and came up with these incorrect trees: {wrong}")
        raise ValueError("Invalid Disambiguation")

