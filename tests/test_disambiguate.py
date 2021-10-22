from gctree import mutation_model, utils
import random
import ete3
import gctree.phylip_parse as pps

newick_tree1 = ("((((12[&&NHX:name=12:sequence=T])4[&&NHX:name=4:sequence=C],"
                "(6[&&NHX:name=6:sequence=C],"
                "7[&&NHX:name=7:sequence=A])5[&&NHX:name=5:sequence=M])"
                "3[&&NHX:name=3:sequence=M],8[&&NHX:name=8:sequence=A],"
                "(11[&&NHX:name=11:sequence=A],10[&&NHX:name=10:sequence=G])"
                "9[&&NHX:name=9:sequence=R])2[&&NHX:name=2:sequence=R])"
                "1[&&NHX:name=1:sequence=G];")

newick_tree2 = ("((((12[&&NHX:name=12:sequence=T])4[&&NHX:name=4:sequence=C],"
                "(6[&&NHX:name=6:sequence=C],"
                "7[&&NHX:name=7:sequence=A])5[&&NHX:name=5:sequence=?])"
                "3[&&NHX:name=3:sequence=?],8[&&NHX:name=8:sequence=A],"
                "(11[&&NHX:name=11:sequence=A],10[&&NHX:name=10:sequence=G])"
                "9[&&NHX:name=9:sequence=?])2[&&NHX:name=2:sequence=?])"
                "1[&&NHX:name=1:sequence=G];")

tree1 = ete3.TreeNode(newick=newick_tree1, format=1)
for node in tree1.traverse():
    node.sequence += node.sequence
tree2 = ete3.TreeNode(newick=newick_tree2, format=1)
tree3 = ete3.TreeNode(newick=newick_tree1, format=1)
for node in tree3.traverse():
    node.sequence += "A" * 5
tree3.sequence = "GAAAAR"


def treeprint(tree: ete3.TreeNode):
    tree = tree.copy()
    for node in tree.traverse():
        node.name = node.sequence
    return(tree.write(format=8))


def sample(tree, n=20, **kwargs):
    print(tree.sequence)
    newickset = set()
    for i in range(n):
        t = tree.copy()
        random.seed(i)
        t = pps.disambiguate(t, random_state=random.getstate(), **kwargs)
        print(treeprint(t))
        newickset.add(treeprint(t))
    return(newickset)


def test_full_ambiguity():
    newickset = sample(tree2)
    correctset = {"((((T)C,(C,A)C)C,A,(A,G)G)G);",
                  "((((T)C,(C,A)A)A,A,(A,G)A)A);",
                  "((((T)C,(C,A)C)C,A,(A,G)A)A);"}
    if not newickset == correctset:
        missing = correctset - newickset
        wrong = newickset - correctset
        print(f"\nDisambiguate function is missing {missing}\n"
              f"and came up with these incorrect trees: {wrong}")
        raise ValueError("Invalid Disambiguation")


def test_restricted_ambiguity():
    newickset = sample(tree1)
    correctset = {'((((TT)CC,(CC,AA)CC)CC,AA,(AA,GG)GG)GG);',
                  '((((TT)CC,(CC,AA)AC)AC,AA,(AA,GG)AG)AG);',
                  '((((TT)CC,(CC,AA)AC)AC,AA,(AA,GG)AA)AA);',
                  '((((TT)CC,(CC,AA)CC)CC,AA,(AA,GG)AA)AA);',
                  '((((TT)CC,(CC,AA)CA)CA,AA,(AA,GG)GA)GA);',
                  '((((TT)CC,(CC,AA)CC)CC,AA,(AA,GG)GA)GA);',
                  '((((TT)CC,(CC,AA)CC)CC,AA,(AA,GG)AG)AG);',
                  '((((TT)CC,(CC,AA)AA)AA,AA,(AA,GG)AA)AA);',
                  '((((TT)CC,(CC,AA)CA)CA,AA,(AA,GG)AA)AA);'}
    if not newickset == correctset:
        missing = correctset - newickset
        wrong = newickset - correctset
        print(f"\nDisambiguate function is missing {missing}\n"
              f"and came up with these incorrect trees: {wrong}")
        raise ValueError("Invalid Disambiguation")


def test_restricted_ambiguity_widewindow_mutability():
    mmodel = mutation_model.MutationModel(mutability_file="S5F/Mutability.csv", substitution_file="S5F/Substitution.csv")
    newickset = sample(tree1, n=5, dist_func=utils.mutability_distance(mmodel), distance_dependence=2)
    correctset = {'((((TT)CC,(CC,AA)CA)CA,AA,(AA,GG)GA)GA);'}
    if not newickset == correctset:
        missing = correctset - newickset
        wrong = newickset - correctset
        print(f"\nDisambiguate function is missing {missing}\n"
              f"and came up with these incorrect trees: {wrong}")
        raise ValueError("Invalid Disambiguation")

def test_restricted_ambiguity_widewindow():
    newickset = sample(tree1, n=100, distance_dependence=-1)
    correctset = {'((((TT)CC,(CC,AA)CC)CC,AA,(AA,GG)GG)GG);',
                  '((((TT)CC,(CC,AA)AC)AC,AA,(AA,GG)AG)AG);',
                  '((((TT)CC,(CC,AA)AC)AC,AA,(AA,GG)AA)AA);',
                  '((((TT)CC,(CC,AA)CC)CC,AA,(AA,GG)AA)AA);',
                  '((((TT)CC,(CC,AA)CA)CA,AA,(AA,GG)GA)GA);',
                  '((((TT)CC,(CC,AA)CC)CC,AA,(AA,GG)GA)GA);',
                  '((((TT)CC,(CC,AA)CC)CC,AA,(AA,GG)AG)AG);',
                  '((((TT)CC,(CC,AA)AA)AA,AA,(AA,GG)AA)AA);',
                  '((((TT)CC,(CC,AA)CA)CA,AA,(AA,GG)AA)AA);'}
    if not newickset == correctset:
        missing = correctset - newickset
        wrong = newickset - correctset
        print(f"\nDisambiguate function is missing {missing}\n"
              f"and came up with these incorrect trees: {wrong}")
        raise ValueError("Invalid Disambiguation")
