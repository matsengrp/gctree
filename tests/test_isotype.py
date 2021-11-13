import gctree.isotype as iso
import ete3

testtrees = [
             ('((((A2)?,(G2,G2)?)?,G2,(A2,A2)?)?)M;', 3.0),
             ('((((A2)?,(G2,G2)?)?,A2,(A2,A2)?)?)M;', 4.0),
             ('((((A2)?,(G2,G3)?)?,G2,(A2,A2)?)?)M;', 5.0),
             ('((((A2)?,(G2,A2)?)?,G2,(A2,A2)?)?)M;', 4.0),
             ]

def test_isotype_disambiguate():
	for newick, weight in testtrees:
		tree = ete3.Tree(newick, format=8)
		for node in tree.traverse():
			node.isotype = iso.Isotype(node.name)
		iso.disambiguate_isotype(tree)
		assert sum(iso.isotype_distance(node.up.isotype, node.isotype) for node in tree.iter_descendants()) == weight
