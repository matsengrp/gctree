from gctree.isotyping import _isotype_dagfuncs, _isotype_annotation_dagfuncs
import gctree.branching_processes as bp
import gctree.isotyping as iso
import gctree.phylip_parse as pp
import ete3

testtrees = [
    ("((((A2)?,(G2,G2)?)?,G2,(A2,A2)?)?)M;", 3.0),
    ("((((A2)?,(G2,G2)?)?,A2,(A2,A2)?)?)M;", 4.0),
    ("((((A2)?,(G2,G3)?)?,G2,(A2,A2)?)?)M;", 5.0),
    ("((((A2)?,(G2,A2)?)?,G2,(A2,A2)?)?)M;", 4.0),
]

trees_seqcounts1 = pp.parse_outfile(
    "tests/example_output/original/small_outfile",
    abundance_file="tests/example_output/original/abundances.csv",
    root="GL",
)

forest = bp.CollapsedForest(trees_seqcounts1)
dag = forest._forest


def test_isotype_disambiguate():
    newisotype = iso.IsotypeTemplate(["M", "G3", "A1", "G2", "G4", "E", "A2"]).new
    for newick, weight in testtrees:
        tree = ete3.Tree(newick, format=8)
        for node in tree.traverse():
            node.isotype = newisotype(node.name)
        iso._disambiguate_isotype(tree)
        assert (
            sum(
                iso.isotype_distance(node.up.isotype, node.isotype)
                for node in tree.iter_descendants()
            )
            == weight
        )


def test_trim_byisotype():
    kwargs = _isotype_annotation_dagfuncs(
        isotypemap_file="example/isotypemap.txt",
        idmap_file="tests/example_output/original/idmap.txt",
    )
    # isotypemap_file=None,
    # idmap=None,
    # idmap_file=None,
    # isotype_names=None,
    tdag = dag.copy()
    tdag.optimal_weight_annotate(**kwargs, optimal_func=lambda l: l[0])  # noqa: E741
    for node in tdag.preorder():
        if node.attr is not None:
            node.attr["isotype"] = node._dp_data
    dag_filter = _isotype_dagfuncs()
    c = tdag.weight_count(**dag_filter)
    key = min(c)
    count = c[key]
    tdag.trim_optimal_weight(**dag_filter)
    assert tdag.weight_count(**dag_filter) == {key: count}
