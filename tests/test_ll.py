import gctree.branching_processes as bp
import gctree.phylip_parse as pp


trees = pp.parse_outfile(
    "tests/small_outfile", root="GL", abundance_file="tests/abundances.csv"
)
forest = bp.CollapsedForest([bp.CollapsedTree(tree) for tree in trees])


def test_ll_genotype_cache():
    p, q = 0.4, 0.6
    c_max, m_max = 10, 10
    bp.CollapsedTree._build_ll_genotype_cache(c_max, m_max, p, q)


def test_forest_likelihoods():
    # make sure ll doesn't have underflow
    p, q = 0.4, 0.6
    forest.ll(p, q, marginal=True)
