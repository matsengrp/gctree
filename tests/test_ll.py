import gctree.branching_processes as bp

def test_ll_genotype_cache():
    p, q = 0.4, 0.6
    c_max, m_max = 10, 10
    bp.CollapsedTree._build_ll_genotype_cache(c_max, m_max, p, q)
