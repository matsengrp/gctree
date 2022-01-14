import gctree.mutation_model as mm
import gctree.isotyping as itp
import gctree.branching_processes as bp
import gctree.utils as utils

from functools import lru_cache
import historydag as hdag
from multiset import FrozenMultiset
import matplotlib.pyplot as plt
from typing import List, Tuple
from numpy import exp, log


def filter_dag(
    dag,
    priority_weights=None,
    verbose: bool = False,
    outbase: str = None,
    mutability_file: str = None,
    substitution_file: str = None,
    isotypemap=None,
    isotypemap_file=None,
    idmap=None,
    idmap_file=None,
    isotype_names=None,
    img_type="svg",
):
    """Trim the provided history DAG, containing sequence labels, by-node abundance attributes,
    and a `parameters` attribute, to minimize a linear combination of
    branching process likelihood, isotype parsimony score, mutability parsimony score, and number of alleles,
    with coefficients provided in `priority_weights`, in that order.

    Args:
        dag: A history DAG object to trim
        priority_weights: A sequence of coefficients for prioritizing tree weights
        verbose: print information about trimming
        outbase: file name stem for a file with information for each tree in the DAG, and rank plots of likelihoods.
            If not provided, no such files will be written.
        mutability_file: A mutability model
        substitution_file: A substitution model
        isotypemap: A mapping of sequences to isotypes
        isotypemap: A dictionary mapping original IDs to observed isotype names
        isotypemap_file: A csv file providing an `isotypemap`
        idmap: A dictionary mapping unique sequence IDs to sets of original IDs of observed sequences
        idmap_file: A csv file providing an `idmap`
        isotype_names: A sequence of isotype names containing values in `isotypemap`, in the correct switching order
        img_type: Format to be used for output plots, if `outbase` string is provided.

    Returns:
        The trimmed history DAG, containing all optimal trees according to specified criteria.
    """
    p, q = dag.attr["parameters"]
    placeholder_dagfuncs = hdag.utils.AddFuncDict(
        {
            "start_func": lambda n: 0,
            "edge_weight_func": lambda n1, n2: 0,
            "accum_func": sum,
        },
        names="PlaceholderWeight",
    )
    try:
        isotype_dagfuncs = make_isotype_dagfuncs(
            isotypemap_file=isotypemap_file,
            isotypemap=isotypemap,
            idmap_file=idmap_file,
            idmap=idmap,
            isotype_names=isotype_names,
        )
        if verbose:
            print("Isotype parsimony will be used as a ranking criterion")
    except ValueError:
        isotype_dagfuncs = placeholder_dagfuncs
    ll_dagfuncs = ll_cmcount_dagfuncs(p, q)
    if mutability_file and substitution_file:
        if verbose:
            print("Mutation model parsimony will be used as a ranking criterion")
        mutability_dagfuncs = make_mutability_dagfuncs(
            mutability_file=mutability_file, substitution_file=substitution_file
        )
    else:
        mutability_dagfuncs = placeholder_dagfuncs
    allele_dagfuncs = hdag.utils.AddFuncDict(
        {
            "start_func": lambda n: 0,
            "edge_weight_func": lambda n1, n2: n1.label != n2.label,
            "accum_func": sum,
        },
        names="NodeCount",
    )

    if priority_weights:
        # a vector of relative weights, for ll, isotype pars, mutability_pars, alleles
        # This is possible because all weights are additive up the tree,
        # and the transformations are strictly increasing/decreasing.
        kwargls = [ll_dagfuncs, isotype_dagfuncs, mutability_dagfuncs, allele_dagfuncs]
        minmaxls = [
            (
                dag.optimal_weight_annotate(**kwargs, optimal_func=min),
                dag.optimal_weight_annotate(**kwargs, optimal_func=max),
            )
            for kwargs in kwargls
        ]
        minmaxls = [
            (mn, mx) if isinstance(kwargs.names, str) else (mn[0], mx[0])
            for (mn, mx), kwargs in zip(minmaxls, kwargls)
        ]

        def transformation(mn, mx):
            if mx - mn < 0.0001:

                def func(n):
                    return n - mn

            else:

                def func(n):
                    return (n - mn) / (mx - mn)

            return func

        transformations = [transformation(*minmax) for minmax in minmaxls]
        # switch first transformation so that we maximize ls, not minimize:
        f = transformations[0]
        transformations[0] = lambda n: -f(n) + 1

        def minfunckey(weighttuple):
            ll, _, isotypepars, _, mutabilitypars, alleles = weighttuple
            weightlist = [ll, isotypepars, mutabilitypars, alleles]
            return sum(
                [
                    priority * transform(weight)
                    for priority, transform, weight in zip(
                        priority_weights, transformations, weightlist
                    )
                ]
            )

    else:

        def minfunckey(weighttuple):
            ll, _, isotypepars, _, mutabilitypars, alleles = weighttuple
            # Sort output by likelihood, then isotype parsimony, then mutability score
            return (-ll, isotypepars, mutabilitypars)

    # Filter by likelihood, isotype parsimony, mutability,
    # and make ctrees, cforest, and render trees
    dagweight_kwargs = (
        ll_dagfuncs + isotype_dagfuncs + mutability_dagfuncs + allele_dagfuncs
    )
    trimdag = dag.copy()
    trimdag.trim_optimal_weight(
        **dagweight_kwargs, optimal_func=lambda l: min(l, key=minfunckey)
    )

    if outbase:
        dag_ls = list(dag.weight_count(**dagweight_kwargs).elements())
        # To clear _dp_data fields of their large cargo
        dag.optimal_weight_annotate(**placeholder_dagfuncs)
        dag_ls.sort(key=minfunckey)

        with open(outbase + ".dag_summary.log", "w") as fh:
            fh.write(f"Parameters: {(p, q)}\n")
            fh.write(
                "tree\talleles\tlogLikelihood\t\tisotype_parsimony\tmutability_parsimony"
                + ("\ttree score" if priority_weights else "")
                + "\n"
            )
            for j, (l, _, isotypepars, _, mutabilitypars, alleles) in enumerate(
                dag_ls, 1
            ):
                if priority_weights:
                    treescore = str(
                        minfunckey((l, 0, isotypepars, 0, mutabilitypars, alleles))
                    )
                else:
                    treescore = ""
                fh.write(
                    f"{j}\t{alleles}\t{l}\t{isotypepars}\t\t\t{mutabilitypars}\t{treescore}\n"
                )

        # For use in later plots
        dag_l = [ls[0] for ls in dag_ls]

        # rank plot of likelihoods
        plt.figure(figsize=(6.5, 2))
        try:
            plt.plot(exp(dag_l), "ko", clip_on=False, markersize=4)
            plt.ylabel("gctree likelihood")
            plt.yscale("log")
            plt.ylim([None, 1.1 * max(exp(dag_l))])
        except FloatingPointError:
            plt.plot(dag_l, "ko", clip_on=False, markersize=4)
            plt.ylabel("gctree log-likelihood")
            plt.ylim([None, 1.1 * max(dag_l)])
        plt.xlabel("parsimony tree")
        plt.xlim([-1, len(dag_l)])
        plt.tick_params(axis="y", direction="out", which="both")
        plt.tick_params(
            axis="x", which="both", bottom="off", top="off", labelbottom="off"
        )
        plt.savefig(outbase + ".inference.likelihood_rank." + img_type)

    if verbose:
        print(f"params: {(p, q)}")
        print("stats for optimal trees:")
        print(
            "alleles\tlogLikelihood\t\tisotype_parsimony\tmutability_parsimony"
            + ("\ttreescore" if priority_weights else "")
        )
        # with format (ll, _, isotypepars, _, mutabilitypars, alleles)
        stats = trimdag.optimal_weight_annotate(
            **dagweight_kwargs, optimal_func=lambda l: min(l, key=minfunckey)
        )
        print(
            f"{stats[5]}\t{stats[0]}\t{stats[2]}\t\t\t{stats[4]}\t"
            + str(minfunckey(stats) if priority_weights else "")
        )

    return trimdag


def ll_cmcount_dagfuncs(p, q):
    """A slower but more numerically stable equivalent to ll_genotype_dagfuncs.
    To use these functions for DAG trimming, use an optimal function like
    `lambda l: max(l, key=lambda ll: ll[0])` for clarity, although min or max should work too."""

    funcdict = bp.cmcounter_dagfuncs()
    cmcount_weight_func = funcdict["edge_weight_func"]
    cmcount_addfunc = funcdict["accum_func"]

    @lru_cache(maxsize=None)
    def ll(cmcounter: FrozenMultiset) -> float:
        if cmcounter:
            if (0, 1) in cmcounter:
                cmcounter = cmcounter - {(0, 1)} + {(1, 1)}
            return bp.lltree(tuple(cmcounter.items()), p, q)[0]
        else:
            return 0.0

    def edge_weight_func(n1, n2):
        st = cmcount_weight_func(n1, n2)
        return (ll(st), st)

    def accum_func(cmsetlist: List[Tuple[float, FrozenMultiset]]):
        st = cmcount_addfunc([t[1] for t in cmsetlist])
        return (ll(st), st)

    return hdag.utils.AddFuncDict(
        {
            "start_func": lambda n: (0, FrozenMultiset()),
            "edge_weight_func": edge_weight_func,
            "accum_func": accum_func,
        },
        names=("LogLikelihood", "cmcounters"),
    )


def ll_genotype_dagfuncs(p, q):
    """Functions for counting tree log likelihood on the history DAG. Although these
    functions are fast, for numerical consistency use :meth:`ll_cmcount_dagfuncs` instead."""

    def edge_weight_ll_genotype(n1, n2):
        """The _ll_genotype weight of the target node, unless it should be collapsed, then 0.
        Expects DAG to have abundances added with :meth:`dag.add_abundances`."""
        if n2.is_leaf() and n2.label == n1.label:
            return 0.0
        else:
            m = len(n2.clades)
            # Check if this edge should be collapsed, and reduce mutant descendants
            if frozenset({n2.label}) in n2.clades:
                m -= 1
            return bp.CollapsedTree._ll_genotype(n2.attr["abundance"], m, p, q)[0]

    return hdag.utils.AddFuncDict(
        {
            "start_func": lambda n: 0,
            "edge_weight_func": edge_weight_ll_genotype,
            "accum_func": sum,
        }
    )


def make_mutability_dagfuncs(*args, **kwargs):
    """Returns a historydag.AddFuncDict containing functions for
    counting tree-wise sums of mutability distances in the history DAG.

    This may not be stable numerically, but we expect every tree to have
    a unique mutability parsimony score, for non-degenerate models, so
    so this shouldn't matter in practice.

    Arguments are passed to :meth:`MutationModel` constructor."""

    mutation_model = mm.MutationModel(*args, **kwargs)
    dist = mutability_distance(mutation_model)

    @hdag.utils.access_field("label")
    @hdag.utils.ignore_ualabel(0)
    @hdag.utils.access_field("sequence")
    def distance(seq1, seq2):
        return dist(seq1, seq2)

    return hdag.utils.AddFuncDict(
        {"start_func": lambda n: 0, "edge_weight_func": distance, "accum_func": sum},
        names="MutabilityParsimony",
    )


def _mutability_distance_precursors(mutation_model):
    # Caching could be moved to the MutationModel class instead.
    context_model = mutation_model.context_model.copy()
    k = mutation_model.k
    h = k // 2
    # Build all sequences with (when k=5) one or two Ns on either end
    templates = [
        ("N" * left, "N" * (k - left - right), "N" * right)
        for left in range(h + 1)
        for right in range(h + 1)
        if left != 0 or right != 0
    ]

    kmers_to_compute = [
        leftns + stub + rightns
        for leftns, ambig_stub, rightns in templates
        for stub in mm.sequence_disambiguations(ambig_stub)
    ]
    # Cache all these mutabilities in context_model also
    context_model.update(
        {kmer: mutation_model.mutability(kmer) for kmer in kmers_to_compute}
    )

    @utils.check_distance_arguments
    def mutpairs(seq1, seq2):
        ns = "N" * h
        seq1N = ns + seq1 + ns
        seq2N = ns + seq2 + ns
        mut_idxs = [
            index
            for index, (base1, base2) in enumerate(zip(seq1N, seq2N))
            if base1 != base2
        ]
        return FrozenMultiset((seq1N[i - h : i + h + 1], seq2N[i]) for i in mut_idxs)

    def sum_minus_logp(pairs):
        # I have chosen to multiply substitution rate for central base
        # with rate of new base. Not sure if this is a good choice.
        if pairs:
            # for floating point behavior
            pairs = sorted(pairs.items())
            p_arr = [
                mult
                * (log(context_model[mer][0]) + log(context_model[mer][1][newbase]))
                for (mer, newbase), mult in pairs
            ]
            return -sum(p_arr)
        else:
            return 0.0

    return (mutpairs, sum_minus_logp)


def mutability_distance(mutation_model):
    """Returns a fast distance function based on mutability_model.
    First, caches computed mutabilities for k-mers with k // 2 N's on either
    end. This is pretty fast for k=5, but the distance function should be created
    once and reused."""
    mutpairs, sum_minus_logp = _mutability_distance_precursors(mutation_model)

    def distance(seq1, seq2):
        return sum_minus_logp(mutpairs(seq1, seq2))

    return distance


def make_isotype_dagfuncs(
    isotypemap=None,
    isotypemap_file=None,
    idmap=None,
    idmap_file=None,
    isotype_names=None,
):
    """Returns a historydag.utils.AddFuncDict which contains functions necessary for filtering
    by isotype parsimony score. These functions return weights which are a tuple containing
    a progressive isotype score, and a frozenset containing isotypes used internally to compute
    isotype scores.

    Args:
        isotypemap: A dictionary mapping original IDs to observed isotype names
        isotypemap_file: A csv file providing an `isotypemap`
        idmap: A dictionary mapping unique sequence IDs to sets of original IDs of observed sequences
        idmap_file: A csv file providing an `idmap`
        isotype_names: A sequence of isotype names containing values in `isotypemap`, in the correct switching order
    """
    if isotypemap_file and isotypemap is None:
        with open(isotypemap_file, "r") as fh:
            isotypemap = dict(map(lambda x: x.strip(), line.split(",")) for line in fh)
    elif isotypemap_file is None:
        raise ValueError("Either isotypemap or isotypemap_file is required")

    if idmap_file and idmap is None:
        with open(idmap_file, "r") as fh:
            idmap = {}
            for line in fh:
                seqid, cell_ids = line.rstrip().split(",")
                cell_idset = {cell_id for cell_id in cell_ids.split(":") if cell_id}
                if len(cell_idset) > 0:
                    idmap[seqid] = cell_idset
        newidmap = itp.explode_idmap(idmap, isotypemap)
        if isotype_names:
            isotype_names = str(isotype_names).split(",")
        else:
            isotype_names = itp.default_isotype_order
        newisotype = itp.IsotypeTemplate(isotype_names, weight_matrix=None).new
    elif idmap is None:
        raise TypeError("Either idmap or idmap_file is required")

    def distance_func(n1, n2):
        if n2.is_leaf():
            seqid = n2.attr["name"]
            if seqid in newidmap:
                isoset = frozenset(newisotype(name) for name in newidmap[seqid].keys())
            else:
                isoset = frozenset()
        else:
            isoset = frozenset()
        return (0, isoset)

    def sumweights(weights):
        ps = [i[0] for i in weights]
        ts = [item for i in weights for item in i[1]]
        if ts:
            newt = frozenset({min(ts)})
        else:
            newt = frozenset()
        return (
            sum(itp.isotype_distance(list(newt)[0], el) for el in ts) + sum(ps),
            newt,
        )

    return hdag.utils.AddFuncDict(
        {
            "start_func": lambda n: (0, frozenset()),
            "edge_weight_func": distance_func,
            "accum_func": sumweights,
        },
        names=("IsotypeParsimony", "Isotype_state"),
    )
