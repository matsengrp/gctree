#! /usr/bin/env python

import gctree.branching_processes as bp
import gctree.phylip_parse as pp
import gctree.mutation_model as mm
import gctree.utils as utils


import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
import pickle
import random
import ete3
import itertools


def test(args):
    """Test subprogram.

    Checks likelihood against a by-hand calculation for a simple tree,
    simulates a forest, computes MLE parameters, and plots some sanity
    check figures to plot_file command line arguments are p, q, number
    of trees to simulate, and plot file name
    """

    import seaborn as sns

    sns.set(style="white", color_codes=True)
    sns.set_style("ticks")

    # compare likelihood to empirical likelihood (might need large n)
    n = 10000
    df = pd.DataFrame(columns=("p", "q", "parameters", "f", "L"))
    ct = 0
    ps = (0.1, 0.2, 0.3, 0.4)
    qs = (0.2, 0.4, 0.6, 0.8)
    np.seterr(all="ignore")
    for p in ps:
        for q in qs:
            forest = bp.CollapsedForest()
            if args.verbose:
                print(f"parameters: p = {p}, q = {q}")
            forest.simulate(p, q, n)
            tree_dict = {}
            for tree in forest:
                tree_hash = tuple(
                    (node.abundance, len(node.children))
                    for node in tree.tree.traverse()
                )
                if tree_hash not in tree_dict:
                    tree_dict[tree_hash] = [tree, tree.ll(p, q)[0], 1]
                else:
                    tree_dict[tree_hash][-1] += 1
            L_empirical, L_theoretical = zip(
                *[
                    (tree_dict[tree_hash][2], np.exp(tree_dict[tree_hash][1]))
                    for tree_hash in tree_dict
                    if tree_dict[tree_hash][2]
                ]
            )
            for tree_hash in tree_dict:
                df.loc[ct] = (
                    p,
                    q,
                    f"p={p}, q={q}",
                    tree_dict[tree_hash][2],
                    np.exp(tree_dict[tree_hash][1]),
                )
                ct += 1
    if args.verbose:
        print()

    plt.figure()
    limx = (1 / n, 1.1)
    limy = (1, 1.1 * n)
    g = sns.lmplot(
        x="L",
        y="f",
        col="p",
        row="q",
        hue="parameters",
        data=df,
        fit_reg=False,
        scatter_kws={"alpha": 0.3},
        height=1.5,
        legend=False,
        row_order=reversed(qs),
    )
    g.set(xscale="log", yscale="log", xlim=limx, ylim=limy)
    for i in range(len(ps)):
        for j in range(len(qs)):
            g.axes[i, j].plot(
                limx, limy, ls="--", c="black", lw=0.5, zorder=0, markeredgewidth=0.1
            )
            g.axes[i, j].set_title(
                f"p={ps[j]}\nq={list(reversed(qs))[i]}",
                x=0.05,
                y=0.7,
                size="x-small",
                ha="left",
            )
    g.set_axis_labels("", "")
    g.fig.text(0.45, 0.02, s="gctree likelihood", multialignment="center")
    g.fig.text(
        0.02,
        0.7,
        s=f"frequency among {n} simulations",
        rotation=90,
        multialignment="center",
    )
    plt.savefig(args.outbase + "." + args.img_type)

    # MLE check
    n = 20
    n2 = 500
    df = pd.DataFrame(columns=("true parameters", "$\\hat{p}$", "$\\hat{q}$"))
    ct = 0
    for p in ps:
        for q in qs:
            for _ in range(n):
                forest = bp.CollapsedForest()
                if args.verbose:
                    print(f"parameters: p = {p}, q = {q}")
                forest.simulate(p, q, n2)
                p_inferred, q_inferred = forest.mle()
                df.loc[ct] = (f"p={p}, q={q}", p_inferred, q_inferred)
                ct += 1

    plt.figure()
    g = sns.lmplot(
        x="$\\hat{p}$",
        y="$\\hat{q}$",
        hue="true parameters",
        data=df,
        fit_reg=False,
        scatter_kws={"alpha": 0.2},
        height=6,
        legend=False,
    )
    g.set(
        xlim=(0.05, 0.45),
        xticks=np.arange(0.0, 0.6, 0.1),
        ylim=(0.1, 0.9),
        yticks=np.arange(0.0, 1.2, 0.2),
    )
    for i in range(len(ps)):
        for j in range(len(qs)):
            plt.scatter([ps[i]], [qs[j]], c="black", marker="+")
    plt.savefig(args.outbase + ".2." + args.img_type)

    return


def infer(args):
    """inference subprogram."""

    def isotype_add(forest):
        forest.add_isotypes(
            isotypemap_file=args.isotype_mapfile,
            idmap_file=args.idmapfile,
            isotype_names=args.isotype_names,
        )

    if len(args.infiles) == 2:
        forest = bp.CollapsedForest(
            *pp.parse_outfile(args.infiles[0], args.infiles[1], args.root)
        )
        if forest.n_trees == 1:
            warnings.warn("only one parsimony tree reported from dnapars")

        if args.verbose:
            print("number of trees with integer branch lengths:", forest.n_trees)
        forest.mle(marginal=True)
        # Add isotypes to forest
        if args.isotype_mapfile:
            isotype_add(forest)
        with open(args.outbase + ".inference.parsimony_forest.p", "wb") as f:
            pickle.dump(forest, f)

    elif len(args.infiles) == 1:
        if args.verbose:
            print(
                "Loading provided parsimony forest. If forest has fit parameters, parameter fitting will be skipped."
            )
        with open(args.infiles[0], "rb") as fh:
            forest = pickle.load(fh)
        # Add isotypes to forest and re-pickle if necessary
        if args.isotype_mapfile:
            isotype_add(forest)
            with open(args.outbase + ".inference.parsimony_forest.p", "wb") as f:
                pickle.dump(forest, f)
    else:
        raise ValueError(
            "The filename of a pickled history DAG object, or a phylipfile and abundance file, are required."
        )

    # Filter the forest according to specified criteria, and along the way,
    # write a log file containing stats for all trees in the forest:
    trimmed_forest, _ = forest.filter_trees(
        ranking_coeffs=args.ranking_coeffs,
        verbose=args.verbose,
        outbase=args.outbase,
        summarize_forest=args.summarize_forest,
        tree_stats=args.tree_stats,
        mutability_file=args.mutability,
        substitution_file=args.substitution,
        chain_split=args.chain_split,
    )

    if args.verbose:
        if trimmed_forest.n_trees > 1:
            n_topologies = trimmed_forest.n_topologies()
            print(
                "Degenerate ranking criteria: trimmed history DAG contains "
                f"{trimmed_forest.n_trees} unique trees, with {n_topologies} unique collapsed topologies."
            )
            if n_topologies > 10:
                print(
                    "A representative from each of the ten most abundant topologies will be sampled randomly for rendering."
                )
            else:
                print(
                    "A representative of each topology will be sampled randomly for rendering."
                )

    random.seed(forest.n_trees)
    ctrees = [
        topoclass.sample_tree()
        for _, topoclass in zip(range(10), trimmed_forest.iter_topology_classes())
    ]

    if args.colormapfile is not None:
        with open(args.colormapfile, "r") as f:
            colormap = {}
            for line in f:
                seqid, color = line.rstrip().split("\t")
                if "," in color:
                    colors = {
                        x.split(":")[0]: int(x.split(":")[1]) for x in color.split(",")
                    }
                    colormap[seqid] = colors
                else:
                    colormap[seqid] = color
    else:
        colormap = None

    # parse position map file(s)
    if args.positionmapfile is not None:
        with open(args.positionmapfile) as f:
            position_map = [int(x) for x in f.read().split()]
    else:
        position_map = None

    if args.positionmapfile2 is not None:
        with open(args.positionmapfile2) as f:
            position_map2 = [int(x) for x in f.read().split()]
    else:
        position_map2 = None

    for j, collapsed_tree in enumerate(ctrees, 1):
        collapsed_tree.render(
            f"{args.outbase}.inference.{j}.svg",
            idlabel=args.idlabel,
            colormap=colormap,
            frame=args.frame,
            position_map=position_map,
            chain_split=args.chain_split,
            frame2=args.frame2,
            position_map2=position_map2,
        )
        collapsed_tree.newick(f"{args.outbase}.inference.{j}.nk")
        with open(f"{args.outbase}.inference.{j}.p", "wb") as f:
            pickle.dump(collapsed_tree, f)

    # rank plot of observed allele frequencies
    y = sorted(
        (node.abundance for node in ctrees[0].tree.traverse() if node.abundance != 0),
        reverse=True,
    )
    plt.figure()
    plt.bar(range(1, len(y) + 1), y, color="black")
    plt.xlabel("genotype")
    plt.ylabel("abundance")
    plt.savefig(args.outbase + ".inference.abundance_rank." + args.img_type)


def simulate(args):
    """Simulation subprogram.

    Simulates a Galton–Watson process, with mutation probabilities
    according to a user defined motif model e.g. S5F
    """
    random.seed(a=args.seed)
    mutation_model = mm.MutationModel(args.mutability, args.substitution)
    if args.lambda0 is None:
        args.lambda0 = [max([1, int(0.01 * len(args.sequence))])]
    args.sequence = args.sequence.upper()
    if args.sequence2 is not None:
        # Use the same mutation rate on both sequences
        if len(args.lambda0) == 1:
            args.lambda0 = [args.lambda0[0], args.lambda0[0]]
        elif len(args.lambda0) != 2:
            raise Exception(
                "Only one or two lambda0 can be defined for a two "
                "sequence simulation."
            )
        # Require both sequences to be in frame 1:
        if args.frame is not None and args.frame != 1:
            if args.verbose:
                print(
                    "Warning: When simulating with two sequences they are "
                    "truncated to be beginning at frame 1."
                )
            args.sequence = args.sequence[
                (args.frame - 1) : (
                    args.frame
                    - 1
                    + (3 * (((len(args.sequence) - (args.frame - 1)) // 3)))
                )
            ]
            args.sequence2 = args.sequence2[
                (args.frame - 1) : (
                    args.frame
                    - 1
                    + (3 * (((len(args.sequence2) - (args.frame - 1)) // 3)))
                )
            ]
        # Extract the bounds between sequence 1 and 2:
        seq_bounds = (
            (0, len(args.sequence)),
            (len(args.sequence), len(args.sequence) + len(args.sequence2)),
        )
        # Merge the two seqeunces to simplify future dealing with the pair:
        args.sequence += args.sequence2
    else:
        seq_bounds = None

    trials = 1000
    # this loop makes us resimulate if size too small, or backmutation
    for trial in range(trials):
        try:
            tree = mutation_model.simulate(
                args.sequence,
                seq_bounds=seq_bounds,
                progeny=lambda seq: args.lambda_,
                lambda0=args.lambda0,
                n=args.n,
                N=args.N,
                T=args.T,
                frame=args.frame,
                verbose=args.verbose,
            )
            # this will fail if backmutations
            collapsed_tree = bp.CollapsedTree(tree=tree)
            tree.ladderize()
            uniques = sum(node.abundance > 0 for node in collapsed_tree.tree.traverse())
            if uniques < 2:
                raise RuntimeError(
                    f"collapsed tree contains {uniques} " "sampled sequences"
                )
            break
        except RuntimeError as e:
            print(f"{e}, trying again")
        else:
            raise
    if trial == trials - 1:
        raise RuntimeError(f"{trials} attempts exceeded")

    # In the case of a sequence pair print them to separate files:
    if args.sequence2 is not None:
        fh1 = open(args.outbase + ".simulation_seq1.fasta", "w")
        fh2 = open(args.outbase + ".simulation_seq2.fasta", "w")
        fh1.write(">root\n")
        fh1.write(args.sequence[seq_bounds[0][0] : seq_bounds[0][1]] + "\n")
        fh2.write(">root\n")
        fh2.write(args.sequence[seq_bounds[1][0] : seq_bounds[1][1]] + "\n")
        for leaf in tree.iter_leaves():
            if leaf.abundance != 0:
                fh1.write(">" + leaf.name + "\n")
                fh1.write(leaf.sequence[seq_bounds[0][0] : seq_bounds[0][1]] + "\n")
                fh2.write(">" + leaf.name + "\n")
                fh2.write(leaf.sequence[seq_bounds[1][0] : seq_bounds[1][1]] + "\n")
    else:
        with open(args.outbase + ".simulation.fasta", "w") as f:
            f.write(">root\n")
            f.write(args.sequence + "\n")
            for leaf in tree.iter_leaves():
                if leaf.abundance != 0:
                    f.write(">" + leaf.name + "\n")
                    f.write(leaf.sequence + "\n")

    # some observable simulation stats to write
    abundance, distance_from_root, degree = zip(
        *[
            (
                node.abundance,
                utils.hamming_distance(node.sequence, args.sequence),
                sum(
                    utils.hamming_distance(node.sequence, node2.sequence) == 1
                    for node2 in collapsed_tree.tree.traverse()
                    if node2.abundance and node2 is not node
                ),
            )
            for node in collapsed_tree.tree.traverse()
            if node.abundance
        ]
    )
    stats = pd.DataFrame(
        {
            "genotype abundance": abundance,
            "Hamming distance to root genotype": distance_from_root,
            "Hamming neighbor genotypes": degree,
        }
    )
    stats.to_csv(args.outbase + ".simulation.stats.tsv", sep="\t", index=False)

    print(
        f"{sum(leaf.abundance for leaf in collapsed_tree.tree.traverse())}"
        " simulated observed sequences"
    )

    # render the full lineage tree
    ts = ete3.TreeStyle()
    ts.rotation = 90
    ts.show_leaf_name = False
    ts.show_scale = False

    colors = {}
    palette = ete3.SVG_COLORS
    palette -= set(["black", "white", "gray"])
    palette = itertools.cycle(list(palette))  # <-- circular iterator

    colors[tree.sequence] = "gray"

    for n in tree.traverse():
        nstyle = ete3.NodeStyle()
        nstyle["size"] = 10
        if args.plotAA:
            if n.AAseq not in colors:
                colors[n.AAseq] = next(palette)
            nstyle["fgcolor"] = colors[n.AAseq]
        else:
            if n.sequence not in colors:
                colors[n.sequence] = next(palette)
            nstyle["fgcolor"] = colors[n.sequence]
        n.set_style(nstyle)

    # this makes the rendered branch lenths correspond to time
    for node in tree.iter_descendants():
        node.dist = node.time - node.up.time
    tree.render(args.outbase + ".simulation.lineage_tree.svg", tree_style=ts)

    # render collapsed tree
    # create an id-wise colormap
    # NOTE: node.name can be a set
    colormap = {
        node.name: colors[node.sequence] for node in collapsed_tree.tree.traverse()
    }
    collapsed_tree.write(args.outbase + ".simulation.collapsed_tree.p")
    collapsed_tree.render(
        args.outbase + ".simulation.collapsed_tree.svg",
        idlabel=args.idlabel,
        colormap=colormap,
        frame=args.frame,
    )
    # print colormap to file
    with open(args.outbase + ".simulation.collapsed_tree.colormap.tsv", "w") as f:
        for name, color in colormap.items():
            f.write(
                (name if isinstance(name, str) else ",".join(name))
                + "\t"
                + color
                + "\n"
            )


def get_parser():
    parser = argparse.ArgumentParser(
        description="genotype collapsed tree inference and simulation"
    )
    subparsers = parser.add_subparsers(
        title="subcommands",
        description="specify one of these",
        required=True,
        help="additional help available for each subcommand",
    )

    # parser for test subprogram
    parser_test = subparsers.add_parser(
        "test",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="run tests on library functions",
    )
    parser_test.set_defaults(func=test)

    # parser for inference subprogram
    parser_infer = subparsers.add_parser(
        "infer",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="likelihood ranking of parsimony trees",
    )
    parser_infer.add_argument(
        "--root",
        type=str,
        default="root",
        help=r'name of root sequence (outgroup root), default ``"root"``',
    )
    parser_infer.add_argument(
        "infiles",
        type=str,
        nargs="+",
        help=(
            "Input files for inference. If two filenames are passed, the first shall be a "
            "dnapars outfile (verbose output with sequences at each site), and the second "
            "shall be an abundance file containing allele frequencies (sequence counts) in "
            "the format: ``SeqID, Nobs``. If a single filename is passed, it shall be the name "
            "of a pickled history DAG object created by gctree. A new pickled forest will be "
            "output only if new annotations (such as isotypes) are added."
        ),
    )
    parser_infer.add_argument(
        "--colormapfile",
        type=str,
        default=None,
        help='File containing color map in the tab-separated format: ``"SeqID\tcolor"``',
    )
    parser_infer.add_argument(
        "--chain_split",
        type=int,
        default=None,
        help=(
            "when using concatenated heavy and light chains, this is the 0-based"
            " index at which the 2nd chain begins, needed for determining coding frame in both chains,"
            " and also to correctly calculate mutability parsimony."
        ),
    )
    parser_infer.add_argument(
        "--frame2",
        type=int,
        default=None,
        choices=(1, 2, 3),
        help=r"codon frame for the second chain when using the ``chain_split`` option",
    )
    parser_infer.add_argument(
        "--positionmapfile",
        type=str,
        default=None,
        help="file containing a list of position numbers (e.g. IMGT) corresponding to indices in sequence",
    )
    parser_infer.add_argument(
        "--positionmapfile2",
        type=str,
        default=None,
        help="positionmapfile for the 2nd chain when using the ``chain_split`` option",
    )
    parser_infer.set_defaults(func=infer)
    parser_infer.add_argument(
        "--idmapfile",
        default=None,
        type=str,
        help=(
            "input filename for a csv file mapping sequence names to original sequence ids. "
            "For use by isotype ranking. "
            "Such a file can be produced by ``deduplicate`` when it is provided the ``--idmapfile`` option. "
        ),
    )
    parser_infer.add_argument(
        "--isotype_mapfile",
        default=None,
        type=str,
        help=(
            "filename for a csv file mapping original sequence ids to observed isotypes"
            ". For example, each line should have the format 'somesequence_id, some_isotype'."
        ),
    )
    parser_infer.add_argument(
        "--isotype_names",
        type=str,
        default=None,
        nargs="+",
        help=(
            "A list of isotype names used in isotype_mapfile, in order of most naive to most differentiated."
            """ Default is equivalent to providing the argument ``--isotype_names IgM IgD IgG3 IgG1 IgG2 IgE IgA``"""
        ),
    )
    parser_infer.add_argument(
        "--mutability",
        type=str,
        default=None,
        help=(
            "Path to mutability model file. If --mutability and --substitution are both provided, "
            "they will be used to rank trees after likelihood and isotype parsimony."
        ),
    )
    parser_infer.add_argument(
        "--substitution",
        type=str,
        default=None,
        help=(
            "Path to substitution model file. If --mutability and --substitution are both provided, "
            "they will be used to rank trees after likelihood and isotype parsimony."
        ),
    )
    parser_infer.add_argument(
        "--ranking_coeffs",
        type=float,
        nargs=3,
        default=None,
        help=(
            "List of coefficients for ranking trees by a linear combination of traits. "
            "Coefficients are in order: isotype parsimony, mutation model parsimony, number of alleles. "
            "A coefficient of -1 will be applied to branching process likelihood. "
            "If not provided, trees will be ranked lexicographically by likelihood, "
            "isotype parsimony, and mutability parsimony in that order."
        ),
    )
    parser_infer.add_argument(
        "--summarize_forest",
        action="store_true",
        help=(
            "write a file `[outbase].forest_summary.log` with a summary of traits for trees in the forest."
        ),
    )
    parser_infer.add_argument(
        "--tree_stats",
        action="store_true",
        help=(
            "write a file `[outbase].tree_stats.log` with stats for all trees in the forest. "
            "For large forests, this is slow and memory intensive."
        ),
    )

    # parser for simulation subprogram
    parser_sim = subparsers.add_parser(
        "simulate",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Neutral, and target selective, model gctree simulation",
    )
    parser_sim.add_argument("sequence", type=str, help="seed root nucleotide sequence")
    parser_sim.add_argument(
        "mutability", type=str, help="path to mutability model file"
    )
    parser_sim.add_argument(
        "substitution", type=str, help="path to substitution model file"
    )
    parser_sim.add_argument(
        "--sequence2",
        type=str,
        default=None,
        help="Second seed root nucleotide sequence. "
        "For simulating heavy/light chain co-evolution.",
    )
    parser_sim.add_argument(
        "--lambda",
        dest="lambda_",
        type=float,
        default=0.9,
        help="poisson branching parameter",
    )
    parser_sim.add_argument(
        "--lambda0",
        type=float,
        default=None,
        nargs="*",
        help="List of one or two elements with the baseline mutation rates. "
        "Space separated input values. First element belonging to seed "
        "sequence one and optionally the next to sequence 2. If only one "
        "rate is provided for two sequences, this rate will be used on "
        "both.",
    )
    parser_sim.add_argument("--n", type=int, default=None, help="cells downsampled")
    parser_sim.add_argument(
        "--N", type=int, default=None, help="target simulation size"
    )
    parser_sim.add_argument(
        "--T",
        type=int,
        nargs="+",
        default=None,
        help="observation time, if None we run until termination and take all"
        " leaves",
    )
    parser_sim.add_argument(
        "--seed", type=int, default=None, help="integer random seed"
    )
    parser_sim.add_argument(
        "--target_dist",
        type=int,
        default=10,
        help="The number of non-synonymous mutations the target should be "
        "away from the root.",
    )
    parser_sim.add_argument(
        "--plotAA",
        type=bool,
        default=False,
        help="Plot trees with collapsing and coloring on amino acid level.",
    )
    parser_sim.set_defaults(func=simulate)

    # shared parameters
    for subparser in [parser_test, parser_infer, parser_sim]:
        subparser.add_argument(
            "--outbase", type=str, default="gctree.out", help="output file base name"
        )
        subparser.add_argument(
            "--img_type", type=str, default="svg", help="output image file type"
        )
        subparser.add_argument(
            "--verbose", action="store_true", help="flag for verbose messaging"
        )

    # common parameters for the inference and simulation subprograms
    for subparser in [parser_infer, parser_sim]:
        subparser.add_argument(
            "--frame", type=int, default=None, choices=(1, 2, 3), help="codon frame"
        )
        subparser.add_argument(
            "--idlabel",
            action="store_true",
            help=r"label nodes with their sequence ids in output tree images, "
            "and write a fasta alignment mapping those sequence ids to sequences. "
            "This is the easiest way to access inferred ancestral sequences.",
        )

    return parser


def main(arg_list=None):
    args = get_parser().parse_args(arg_list)
    args.func(args)
