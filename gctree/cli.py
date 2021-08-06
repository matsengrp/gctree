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
import scipy.stats
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
            for tree in forest.forest:
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
    outfiles = [pp.parse_outfile(args.phylipfile, args.abundance_file, args.root)]
    if args.bootstrap_phylipfile is not None:
        outfiles.extend(
            pp.parse_outfile(args.bootstrap_phylipfile, args.abundance_file, args.root)
        )
    bootstrap = len(outfiles) > 1
    if bootstrap:
        # store bootstrap parameter data
        df = pd.DataFrame(columns=("$\\hat{p}$", "$\\hat{q}$"))
        # we'll store the mle gctrees here for computing support later
        gctrees = []

    for i, content in enumerate(outfiles):
        if i > 0:
            if args.verbose:
                print(f"bootstrap sample {i}")
                print("----------------------")
            outbase = args.outbase + f".bootstrap_{i}"
        else:
            outbase = args.outbase
        phylip_collapsed = [
            bp.CollapsedTree(tree=tree, allow_repeats=(i > 0)) for tree in content
        ]
        phylip_collapsed_unique = []
        for tree in phylip_collapsed:
            if (
                sum(
                    tree.compare(tree2, method="identity")
                    for tree2 in phylip_collapsed_unique
                )
                == 0
            ):
                phylip_collapsed_unique.append(tree)

        parsimony_forest = bp.CollapsedForest(forest=phylip_collapsed_unique)

        if parsimony_forest.n_trees == 1:
            warnings.warn("only one parsimony tree reported from dnapars")

        if args.verbose:
            print(
                "number of trees with integer branch lengths:", parsimony_forest.n_trees
            )

        # check for unifurcations at root
        unifurcations = sum(
            tree.tree.abundance == 0 and len(tree.tree.children) == 1
            for tree in parsimony_forest.forest
        )
        if unifurcations and args.verbose:
            print(
                f"{unifurcations} trees exhibit unobserved unifurcation from"
                " root. Adding psuedocounts to these roots"
            )

        # fit p and q using all trees
        # if we get floating point errors, try a few more times
        # (starting params are random)
        max_tries = 10
        for tries in range(max_tries):
            try:
                p, q = parsimony_forest.mle(marginal=True)
                break
            except FloatingPointError as e:
                if tries + 1 < max_tries and args.verbose:
                    print(
                        f"floating point error in MLE: {e}. "
                        f"Attempt {tries + 1} of {max_tries}. "
                        "Rerunning with new random start."
                    )
                else:
                    raise
            else:
                raise

        if args.verbose:
            print(f"params: {(p, q)}")

        if i > 0:
            df.loc[i - 1] = p, q

        # get likelihoods and sort by them
        ls = [tree.ll(p, q)[0] for tree in parsimony_forest.forest]
        ls, parsimony_forest.forest = zip(
            *sorted(zip(ls, parsimony_forest.forest), key=lambda x: x[0], reverse=True)
        )

        if bootstrap:
            gctrees.append(parsimony_forest.forest[0])

        with open(outbase + ".inference.parsimony_forest.p", "wb") as f:
            pickle.dump(parsimony_forest, f)

        if args.colormapfile is not None:
            with open(args.colormapfile, "r") as f:
                colormap = {}
                for line in f:
                    seqid, color = line.rstrip().split("\t")
                    if "," in color:
                        colors = {
                            x.split(":")[0]: int(x.split(":")[1])
                            for x in color.split(",")
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

        if args.verbose:
            print("tree\talleles\tlogLikelihood")
        for j, (l, collapsed_tree) in enumerate(zip(ls, parsimony_forest.forest), 1):
            alleles = sum(1 for _ in collapsed_tree.tree.traverse())
            if args.verbose:
                print(f"{j}\t{alleles}\t{l}")
            collapsed_tree.render(
                f"{outbase}.inference.{j}.svg",
                idlabel=args.idlabel,
                colormap=colormap,
                frame=args.frame,
                position_map=position_map,
                chain_split=args.chain_split,
                frame2=args.frame2,
                position_map2=position_map2,
            )
            collapsed_tree.newick(f"{outbase}.inference.{j}.nk")
            with open(f"{outbase}.inference.{j}.p", "wb") as f:
                pickle.dump(collapsed_tree, f)

        # rank plot of likelihoods
        plt.figure(figsize=(6.5, 2))
        try:
            plt.plot(np.exp(ls), "ko", clip_on=False, markersize=4)
            plt.ylabel("gctree likelihood")
            plt.yscale("log")
            plt.ylim([None, 1.1 * max(np.exp(ls))])
        except FloatingPointError:
            plt.plot(ls, "ko", clip_on=False, markersize=4)
            plt.ylabel("gctree log-likelihood")
            plt.ylim([None, 1.1 * max(ls)])
        plt.xlabel("parsimony tree")
        plt.xlim([-1, len(ls)])
        plt.tick_params(axis="y", direction="out", which="both")
        plt.tick_params(
            axis="x", which="both", bottom="off", top="off", labelbottom="off"
        )
        plt.savefig(outbase + ".inference.likelihood_rank." + args.img_type)

        # rank plot of observed allele frequencies
        y = sorted(
            (
                node.abundance
                for node in parsimony_forest.forest[0].tree.traverse()
                if node.abundance != 0
            ),
            reverse=True,
        )
        plt.figure()
        plt.bar(range(1, len(y) + 1), y, color="black")
        plt.xlabel("genotype")
        plt.ylabel("abundance")
        plt.savefig(outbase + ".inference.abundance_rank." + args.img_type)

    if bootstrap:
        import seaborn as sns

        sns.set(style="white", color_codes=True)
        sns.set_style("ticks")
        plt.figure()
        np.seterr(all="ignore")
        sns.jointplot(
            "$\\hat{p}$",
            "$\\hat{q}$",
            data=df,
            joint_kws={"alpha": 0.1, "s": 0.1},
            space=0,
            color="k",
            stat_func=None,
            xlim=(0, 1),
            ylim=(0, 1),
            height=3,
        )
        plt.savefig(args.outbase + ".inference.bootstrap_theta." + args.img_type)
        gctrees[0].support(gctrees[1:])
        gctrees[0].render(
            args.outbase + ".inference.bootstrap_support.svg",
            colormap=colormap,
            idlabel=args.idlabel,
            frame=args.frame,
            position_map=position_map,
            show_support=True,
            chain_split=args.chain_split,
            frame2=args.frame2,
            position_map2=position_map2,
        )
        gctrees[0].support(gctrees[1:], compatibility=True)
        gctrees[0].render(
            args.outbase + ".inference.bootstrap_compatibility.svg",
            colormap=colormap,
            idlabel=args.idlabel,
            frame=args.frame,
            position_map=position_map,
            show_support=True,
            chain_split=args.chain_split,
            frame2=args.frame2,
            position_map2=position_map2,
        )


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
                progeny=scipy.stats.poisson(args.lambda_),
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
        "phylipfile",
        type=str,
        help="dnapars outfile (verbose output with sequences at each site)",
    )
    parser_infer.add_argument(
        "abundance_file",
        type=str,
        help="File containing allele frequencies (sequence counts) in the "
        'format: "SeqID,Nobs"',
    )
    parser_infer.add_argument(
        "--bootstrap_phylipfile",
        type=str,
        help="dnapars outfile from seqboot (multiple data sets)",
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
        help="when using concatenated heavy and light chains, this is the 0-based index at which the 2nd chain begins, needed for determining coding frame in both chains",
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
            help=r"flag for labeling the sequence ids of the nodes in the "
            "output tree images, also write associated fasta alignment "
            "if ``True``",
        )

    return parser


def main(arg_list=None):
    args = get_parser().parse_args(arg_list)
    args.func(args)
