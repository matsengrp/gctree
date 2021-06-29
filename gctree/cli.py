#! /usr/bin/env python

import gctree.branching_processes as bp
import gctree.phylip_parse as pp
import gctree.mutation_model as mm
import gctree.selection_utils as su
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
            forest = bp.CollapsedForest((p, q), n)
            print(f"parameters: p = {p}, q = {q}")
            forest.simulate()
            tree_dict = {}
            for tree in forest.forest:
                tree_hash = tuple(
                    (node.frequency, len(node.children))
                    for node in tree.tree.traverse()
                )
                if tree_hash not in tree_dict:
                    tree_dict[tree_hash] = [tree, tree.ll((p, q))[0], 1]
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
    g.fig.text(0.45, 0.02, s="GCtree likelihood", multialignment="center")
    g.fig.text(
        0.02,
        0.7,
        s=f"frequency among {n} simulations",
        rotation=90,
        multialignment="center",
    )
    plt.savefig(args.outbase + ".pdf")

    # MLE check
    n = 20
    n2 = 500
    df = pd.DataFrame(columns=("true parameters", "$\\hat{p}$", "$\\hat{q}$"))
    ct = 0
    for p in ps:
        for q in qs:
            for _ in range(n):
                forest = bp.CollapsedForest((p, q), n2)
                print(f"parameters: p = {p}, q = {q}")
                forest.simulate()
                result = forest.mle()
                df.loc[ct] = (f"p={p}, q={q}", result.x[0], result.x[1])
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
    plt.savefig(args.outbase + ".2.pdf")

    return


def infer(args):
    """inference subprogram."""
    outfiles = [pp.parse_outfile(args.phylipfile, args.countfile, args.naive)]
    if args.bootstrap_phylipfile is not None:
        outfiles.extend(
            pp.parse_outfile(args.bootstrap_phylipfile, args.countfile, args.naive)
        )
    bootstrap = len(outfiles) > 1
    if bootstrap:
        # store bootstrap parameter data
        df = pd.DataFrame(columns=("$\\hat{p}$", "$\\hat{q}$"))
        # we'll store the mle gctrees here for computing support later
        gctrees = []

    for i, content in enumerate(outfiles):
        if i > 0:
            print(f"bootstrap sample {i}")
            print("----------------------")
            outbase = args.outbase + f".bootstrap_{i}"
        else:
            outbase = args.outbase
        phylip_collapsed = [
            bp.CollapsedTree(tree=tree, frame=args.frame, allow_repeats=(i > 0))
            for tree in content
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

        print("number of trees with integer branch lengths:", parsimony_forest.n_trees)

        # check for unifurcations at root
        unifurcations = sum(
            tree.tree.frequency == 0 and len(tree.tree.children) == 1
            for tree in parsimony_forest.forest
        )
        if unifurcations:
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
                parsimony_forest.mle(empirical_bayes_sum=True)
                break
            except FloatingPointError as e:
                if tries + 1 < max_tries:
                    print(
                        f"floating point error in MLE: {e}. "
                        f"Attempt {tries + 1} of {max_tries}. "
                        "Rerunning with new random start."
                    )
                else:
                    raise
            else:
                raise

        print(f"params = {parsimony_forest.params}")

        if i > 0:
            df.loc[i - 1] = parsimony_forest.params

        # get likelihoods and sort by them
        ls = [
            tree.ll(parsimony_forest.params, build_cache=False)[0]
            for tree in parsimony_forest.forest
        ]
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

        print("tree\talleles\tlogLikelihood")
        for j, (l, collapsed_tree) in enumerate(zip(ls, parsimony_forest.forest), 1):
            alleles = sum(1 for _ in collapsed_tree.tree.traverse())
            print(f"{j}\t{alleles}\t{l}")
            collapsed_tree.render(
                f"{outbase}.inference.{j}.svg",
                idlabel=args.idlabel,
                colormap=colormap,
                chain_split=args.chain_split,
            )
            collapsed_tree.newick(f"{outbase}.inference.{j}.nk")

        # rank plot of likelihoods
        plt.figure(figsize=(6.5, 2))
        try:
            plt.plot(np.exp(ls), "ko", clip_on=False, markersize=4)
            plt.ylabel("GCtree likelihood")
            plt.yscale("log")
            plt.ylim([None, 1.1 * max(np.exp(ls))])
        except FloatingPointError:
            plt.plot(ls, "ko", clip_on=False, markersize=4)
            plt.ylabel("GCtree log-likelihood")
            plt.ylim([None, 1.1 * max(ls)])
        plt.xlabel("parsimony tree")
        plt.xlim([-1, len(ls)])
        plt.tick_params(axis="y", direction="out", which="both")
        plt.tick_params(
            axis="x", which="both", bottom="off", top="off", labelbottom="off"
        )
        plt.savefig(outbase + ".inference.likelihood_rank.pdf")

        # rank plot of observed allele frequencies
        y = sorted(
            (
                node.frequency
                for node in parsimony_forest.forest[0].tree.traverse()
                if node.frequency != 0
            ),
            reverse=True,
        )
        plt.figure()
        plt.bar(range(1, len(y) + 1), y, color="black")
        plt.xlabel("genotype")
        plt.ylabel("abundance")
        plt.savefig(outbase + ".inference.abundance_rank.pdf")

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
        plt.savefig(args.outbase + ".inference.bootstrap_theta.pdf")
        gctrees[0].support(gctrees[1:])
        gctrees[0].render(
            args.outbase + ".inference.bootstrap_support.svg",
            colormap=colormap,
            idlabel=args.idlabel,
            chain_split=args.chain_split,
            show_support=True,
        )
        gctrees[0].support(gctrees[1:], compatibility=True)
        gctrees[0].render(
            args.outbase + ".inference.bootstrap_compatibility.svg",
            colormap=colormap,
            idlabel=args.idlabel,
            chain_split=args.chain_split,
            show_support=True,
        )


def simulate(args):
    """Simulation subprogram.

    Can simulate in two modes. a) Neutral mode. A Galton–Watson process,
    with mutation probabilities according to a user defined motif model
    e.g. S5F b) Selection mode. Using the same mutation process as in
    a), but in selection mode the poisson progeny distribution's lambda
    is variable according to the hamming distance to a list of target
    sequences. The closer a sequence gets to one of the targets the
    higher fitness and the closer lambda will approach 2, vice versa
    when the sequence is far away lambda approaches 0.
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
        args.sequence += args.sequence2.upper()
    else:
        seq_bounds = None
    if args.selection:
        if args.frame is None:
            raise Exception("Frame must be defined when simulating with " "selection.")
        # the fully activating fraction on BA must be possible to reach within
        # B_total
        assert args.B_total >= args.f_full
        # Make a list of target sequences:
        targetAAseqs = [
            mutation_model.one_mutant(args.sequence, args.target_dist, frame=args.frame)
            for i in range(args.target_count)
        ]
        # Find the total amount of A necessary for sustaining the inputted
        # carrying capacity:
        print((args.carry_cap, args.B_total, args.f_full, args.mature_affy))
        A_total = su.find_A_total(
            args.carry_cap, args.B_total, args.f_full, args.mature_affy, args.U
        )
        # Calculate the parameters for the logistic function:
        Lp = su.find_Lp(args.f_full, args.U)
        selection_params = [
            args.stop_dist,
            args.mature_affy,
            args.naive_affy,
            args.target_dist,
            args.skip_update,
            targetAAseqs,
            A_total,
            args.B_total,
            Lp,
            args.k,
            args.outbase,
        ]
    else:
        selection_params = None

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
                selection_params=selection_params,
            )
            if args.selection:
                collapsed_tree = bp.CollapsedTree(
                    tree=tree, frame=args.frame, collapse_syn=False, allow_repeats=True
                )
            else:
                # this will fail if backmutations
                collapsed_tree = bp.CollapsedTree(tree=tree, frame=args.frame)
            tree.ladderize()
            uniques = sum(node.frequency > 0 for node in collapsed_tree.tree.traverse())
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
        fh1.write(">naive\n")
        fh1.write(args.sequence[seq_bounds[0][0] : seq_bounds[0][1]] + "\n")
        fh2.write(">naive\n")
        fh2.write(args.sequence[seq_bounds[1][0] : seq_bounds[1][1]] + "\n")
        for leaf in tree.iter_leaves():
            if leaf.frequency != 0:
                fh1.write(">" + leaf.name + "\n")
                fh1.write(leaf.sequence[seq_bounds[0][0] : seq_bounds[0][1]] + "\n")
                fh2.write(">" + leaf.name + "\n")
                fh2.write(leaf.sequence[seq_bounds[1][0] : seq_bounds[1][1]] + "\n")
    else:
        with open(args.outbase + ".simulation.fasta", "w") as f:
            f.write(">naive\n")
            f.write(args.sequence + "\n")
            for leaf in tree.iter_leaves():
                if leaf.frequency != 0:
                    f.write(">" + leaf.name + "\n")
                    f.write(leaf.sequence + "\n")

    # some observable simulation stats to write
    frequency, distance_from_naive, degree = zip(
        *[
            (
                node.frequency,
                utils.hamming_distance(node.sequence, args.sequence),
                sum(
                    utils.hamming_distance(node.sequence, node2.sequence) == 1
                    for node2 in collapsed_tree.tree.traverse()
                    if node2.frequency and node2 is not node
                ),
            )
            for node in collapsed_tree.tree.traverse()
            if node.frequency
        ]
    )
    stats = pd.DataFrame(
        {
            "genotype abundance": frequency,
            "Hamming distance to root genotype": distance_from_naive,
            "Hamming neighbor genotypes": degree,
        }
    )
    stats.to_csv(args.outbase + ".simulation.stats.tsv", sep="\t", index=False)

    print(
        f"{sum(leaf.frequency for leaf in collapsed_tree.tree.traverse())}"
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

    # Either plot by DNA sequence or amino acid sequence:
    if args.plotAA and args.selection:
        colors[tree.AAseq] = "gray"
    else:
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
    if args.plotAA and args.selection:
        colormap = {
            node.name: colors[node.AAseq] for node in collapsed_tree.tree.traverse()
        }
    else:
        colormap = {
            node.name: colors[node.sequence] for node in collapsed_tree.tree.traverse()
        }
    collapsed_tree.write(args.outbase + ".simulation.collapsed_tree.p")
    collapsed_tree.render(
        args.outbase + ".simulation.collapsed_tree.svg",
        idlabel=args.idlabel,
        colormap=colormap,
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

    if args.selection:
        # Define a list a suitable colors that are easy to distinguish:
        palette = [
            "crimson",
            "purple",
            "hotpink",
            "limegreen",
            "darkorange",
            "darkkhaki",
            "brown",
            "lightsalmon",
            "darkgreen",
            "darkseagreen",
            "darkslateblue",
            "teal",
            "olive",
            "wheat",
            "magenta",
            "lightsteelblue",
            "plum",
            "gold",
        ]
        # circular iterator
        palette = itertools.cycle(list(palette))
        colors = {i: next(palette) for i in range(int(len(args.sequence) // 3))}
        # The minimum distance to the target is colored:
        colormap = {
            node.name: colors[node.target_dist]
            for node in collapsed_tree.tree.traverse()
        }
        collapsed_tree.write(
            args.outbase + ".simulation.collapsed_runstat_color_tree.p"
        )
        collapsed_tree.render(
            args.outbase + ".simulation.collapsed_runstat_color_tree.svg",
            idlabel=args.idlabel,
            colormap=colormap,
        )
        # Write a file with the selection run stats. These are also plotted:
        with open(args.outbase + "selection_sim.runstats.p", "rb") as fh:
            runstats = pickle.load(fh)
            su.plot_runstats(runstats, args.outbase, colors)


def get_parser():
    parser = argparse.ArgumentParser(
        description="germinal center tree inference and simulation"
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
        "--naive", type=str, default=None, help="name of naive sequence (outgroup root)"
    )
    parser_infer.add_argument(
        "phylipfile",
        type=str,
        help="dnapars outfile (verbose output with sequences at each site)",
    )
    parser_infer.add_argument(
        "countfile",
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
        help='File containing color map in the format: "SeqID\tcolor"',
    )
    parser_infer.add_argument(
        "--chain_split",
        type=int,
        default=None,
        help="split between heavy and light for combined seqs",
    )
    parser_infer.set_defaults(func=infer)

    # parser for simulation subprogram
    parser_sim = subparsers.add_parser(
        "simulate",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        help="Neutral, and target selective, model gctree simulation",
    )
    parser_sim.add_argument("sequence", type=str, help="seed naive nucleotide sequence")
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
        help="Second seed naive nucleotide sequence. "
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
        "--selection",
        type=bool,
        default=False,
        help="Simulation with selection? true/false. When doing simulation "
        "with selection an observation time cut must be set.",
    )
    parser_sim.add_argument(
        "--stop_dist",
        type=int,
        default=None,
        help="Stop when this distance has been reached in the selection " "model.",
    )
    parser_sim.add_argument(
        "--carry_cap",
        type=int,
        default=1000,
        help="The carrying capacity of the simulation with selection. This "
        "number affects the fixation time of a new mutation. Fixation "
        "time is approx. log2(carry_cap), e.g. log2(1000) ~= 10.",
    )
    parser_sim.add_argument(
        "--target_count",
        type=int,
        default=10,
        help="The number of targets to generate.",
    )
    parser_sim.add_argument(
        "--target_dist",
        type=int,
        default=10,
        help="The number of non-synonymous mutations the target should be "
        "away from the naive.",
    )
    parser_sim.add_argument(
        "--naive_affy",
        type=float,
        default=100,
        help="Affinity of the naive sequence in nano molar.",
    )
    parser_sim.add_argument(
        "--mature_affy",
        type=float,
        default=1,
        help="Affinity of the mature sequences in nano molar.",
    )
    parser_sim.add_argument(
        "--skip_update",
        type=int,
        default=100,
        help="When iterating through the leafs the B:A fraction is "
        "recalculated every time. It is possible though to update less "
        "often and get the same approximate results. This parameter sets "
        "the number of iterations to skip, before updating the B:A "
        "results. skip_update < carry_cap/10 recommended.",
    )
    parser_sim.add_argument(
        "--B_total",
        type=float,
        default=1,
        help="Total number of BCRs per B cell normalized to 10e4. So 1 equals "
        "10e4, 100 equals 10e6 etc. It is recommended to keep this as "
        "the default.",
    )
    parser_sim.add_argument(
        "--U",
        type=float,
        default=5,
        help="Controls the fraction of BCRs binding antigen necessary to "
        "only sustain the life of the B cell It is recommended to keep "
        "this as the default.",
    )
    parser_sim.add_argument(
        "--f_full",
        type=float,
        default=1,
        help="The fraction of antigen bound BCRs on a B cell that is needed "
        "to elicit close to maximum reponse. Cannot be smaller than "
        "B_total. It is recommended to keep this as the default.",
    )
    parser_sim.add_argument(
        "--k",
        type=float,
        default=2,
        help="The exponent in the function to map hamming distance to "
        "affinity. It is recommended to keep this as the default.",
    )
    parser_sim.add_argument(
        "--plotAA",
        type=bool,
        default=False,
        help="Plot trees with collapsing and coloring on amino acid level.",
    )
    parser_sim.add_argument(
        "--verbose",
        type=bool,
        default=False,
        help="Print progress during simulation. Mostly useful for simulation "
        "with selection since this can take a while.",
    )
    parser_sim.set_defaults(func=simulate)

    # a common outbase parameter
    for subparser in [parser_test, parser_infer, parser_sim]:
        subparser.add_argument(
            "--outbase", type=str, default="gctree.out", help="output file base name"
        )

    # common parameters for the inference and simulation subprograms
    for subparser in [parser_infer, parser_sim]:
        subparser.add_argument(
            "--frame", type=int, default=None, choices=(1, 2, 3), help="codon frame"
        )
        subparser.add_argument(
            "--idlabel",
            action="store_true",
            help="flag for labeling the sequence ids of the nodes in the "
            "output tree images, also write associated fasta alignment "
            "if True",
        )

    return parser


def main(arg_list=None):
    args = get_parser().parse_args(arg_list)
    args.func(args)
