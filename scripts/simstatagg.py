#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""aggregation plots of metrics from several simulation runs with same
parameters."""

from matplotlib import pyplot as plt
import pandas as pd
import argparse
from gctree.deduplicate import fasta_parse
from gctree.utils import hamming_distance
import seaborn as sns

sns.set(style="white", color_codes=True)


parser = argparse.ArgumentParser(
    description="aggregate validation of repeated runs with same parameters"
)
parser.add_argument(
    "input", type=str, nargs="+", help="gctree.simulation.stats.tsv files"
)
parser.add_argument("--outbase", type=str, help="output file base name")
parser.add_argument(
    "--experimental",
    type=str,
    help="experimental .fasta (like Tas et al. data)",
    default=None,
)
args = parser.parse_args()

for i, fname in enumerate(args.input):
    df = pd.read_csv(fname, sep="\t")
    df["simulation"] = i + 1
    if i == 0:
        aggdat = df
    else:
        aggdat = aggdat.append(df, ignore_index=True)

sims = set(aggdat["simulation"])
nsims = len(sims)

if args.experimental is not None:
    new_aln, counts = fasta_parse(args.experimental, naive="GL", id_abundances=True)[:2]
    exp_dict = {seq.id: str(seq.seq) for seq in new_aln}
    naive_id = [seq for seq in exp_dict if "gl" in seq][0]
    frequency, distance_from_naive, degree = zip(
        *[
            (
                counts[seq],
                hamming_distance(exp_dict[seq], exp_dict[naive_id]),
                sum(
                    hamming_distance(exp_dict[seq], exp_dict[seq2]) == 1
                    for seq2 in exp_dict
                    if seq2 is not seq and counts[seq2] != 0
                ),
            )
            for seq in exp_dict
            if counts[seq] != 0
        ]
    )
    exp_stats = pd.DataFrame(
        {
            "genotype abundance": frequency,
            "Hamming distance to root genotype": distance_from_naive,
            "Hamming neighbor genotypes": degree,
        }
    )

# bw = .3
alpha = min([0.9, 20 / nsims])
bins = range(
    max(
        aggdat["Hamming distance to root genotype"].max(),
        exp_stats["Hamming distance to root genotype"].max()
        if args.experimental is not None
        else 0,
    )
    + 2
)

plt.figure(figsize=(6, 3))
plt.subplot(1, 2, 1)
for simulation, simulation_aggdat in aggdat.groupby("simulation"):
    sns.distplot(
        simulation_aggdat["Hamming distance to root genotype"],
        bins=bins,
        kde=False,
        hist_kws={"histtype": "step", "cumulative": True, "alpha": alpha, "lw": 1},
    )
if args.experimental is not None:
    sns.distplot(
        exp_stats["Hamming distance to root genotype"],
        bins=bins,
        kde=False,
        hist_kws={
            "histtype": "step",
            "cumulative": True,
            "color": "k",
            "lw": 3,
            "alpha": 0.8,
        },
    )
plt.xlabel("Hamming distance to root genotype")
plt.xlim([0, bins[-1]])
plt.ylabel("observed genotypes")
plt.tight_layout()

plt.subplot(1, 2, 2)
xbins = range(
    max(
        aggdat["genotype abundance"].max(),
        exp_stats["genotype abundance"].max() if args.experimental is not None else 0,
    )
    + 2
)
ybins = range(
    max(
        aggdat["Hamming neighbor genotypes"].max(),
        exp_stats["Hamming neighbor genotypes"].max()
        if args.experimental is not None
        else 0,
    )
    + 2
)
for simulation, simulation_aggdat in aggdat.groupby("simulation"):
    plt.plot(
        simulation_aggdat["genotype abundance"],
        simulation_aggdat["Hamming neighbor genotypes"],
        "+",
        mew=1,
        alpha=alpha / 2,
    )
if args.experimental is not None:
    plt.plot(
        exp_stats["genotype abundance"],
        exp_stats["Hamming neighbor genotypes"],
        "o",
        mew=2,
        alpha=0.8,
        markerfacecolor="none",
        color="k",
    )
plt.xlabel("genotype abundance")
plt.ylabel("Hamming neighbor genotypes")
plt.xscale("symlog")
plt.yscale("symlog")
plt.xlim([0.9, None])
plt.ylim([-0.1, None])
plt.tight_layout()
plt.savefig(args.outbase + ".pdf")
