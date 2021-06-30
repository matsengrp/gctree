#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""summary stats of the seq data going into tree inference for CFT."""

from gctree.deduplicate import fasta_parse
from matplotlib import pyplot as plt
import pandas as pd
import argparse
import seaborn as sns

sns.set(style="white", color_codes=True)
from gctree.utils import hamming_distance


def main():
    parser = argparse.ArgumentParser(description="summary statistics of pre-tree data")
    parser.add_argument("input", type=str, nargs="+", help="simulated fasta files")
    parser.add_argument(
        "--experimental", type=str, nargs="*", help="experimental fasta files"
    )
    parser.add_argument("--outbase", type=str, help="output file base name")
    parser.add_argument(
        "--root_idexp", type=str, default="root0", help="root sequence ID"
    )
    args = parser.parse_args()

    # simulations
    root_id = "root"
    for i, fname in enumerate(args.input):
        print(fname)
        seqs = {seq.id: str(seq.seq) for seq in fasta_parse(fname, "root")[0]}
        nseqs = len(seqs)
        if nseqs <= 2:
            continue

        distance_from_root, degree = zip(
            *[
                (
                    hamming_distance(seqs[seqid], seqs[root_id]),
                    min(
                        hamming_distance(seqs[seqid], seqs[seqid2])
                        for seqid2 in seqs
                        if seqid2 != root_id and seqid2 != seqid
                    ),
                )
                for seqid in seqs
                if seqid != root_id
            ]
        )
        df = pd.DataFrame(
            {
                "distance to root sequence": distance_from_root,
                "nearest neighbor distance": degree,
            }
        )
        df["data set"] = i + 1
        if i == 0:
            aggdat = df
        else:
            aggdat = aggdat.append(df, ignore_index=True)

    ndatasets = len(set(aggdat["data set"]))

    # experimental
    if args.experimental is not None:
        for i, fname in enumerate(args.experimental):
            print(fname)
            seqs = {
                seq.id: str(seq.seq) for seq in fasta_parse(fname, args.root_idexp)[0]
            }
            nseqs = len(seqs)
            if nseqs <= 2:
                continue

            distance_from_root, degree = zip(
                *[
                    (
                        hamming_distance(seqs[seqid], seqs[args.root_idexp]),
                        min(
                            hamming_distance(seqs[seqid], seqs[seqid2])
                            for seqid2 in seqs
                            if seqid2 != args.root_idexp and seqid2 != seqid
                        ),
                    )
                    for seqid in seqs
                    if seqid != args.root_idexp
                ]
            )
            df = pd.DataFrame(
                {
                    "distance to root sequence": distance_from_root,
                    "nearest neighbor distance": degree,
                }
            )
            df["data set"] = i + 1
            if i == 0:
                aggdat_exp = df
            else:
                aggdat_exp = aggdat_exp.append(df, ignore_index=True)

        ndatasets += len(set(aggdat_exp["data set"]))

    # bw = .3
    alpha = min([0.9, 20 / ndatasets])
    bins = range(
        max(
            aggdat["distance to root sequence"].max(),
            aggdat_exp["distance to root sequence"].max()
            if args.experimental is not None
            else 0,
        )
        + 2
    )

    plt.figure(figsize=(6, 3))
    plt.subplot(1, 2, 1)
    ct = 0
    for dataset, dataset_aggdat in aggdat.groupby("data set"):
        ct += 1
        sns.distplot(
            dataset_aggdat["distance to root sequence"],
            bins=bins,
            kde=False,
            color="gray",
            hist_kws={"histtype": "step", "cumulative": True, "alpha": alpha, "lw": 1},
        )

    if args.experimental is not None:
        for dataset, dataset_aggdat in aggdat_exp.groupby("data set"):
            ct += 1
            sns.distplot(
                dataset_aggdat["distance to root sequence"],
                bins=bins,
                kde=False,
                color="black",
                hist_kws={
                    "histtype": "step",
                    "cumulative": True,
                    "alpha": 0.5,
                    "lw": 3,
                },
            )
    plt.xlabel("distance to root sequence")
    plt.xlim([0, bins[-1]])
    plt.ylabel("observed sequences")
    plt.tight_layout()

    bins = range(
        max(
            aggdat["nearest neighbor distance"].max(),
            aggdat_exp["nearest neighbor distance"].max()
            if args.experimental is not None
            else 0,
        )
        + 2
    )

    plt.subplot(1, 2, 2)
    ct = 0
    for dataset, dataset_aggdat in aggdat.groupby("data set"):
        ct += 1
        sns.distplot(
            dataset_aggdat["nearest neighbor distance"],
            bins=bins,
            kde=False,
            color="gray",
            hist_kws={"histtype": "step", "cumulative": True, "alpha": alpha, "lw": 1},
        )
    if args.experimental is not None:
        for dataset, dataset_aggdat in aggdat_exp.groupby("data set"):
            ct += 1
            sns.distplot(
                dataset_aggdat["nearest neighbor distance"],
                bins=bins,
                kde=False,
                color="black",
                hist_kws={
                    "histtype": "step",
                    "cumulative": True,
                    "alpha": 0.5,
                    "lw": 3,
                },
            )
    plt.xlabel("nearest neighbor distance")
    plt.xlim([0, bins[-1]])
    plt.ylabel("")
    plt.tight_layout()

    plt.savefig(args.outbase + ".pdf")


if __name__ == "__main__":
    main()
