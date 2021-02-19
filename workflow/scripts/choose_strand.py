"""
This module has functions to choose strand for each splice junction.
"""

import argparse
import gzip

from typing import Dict, DefaultDict, List, Tuple
from collections import defaultdict, namedtuple

import pandas as pd


JP = namedtuple(
    "JunctionParams",
    ["strand", "total_count", "staggered_count", "entropy", "annotation_status", "sequence"]
)


def parse_cli_args() -> Dict:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Choose strand in junctions file"
    )
    parser.add_argument(
        "-i", "--input", type=str, metavar="FILE",
        required=True, help="input annotated junctions file (step 3)"
    )
    parser.add_argument(
        "-r", "--ranked", type=str, metavar="FILE",
        required=True, help="ranked splice site sequences file"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE",
        required=True, help="output file name"
    )
    parser.add_argument(
        "-s", "--stats", type=str, metavar="FILE",
        required=True, help="junction stats file"
    )
    args = parser.parse_args()
    return vars(args)


def read_junctions(filename: str) -> DefaultDict[Tuple, List[JP]]:
    """
    Load annotated junctions from step 3 to dictionary.
    """
    dd = defaultdict(list)

    with gzip.open(filename, "rt") as junctions:

        for line in junctions:
            junction_id, total_count, staggered_count, \
                entropy, annotation_status, sequence = line.strip().split("\t")

            ref_name, left, right, strand = junction_id.split("_")

            dd[(ref_name, left, right)].append((
                JP(strand=strand, total_count=int(total_count),
                   staggered_count=int(staggered_count), entropy=entropy,
                   annotation_status=int(annotation_status), sequence=sequence)
            ))

    return dd


def choose_strand(dd: DefaultDict[Tuple, List[JP]], filename: str) -> pd.DataFrame:
    """
    Choose strand for each junction using annotation status of each strand and
    ranked list of splice site sequences and return dataframe.
    """
    ranked = dict()  # prepare ranked sequences
    with open(filename, "r") as f:
        for line in f:
            seq, rank = line.strip().split("\t")
            ranked[seq] = int(rank)

    rows_list = []  # store rows of dataframe

    for junction, candidates in dd.items():
        chosen = max(
            candidates,
            key=lambda x: (x.annotation_status, ranked.get(x.sequence, 0))
        )

        ref_name, left, right = junction
        rows_list.append((
            "_".join([ref_name, left, right, chosen.strand]),
            chosen.total_count, chosen.staggered_count,
            chosen.entropy, chosen.annotation_status, chosen.sequence
        ))

    df = pd.DataFrame(rows_list,
                      columns=["junction_id", "total_count", "staggered_count",
                               "entropy", "ann", "seq"])
    return df


def compute_junction_stats(df: pd.DataFrame) -> Dict[str, int]:
    """
    Compute how many junctions correspond to GTAG and non-GTAG splice sites
    and how many junctions have at least one annotated end.

    :param df: DataFrame with junctions
    :return: dictionary with computed stats
    """
    stats = dict()

    gtag_rows = df["seq"] == "GTAG"
    ctac_rows = df["seq"] == "CTAC"
    stats["GTAG"] = df["total_count"][gtag_rows].sum()
    stats["CTAC"] = df["total_count"][ctac_rows].sum()
    stats["Non-GTAG"] = df["total_count"][~gtag_rows].sum()
    stats["GTAG / non-GTAG ratio"] = round(stats["GTAG"] / stats["Non-GTAG"], 2)

    annotated_rows = df["ann"] >= 1
    stats["Annotated"] = df["total_count"][annotated_rows].sum()
    stats["Not annotated"] = df["total_count"][~annotated_rows].sum()
    stats["Annotated / not annotated ratio"] = round(stats["Annotated"] / stats["Not annotated"], 2)

    return stats


def main():
    args = parse_cli_args()
    dd = read_junctions(args["input"])
    df = choose_strand(dd=dd, filename=args["ranked"])
    df.sort_values(by=["junction_id"], inplace=True)
    junction_stats = compute_junction_stats(df=df)
    with open(args["stats"], "w") as f:
        for key, value in junction_stats.items():
            f.write(f"{key}: {value}\n")
    df.to_csv(args["output"], sep="\t", index=False, header=False, compression="gzip")


if __name__ == '__main__':
    main()
