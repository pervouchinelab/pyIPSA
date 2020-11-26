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
    ["strand", "total_count", "staggered_count", "entropy",
     "annotation_status", "sequence"]
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
    args = parser.parse_args()
    return vars(args)


def read_junctions(filename: str) -> DefaultDict[Tuple, List[JP]]:
    """Read annotated junctions file from step 3 and return defaultdict."""
    dd = defaultdict(list)

    with gzip.open(filename, "rt") as junctions:

        for line in junctions:
            junction_id, total_count, staggered_count, \
                entropy, annotation_status, sequence = line.strip().split("\t")

            ref_name, left, right, strand = junction_id.split("_")

            dd[(ref_name, left, right)].append((
                JP(strand=strand, total_count=total_count,
                   staggered_count=staggered_count, entropy=entropy,
                   annotation_status=annotation_status, sequence=sequence)
            ))

    return dd


def choose_strand(dd: DefaultDict[Tuple, List[JP]], filename: str) -> pd.DataFrame:
    """Choose strand for each junction using annotation status of each strand and
    ranked list of splice site sequences and return dataframe."""
    ranked = dict()  # prepare ranked sequences
    with open(filename, "r") as f:
        for line in f:
            seq, rank = line.strip().split("\t")
            ranked[seq] = int(rank)

    rows_list = []  # store rows of dataframe

    for junction, candidates in dd.items():
        chosen = max(
            candidates,
            key=lambda x: (int(x.annotation_status), ranked.get(x.sequence, 0))
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


def main():
    args = parse_cli_args()
    dd = read_junctions(args["input"])
    df = choose_strand(dd=dd, filename=args["ranked"])
    df.sort_values(by=["junction_id"], inplace=True)
    df.to_csv(args["output"], sep="\t", index=False, header=False, compression="gzip")


if __name__ == '__main__':
    main()
