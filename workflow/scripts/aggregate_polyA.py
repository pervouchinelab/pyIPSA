import argparse
from collections import namedtuple

from typing import Dict, Tuple

import numpy as np
import pandas as pd

NewRow = namedtuple("Pandas", ["polyA_id", "offset", "count"])


def parse_cli_args() -> Dict:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Aggregate polyA sites"
    )
    parser.add_argument(
        "-i", "--input", type=str, metavar="FILE",
        required=True, help="input counts file (step 1)"
    )
    parser.add_argument(
        "-s", "--stats", type=str, metavar="FILE",
        required=True, help="input stats file (from step 1)"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE",
        required=True, help="output file name (step 2)"
    )
    parser.add_argument(
        "--min_overhang", type=int, default=1,
        metavar="INT", help="minimal overhang"
    )
    args = parser.parse_args()
    return vars(args)


def read_stats(filename: str) -> Tuple[int, str, bool, bool, str]:
    """Read sample stats file and return read length,
    genome and its version, and if reads are paired or stranded."""
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("-"):
                break
            left, right = line.strip().split(": ")
            if left[0] == "r":
                read_length = int(right)
            elif left[0] == "g":
                genome = right
            elif left[0] == "p":
                paired = (right == "True")
            elif left[0] == "s":
                stranded = (right == "True")
            elif left[0] == "l":
                library_type = right

    return read_length, genome, paired, stranded, library_type


def read_counts(filename: str) -> pd.DataFrame:
    """Read junctions with counts file from step 1 and return dataframe."""
    df = pd.read_table(filename, header=None)
    df.columns = ["polyA_id", "overhang", "F1", "R1", "F2", "R2"]
    return df


def filter_strand(df: pd.DataFrame, stranded: bool, strand_mode: str) -> pd.DataFrame:
    """Update dataframe with strand information (only for junctions)."""
    rows_list = []

    for row in df.itertuples():  # iterate all junctions

        strand = row.polyA_id.split('_')[2]
        forward_count = row.F2 + row.R1
        reverse_count = row.F1 + row.R2

        if stranded:

            if strand_mode == "F1R2":
                reverse_count, forward_count = forward_count, reverse_count

        else:
            forward_count = row.F1 + row.R1 + row.F2 + row.R2
            reverse_count = row.F1 + row.R1 + row.F2 + row.R2

        if strand == "+" and forward_count > 0:
            rows_list.append(NewRow(
                polyA_id=row.polyA_id, offset=row.overhang, count=forward_count
            ))

        if strand == "-" and reverse_count > 0:
            rows_list.append(NewRow(
                polyA_id=row.polyA_id, offset=row.overhang, count=reverse_count
            ))

    return pd.DataFrame(rows_list, columns=["polyA_id", "overhang", "count"])


def compute_entropy(x: pd.Series) -> pd.Series:
    """Compute entropy."""
    s = x.sum()
    result = np.log2(s) - (np.sum(x * np.log2(x)) / s)
    return np.round(result, 2)


def aggregate(df: pd.DataFrame, min_overhang: int) -> pd.DataFrame:
    """Aggregate counts by sites/junction id."""
    df = df[(df["overhang"] >= min_overhang)]
    df = df.groupby(["polyA_id"], as_index=False).agg(
        total_count=pd.NamedAgg(column="count", aggfunc=sum),
        staggered_count=pd.NamedAgg(column="count", aggfunc="count"),
        entropy=pd.NamedAgg(
            column="count",
            aggfunc=compute_entropy
        )
    )
    return df


def main():
    args = parse_cli_args()
    counts = read_counts(filename=args["input"])
    read_length, genome, paired, stranded, library_type = read_stats(args["stats"])
    counts = filter_strand(df=counts, stranded=stranded, strand_mode=library_type)
    counts = aggregate(df=counts, min_overhang=args["min_overhang"])
    counts.sort_values(by=["polyA_id"], inplace=True)
    counts.to_csv(args["output"], sep="\t", index=False, header=False, compression="gzip")


if __name__ == '__main__':
    main()
