"""
This module contains functions that:
- aggregate read counts
"""
import argparse
from typing import Dict, Tuple

import numpy as np
import pandas as pd


def parse_cli_args() -> Dict:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Aggregate splice site junctions counts"
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
        "-m", "--min_offset", type=int, default=1,
        metavar="INT", help="minimal offset"
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
                paired = bool(right)
            elif left[0] == "s":
                stranded = bool(right)
            elif left[0] == "l":
                library_type = right

    return read_length, genome, paired, stranded, library_type


def read_counts(filename: str) -> pd.DataFrame:
    """Read sites with counts file from step 1 and return dataframe."""
    df = pd.read_table(filename, header=None)
    df.columns = ["site_id", "offset", "count"]
    return df


def compute_entropy(x: pd.Series) -> pd.Series:
    """Compute entropy."""
    s = x.sum()
    result = np.log2(s) - (np.sum(x * np.log2(x)) / s)
    return np.round(result, 2)


def aggregate(df: pd.DataFrame, min_offset: int, read_length: int) -> pd.DataFrame:
    """Aggregate counts by sites/junction id."""
    max_offset = read_length - min_offset
    df = df[(df["offset"] >= min_offset) & (df["offset"] <= max_offset)]
    df = df.groupby(["site_id"], as_index=False).agg(
        total_count=pd.NamedAgg(column="count", aggfunc=sum),
        staggered_count=pd.NamedAgg(column="count", aggfunc="count"),
        entropy=pd.NamedAgg(column="count", aggfunc=compute_entropy)
    )
    return df


def main():
    args = parse_cli_args()
    counts = read_counts(filename=args["input"])
    read_length, genome, paired, stranded, library_type = read_stats(args["stats"])
    counts = aggregate(df=counts, min_offset=args["min_offset"], read_length=read_length)
    counts.sort_values(by=["site_id"], inplace=True)
    counts.to_csv(args["output"], sep="\t", index=False, header=False, compression="gzip")


if __name__ == '__main__':
    main()
