"""
This module contains functions that:
- aggregate read counts
- add strand info to junctions
- filter out too small or too large introns
"""
import argparse
from collections import namedtuple

from typing import Dict, Tuple

import numpy as np
import pandas as pd

NewRow = namedtuple("Pandas", ["sj_id", "offset", "count"])


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
        "--strand_mode", type=str, metavar="STR",  choices=["F2R1", "F1R2"],
        help="assign F2R1 to forward strand F1R2 (only for stranded data)"
    )
    parser.add_argument(
        "--min_offset", type=int, default=1,
        metavar="INT", help="minimal offset"
    )
    parser.add_argument(
        "--min_intron_length", type=int, default=40,
        metavar="INT", help="minimal intron length"
    )
    parser.add_argument(
        "--max_intron_length", type=int, default=100000,
        metavar="INT", help="maximal intron length"
    )
    args = parser.parse_args()
    return vars(args)


def read_stats(filename: str) -> Tuple[str, int, bool, bool]:
    """Read sample stats file and return read length,
    genome and its version, and if reads are paired or stranded."""
    with open(filename, "r") as f:
        for line in f:
            if " is " not in line:
                continue
            left, right = line.strip().split(" is ")
            if "genome" in left:
                genome = right
            elif "Read" in left:
                read_length = int(right)
            elif right.endswith("-end"):
                paired = True if right[:-4] == "pair" else False
            elif right.endswith("stranded"):
                stranded = True if len(right) == 8 else False
    return genome, read_length, paired, stranded


def read_counts(filename: str) -> pd.DataFrame:
    """Read junctions with counts file from step 1 and return dataframe."""
    df = pd.read_table(filename, header=None)
    df.columns = ["sj_id", "offset", "F1", "R1", "F2", "R2"]
    return df


def add_strand(df: pd.DataFrame, stranded: bool, strand_mode: str) -> pd.DataFrame:
    """Update dataframe with strand information (only for junctions)."""
    rows_list = []

    for row in df.itertuples():  # iterate all junctions

        if stranded:
            forward_count = row.F2 + row.R1
            reverse_count = row.F1 + row.R2

            if strand_mode == "F1R2":
                forward_count, reverse_count = reverse_count, forward_count

        else:
            forward_count = row.F1 + row.R1 + row.F2 + row.R2
            reverse_count = row.F1 + row.R1 + row.F2 + row.R2

        if forward_count > 0:
            rows_list.append(NewRow(
                sj_id=row.sj_id+"_+", offset=row.offset, count=forward_count
            ))

        if reverse_count > 0:
            rows_list.append(NewRow(
                sj_id=row.sj_id+"_-", offset=row.offset, count=reverse_count
            ))

    return pd.DataFrame(rows_list, columns=["sj_id", "offset", "count"])


def compute_entropy(x: pd.Series) -> pd.Series:
    """Compute entropy."""
    s = x.sum()
    result = np.log2(s) - (np.sum(x * np.log2(x)) / s)
    return np.round(result, 2)


def aggregate(df: pd.DataFrame, min_offset: int, read_length: int) -> pd.DataFrame:
    """Aggregate counts by sites/junction id."""
    max_offset = read_length - min_offset
    df = df[(df["offset"] >= min_offset) & (df["offset"] <= max_offset)]
    df = df.groupby(["sj_id"]).agg(
        total_count=pd.NamedAgg(column="count", aggfunc=sum),
        staggered_count=pd.NamedAgg(column="count", aggfunc="count"),
        entropy=pd.NamedAgg(
            column="count",
            aggfunc=compute_entropy
        )
    )
    return df


def row_filter(row: pd.Series,
               min_intron_length: int,
               max_intron_length: int) -> bool:
    """Row filter to apply on dataframe."""
    ref_name, left, right, strand = row["sj_id"].split("_")
    length = int(right) - int(left) - 1
    if min_intron_length <= length <= max_intron_length:
        return True
    else:
        return False


def filter_introns(df: pd.DataFrame,
                   min_intron_length: int,
                   max_intron_length: int) -> pd.DataFrame:
    """Filter out too small or too large introns."""
    df = df.reset_index()
    keep = df.apply(func=row_filter, axis=1, args=(min_intron_length, max_intron_length))
    return df[keep]


def main():
    args = parse_cli_args()
    counts = read_counts(filename=args["input"])
    genome, read_length, paired, stranded = read_stats(args["stats"])
    counts = add_strand(df=counts, stranded=stranded, strand_mode=args["strand_mode"])
    counts = aggregate(df=counts, min_offset=args["min_offset"], read_length=read_length)
    counts = filter_introns(df=counts,
                            min_intron_length=args["min_intron_length"],
                            max_intron_length=args["max_intron_length"])
    counts.sort_values(by=["sj_id"], inplace=True)
    counts.to_csv(args["output"], sep="\t", index=False, header=False, compression="gzip")


if __name__ == '__main__':
    main()
