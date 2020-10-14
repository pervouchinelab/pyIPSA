import pandas as pd
import numpy as np
import argparse
from collections import namedtuple

NewRow = namedtuple("Pandas", ["sj_id", "offset", "count"])


def parse_cli_args():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(
        description="Aggregates counts by the offset with three different functions"
    )
    parser.add_argument(
        "-i", "--input_tsv", type=str, metavar="FILE",
        required=True, help="input counts file (.gz)", dest="tsv"
    )
    parser.add_argument(
        "-l", "--log", type=str, metavar="FILE", required=True,
        help="input log file"
    )
    parser.add_argument(
        "-o", "--output_tsv", type=str, metavar="FILE", required=True,
        help="output file name")
    parser.add_argument(
        "-m", "--min_offset", type=int, default=0,
        metavar="", help="minimal offset")
    # TODO: add min and max intron lengths
    args = parser.parse_args()
    return vars(args)


def read_stats(filename: str) -> (str, int, bool, bool):
    """
    This function reads a library analysis file and returns read length
    inferred from the alignment, and two booleans: is the data paired-end
    """
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


def read_sj(filename):
    df = pd.read_table(filename, header=None)
    df.columns = ["sj_id", "offset", "F1", "R1", "F2", "R2"]
    return df


def update_sj(df, read_length, min_offset, stranded):
    max_offset = read_length - min_offset
    rows_list = []
    i = 0
    if stranded:
        for row in df.itertuples():
            if row.offset < min_offset or row.offset > max_offset:
                continue
            # TODO: add strand change mod
            sm = row.F1 + row.R2
            sp = row.F2 + row.R1
            if sp > 0:
                rows_list.append(NewRow(sj_id=row.sj_id+"_+", offset=row.offset, count=sp))
                i += 1
            if sm > 0:
                rows_list.append(NewRow(sj_id=row.sj_id+"_-", offset=row.offset, count=sm))
                i += 1
    else:
        for row in df.itertuples():
            if row.offset < min_offset or row.offset > max_offset:
                continue
            s = row.F1 + row.R1 + row.F2 + row.R2
            rows_list.append(NewRow(sj_id=row.sj_id+"_+", offset=row.offset, count=s))
            i += 1
            rows_list.append(NewRow(sj_id=row.sj_id+"_-", offset=row.offset, count=s))
            i += 1
    return pd.DataFrame(rows_list, columns=["sj_id", "offset", "count"])


def aggregate_sj(df):
    agg_df = df.groupby(["sj_id"]).agg(
        total_count=pd.NamedAgg(column="count", aggfunc=sum),
        staggered_count=pd.NamedAgg(column="count", aggfunc="count"),
        entropy=pd.NamedAgg(column="count", aggfunc=lambda x: np.log2(x.sum()) - (np.sum(x * np.log2(x)) / x.sum()))
    )
    agg_df["entropy"] = np.round(agg_df["entropy"], 2)
    return agg_df


def main():
    args = parse_cli_args()
    genome, read_length, paired, stranded = read_stats(args["log"])
    sj_df = read_sj(args["tsv"])
    ssj = update_sj(sj_df, read_length=read_length, min_offset=args["min_offset"], stranded=stranded)
    agg_ssj = aggregate_sj(ssj)
    agg_ssj.to_csv(args["output_tsv"], sep="\t", header=None, compression="gzip")


if __name__ == '__main__':
    main()
