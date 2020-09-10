import pandas as pd
import numpy as np
import argparse
from collections import namedtuple

NewRow = namedtuple("Pandas", ["junction_id", "offset", "count"])


def library_params(library_file):
    with open(library_file, "r") as lf:
        read_length = int(lf.readline().split("=")[1])
        [lf.readline() for _ in range(5)]
        paired = "pair" in lf.readline()
        [lf.readline() for _ in range(11)]
        stranded = "unstranded" not in lf.readline()
        return read_length, paired, stranded


def read_ssj(filename):
    ssj = pd.read_table(filename, header=None)
    ssj.columns = ["junction_id", "offset", "F1", "R1", "F2", "R2"]
    return ssj


def update_ssj(df, read_length, min_offset, stranded):
    max_offset = read_length - min_offset
    rows_list = []
    i = 0
    if stranded:
        for row in df.itertuples():
            if row.offset < min_offset or row.offset > max_offset:
                continue
            sm = row.F1 + row.R2
            sp = row.F2 + row.R1
            if sp > 0:
                rows_list.append(NewRow(junction_id=row.junction_id+"_+", offset=row.offset, count=sp))
                i += 1
            if sm > 0:
                rows_list.append(NewRow(junction_id=row.junction_id+"_-", offset=row.offset, count=sm))
                i += 1
    else:
        for row in df.itertuples():
            if row.offset < min_offset or row.offset > max_offset:
                continue
            s = row.F1 + row.R1 + row.F2 + row.R2
            rows_list.append(NewRow(junction_id=row.junction_id+"_+", offset=row.offset, count=s))
            i += 1
            rows_list.append(NewRow(junction_id=row.junction_id+"_-", offset=row.offset, count=s))
            i += 1
    return pd.DataFrame(rows_list, columns=["junction_id", "offset", "count"])


def aggregate_ssj(ssj):
    agg_ssj = ssj.groupby(["junction_id"]).agg(
        total_count=pd.NamedAgg(column="count", aggfunc=sum),
        staggered_count=pd.NamedAgg(column="count", aggfunc="count"),
        entropy=pd.NamedAgg(column="count", aggfunc=lambda x: np.log2(x.sum()) - (np.sum(x * np.log2(x)) / x.sum()))
    )
    agg_ssj["entropy"] = np.round(agg_ssj["entropy"], 6)
    return agg_ssj


def main():
    parser = argparse.ArgumentParser(
        description="Aggregates counts by the offset with three different functions"
    )
    parser.add_argument("-i", "--input_tsv", type=str,
                        metavar="", required=True, help="Input counts file (TSV)")
    parser.add_argument("-o", "--output_tsv", type=str,
                        metavar="", required=True, help="Output aggregated counts file (TSV")
    parser.add_argument("-l", "--library_log", type=str,
                        metavar="", required=True, help="Library analysis file")
    parser.add_argument("-m", "--min_offset", type=int, default=0, metavar="", help="Minimal offset")
    cli_args = parser.parse_args()
    read_length, paired, stranded = library_params(cli_args.library_log)
    ssj = read_ssj(cli_args.input_tsv)
    ssj = update_ssj(ssj, read_length=read_length, min_offset=cli_args.min_offset, stranded=stranded)
    agg_ssj = aggregate_ssj(ssj)
    agg_ssj.to_csv(cli_args.output_tsv, sep="\t", header=None, compression="gzip")


if __name__ == '__main__':
    main()
