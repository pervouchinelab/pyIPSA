import argparse
import gzip
from collections import defaultdict
import pandas as pd


def parse_cli_args():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(
        description="Chooses strand in junctions file"
    )
    parser.add_argument(
        "-j", "--junctions", type=str, metavar="FILE", required=True,
        help="input junctions file (.gz)", dest="input"
    )
    parser.add_argument(
        "-r", "--ranked", type=str, metavar="FILE", required=True,
        help="input ranked list of splice sites"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="output file name"
    )
    args = parser.parse_args()
    return vars(args)


def read_junctions(filename: str):
    dd = defaultdict(list)
    with gzip.open(filename, "rt") as j:
        for line in j:
            junction_id, total_counts, staggered_count, \
                entropy, ann, seq = line.strip().split("\t")
            ref_name, left, right, strand = junction_id.split("_")

            dd[(ref_name, left, right)].append((
                strand, total_counts, staggered_count, entropy, ann, seq
            ))
    return dd


def choose(dd, filename):
    ranked = dict()
    with open(filename, "r") as f:
        for line in f:
            seq, rank = line.strip().split("\t")
            ranked[seq] = int(rank)

    rows_list = []
    for junction, candidates in dd.items():
        candidate = max(candidates, key=lambda x: (int(x[4]), ranked.get(x[5], 0)))
        strand, total_counts, staggered_count, entropy, ann, seq = candidate
        ref_name, left, right = junction
        rows_list.append((
            "_".join([ref_name, left, right, strand]),
            total_counts, staggered_count, entropy, ann, seq
        ))
    df = pd.DataFrame(rows_list,
                      columns=["junction_id", "total_counts", "staggered_count",
                               "entropy", "ann", "seq"])
    df.sort_values(by="junction_id")
    return df


def main():
    args = parse_cli_args()
    dd = read_junctions(args["input"])
    df = choose(dd, args["ranked"])
    df.to_csv(args["output"], sep="\t", index=False, header=None, compression="gzip")


if __name__ == '__main__':
    main()
