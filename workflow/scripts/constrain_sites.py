import argparse
import gzip
import pandas as pd


def parse_cli_args():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(
        description="Chooses strand in junctions file"
    )
    parser.add_argument(
        "-j", "--junctions", type=str, metavar="FILE", required=True,
        help="input junctions file (.gz)"
    )
    parser.add_argument(
        "-s", "--sites", type=str, metavar="FILE", required=True,
        help="input sites file (.gz)"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="output file name"
    )
    args = parser.parse_args()
    return vars(args)


def read_junctions_to_singletons(filename: str):
    singletons = set()
    with gzip.open(filename, "rt") as j:
        for line in j:
            junction_id, *_ = line.strip().split("\t")
            ref_name, left, right, strand = junction_id.split("_")
            singletons.add((ref_name, left, strand))
            singletons.add((ref_name, right, strand))
    return singletons


def filter_sites(filename: str, singletons):
    rows_list = []
    with gzip.open(filename, "rt") as s:
        for line in s:
            site_id, total_count, staggered_count, entropy = line.strip().split("\t")
            if tuple(site_id.split("_")) in singletons:
                rows_list.append((site_id, total_count, staggered_count, entropy))
    df = pd.DataFrame(rows_list, columns=["site_id", "total_count", "staggered_count", "entropy"])
    df.sort_values(by="site_id")
    return df


def main():
    args = parse_cli_args()
    singletons = read_junctions_to_singletons(args["junctions"])
    df = filter_sites(args["sites"], singletons)
    df.to_csv(args["output"], sep="\t", index=False, header=None, compression="gzip")


if __name__ == '__main__':
    main()
