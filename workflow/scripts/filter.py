import argparse

import pandas as pd


def parse_cli_args():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(
        description="Filters splice junctions and sites by different criteria"
    )
    parser.add_argument(
        "-i", "--input", type=str, metavar="FILE", required=True,
        help="input junctions file (step >= 3)"
    )
    parser.add_argument(
        "-e", "--entropy", type=float, metavar="FLOAT", default=1.5,
        help="entropy threshold"
    )
    parser.add_argument(
        "-c", "--count", type=int, metavar="INT", default=1,
        help="total count threshold"
    )
    parser.add_argument(
        "-g", "--gtag", action="store_true",
        help="use only GTAG splice sites"
    )
    args = parser.parse_args()
    return vars(args)


def read_and_filter(filename: str, entropy: float,
                    count: int, gtag: bool) -> pd.DataFrame:
    """Reads junctions, filters them and outputs sorted DataFrame."""
    df = pd.read_table(filename, header=None)
    df.columns = ["sj_id", "total_count", "staggered_count", "entropy",
                  "status", "nucleotides"]
    # TODO: add filtration
    return df


def main():
    args = parse_cli_args()
    print(args)


if __name__ == '__main__':
    main()