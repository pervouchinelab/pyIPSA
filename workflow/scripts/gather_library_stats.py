import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd


def parse_cli_args() -> Dict:
    """Parses command line arguments."""
    parser = argparse.ArgumentParser(
        description="Gather library stats from given replicates (pivot table).",
        usage="%(prog)s [-h] replicate1 [replicate2 ...] -o FILE"
    )
    parser.add_argument(
        "stats", type=str, nargs="+", metavar="FILES",
        help="input stats files (step = 1)"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="output file (table describing replicates)"
    )
    args = parser.parse_args()
    return vars(args)


def stats2df(replicates: List[str]) -> pd.DataFrame:
    """Read replicates' stats and gather all the info into one table."""
    records = []

    for replicate in replicates:
        p = Path(replicate)
        name = Path(p.stem).stem

        with p.open("r") as f:
            for line in f:
                if line.startswith("-"):
                    break
                left, right = line.strip().split(": ")
                if left[0] == "r":
                    read_length = right
                elif left[0] == "g":
                    genome = right
                elif left[0] == "p":
                    paired = right
                elif left[0] == "s":
                    stranded = right
                elif left[0] == "l":
                    library_type = right

        records.append((name, read_length, genome, paired, stranded, library_type))

    df = pd.DataFrame(records)
    df.columns = ["replicate", "read length", "genome", "paired", "stranded", "library type"]
    df.sort_values(by=["replicate"])
    return df


def main():
    args = parse_cli_args()
    df = stats2df(args["stats"])
    df.to_csv(args["output"], index=False, sep="\t")


if __name__ == '__main__':
    main()
