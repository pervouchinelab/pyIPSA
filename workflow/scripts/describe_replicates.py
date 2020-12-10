import argparse
from pathlib import Path
from typing import Dict, List

import pandas as pd


def parse_cli_args() -> Dict:
    """Parses command line arguments."""
    parser = argparse.ArgumentParser(
        description="Describe replicates (pivot table).",
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
        name = p.stem

        with p.open("r") as f:
            for line in f:
                if " is " not in line:
                    continue
                left, right = line.strip().split(" is ")
                if "genome" in left:
                    genome = right
                if "Read" in left:
                    read_length = int(right)
                if right.endswith("-end"):
                    paired = True if right[:-4] == "pair" else False
                if right.endswith("stranded"):
                    stranded = True if len(right) == 8 else False

        records.append((name, genome, read_length, paired, stranded))

    df = pd.DataFrame(records)
    df.columns = ["replicate", "genome", "read_length", "paired", "stranded"]
    df.sort_values(by=["replicate"])
    return df


def main():
    args = parse_cli_args()
    df = stats2df(args["stats"])
    df.to_csv(args["output"], index=False, sep="\t")


if __name__ == '__main__':
    main()
