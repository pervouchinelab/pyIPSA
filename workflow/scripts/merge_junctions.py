import argparse
import gzip
from typing import List, Dict


def parse_cli_args() -> Dict:
    """Parses command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge junctions by reference, start, stop and strand.",
        usage="%(prog)s [-h] file1 [file2 ...] -o FILE"
    )
    parser.add_argument(
        "files", type=str, nargs="+", metavar="FILES",
        help="input junction files (step >= 2)"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="output file with merged junctions"
    )
    args = parser.parse_args()
    return vars(args)


def get_unique_junctions(files: List[str]) -> List[str]:
    """
    Extract splice junctions from files and return sorted list
    of unique junctions.

    Parameters
    ----------
    files : list of str
        List of input files with junctions

    Returns
    -------
    list
        Sorted list of unique junctions
    """
    unique_junctions = set()
    for file in files:
        with gzip.open(file, "rt") as lines:
            for line in lines:
                junction_id, *_ = line.strip().split()
                unique_junctions.add(junction_id)
    unique_junctions = [junction_id + "\n" for junction_id in unique_junctions]
    unique_junctions.sort()
    return unique_junctions


def main():
    args = parse_cli_args()
    unique_junctions = get_unique_junctions(args["files"])
    with gzip.open(args["output"], "wt") as output:
        output.writelines(unique_junctions)


if __name__ == '__main__':
    main()
