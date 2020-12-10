import argparse
import gzip
from collections import defaultdict, namedtuple
from typing import Dict, DefaultDict, List, Tuple

import pandas as pd

Row = namedtuple("Row", ["sj_id", "inclusion", "exclusion", "retention"])


def parse_cli_args() -> Dict:
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(
        description="Computes inclusion and processing rates of splice junctions"
    )
    parser.add_argument(
        "-j", "--junctions", type=str, metavar="FILE", required=True,
        help="input junctions file (step >= 4)"
    )
    parser.add_argument(
        "-s", "--sites", type=str, metavar="FILE", required=True,
        help="input sites file (step >= 4)"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="output file name (step 7)"
    )
    args = parser.parse_args()
    return vars(args)


def read_junctions(filename: str) -> Tuple[Dict, List]:
    sites = defaultdict(lambda: {"+": set(), "-": set()})
    with gzip.open(filename, "rt") as j:
        for line in j:
            junction_id, *_ = line.strip().split("\t")
            ref_name, start, stop, strand = junction_id.split("_")
            start, stop = int(start), int(stop)
            sites[ref_name][strand].add(start)
            sites[ref_name][strand].add(stop)

    index = {}
    back = []
    i = 0
    for ref_name in sorted(sites.keys()):
        for strand in sorted(sites[ref_name].keys()):
            for pos in sorted(sites[ref_name][strand]):
                index[(ref_name, strand, pos)] = i
                back.append((ref_name, strand, pos))
                i += 1
    return index, back


def compute_rates(junctions_filename: str,
                  sites_filename: str,
                  index: dict,
                  back: list) -> Tuple[DefaultDict, DefaultDict, DefaultDict]:
    inclusion = defaultdict(int)
    exclusion = defaultdict(int)
    retention = defaultdict(int)
    with gzip.open(junctions_filename, "rt") as j:
        for line in j:
            junction_id, total_count, *_ = line.strip().split("\t")
            ref_name, start, stop, strand = junction_id.split("_")
            start, stop, total_count = int(start), int(stop), int(total_count)

            left = index.get((ref_name, strand, start))
            right = index.get((ref_name, strand, stop))
            if left is None or right is None or left >= right:
                continue
            inclusion[(ref_name, strand, start, "D" if strand == "+" else "A")] += total_count
            inclusion[(ref_name, strand, stop, "A" if strand == "+" else "D")] += total_count

            for pos in range(left + 1, right):
                (ref_name, strand, pos) = back[pos]
                exclusion[(ref_name, strand, pos)] += total_count

    with gzip.open(sites_filename, "rt") as s:
        for line in s:
            site_id, total_count, *_ = line.strip().split("\t")
            ref_name, pos, strand = site_id.split("_")
            pos, total_count = int(pos), int(total_count)
            retention[(ref_name, strand, pos)] += total_count

    return inclusion, exclusion, retention


def output(inclusion: defaultdict,
           exclusion: defaultdict,
           retention: defaultdict,
           back: list) -> pd.DataFrame:
    rows_list = []
    for ref_name, strand, pos in back:
        key = (ref_name, strand, pos, "A")
        if key in inclusion:
            rows_list.append(Row(
                sj_id="_".join([ref_name, str(pos), strand, key[3]]),
                inclusion=inclusion[key],
                exclusion=exclusion[(ref_name, strand, pos)],
                retention=retention[(ref_name, strand, pos)]
            ))
        key = (ref_name, strand, pos, "D")
        if key in inclusion:
            rows_list.append(Row(
                sj_id="_".join([ref_name, str(pos), strand, key[3]]),
                inclusion=inclusion[key],
                exclusion=exclusion[(ref_name, strand, pos)],
                retention=retention[(ref_name, strand, pos)]
            ))
    df = pd.DataFrame(rows_list)
    df.sort_values(by=["sj_id"])
    return df


def main():
    args = parse_cli_args()
    index, back = read_junctions(filename=args["junctions"])
    inclusion, exclusion, retention = compute_rates(
        junctions_filename=args["junctions"],
        sites_filename=args["sites"],
        index=index, back=back
    )
    output(inclusion, exclusion, retention, back).to_csv(
        args["output"], sep="\t", index=False,
        header=False, compression="gzip"
    )


if __name__ == '__main__':
    main()
