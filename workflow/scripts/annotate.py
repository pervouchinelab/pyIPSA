import argparse
import gzip
import pandas as pd
import pysam


def parse_cli_args():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(description="Annotate splice site junctions")
    parser.add_argument(
        "-j", "--junctions", type=str, metavar="FILE",
        required=True, help="input file with junctions (.gz)"
    )
    parser.add_argument(
        "-k", "--known_sj", type=str, metavar="FILE", required=True,
        help="input file with introns (.gz)")
    parser.add_argument(
        "-f", "--fasta", type=str, metavar="FILE", required=True,
        help="genome file (FASTA)"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="output file name"
    )
    args = parser.parse_args()
    return vars(args)


def read_sj(filename):
    df = pd.read_table(filename, header=None)
    df.columns = ["sj_id", "total_count", "staggered_count", "entropy"]
    return df


def annotate_splice_sites(df, filename):
    pairs, singletons = set(), set()
    with gzip.open(filename, "rt") as f:
        for line in f:
            ref, left, right, strand = line.strip().split("\t")
            left, right = int(left), int(right)
            pairs.add((ref, left, right, strand))
            singletons.add((ref, left, strand))
            singletons.add((ref, right, strand))

    anno_status = []
    for junction_id in df["sj_id"]:
        ref, start, stop, strand = junction_id.split("_")
        start, stop = int(start), int(stop)
        status = 0
        if (ref, start, stop, strand) in pairs:
            status += 1

        if (ref, start, strand) in singletons:
            status += 1

        if (ref, stop, strand) in singletons:
            status += 1

        anno_status.append(status)
    df["annotation_status"] = anno_status
    return df


def add_splice_sites(df, genome):
    genome = pysam.FastaFile(genome)
    complement = {"A": "T", "G": "C", "C": "G", "T": "A", "N": "N"}
    splice_sites = []
    for junction_id in df["sj_id"]:
        ref, start, stop, strand = junction_id.split("_")
        start, stop = int(start), int(stop)
        c = strand == "-"
        n = genome.fetch(ref, start, start + 2).upper() + genome.fetch(ref, stop - 3, stop - 1).upper()
        if c:
            n = "".join(complement[n] for n in n)[::-1]
        splice_sites.append(n)
    df["splice_sites"] = splice_sites
    return df


def main():
    args = parse_cli_args()
    agg_j = read_sj(args["junctions"])
    ann_j = annotate_splice_sites(agg_j, args["known_sj"])
    ann_j = add_splice_sites(ann_j, args["fasta"])
    ann_j.to_csv(args["output"], sep="\t", index=False, header=None, compression="gzip")


if __name__ == '__main__':
    main()
