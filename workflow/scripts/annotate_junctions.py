"""
This module contains functions that:
- compute annotation status of splice sites of junctions
- retrieve splice junctions sequences
"""
import argparse
import gzip

from typing import Dict

import pandas as pd
import pysam


def parse_cli_args() -> Dict:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Annotate splice site junctions"
    )
    parser.add_argument(
        "-i", "--input", type=str, metavar="FILE",
        required=True, help="input aggregated counts file (step 2)"
    )
    parser.add_argument(
        "-k", "--known", type=str, metavar="FILE",
        required=True, help="annotated introns file")
    parser.add_argument(
        "-f", "--fasta", type=str, metavar="FILE",
        required=True, help="genome sequence file (FASTA)"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE",
        required=True, help="output file name"
    )
    args = parser.parse_args()
    return vars(args)


def read_aggregated_counts(filename: str) -> pd.DataFrame:
    """Read aggregated counts file from step 2 and return dataframe."""
    df = pd.read_table(filename, header=None)
    df.columns = ["sj_id", "total_count", "staggered_count", "entropy"]
    return df


def annotate_splice_sites(df: pd.DataFrame, filename: str) -> pd.DataFrame:
    """Update dataframe with annotation status using annotated introns file."""
    pairs, singletons = set(), set()  # prepare annotated splice sites from introns
    with gzip.open(filename, "rt") as f:
        for line in f:
            ref_name, left, right, strand = line.strip().split("\t")
            left, right = int(left), int(right)
            pairs.add((ref_name, left, right, strand))
            singletons.add((ref_name, left, strand))
            singletons.add((ref_name, right, strand))

    annotation_status = []  # store annotation status column
    for junction_id in df["sj_id"]:

        ref_name, left, right, strand = junction_id.split("_")
        left, right = int(left), int(right)

        status = 0
        if (ref_name, left, right, strand) in pairs:
            status += 1

        if (ref_name, left, strand) in singletons:
            status += 1

        if (ref_name, right, strand) in singletons:
            status += 1

        annotation_status.append(status)

    df["annotation_status"] = annotation_status
    return df


def add_sequences(df: pd.DataFrame, genome: str) -> pd.DataFrame:
    """Update dataframe with splice sites sequences using genome sequence file."""
    genome = pysam.FastaFile(genome)
    complement = {"A": "T", "G": "C", "C": "G", "T": "A", "N": "N"}

    sequences = []  # store sequence column
    for junction_id in df["sj_id"]:
        ref_name, left, right, strand = junction_id.split("_")
        left, right = int(left), int(right)
        c = (strand == "-")  # flag for reverse complement

        seq = genome.fetch(ref_name, left, left + 2).upper() + \
            genome.fetch(ref_name, right - 3, right - 1).upper()

        if c:
            seq = "".join(complement[n] for n in seq)[::-1]

        sequences.append(seq)

    df["sequence"] = sequences
    return df


def main():
    args = parse_cli_args()
    df = read_aggregated_counts(args["input"])
    df = annotate_splice_sites(df=df, filename=args["known"])
    df = add_sequences(df=df, genome=args["fasta"])
    df.sort_values(by=["sj_id"], inplace=True)
    df.to_csv(args["output"], sep="\t", index=False, header=False, compression="gzip")


if __name__ == '__main__':
    main()
