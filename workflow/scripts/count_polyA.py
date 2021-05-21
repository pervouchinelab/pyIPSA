import argparse
from collections import defaultdict, namedtuple
from typing import Dict

import pandas as pd
import pysam

from .ipsa_config import *

MIN_TAIL = 8
MIN_MATCH = 10
MIN_PERCENT = 0.8

RawPolyASite = namedtuple("RawPolyA", ["ref_name", "start", "strand", "overhang"])
PolyASite = namedtuple("PolyASite", ["site_id", "overhang", "F1", "R1", "F2", "R2"])


def parse_cli_args() -> Dict:
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Extract polyA sites from BAM"
    )
    parser.add_argument(
        "-i", "--input", type=str, metavar="FILE", required=True,
        help="input alignment file (BAM)"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="output file name"
    )
    parser.add_argument(
        "-u", "--unique", action="store_true",
        help="only use uniquely mapped reads"
    )
    parser.add_argument(
        "-p", "--primary", action="store_true",
        help="only use primary alignments of multimapped reads"
    )
    parser.add_argument(
        "-t", "--threads", type=int, metavar="INT", default=1,
        help="number of threads to use"
    )
    args = parser.parse_args()
    return vars(args)


def segment2polya(segment: pysam.AlignedSegment, polya_with_counts: defaultdict):
    """
    Take a read (AlignedSegment) and add new counts to junctions.

    :param segment: aligned segment (read)
    :param polya_with_counts: dictionary storing junctions and their counts
    """

    # name of reference sequence where read is aligned
    ref_name = segment.reference_name

    # segment (read) type (0 - read1+, 1 - read1-, 2 - read2+, 3 - read2-)
    seg_type = 2 * segment.is_read2 + segment.is_reverse

    # initial positions in segment and reference
    seg_pos = 0
    ref_pos = segment.reference_start + BASE  # from 0-based to 1-based

    if segment.cigartuples[0][0] == 4:
        soft_clip_length = segment.cigartuples[0][1]
        next_block_length = segment.cigartuples[1][1]

        if soft_clip_length >= MIN_TAIL and next_block_length >= MIN_MATCH:
            soft_clip_seq = segment.query_sequence[:soft_clip_length]
            if soft_clip_seq.count("T") >= MIN_PERCENT * soft_clip_length:
                raw_polya = RawPolyASite(
                    ref_name=ref_name,
                    start=ref_pos,
                    strand="-",
                    overhang=soft_clip_length
                )
                polya_with_counts[raw_polya][seg_type] += 1

    if segment.cigartuples[-1][0] == 4:
        soft_clip_length = segment.cigartuples[-1][1]
        prev_block_length = segment.cigartuples[-2][1]

        if soft_clip_length >= MIN_TAIL and prev_block_length >= MIN_MATCH:
            soft_clip_seq = segment.query_sequence[-soft_clip_length:]
            if soft_clip_seq.count("A") >= MIN_PERCENT * soft_clip_length:
                raw_polya = RawPolyASite(
                    ref_name=ref_name,
                    start=segment.reference_end,
                    strand="+",
                    overhang=soft_clip_length
                )
                polya_with_counts[raw_polya][seg_type] += 1

    return


def alignment2polya(
        alignment: pysam.AlignmentFile,
        unique: bool,
        primary: bool
):
    """
    Find all junctions in alignment file and compute their counts.
    Also determine most common read length.

    :param alignment: alignment file
    :param unique: account only uniquely mapped (aligned) segments
    :param primary: account only primary alignment if segment is a multi-mapper
    :return: dictionary mapping junction to its counts in each offset and most common read length
    """
    segment: pysam.AlignedSegment  # just annotation line for convenience
    polya_with_counts = defaultdict(lambda: [0] * 4)

    # iterating through the alignment file
    for segment in alignment.fetch():

        # if read is not mapped
        if segment.is_unmapped:
            continue

        # if read is a multi-mapper
        if unique and segment.has_tag("NH") and segment.get_tag("NH") > 1:
            continue

        # if read is a multi-mapper consider only primary alignment
        if primary and segment.is_secondary:
            continue

        # if read is a multi-mapper consider only primary alignment
        if primary and segment.is_supplementary:
            continue

        segment2polya(segment=segment, polya_with_counts=polya_with_counts)

    return polya_with_counts


def polya2dataframe(polya_with_counts) -> pd.DataFrame:
    """
    Collect all junctions into sorted DataFrame.

    :param polya_with_counts: dictionary mapping junction to its counts in each offset
    :return: DataFrame with junctions and their counts in each offset
    """
    rows_list = []

    for raw_polya, counts in polya_with_counts.items():
        ref_name, start, strand, overhang = raw_polya

        # filtration of reference names
        if ref_name not in ALLOWED_REFERENCE_NAMES:
            lowered = ref_name.lower()
            if lowered in TRANSLATION:
                ref_name = TRANSLATION[lowered]
            else:
                continue

        polya = PolyASite(
            site_id="_".join([ref_name, str(start), strand]),
            overhang=overhang,
            F1=counts[0], R1=counts[1], F2=counts[2], R2=counts[3]
        )
        rows_list.append(polya)

    df = pd.DataFrame(rows_list, columns=["site_id", "overhang", "F1", "R1", "F2", "R2"])
    df = df.sort_values(by=["site_id", "overhang"])

    return df


def main():
    args = parse_cli_args()
    # read BAM and get junctions
    bam = pysam.AlignmentFile(args["input"], threads=args["threads"])
    polya = alignment2polya(
        alignment=bam, primary=args["primary"], unique=args["unique"]
    )
    df_polya = polya2dataframe(polya_with_counts=polya)
    # write outputs
    df_polya.to_csv(args["output"], sep="\t", index=False, header=False)


if __name__ == "__main__":
    main()
