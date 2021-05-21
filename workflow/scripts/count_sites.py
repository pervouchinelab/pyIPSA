"""
This module contains functions that:
- extract splice sites from BAM
-
"""
import argparse
import gzip
from collections import Counter, defaultdict, namedtuple
from typing import Dict, DefaultDict, Set, Tuple, Counter

import pandas as pd
import pysam

from .ipsa_config import *

SiteWithOffset = namedtuple("SiteWithOffset", ["site_id", "offset"])
SiteWithCount = namedtuple("Site", ["site_id", "offset", "total_count"])


def parse_cli_args():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(
        description="Extracts splice sites from BAM"
    )
    parser.add_argument(
        "-i", "--input", type=str, metavar="FILE", required=True,
        help="input alignment file (BAM)"
    )
    parser.add_argument(
        "-j", "--junctions", type=str, metavar="FILE", required=True,
        help="input junctions file (step >= 4)"
    )
    parser.add_argument(
        "-s", "--stats", type=str, metavar="FILE",
        required=True, help="input stats file (from step 1)"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="output file name (step 1)"
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


def read_stats(filename: str) -> Tuple[int, str, bool, bool, str]:
    """Read sample stats file and return read length,
    genome and its version, and if reads are paired or stranded."""
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("-"):
                break
            left, right = line.strip().split(": ")
            if left[0] == "r":
                read_length = int(right)
            elif left[0] == "g":
                genome = right
            elif left[0] == "p":
                paired = (right == "True")
            elif left[0] == "s":
                stranded = (right == "True")
            elif left[0] == "l":
                library_type = right

    return read_length, genome, paired, stranded, library_type


def junctions2sites(filename: str) -> DefaultDict[Tuple[str, int], Set]:
    """
    Extract sites from junctions found in J pipeline.

    :param filename: input junctions file with strand info (step >= 4)
    :return: dictionary with splice sites (flanking nucleotides)
    """
    sites = defaultdict(set)

    # read junctions file
    with gzip.open(filename, 'rt') as jf:

        for line in jf:

            junction_id = line.strip().split('\t')[0]
            ref_name, left, right, strand = junction_id.split('_')
            sites[(ref_name, int(left))].add(strand)
            sites[(ref_name, int(right))].add(strand)

    return sites


def prepare_ref_names(alignment: pysam.AlignmentFile) -> Dict[str, str]:
    """
    Fix reference names if needed.

    :param alignment: alignment file (pysam object)
    :return: mapping from incorrect names to correct ones
    """
    trans = dict()

    for ref_name in alignment.references:

        if ref_name not in ALLOWED_REFERENCE_NAMES:

            lowered = ref_name.lower()

            if lowered in TRANSLATION:
                trans[ref_name] = TRANSLATION[lowered]
        else:

            trans[ref_name] = ref_name

    return trans


def segment2counts(
        segment: pysam.AlignedSegment,
        ref_name: str,
        sites: DefaultDict[Tuple[str, int], Set],
        counts: Counter[SiteWithOffset],
        stranded: bool,
        strand_mode: str
) -> None:
    """
    Take a read from alignment and add counts for existing sites present in read.

    :param segment: read from alignment
    :param ref_name: corrected reference name
    :param sites: dictionary with sites extracted from junctions
    :param counts: dictionary with sites and counts
    :param stranded: true if data is stranded
    :param strand_mode: F1R2 or F2R1 - determines strand
    :return: does not return anything
    """

    # segment (read) type (0 - read1+, 1 - read1-, 2 - read2+, 3 - read2-)
    seg_type = 2 * segment.is_read2 + segment.is_reverse
    if stranded:
        seg_strand = "+" if seg_type in (1, 2) else '-'  # choose strand
        if strand_mode == "F1R2":  # other forward pair
            seg_strand = ("+", "-")[seg_strand == "+"]
    else:
        seg_strand = "*"

    # initial positions in segment and reference
    seg_pos = 0
    ref_pos = segment.reference_start + BASE

    for cigar_tuple in segment.cigartuples:

        if cigar_tuple[0] in (0, 7, 8):  # alignment match (M), sequence match (=), sequence mismatch (X)

            for pos in range(ref_pos + 1, ref_pos + cigar_tuple[1] - 1):

                if (ref_name, pos) in sites:
                    site = None

                    if stranded:
                        if seg_strand in sites[(ref_name, pos)]:
                            site = SiteWithOffset(
                                site_id="_".join((ref_name, str(pos), seg_strand)),
                                offset=seg_pos + (pos - ref_pos)
                            )
                        else:
                            continue
                    else:
                        for strand in sites[(ref_name, pos)]:
                            site = SiteWithOffset(
                                site_id="_".join((ref_name, str(pos), strand)),
                                offset=seg_pos + (pos - ref_pos)
                            )

                    if site is not None:
                        counts[site] += 1

            seg_pos += cigar_tuple[1]
            ref_pos += cigar_tuple[1]

        elif cigar_tuple[0] == 1:  # insertion to the reference (I)
            seg_pos += cigar_tuple[1]

        elif cigar_tuple[0] in (2, 3):  # deletion from the reference (D), skipped region from the reference (N)
            ref_pos += cigar_tuple[1]

    return None


def alignment2counts(
        alignment: pysam.AlignmentFile,
        sites: DefaultDict[Tuple[str, int], Set],
        unique: bool,
        primary: bool,
        stranded: bool,
        strand_mode: str
) -> Counter[SiteWithOffset]:
    """
    Take an alignment and dictionary of splice sites as input
    and returns dictionary of splice sites with their counts
    """

    trans = prepare_ref_names(alignment)

    counts = Counter()
    segment: pysam.AlignedSegment  # just annotation line

    # iterating through the alignment file
    for segment in alignment.fetch():

        # check if read is mapped
        if segment.is_unmapped:
            continue

        # check if read is not multimapped
        if unique and segment.has_tag("NH") and segment.get_tag("NH") > 1:
            continue

        # if read is multimapped only primary alignment will be considered
        if primary and segment.is_secondary:
            continue

        # if read is multimapped only primary alignment will be considered
        if primary and segment.is_supplementary:
            continue

        ref_name = trans.get(segment.reference_name)
        if ref_name is None:
            continue

        # adding new counts to splice site
        segment2counts(
            segment=segment,
            ref_name=ref_name,
            sites=sites,
            counts=counts,
            stranded=stranded,
            strand_mode=strand_mode
        )

    return counts


def counts2dataframe(counts: Counter[SiteWithOffset]) -> pd.DataFrame:
    """Takes a dictionary of sites with their counts
     and creates a pandas DataFrame, then sorts it by site_id and offset"""

    # list of rows for pandas DataFrame
    rows = []

    # filling this list of rows
    for site, count in counts.items():
        rows.append(
            SiteWithCount(
                site_id=site.site_id,
                offset=site.offset,
                total_count=count
            )
        )

    # creating dataframe
    df = pd.DataFrame(rows, columns=["site_id", "offset", "total_count"])
    # sorting it
    df = df.sort_values(by=["site_id", "offset"])

    return df


def main():
    args = parse_cli_args()
    bam = pysam.AlignmentFile(args["input"], threads=args["threads"])
    sites = junctions2sites(filename=args["junctions"])
    read_length, genome, paired, stranded, library_type = read_stats(args["stats"])
    counts = alignment2counts(
        alignment=bam, sites=sites,
        unique=args["unique"], primary=args["primary"],
        stranded=stranded, strand_mode=library_type
    )
    df = counts2dataframe(counts=counts)
    df.to_csv(args["output"], sep="\t", index=False, header=False, compression="gzip")


if __name__ == '__main__':
    main()
