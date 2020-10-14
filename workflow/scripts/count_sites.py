"""
This module contains functions for extraction of
splice sites from the alignment file (BAM)
"""
import gzip
import argparse
from collections import defaultdict, namedtuple
from typing import Tuple, List, Dict, DefaultDict

import pysam
import pandas as pd

from .ipsa_config import *

Site = namedtuple("Site", ["site_id", "offset"])
SiteWithCounts = namedtuple("Site", ["site_id", "offset", "F1", "R1", "F2", "R2"])


def parse_cli_args():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(
        description="Extracts splice sites from BAM"
    )
    parser.add_argument(
        "-i", "--input_bam", type=str, metavar="FILE", required=True,
        help="input alignment file (BAM)", dest="bam"
    )
    parser.add_argument(
        "-j", "--junctions", type=str, metavar="FILE", required=True,
        help="input file with junctions (.gz)"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="output file name", dest="tsv"
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


def junctions_to_splice_sites(filename: str) -> DefaultDict[str, set]:
    """Takes a file with splice junctions as input and returns
    a dictionary of all splice sites
    """
    sites = defaultdict(set)

    # read junctions file
    with gzip.open(filename, 'rt') as junctions_file:
        for line in junctions_file:
            junction_id = line.strip().split('\t')[0]
            ref_name, left, right = junction_id.split('_')
            # new keys for sites in dictionary
            sites[ref_name].update({int(left), int(right)})

    return sites


def segment_to_sites(segment: pysam.AlignedSegment,
                     ref_name: str,
                     sites: DefaultDict[str, set],
                     sites_with_counts: defaultdict):
    """Takes a read (AlignedSegment) and dictionary of splice sites
    and updates dictionary of splice sites with counts"""

    # segment (read) type (0 - read1+, 1 - read1-, 2 - read2+, 3 - read2-)
    seg_type = 2 * segment.is_read2 + segment.is_reverse

    # initial positions in segment and reference
    seg_pos = 0
    ref_pos = segment.reference_start + BASE

    # extracting
    for cigar_tuple in segment.cigartuples:

        if cigar_tuple[0] in (0, 7, 8):  # alignment match (M), sequence match (=), sequence mismatch (X)

            # # iterating through positions where read is matched
            # for inc in range(1, cigar_tuple[1] - 1):
            #     pos = ref_pos + inc
            #
            #     # check if such position among splice sites
            #     if pos in sites[ref_name]:
            #         # add new count to splice site
            #         offset = seg_pos + inc
            #         site = Site(site_id="_".join((ref_name, str(pos))), offset=offset)
            #         sites_with_counts[site][seg_type] += 1

            positions = {ref_pos + inc for inc in range(1, cigar_tuple[1] - 1)}
            intersection = positions & sites[ref_name]
            for pos in intersection:
                inc = pos - ref_pos
                offset = seg_pos + inc
                site = Site(site_id="_".join((ref_name, str(pos))), offset=offset)
                sites_with_counts[site][seg_type] += 1

            seg_pos += cigar_tuple[1]
            ref_pos += cigar_tuple[1]

        elif cigar_tuple[0] == 1:  # insertion to the reference (I)
            seg_pos += cigar_tuple[1]

        elif cigar_tuple[0] in (2, 3):  # deletion from the reference (D), skipped region from the reference (N)
            ref_pos += cigar_tuple[1]

    return None


def prepare_ref_names(alignment: pysam.AlignmentFile) -> Dict[str, str]:
    # TODO: document the function
    trans = dict()
    for ref_name in alignment.references:
        if ref_name not in ALLOWED_REFERENCE_NAMES:
            lowered = ref_name.lower()
            if lowered in TRANSLATION:
                trans[ref_name] = TRANSLATION[lowered]
        else:
            trans[ref_name] = ref_name
    return trans


def alignment_to_sites(alignment: pysam.AlignmentFile, sites: defaultdict, unique: bool, primary: bool, ):
    """Takes an alignment and dictionary of splice sites as input
     and returns dictionary of splice sites with their counts"""

    # TODO: update reference names if needed
    trans = prepare_ref_names(alignment)

    sites_with_counts = defaultdict(lambda: [0] * 4)
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

        # TODO: update reference name
        ref_name = trans.get(segment.reference_name)
        if ref_name is None:
            continue

        # adding new counts to splice site
        segment_to_sites(segment=segment, ref_name=ref_name,
                         sites=sites, sites_with_counts=sites_with_counts)

    return sites_with_counts


def sites_to_dataframe(sites_with_counts: defaultdict):
    """Takes a dictionary of sites with their counts
     and creates a pandas DataFrame, then sorts it by site_id and offset"""

    # list of rows for pandas DataFrame
    rows_list = []
    site: Site  # just annotation line

    # filling this list of rows
    for site, counts in sites_with_counts.items():
        site_with_counts = SiteWithCounts(
            site_id=site.site_id,
            offset=site.offset,
            F1=counts[0], R1=counts[1], F2=counts[2], R2=counts[3]
        )
        rows_list.append(site_with_counts)

    # creating dataframe
    df = pd.DataFrame(rows_list, columns=["site_id", "offset", "F1", "R1", "F2", "R2"])
    # sorting it
    df = df.sort_values(by=["site_id", "offset"])

    return df


def main():
    # TODO: document
    args = parse_cli_args()
    bam = pysam.AlignmentFile(args["bam"], threads=args["threads"])
    sites = junctions_to_splice_sites(args["junctions"])
    sites_with_counts = alignment_to_sites(bam, primary=args["primary"],
                                           unique=args["unique"], sites=sites)
    df_sites = sites_to_dataframe(sites_with_counts)
    df_sites.to_csv(args["tsv"], sep="\t", index=False,
                    header=None, compression="gzip")


if __name__ == '__main__':
    main()
