"""
This module contains functions for extracting
splice junctions from the alignment file (BAM) and
counting reads supporting them.
"""
import pysam
import pandas as pd
import argparse
from collections import namedtuple, defaultdict

BASE = 1

JunctionWithReadType = namedtuple("JunctionWithReadType", ["junction_id", "offset", "read_type"])
Junction = namedtuple("Junction", ["junction_id", "offset"])
JunctionWithCounts = namedtuple("JunctionWithCounts", ["junction_id", "offset", "F1", "R1", "F2", "R2"])


def parse_cli_args():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(
        description="Extracts splice junctions from the alignment file \
        and counts all the reads supporting these splicing junctions"
    )
    parser.add_argument(
        "-i", "--input_bam", type=str, metavar="FILE", required=True,
        help="Input alignment file (BAM)", dest="bam"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="Output gzipped file with junctions and their counts (TSV)", dest="tsv"
    )
    parser.add_argument(
        "-u", "--unique", action="store_true", help="Consider only uniquely mapping reads"
    )
    parser.add_argument(
        "-p", "--primary", action="store_true",
        help="Consider only primary alignments of multimapped reads"
    )
    parser.add_argument(
        "-t", "--threads", type=int, metavar="INT", default=1,
        help="Number of threads used to read alignment file"
    )
    args = parser.parse_args()
    return vars(args)


def segment_to_junctions(segment: pysam.AlignedSegment):
    """Takes a read (AlignedSegment) as input and returns
     a list of junctions found in the read"""

    # name of reference sequence where read is aligned
    ref_name = segment.reference_name

    # segment (read) type (0 - read1+, 1 - read1-, 2 - read2+, 3 - read2-)
    seg_type = 2 * segment.is_read2 + segment.is_reverse

    # initial positions in segment and reference
    seg_pos = 0
    ref_pos = segment.reference_start + BASE  # from 0-based to 1-based

    # list of junctions supported by the segment
    seg_junctions = []

    # extracting junctions from cigar string:
    for cigar_tuple in segment.cigartuples:

        if cigar_tuple[0] in (0, 7, 8):  # alignment match (M), sequence match (=), sequence mismatch (X)
            seg_pos += cigar_tuple[1]
            ref_pos += cigar_tuple[1]

        elif cigar_tuple[0] == 1:  # insertion to the reference (I)
            seg_pos += cigar_tuple[1]

        elif cigar_tuple[0] == 2:  # deletion from the reference (D)
            ref_pos += cigar_tuple[1]

        elif cigar_tuple[0] == 3:  # skipped region from the reference (N)
            # a junction is saved as a named tuple with fields junction_id, offset and read type
            junction_id = "_".join((
                ref_name,
                str(ref_pos - 1),  # last position of the left exon (offset)
                str(ref_pos + cigar_tuple[1])  # first position of the right exon
            ))
            # TODO: update dictionary with junctions right here without saving to list etc.
            junction_with_read_type = JunctionWithReadType(junction_id=junction_id,
                                                           offset=seg_pos,
                                                           read_type=seg_type)
            seg_junctions.append(junction_with_read_type)
            ref_pos += cigar_tuple[1]

    return seg_junctions


def alignment_to_junctions(alignment: pysam.AlignmentFile, unique: bool, primary: bool):
    """Takes an alignment as input and returns a dictionary of all found junctions"""
    all_junctions = defaultdict(lambda: [0] * 4)
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

        # junctions found in the segment (read)
        seg_junctions = segment_to_junctions(segment)

        # adding these junctions to the dictionary
        for junction_with_read_type in seg_junctions:
            # junction without read type is used as key in dictionary
            junction = Junction(
                junction_id=junction_with_read_type.junction_id,
                offset=junction_with_read_type.offset
            )
            read_type = junction_with_read_type.read_type
            all_junctions[junction][read_type] += 1

    return all_junctions


def junctions_to_dataframe(junctions: defaultdict):
    """Takes a dictionary of junctions with their counts
     and creates a pandas DataFrame, then sorts it by junction_id and offset"""

    # list of rows for pandas DataFrame
    rows_list = []
    junction: Junction  # just annotation line

    # filling this list of rows
    for junction, counts in junctions.items():

        junction_with_counts = JunctionWithCounts(
            junction_id=junction.junction_id,
            offset=junction.offset,
            F1=counts[0], R1=counts[1], F2=counts[2], R2=counts[3]
        )
        rows_list.append(junction_with_counts)

    # creating dataframe
    df = pd.DataFrame(rows_list, columns=["junction_id", "offset", "F1", "R1", "F2", "R2"])
    # sorting it
    df = df.sort_values(by=["junction_id", "offset"])

    return df


def main():
    args = parse_cli_args()
    bam = pysam.AlignmentFile(args["bam"], threads=args["threads"])
    junctions = alignment_to_junctions(bam, primary=args["primary"], unique=args["unique"])
    df = junctions_to_dataframe(junctions)
    df.to_csv(args["tsv"], sep="\t", index=False, header=None, compression="gzip")


if __name__ == '__main__':
    main()
