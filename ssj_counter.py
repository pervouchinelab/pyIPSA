import pysam
import argparse
import sys

BASE = 1


def segment_to_junctions(segment):

    # if not segment.is_proper_pair or segment.is_secondary or segment.is_supplementary:
    #     return []

    # reference name
    reference_name = segment.reference_name
    # reference_name = 'chr' + reference_name if reference_name[:3] != 'chr' else reference_name

    # segment type (0 - read1+, 1 - read1-, 2 - read2+, 3 - read2-)
    segment_type = 2 * segment.is_read2 + segment.is_reverse

    # initial positions in segment and reference
    segment_pos = 0
    reference_pos = segment.reference_start + BASE

    # list of junctions supported by the segment
    segment_junctions = []

    # extracting junctions from cigartuple
    for cigartuple in segment.cigartuples:

        if cigartuple[0] in (0, 7, 8):  # alignment match (M), sequence match (=), sequence mismatch (X)
            segment_pos += cigartuple[1]
            reference_pos += cigartuple[1]
        elif cigartuple[0] == 1:  # insertion to the reference (I)
            segment_pos += cigartuple[1]
        elif cigartuple[0] == 2:  # deletion from the reference (D)
            reference_pos += cigartuple[1]
        elif cigartuple[0] == 3:  # skipped region from the reference (N)
            # a junction is saved as a tuple
            # (chr1_114559_115151, 22, 2)
            # the first position is junction id
            # the second is offset
            # the last one is type of segment
            junction = (
                '_'.join((
                    reference_name,
                    str(reference_pos - 1),  # last exonic position on the left
                    str(reference_pos + cigartuple[1])  # first exonic position on the right
                )),
                segment_pos,
                segment_type
            )
            segment_junctions.append(junction)
            reference_pos += cigartuple[1]

    return segment_junctions


def alignment_to_junctions(alignment):
    # dictionary to store different splice junctions with their offsets
    # storing scheme
    # (chr1_114559_115151, 22): [0, 5, 4, 0]
    all_junctions = dict()
    # iterating through alignments
    for segment in alignment.fetch():

        extracted_junctions = segment_to_junctions(segment)
        # skipping segment alignment without junctions
        if not extracted_junctions:
            continue

        for junction in extracted_junctions:
            # adding newly discovered junction
            if not (junction[:2] in all_junctions):
                all_junctions[junction[:2]] = [0] * 4
            # incrementing number of supporting reads
            all_junctions[junction[:2]][junction[2]] += 1

    return all_junctions


def output_junctions_to_stdout(junctions):
    # yield '\t'.join(('junction_id', 'offset', 'F1', 'R1', 'F2', 'R2')) + '\n'
    for j in junctions:
        line = (j[0], str(j[1]), *map(str, junctions[j]))
        yield '\t'.join(line) + '\n'


def main():
    parser = argparse.ArgumentParser(description="Counting reads supporting each splice site junction")
    parser.add_argument("-i", "--input_bam", type=str, metavar="", required=True, help="Input alignment file (BAM)")
    cli_args = parser.parse_args()
    sam_file = pysam.AlignmentFile(cli_args.input_bam, 'rb')
    junctions_with_counts = alignment_to_junctions(sam_file)
    sys.stdout.writelines(output_junctions_to_stdout(junctions_with_counts))


if __name__ == '__main__':
    main()
