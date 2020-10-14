"""
This module contains functions for extraction of
splice junctions from the alignment file (BAM) and
functions that count reads which support these splice junctions.
"""
import pysam
import pandas as pd
import argparse
from collections import namedtuple, defaultdict, Counter
from typing import Tuple, List, Dict, DefaultDict
from .ipsa_config import *
from .ipsa_utils import load_splice_sites

H = 60  # Log file horizontal rule length
RawJunction = namedtuple("RawJunction", ["ref_name", "start", "stop", "offset"])
Junction = namedtuple("Junction", ["junction_id", "offset", "F1", "R1", "F2", "R2"])


def parse_cli_args():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(
        description="Extract splice junctions from BAM"
    )
    parser.add_argument(
        "-i", "--input_bam", type=str, metavar="FILE", required=True,
        help="input alignment file (BAM)", dest="bam"
    )
    parser.add_argument(
        "-s", "--splice_sites", type=str, metavar="DIR", required=True,
        help="directory with annotated splice sites"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="output file name", dest="tsv"
    )
    parser.add_argument(
        "-l", "--log", type=str, metavar="FILE", required=True,
        help="output log name"
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


def segment_to_junctions(segment: pysam.AlignedSegment, junctions_with_counts: defaultdict):
    """Takes a read (AlignedSegment) and dictionary of junctions with counts
     as input. Then updates the dictionary with new data extracted from the read"""

    # name of reference sequence where read is aligned
    ref_name = segment.reference_name

    # segment (read) type (0 - read1+, 1 - read1-, 2 - read2+, 3 - read2-)
    seg_type = 2 * segment.is_read2 + segment.is_reverse

    # initial positions in segment and reference
    seg_pos = 0
    ref_pos = segment.reference_start + BASE  # from 0-based to 1-based

    # extracting junctions from cigar string:
    for cigar_tuple in segment.cigartuples:

        if cigar_tuple[0] in (0, 7, 8):
            # alignment match (M), sequence match (=) and mismatch (X)
            seg_pos += cigar_tuple[1]
            ref_pos += cigar_tuple[1]

        elif cigar_tuple[0] == 1:
            # insertion to the reference (I)
            seg_pos += cigar_tuple[1]

        elif cigar_tuple[0] == 2:
            # deletion from the reference (D)
            ref_pos += cigar_tuple[1]

        elif cigar_tuple[0] == 3:
            # skipped region from the reference (N)
            raw_junction = RawJunction(
                ref_name=ref_name,
                start=str(ref_pos - 1),  # last position of the left exon (offset)
                stop=str(ref_pos + cigar_tuple[1]),  # first position of the right exon,
                offset=seg_pos
            )
            # updating dictionary of junctions
            junctions_with_counts[raw_junction][seg_type] += 1
            ref_pos += cigar_tuple[1]

    return None


def alignment_to_junctions(alignment: pysam.AlignmentFile, unique: bool, primary: bool) -> \
        Tuple[DefaultDict[RawJunction, List[int]], int]:
    """Takes an alignment file as input.
    Returns a dictionary of all junctions with counts and inferred read length"""
    segment: pysam.AlignedSegment  # just annotation line
    junctions_with_counts = defaultdict(lambda: [0] * 4)
    read_lengths_counter = Counter()

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

        read_lengths_counter[segment.infer_read_length()] += 1

        # adding new junctions if present and updating counts
        segment_to_junctions(segment=segment, junctions_with_counts=junctions_with_counts)

    return junctions_with_counts, read_lengths_counter.most_common()[0][0]


def junctions_to_dataframe(junctions_with_counts: DefaultDict[RawJunction, List[int]]) -> pd.DataFrame:
    """Takes a dictionary of junctions with their counts, filters them
     and creates a pandas DataFrame, then sorts it by junction_id and offset"""

    # annotation lines
    junction: Junction
    raw_junction: RawJunction

    rows_list = []

    for raw_junction, counts in junctions_with_counts.items():
        ref_name, start, stop, offset = raw_junction

        # filtration of reference names
        if ref_name not in ALLOWED_REFERENCE_NAMES:
            lowered = ref_name.lower()
            if lowered in TRANSLATION:
                ref_name = TRANSLATION[lowered]
            else:
                continue

        junction = Junction(
            junction_id="_".join([ref_name, start, stop]),
            offset=offset,
            F1=counts[0], R1=counts[1], F2=counts[2], R2=counts[3]
        )
        rows_list.append(junction)

    df = pd.DataFrame(rows_list, columns=["junction_id", "offset", "F1", "R1", "F2", "R2"])
    df = df.sort_values(by=["junction_id", "offset"])

    return df


def guess_genome(df: pd.DataFrame, dir_with_splice_sites: str) -> Dict[str, int]:
    """Takes DataFrame of junctions and path to directory with annotated splice sites.
    Tries to guess genome (organism) by intersecting novel junctions with annotated ones."""
    # preparing novel start-stop pairs
    pairs = df["junction_id"].apply(lambda jid: tuple(map(int, jid.split("_")[1:])))
    pairs = set(pairs)
    # pairs from annotations
    pairs_by_genome = load_splice_sites(dir_with_splice_sites)

    hits_by_genome = dict()
    for genome in pairs_by_genome:
        hits_by_genome[genome] = len(pairs & pairs_by_genome[genome])
    return hits_by_genome


def get_stats(df: pd.DataFrame) -> List[str]:
    """Takes DataFrame of junctions and outputs list of alignment statistics"""
    a = []
    dfi = df.set_index(['junction_id', 'offset'])

    # How many different read pairs?
    col_sums = dfi.sum(axis=0)
    a += ["Sum of counts for each read type:\n", f"{col_sums.to_string()}\n\n"]

    # Are reads paired?
    paired = all(col_sums > 0)
    a += [f"Guess: the data is {('single', 'pair')[paired]}-end\n", "-" * H + "\n"]

    # Correlation analysis
    dfi['F1+R2'] = dfi['F1'] + dfi['R2']
    dfi['F2+R1'] = dfi['F2'] + dfi['R1']
    corr_mat = dfi.groupby("junction_id").sum().corr()
    a += ["Correlation matrix:\n\n", f"{corr_mat.round(2).to_string()}\n\n"]

    # Is data stranded?
    if paired:
        stranded = bool(corr_mat.loc['F1+R2', 'F2+R1'] <= 0.5)
    else:
        if all(col_sums[:2]):
            stranded = bool(corr_mat.loc['F1', 'R1'] <= 0.5)
        else:
            stranded = bool(corr_mat.loc['F2', 'R2'] <= 0.5)
    a += [f"Guess: the data is {('un', '')[stranded]}stranded\n", "-" * H + "\n"]

    # Distribution of offsets
    offsets = dfi.groupby(by="offset").sum()
    a += ["Offsets distribution:\n\n", offsets.reset_index().to_string(index=False)]

    return a


def write_analysis(filename: str, read_length: int,
                   hits_by_genome: Dict[str, int], stats: List[str]):
    tmp = [f"Guess: genome is {max(hits_by_genome, key=hits_by_genome.get)}\n"]
    tmp += [f"Number of hits with {genome} genome: {hits}\n" for genome, hits in hits_by_genome.items()]
    tmp += ["-" * H + "\n", f"Read length is {read_length}\n", "-" * H + "\n"]
    with open(filename, "w") as txt:
        txt.writelines(tmp + stats)


def main():
    args = parse_cli_args()
    # read BAM and get junctions
    bam = pysam.AlignmentFile(args["bam"], threads=args["threads"])
    junctions, read_length = alignment_to_junctions(alignment=bam,
                                                    primary=args["primary"],
                                                    unique=args["unique"])
    df_junctions = junctions_to_dataframe(junctions_with_counts=junctions)
    # guess the genome
    hits_by_genome = guess_genome(df=df_junctions,
                                  dir_with_splice_sites=args["splice_sites"])
    # get alignment stats
    stats = get_stats(df=df_junctions)
    # write outputs
    write_analysis(filename=args["log"], read_length=read_length,
                   hits_by_genome=hits_by_genome, stats=stats)
    df_junctions.to_csv(args["tsv"], sep="\t", index=False,
                        header=None, compression="gzip")


if __name__ == '__main__':
    main()
