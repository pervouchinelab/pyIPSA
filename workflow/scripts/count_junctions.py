"""
This module contains functions for extraction of
splice junctions from the alignment file (BAM) and
functions that count reads which support these splice junctions.
"""
import argparse
from collections import namedtuple, defaultdict, Counter
from pathlib import Path
from typing import Any, DefaultDict, Dict, List, Tuple


import pandas as pd
import pysam

from .ipsa_config import *
from .ipsa_utils import load_pairs

H = 60  # Log file horizontal rule length
RawJunction = namedtuple("RawJunction", ["ref_name", "start", "stop", "offset"])
Junction = namedtuple("Junction", ["junction_id", "offset", "F1", "R1", "F2", "R2"])


def parse_cli_args() -> Dict:
    """
    Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Extract splice junctions from BAM"
    )
    parser.add_argument(
        "-i", "--input", type=str, metavar="FILE", required=True,
        help="input alignment file (BAM)"
    )
    parser.add_argument(
        "-k", "--known", type=str, metavar="DIR", required=True,
        help="directory with annotated junctions"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="output file name (step 1)"
    )
    parser.add_argument(
        "-l", "--lib_stats", type=str, metavar="FILE", required=True,
        help="filename to save library parameters and some stats"
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


def segment2junctions(segment: pysam.AlignedSegment, junctions_with_counts: defaultdict):
    """
    Take a read (AlignedSegment) and add new counts to junctions.

    :param segment: aligned segment (read)
    :param junctions_with_counts: dictionary storing junctions and their counts
    """

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

    return


def alignment2junctions(
        alignment: pysam.AlignmentFile,
        unique: bool,
        primary: bool
) -> Tuple[DefaultDict[RawJunction, List[int]], int]:
    """
    Find all junctions in alignment file and compute their counts.
    Also determine most common read length.

    :param alignment: alignment file
    :param unique: account only uniquely mapped (aligned) segments
    :param primary: account only primary alignment if segment is a multi-mapper
    :return: dictionary mapping junction to its counts in each offset and most common read length
    """
    segment: pysam.AlignedSegment  # just annotation line for convenience
    junctions_with_counts = defaultdict(lambda: [0] * 4)
    read_lengths_counter = Counter()

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

        read_lengths_counter[segment.infer_read_length()] += 1

        # adding new junctions if present and updating counts
        segment2junctions(segment=segment, junctions_with_counts=junctions_with_counts)

    return junctions_with_counts, read_lengths_counter.most_common()[0][0]


def junctions2dataframe(junctions_with_counts: DefaultDict[RawJunction, List[int]]) -> pd.DataFrame:
    """
    Collect all junctions into sorted DataFrame.

    :param junctions_with_counts: dictionary mapping junction to its counts in each offset
    :return: DataFrame with junctions and their counts in each offset
    """
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


def guess_genome(
        df: pd.DataFrame,
        dir_path: str
) -> Dict[str, int]:
    """
    Compute how many newly discovered junctions intersect with known ones in all available genomes.
    Genome with most hits is the most probable source of newly discovered junctions.

    :param df: DataFrame with junctions and their counts in each offset
    :param dir_path: directory storing known junctions from available genomes
    :return: dictionary mapping genome to number of hits
    """
    # prepare start-stop pairs from novel junctions
    pairs = df["junction_id"].apply(lambda junction_id: tuple(map(int, junction_id.split("_")[1:])))
    pairs = set(pairs)
    # prepare start-stop from known junctions
    pairs_by_genome = load_pairs(dir_path)

    hits_by_genome = dict()
    for genome in pairs_by_genome:
        hits_by_genome[genome] = len(pairs & pairs_by_genome[genome])

    return hits_by_genome


def guess_lib_type(
        df: pd.DataFrame,
        dir_path: str,
        genome: str
) -> str:
    """
    Guess library type by intersecting novel junctions with known ones.

    :param df: DataFrame with junctions and their counts in each offset
    :param dir_path: directory storing known junctions from available genomes
    :param genome: genome string
    :return: library type - F2R1 or F1R2 or unknown
    """
    known = pd.read_table(Path(dir_path, f"{genome}.ss.tsv.gz"), header=None)
    known["junction_id"] = known[0] + "_" + known[1].astype(str) + "_" + known[2].astype(str)
    known = known.set_index("junction_id").drop([0, 1, 2], axis=1)
    known.columns = ["strand"]

    novel = df.drop(["offset", "F1", "R1", "F2", "R2"], axis=1).groupby(["junction_id"]).sum()

    table = known.join(novel, on=["junction_id"], how="inner").groupby("strand").sum()

    if table.loc["+", "F2+R1"] > 10 * table.loc["+", "F1+R2"]:
        return "F2R1"
    elif table.loc["+", "F1+R2"] > 10 * table.loc["+", "F2+R1"]:
        return "F1R2"
    else:
        return "unknown"


def compute_lib_params(
        df: pd.DataFrame,
        dir_path: str,
        read_length: int
) -> Dict[str, Any]:
    """
    Compute library parameters and some stats from junctions DataFrame.

    :param df: DataFrame with junctions and their counts in each offset
    :param dir_path: directory storing known junctions from available genomes
    :param read_length: precomputed read length
    :return: dictionary of library parameters
    """
    params = {"read length": read_length}
    df_indexed = df.set_index(["junction_id", "offset"])

    # Guess organism and its genome version
    params["hits by genome"] = guess_genome(df=df, dir_path=dir_path)
    params["genome"] = max(params["hits by genome"], key=params["hits by genome"].get)

    # Total counts by read type
    params["counts"] = df_indexed.sum(axis=0)

    # Are reads paired?
    params["paired"] = (params["counts"] > 0).all()

    # Correlation matrix for read pairs
    df_indexed["F1+R2"] = df_indexed["F1"] + df_indexed["R2"]
    df_indexed["F2+R1"] = df_indexed["F2"] + df_indexed["R1"]
    correlation_matrix = df_indexed.groupby("junction_id").sum().corr()
    params["correlation matrix"] = correlation_matrix

    # Is data stranded?
    if params["paired"]:
        params["stranded"] = bool(correlation_matrix.loc["F1+R2", "F2+R1"] <= 0.5)
    else:
        if (params["counts"][:2] > 0).all() and (params["counts"][2:] == 0).all():
            params["stranded"] = bool(correlation_matrix.loc["F1", "R1"] <= 0.5)
        else:
            params["stranded"] = bool(correlation_matrix.loc["F2", "R2"] <= 0.5)

    # Library - F1R2 or F2R1?
    if params["stranded"]:
        params["library type"] = guess_lib_type(
            df=df,
            dir_path=dir_path,
            genome=params["genome"]
        )
    else:
        params["library type"] = "Unstranded"

    # Offsets distribution
    params["offsets"] = df_indexed.groupby(by="offset").sum().reset_index()

    return params


def write_lib_params(
        params: Dict[str, Any],
        filename: str
):
    with open(filename, "w") as f:
        for par in ("read length", "genome", "paired", "stranded", "library type"):
            f.write(f"{par}: {params[par]}\n")
        sep = "-" * H + "\n"
        f.write(sep)
        f.write("Number of hits for each genome:\n")
        f.writelines(f"{genome}: {hits}\n" for genome, hits in params["hits by genome"].items())
        f.write(sep)
        f.write("Number of reads of each type:\n")
        f.write(params["counts"].to_string() + "\n")
        f.write(sep)
        f.write("Correlation matrix:\n\n")
        f.write(params["correlation matrix"].to_string() + "\n")
        f.write(sep)
        f.write("Offsets distribution:\n\n")
        f.write(params["offsets"].to_string())


def main():
    args = parse_cli_args()
    # read BAM and get junctions
    bam = pysam.AlignmentFile(args["input"], threads=args["threads"])
    junctions, read_length = alignment2junctions(
        alignment=bam, primary=args["primary"], unique=args["unique"]
    )
    df_junctions = junctions2dataframe(junctions_with_counts=junctions)
    # get library parameters
    lib_params = compute_lib_params(
        df=df_junctions,
        dir_path=args["known"],
        read_length=read_length
    )
    # write outputs
    df_junctions.to_csv(args["output"], sep="\t", index=False, header=False, compression="gzip")
    write_lib_params(params=lib_params, filename=args["lib_stats"])


if __name__ == "__main__":
    main()
