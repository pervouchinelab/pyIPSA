"""
This module contains functions for extracting splice junctions (introns)
from GTF/GFF annotation files
"""
import gzip
import argparse
import re
from collections import namedtuple
from typing import Iterator


def parse_cli_args():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(
        description="""
        Extracts splice junctions (introns) from the annotation file
        They are stored as 1-based positions of flanking nucleotides
        """
    )
    parser.add_argument(
        "-i", "--input", type=str, metavar="FILE", required=True,
        help="input annotation file (GTF/GFF)", dest="gtf"
    )
    parser.add_argument(
        "-m", "--min_intron_length", type=int, metavar="INT", default=40,
        help="Minimal intron length", dest="min"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="Output gzipped file with splice junctions", dest="tsv"
    )
    args = parser.parse_args()
    return vars(args)


Exon = namedtuple("Exon", ["ref_name", "start", "end", "strand", "gene_id", "transcript_id"])


def annotation_to_exons(filename: str) -> Iterator[Exon]:
    """This generator takes gzipped GTF/GFF file as input and yields
    only lines corresponding to exons without redundant info"""
    with gzip.open(filename, 'rt') as gtf:
        for line in gtf:

            # filtering out comment lines
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            if "#" in line:
                line = line.split("#")[0].strip()

            # splitting to fields
            ref_name, source, feature, start, end, \
                score, strand, frame, attributes = line.split("\t")
            start, end = int(start), int(end)

            # filtering out non-exon lines
            if feature != "exon" or start >= end:
                continue

            # extracting gene_id and transcript_id from attributes
            gene_id_finder = re.compile('(?<=gene_id) "([\\w.]+)";')
            gene_id = gene_id_finder.search(attributes).group(1)
            transcript_id_finder = re.compile('(?<=transcript_id) "([\\w.]+)";')
            transcript_id = transcript_id_finder.search(attributes).group(1)

            # check them
            if not all((gene_id, transcript_id)):
                continue

            yield Exon(ref_name=ref_name, start=start, end=end, strand=strand,
                       gene_id=gene_id, transcript_id=transcript_id)


def exons_to_junctions(exons: Iterator[Exon], min_intron_length: int) -> Iterator[str]:
    """This generator takes iterable of exons as input
     and yields unique splice junctions in sorted way"""
    transcripts = {}
    # genes = {}

    # iterating through exons
    for exon in exons:
        # aggregating exons to transcripts
        if exon.transcript_id not in transcripts:
            transcripts[exon.transcript_id] = [
                exon.ref_name,
                exon.strand,
                [[exon.start, exon.end]]
            ]
        else:
            transcripts[exon.transcript_id][2].append([exon.start, exon.end])
        # aggregating transcript to genes
        # if exon.gene_id not in genes:
        #     genes[exon.gene_id] = [exon.transcript_id]
        # else:
        #     genes[exon.gene_id].append(exon.transcript_id)

    # sorting, filtering and merging exons
    for transcript_id, [ref_name, strand, exons] in transcripts.items():
        exons.sort()
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] < min_intron_length:
                tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        transcripts[transcript_id] = [ref_name, strand, tmp_exons]

    # unique splice junctions
    junctions = set()
    for ref_name, strand, exons in transcripts.values():
        for i in range(1, len(exons)):
            junctions.add((ref_name, str(exons[i-1][1]), str(exons[i][0]), strand))

    # sorting them
    junctions = sorted(junctions)

    for ref_name, start, end, strand in junctions:
        yield "\t".join((ref_name, str(start), str(end), strand)) + "\n"


def main():
    args = parse_cli_args()
    exons = annotation_to_exons(args["gtf"])
    junctions = exons_to_junctions(exons=exons, min_intron_length=args["min"])
    with gzip.open(args["tsv"], "wt") as f:
        f.writelines(junctions)


if __name__ == '__main__':
    main()
