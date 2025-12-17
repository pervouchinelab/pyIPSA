import pysam 

import click
import polars as pl

from typing import Any


def annotate_single_junction(row: dict[str, Any], fa: pysam.FastaFile):
    complement = {"A": "T", "G": "C", "C": "G", "T": "A", "N": "N"}
    c = (row['strand'] == "-")  # flag for reverse complement

    try:
        seq = fa.fetch(row['chrom'], row['start'], row['start'] + 2).upper() + \
            fa.fetch(row['chrom'], row['end'] - 3, row['end'] - 1).upper()
    except KeyError:
        return None

    if c:
        seq = "".join(complement[n] for n in seq)[::-1]

    return seq


def annotate_junctions(dfj: pl.DataFrame, fa: pysam.FastaFile) -> pl.DataFrame:
    seql = [annotate_single_junction(row, fa) for row in dfj.iter_rows(named=True)]
    return dfj.with_columns(
        seq = pl.Series(seql, dtype=pl.String)
    )


def rank_junction_seq(dfj) -> pl.DataFrame:
    return dfj.group_by('seq').agg(pl.len()).filter(~pl.col('seq').is_null()).sort(by='len', descending=True)\


@click.command()
@click.option("--input-sj", required=True)
@click.option("--input-fasta", required=True)
@click.option("--output-sj-ranked", required=True)
def main(input_sj, input_fasta, output_sj_ranked):
    dfj0 = pl.read_csv(input_sj, separator='\t', has_header=False, new_columns=['chrom', 'start', 'end', 'strand'])
    genome_fa = pysam.FastaFile(input_fasta)

    dfj1 = annotate_junctions(dfj0, genome_fa)
    dfj2 = rank_junction_seq(dfj1)

    dfj2.write_csv(output_sj_ranked, include_header=False, separator='\t')


if __name__ == "__main__":
    main()