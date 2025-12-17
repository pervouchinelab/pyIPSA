import click
import polars as pl

from gtfparse import read_gtf
from .gtf_utils import get_exons_from_gtf


def get_introns_from_gtf(df: pl.DataFrame) -> pl.DataFrame:
    dfe1 = get_exons_from_gtf(df)
    dfi1 = dfe1\
        .select(['seqname', 'end', 'coord_next', 'strand'])\
        .filter(~pl.col('coord_next').is_null())\
        .unique()\
        .rename({'end': 'start', 'coord_next': 'end'})\
        .sort(by=['seqname', 'start', 'end'])
    return dfi1


def filter_intron_rows(df: pl.DataFrame) -> pl.DataFrame:
    return df.filter(~pl.col('seqname').cast(pl.String).str.contains('_'))


@click.command()
@click.option("--input-gtf", required=True)
@click.option("--output-sj", required=True)
def main(input_gtf, output_sj):
    dfa1: pl.DataFrame = read_gtf(input_gtf)  # pyright: ignore[reportAssignmentType]

    dfi1 = get_introns_from_gtf(dfa1)
    dfi2 = filter_intron_rows(dfi1)

    dfi2.write_csv(output_sj, include_header=False, separator='\t')


if __name__ == "__main__":
    main()