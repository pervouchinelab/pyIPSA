import polars as pl


def get_exons_from_gtf(df1: pl.DataFrame, feature: str = "exon") -> pl.DataFrame:
	# get only exons; for each exon find left and right junctions from transcripts; remove terminal
	df2 = df1.filter(pl.col('feature') == feature)
	df2 = df2\
        .sort(by=['transcript_id', 'start'], descending=False)\
        .with_columns(
            pl.col('end').shift(1).over('transcript_id').alias('coord_prev'),
            pl.col('start').shift(-1).over('transcript_id').alias('coord_next')
	    )\
        .rename({'exon_id': 'exon_id_orig'})\
		.with_columns(
        (
            pl.col("seqname")
            + "_"
            + pl.col("coord_prev").cast(str)
            + "_"
            + pl.col("start").cast(str)
            + "_"
            + pl.col("strand")
        ).alias("junction_id_l"),
        (
            pl.col("seqname")
            + "_"
            + pl.col("end").cast(str)
            + "_"
            + pl.col("coord_next").cast(str)
            + "_"
            + pl.col("strand")
        ).alias("junction_id_r"),
        (
            pl.col("seqname")
            + "_"
            + pl.col("coord_prev").cast(str)
            + "_"
            + pl.col("coord_next").cast(str)
            + "_"
            + pl.col("strand")
        ).alias("junction_id_o"),
        (
            pl.col("seqname")
            + "_"
            + pl.col("start").cast(str)
            + "_"
            + pl.col("end").cast(str)
            + "_"
            + pl.col("strand")
        ).alias("exon_id"),
        (
            pl.col("seqname")
            + "_"
            + pl.col("coord_prev").cast(str)
            + "_"
            + pl.col("start").cast(str)
            + "_"
            + pl.col("end").cast(str)
            + "_"
            + pl.col("coord_next").cast(str)
            + "_"
            + pl.col("strand")
        ).alias("exon_id_full")
        )
	return df2