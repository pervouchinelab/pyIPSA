import argparse
import pandas as pd
import pysam


def read_ssj(filename):
    ssj = pd.read_table(filename, header=None)
    ssj.columns = ["junction_id", "total_count", "staggered_count", "entropy"]
    return ssj


def add_splice_sites(df, genome):
    genome = pysam.FastaFile(genome)
    complement = {"A": "T", "G": "C", "C": "G", "T": "A", "N": "N"}
    splice_sites = []
    for junction_id in df["junction_id"]:
        ref, start, stop, strand = junction_id.split("_")
        ref, start, stop = "chr" + ref, int(start), int(stop)
        c = strand == "-"
        n = genome.fetch(ref, start, start + 2).upper() + genome.fetch(ref, stop - 3, stop - 1).upper()
        if c:
            n = "".join(complement[n] for n in n)[::-1]
        splice_sites.append(n)
    df["splice_sites"] = splice_sites
    return df


def main():
    parser = argparse.ArgumentParser(description="Annotating splice site junctions")
    parser.add_argument("-ssj", "--splice_site_junctions", type=str,
                        metavar="", required=True, help="Splice junctions file (TSV)")
    # parser.add_argument("-a", "--annotated_splice_sites", type=argparse.FileType("r"),
    #                     metavar="", required=True, help="Annotated splice sites file")
    parser.add_argument("-fa", "--genome_fa", type=str,
                        metavar="", required=True, help="Genome file (FASTA)")
    parser.add_argument("-o", "--output_tsv", type=str,
                        metavar="", required=True, help="Output annotated counts file (TSV")
    cli_args = parser.parse_args()
    agg_ssj = read_ssj(cli_args.splice_site_junctions)
    ann_ssj = add_splice_sites(agg_ssj, cli_args.genome_fa)
    ann_ssj.to_csv(cli_args.output_tsv, sep="\t", index=False, header=None, compression="gzip")


if __name__ == '__main__':
    main()
