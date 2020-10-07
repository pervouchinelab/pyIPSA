import argparse
import pandas as pd
import pysam
from collections import Counter


def parse_cli_args():
    """Parses command line arguments"""
    parser = argparse.ArgumentParser(
        description="""
            Uses file with junctions and alignment file to guess organism 
            (genome and its version) and sequencing library parameters:
            read length, if the data is stranded or not, etc.
            """
    )
    parser.add_argument(
        "-j", "--junctions", type=str, metavar="FILE", required=True,
        help="Input junctions file (TSV)"
    )
    parser.add_argument(
        "-i", "--input_bam", type=str, metavar="FILE", required=True,
        help="Input alignment file (BAM)", dest="bam"
    )
    parser.add_argument(
        "-o", "--output", type=str, metavar="FILE", required=True,
        help="Output file file with analysis (TXT)"
    )
    args = parser.parse_args()
    return vars(args)


def infer_read_length(filename: str) -> int:
    """This function takes alignment file (BAM)
     and returns its most common read length"""
    c = Counter()
    with pysam.AlignmentFile(filename, "rb") as bam:
        i = 0
        for segment in bam.fetch():
            c[segment.infer_read_length()] += 1
            i += 1
            if i == 10**5:
                break
    return c.most_common()[0][0]


def print_stats(sam_file, ssj_file):
    # print(f"Read length = {Counter(read_lengths).most_common()[0][0]}")
    junctions = pd.read_table(ssj_file, header=None)
    junctions.columns = ['junction_id', 'offset', 'F1', 'R1', 'F2', 'R2']
    junctions = junctions.set_index(['junction_id', 'offset'])
    # column sums
    col_sums = junctions.sum(axis=0)
    print(f'Sum of counts for each read type:\n{col_sums.to_string()}')
    paired = all(col_sums > 0)
    print(f'Guessing the data is {("single", "pair")[paired]}-end')
    # computing sums of pairs
    junctions['F1+R2'] = junctions['F1'] + junctions['R2']
    junctions['F2+R1'] = junctions['F2'] + junctions['R1']
    # correlation
    print(f'\n{"-" * 60}\n')
    correlation_matrix = junctions.corr()
    print(f'Correlation matrix:\n{correlation_matrix.to_string()}')
    if paired:
        stranded = bool(correlation_matrix.loc['F1+R2', 'F2+R1'] <= 0.5)
    else:
        if all(col_sums[:2]):
            stranded = bool(correlation_matrix.loc['F1', 'R1'] <= 0.5)
        else:
            stranded = bool(correlation_matrix.loc['F2', 'R2'] <= 0.5)
    print(f'Guessing the data is {("un", "")[stranded]}stranded')
    print(f'\n{"-" * 60}\n')
    # offset distribution
    print('Offset distribution:')
    print(junctions.groupby(by='offset').sum().iloc[:, :4].to_string())


def main():
    pass


if __name__ == '__main__':
    main()
