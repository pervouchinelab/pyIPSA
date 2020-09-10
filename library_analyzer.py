import argparse
import pandas as pd
import pysam
from collections import Counter


def print_stats(sam_file, ssj_file):
    read_lengths = []
    i = 0
    with pysam.AlignmentFile(sam_file, "rb") as alignment:
        for segment in alignment.fetch():
            read_lengths.append(segment.infer_read_length())
            i += 1
            if i == 100000:
                break
    print(f"Read length = {Counter(read_lengths).most_common()[0][0]}")
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
    parser = argparse.ArgumentParser(description="Analyzing the RNA-Seq library type")
    parser.add_argument("-i", "--input_bam", type=str,
                        metavar="", required=True, help="Input alignment file (BAM)")
    parser.add_argument("-ssj", "--splice_site_junctions", type=str,
                        metavar="", required=True, help="Splice site junctions file")
    cli_args = parser.parse_args()
    print_stats(cli_args.input_bam, cli_args.splice_site_junctions)


if __name__ == '__main__':
    main()
