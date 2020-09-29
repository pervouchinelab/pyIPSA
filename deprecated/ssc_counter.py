import pysam
import gzip
import argparse
import sys

BASE = 1


def junctions_to_splice_sites(filename):
    sites = {}

    with gzip.open(filename, 'rt') as junctions_file:
        for line in junctions_file:
            junction_id = line.strip().split('\t')[0]
            reference_name, left, right = junction_id.split('_')
            sites[(reference_name, int(left))] = True
            sites[(reference_name, int(right))] = True

    return sites


def sites_coverage(alignment: pysam.AlignmentFile, sites: dict):

    coverage = dict()

    for segment in alignment.fetch():

        reference_name = segment.reference_name
        segment_type = 2 * segment.is_read2 + segment.is_reverse
        segment_pos = 0
        reference_pos = segment.reference_start + BASE

        for cigartuple in segment.cigartuples:

            if cigartuple[0] == 0:
                for inc in range(1, cigartuple[1] - 1):
                    matched_position = reference_pos + inc
                    if sites.get((reference_name, matched_position), False):
                        offset = segment_pos + inc
                        key = (reference_name + '_' + str(matched_position), offset)
                        if key not in coverage:
                            coverage[key] = [0] * 4
                        coverage[key][segment_type] += 1

                segment_pos += cigartuple[1]
                reference_pos += cigartuple[1]

            elif cigartuple[0] == 1:
                segment_pos += cigartuple[1]
            elif cigartuple[0] in (2, 3):
                reference_pos += cigartuple[1]

    return coverage


def output_sites_to_stdout(sites):
    for key in sites:
        line = (key[0], str(key[1]), *map(str, sites[key]))
        yield '\t'.join(line) + '\n'


def main():
    parser = argparse.ArgumentParser(description="Counting continuous reads for each splice site")
    parser.add_argument("-i", "--input_bam", type=str,
                        metavar="", required=True, help="Input alignment file (BAM)")
    parser.add_argument("-ssj", "--splice_site_junctions", type=str,
                        metavar="", required=True, help="Splice site junctions file")
    cli_args = parser.parse_args()
    splice_sites = junctions_to_splice_sites(cli_args.splice_site_junctions)
    samfile = pysam.AlignmentFile(cli_args.input_bam, 'rb')
    sys.stdout.writelines(output_sites_to_stdout(sites_coverage(samfile, splice_sites)))


if __name__ == '__main__':
    main()
