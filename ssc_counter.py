import pysam
import sys
import time

BASE = 1


def timer(func):

    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        sys.stderr.write(f'Execution time of "{func.__name__}" function: '
                         f'{round(time.time() - start, 2)} sec\n')
        return result

    return wrapper


@timer
def junctions_to_sites(filename):
    sites = {}

    with open(filename, 'r') as junctions_file:
        for line in junctions_file:
            junction_id = line.strip().split('\t')[0]
            reference_name, left, right = junction_id.split('_')
            sites[(reference_name, int(left))] = True
            sites[(reference_name, int(right))] = True

    return sites


@timer
def sites_coverage(alignment, sites):

    coverage = dict()

    for segment in alignment.fetch():

        reference_name = segment.reference_name
        segment_pos = 0
        reference_pos = segment.reference_start + BASE

        for cigartuple in segment.cigartuples:

            if cigartuple[0] == 0:
                for inc in range(1, cigartuple[1] - 1):
                    matched_position = reference_pos + inc
                    if sites.get((reference_name, matched_position), False):
                        offset = segment_pos + inc
                        key = (reference_name + '_' + str(matched_position), offset)
                        if key in coverage:
                            coverage[key] += 1
                        else:
                            coverage[key] = 1

                segment_pos += cigartuple[1]
                reference_pos += cigartuple[1]

            elif cigartuple[0] == 1:
                segment_pos += cigartuple[1]
            elif cigartuple[0] in (2, 3):
                reference_pos += cigartuple[1]

    return coverage


@timer
def output(sites):
    for key in sites:
        line = (key[0], str(key[1]), str(sites[key]))
        yield '\t'.join(line) + '\n'


@timer
def main():
    splicing_sites = junctions_to_sites(sys.argv[1])
    samfile = pysam.AlignmentFile(sys.argv[2], 'rb')
    sys.stdout.writelines(output(sites_coverage(samfile, splicing_sites)))


if __name__ == '__main__':
    main()
