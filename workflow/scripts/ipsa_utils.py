import gzip
import sys
import time
from typing import List, Set, Tuple, Dict

from .ipsa_config import *


def load_pairs(dir_path: str) -> Dict[str, Set[Tuple[int, int]]]:
    """
    Load start-stop pairs from all known junctions.
    
    :param dir_path: directory storing known junctions
    :return: dictionary mapping genome to all its junctions' start-stop pairs
    """
    junctions_by_genome = dict()

    for genome in AVAILABLE_GENOMES:
        pairs = set()
        with gzip.open(dir_path + f"/{genome}.ss.tsv.gz", "rt") as f:
            for line in f:
                _, left, right, _ = line.strip().split("\t")
                left, right = int(left), int(right)
                pairs.add((left, right))
        junctions_by_genome[genome] = pairs

    return junctions_by_genome


def ref_names_from_file(filename: str) -> List[str]:
    """
    Extract reference names from annotation file or file with known junctions.

    :param filename: file with known junctions or any tab-delimited file with reference names in the first column
    :return: sorted list of reference names
    """
    a = set()
    with gzip.open(filename, 'rt') as gtf:
        for line in gtf:

            # filtering out comment lines
            line = line.strip()
            if line.startswith("#") or not line:
                continue
            if "#" in line:
                line = line.split("#")[0].strip()

            # splitting to fields
            ref_name, *_ = line.split("\t")

            a.add(ref_name)

    a = list(a)
    a.sort()
    return a


def timer(func):
    """Execution time profiling decorator."""
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        sys.stderr.write(f'Execution time of "{func.__name__}" function: '
                         f'{round(time.time() - start, 2)} sec\n')
        return result

    return wrapper
