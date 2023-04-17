import argparse
from pathlib import Path
from typing import Dict, List
import pandas as pd
import gzip


OUTPUT_DIR = config["output_dir"]
THREADS = config["threads"]
POOLING_MODE = "P" if config["pooled"] else ""

SAMPLES_TABLE = pd.read_table(config["input_files"])
SAMPLES = SAMPLES_TABLE["name"]
SAMPLES_DICT = {r["name"] : r["path"] for i, r in SAMPLES_TABLE.iterrows()}


def get_bam_by_sample(wildcards):
	return SAMPLES_DICT[wildcards["sample"]]


def get_bai_by_sample(wildcards):
	return SAMPLES_DICT[wildcards["sample"]] + ".bai"


def gather_library_stats2df(fnames, names) -> pd.DataFrame:
    """Read replicates' stats and gather all the info into one table."""
    records = []

    for p, name in zip(map(Path, fnames), names):

        with p.open("r") as f:
            for line in f:
                if line.startswith("-"):
                    break
                left, right = line.strip().split(": ")
                if left[0] == "r":
                    read_length = right
                elif left[0] == "g":
                    genome = right
                elif left[0] == "p":
                    paired = right
                elif left[0] == "s":
                    stranded = right
                elif left[0] == "l":
                    library_type = right

        records.append((name, read_length, genome, paired, stranded, library_type))

    df = pd.DataFrame(records)
    df.columns = ["replicate", "read length", "genome", "paired", "stranded", "library type"]
    df.sort_values(by=["replicate"])
    return df


def gather_library_stats_fun(in_files, sample_names, out_file):
    df = gather_library_stats2df(in_files, sample_names)
    df.to_csv(out_file, index=False, sep="\t")
    