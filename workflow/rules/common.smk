import pandas as pd

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