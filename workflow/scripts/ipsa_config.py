"""This modules contains some globally available variables"""

BASE = 1

AVAILABLE_GENOMES = [
    "dm3", "dm6", "mm9", "mm10", "hg19", "hg38"
]

# This set contains only relevant reference names.
# It is used to filter out redundant names
# which are sometimes present in alignment files
ALLOWED_REFERENCE_NAMES = {
    # full reference names from dm3 annotation
    "chr2L", "chr2LHet", "chr2R", "chr2RHet", "chr3L", "chr3LHet", "chr3R", "chr3RHet",
    "chr4", "chrM", "chrU", "chrUextra", "chrX", "chrXHet", "chrYHet",
    # full reference names from dm6 annotation
    "chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrM", "chrX", "chrY",
    # full reference names from mm9 and mm10 annotations
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
    "chr17", "chr18", "chr19", "chrM", "chrX", "chrY",
    # full reference names from hg19 and hg38 annotations
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
    "chr17", "chr18", "chr19",  "chr20", "chr21", "chr22", "chrM", "chrX", "chrY"
}

# This dictionary is for translation of incomplete reference names to complete ones.
# Incomplete reference names sometimes occur in alignment files.
# Examples:
# 1 -> chr1
# MT -> chrM
TRANSLATION = {
    ref_name.lower()[3:]: ref_name for ref_name in sorted(ALLOWED_REFERENCE_NAMES)
}
TRANSLATION["mt"] = "chrM"


if __name__ == '__main__':
    print("ALLOWED_REFERENCE_NAMES =", ALLOWED_REFERENCE_NAMES)
    print("TRANSLATION =", TRANSLATION)
