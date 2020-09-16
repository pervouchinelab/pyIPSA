# pyIPSA
Integrative Pipeline for Splicing Analysis

To use this program you must have python environment
with the following programs and libraries installed:
- `snakemake`
- `pysam`
- `pandas`
- `numpy`

Such environment may be created via conda:  
1. `conda create -n ipsa`  
2. `conda activate ipsa`
3. `conda install -c conda-forge -c bioconda
snakemake pysam pandas`

To install the program:  
`git clone https://github.com/Leoberium/pyIPSA.git`

To run (test run if samples folder is empty):
1. Python environment with required libraries
must be active
2. Make directory with Snakefile active
3. Run `snakemake` command

Options for `snakemake` command are available in 
[snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html).
