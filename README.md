# pyIPSA
Integrative Pipeline for Splicing Analysis

To use this program you must have python environment
with the following programs and libraries installed:
- `snakemake`
- `pysam`
- `pandas`
- `numpy`

To install:  
`git clone https://github.com/Leoberium/pyIPSA.git`

To run:
1. Python environment with required libraries
must be active
2. Make directory with Snakefile active
3. Run `snakemake` command

To test:  
`snakemake test`

Options for `snakemake` command are available in 
[snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

Such environment may be created via conda:  
1. `conda create -n ipsa`  
2. `conda activate ipsa`
3. `conda install -c conda-forge -c bioconda
snakemake pysam pandas`

