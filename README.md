# pyIPSA
Integrative Pipeline for Splicing Analysis

### Installation & Run

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

### Working folders

1. `input` - the folder where all BAMs are located
(1 sample - 1 BAM)
2. `output` - the folder where all results are saved.
Each sample has its own directory in this folder.
3. `genomes` - the folder with FASTA files of genomes
4. `annotations` - the folder with GTF/GFF files of annotations
5. `introns` - the folder with annotated introns