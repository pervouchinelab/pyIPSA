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
3. `conda install -c conda-forge -c bioconda snakemake pysam pandas`

To install the pipeline:  
`git clone https://github.com/Leoberium/pyIPSA.git`

To run (test run if samples folder is empty):
1. Python environment with required libraries
must be active
2. Make directory with Snakefile active
3. Run `snakemake` command

Options for `snakemake` command are available in 
[snakemake documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

#### For Arcuda users

Just load module with Python and install libraries:
1. `module load ScriptLang/python/3.8.3`
2. `pip3 install --user --upgrade snakemake pandas pysam`

After that you can run pipeline using cluster engine:
`snakemake --cluster qsub --j <number of jobs>`

### Working folders

Folder in repository:
1. `config` - the folder with config file, where you set up your pipeline
2. `deprecated` - the folder with old scripts not used in workflow
3. `known_SJ` - the folder with annotated splice junctions
4. `workflow` - the folder with working scripts of the pipeline

Additional directories created
1. `genomes` - the folder which stores all downloaded genomes
2. `annotations` - the folder which stores all downloaded annotations
