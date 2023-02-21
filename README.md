# pyIPSA
Integrative Pipeline for Splicing Analysis

[![Documentation Status](https://readthedocs.org/projects/pyipsa/badge/?version=latest)](https://pyipsa.readthedocs.io/en/latest/?badge=latest)


## Installation & Run

### Step 1: Obtain a copy of this workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system, into the place where you want to perform the data analysis.

    git clone https://github.com/pervouchinelab/pyIPSA.git
    cd pyIPSA

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    # install mamba package manager if you don't have it
    conda install -n base -c conda-forge mamba
    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

#### For Arkuda users

Install a profile for Arkuda (adapted from [Snakemake-Profiles/generic](https://github.com/Snakemake-Profiles/generic)):

    mkdir -p ~/.config/snakemake
    cp -R workflow/profiles/arkuda ~/.config/snakemake
    chmod +x ~/.config/snakemake/arkuda/*

Arkuda creates conda environments without execute permissions by default, 
so you have to create the environments and add the permissions manually. 

    snakemake --use-conda -c1 --conda-create-envs-only
    chmod +x .snakemake/conda/**/bin/*

Execute the workflow:

    snakemake --use-conda --default-resources 'mem_mb=5000' 'time_min=60' --profile arkuda  -j<number of jobs>


## Working folders

Folder in repository:
1. `config` - the folder with config file, where you set up your pipeline
2. `known_SJ` - the folder with annotated splice junctions
3. `workflow` - the folder with working scripts of the pipeline

Additional directories created
1. `genomes` - the folder which stores all downloaded genomes
2. `annotations` - the folder which stores all downloaded annotations
