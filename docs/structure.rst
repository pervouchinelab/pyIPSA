Pipeline
=========

Input
-----

``input_dir`` in configuration file. Must contain at least one alignment file (BAM).

Outputs
-------

``output_dir`` in configuration file. If the run was fully successful, the output has the following structures:

* ``J1`` - **Step 1** gzipped files with junctions and their counts. One file per each sample.
* ``J2`` - **Step 2** gzipped files with aggregated junctions.
* ``J3`` - **Step 3** gzipped files with annotated junctions.
* ``J4`` - **Step 4** gzipped files with annotated junctions after choosing strand.
* ``J6`` - **Step 6** gzipped files with filtered junctions.
* ``stats`` - files with some statistics from **Step 1** for each sample.
* ``overview.tsv`` - aggregates some stats from all samples.
* ``S1`` - **Step 1** gzipped files with sites and their counts. ``PS1`` if junctions were pooled.
* ``S2`` - **Step 2** gzipped files with aggregated sites. ``PS2`` if junctions were pooled.
* ``S6`` - **Step 6** gzipped files with filtered sites. ``PS6`` if junctions were pooled.
* ``R`` - final step files with inclusion, exclusion and retention rates for splice sites from each sample.


Workflow
--------

The root directory of **pyIPSA** has several folders:

* ``config`` - contains configuration files in YAML format
* ``deprecated`` - contains obsolete scripts
* ``docs`` - documentation source
* ``known_SJ`` - has 2 files for each genome:

    * ``*.ranked.txt`` - splice site and its rate of usage in that genome
    * ``*.ss.tsv.gz`` - all introns from corresponding annotation

* ``workflow`` - the workflow itself, consists of:

    * ``rules`` - snakemake rules
    * ``scripts`` - python scripts which process the data
    * ``Snakefile`` - main snakemake file

After usage additional folders may appear:

* ``genomes`` - stores genome files used in analysis
* ``annotations`` - stores annotation files used in analysis