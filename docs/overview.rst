===============
pyIPSA Workflow
===============

IPSA (Integrative Pipeline for Splicing Analysis) is a workflow and a set of tools for splicing analysis.
It extracts and processes local splicing estimates.

As input, IPSA takes a set of BAM files. IPSA extracts split reads and continuous reads from BAMs and processes them in
a series of steps.

1. Extract split reads with offset. The input of this step is a BAM file.
The output is a tab-delimited file that is saved in J1 directory and has the following columns:
    # splice junction id