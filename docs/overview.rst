===============
pyIPSA Workflow
===============

IPSA (Integrative Pipeline for Splicing Analysis) is a workflow and a set of tools for splicing analysis.
It extracts and processes local splicing estimates.

As input, IPSA takes a set of BAM files. IPSA extracts split reads and continuous reads from BAMs and processes
them in several steps.

Step 1: count
=============

Extract split reads with offset. The input of this step is a BAM file.
The output is a gzipped tab-delimited file named by sample ``<sample>.J1.gz``.
It is saved in J1 directory and has the following content::
  chr22_17629450_17630432    63    0    5    2    1
  chr22_17629450_17630432    64    0    0    0    2
  chr22_17629450_17630432    68    2    0    1    0
  chr22_17629450_17630432    69    2    0    4    4

The columns are:
    * ``junction_id`` - reference sequence (chromosome), 1-based start position and 1-based end position joined by _
    * ``offset`` - distance from start of the read to the start of the junctions
    * ``F1`` - read count for the first read on the forward strand
    * ``R1`` - read count for the first read on the reverse strand
    * ``F2`` - read count for the second read on the forward strand
    * ``R2`` - read count for the second read on the reverse strand

Step 2: aggregate
=================