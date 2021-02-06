pyIPSA Workflow
===============

pyIPSA (Integrative Pipeline for Splicing Analysis) is a workflow and a set of tools for splicing analysis.
It extracts and processes local splicing estimates.

As input, **pyIPSA** takes a set of BAM files. **pyIPSA** extracts split reads and continuous reads from BAMs and processes
them in several steps.

Step 1: count
-------------

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
-----------------

Aggregate all offsets for each junction/site.
The input of this step is ``<sample>.J1.gz`` file from **Step1**.
The output is a gzipped tab-delimited file named by sample ``<sample>.J2.gz``.
It is saved in J2 directory and has the following content::

    chr22_16151821_16162397_+    10    8    2.92
    chr22_16151821_16162397_-    10    8    2.92
    chr22_16159389_16162397_+    2     2    1.0


The columns are:

* ``junction_id`` - reference sequence (chromosome), 1-based start position, 1-based end position and strand joined by _
* ``total count`` - total count of this junction from all reads in all offsets
* ``staggered count`` - number of offsets that have read with count for this junction
* ``entropy`` - entropy for offset distribution

Step 3: annotate
----------------

Annotate all junctions: add annotation status and splice site nucleotides.
It is saved in J3 directory and has the following content::

    chr22_16162487_16186811_-    9    9    3.17    3    GTAG
    chr22_16186946_16187165_+    2    2    1.0     0    CTAC


The columns are:

* ``junction_id`` - reference sequence (chromosome), 1-based start position, 1-based end position and strand joined by _
* ``total count`` - total count of this junction from all reads in all offsets
* ``staggered count`` - number of offsets that have read with count for this junction
* ``entropy`` - entropy for offset distribution
* ``annotation status`` - possible values are:

    * ``0`` - both junction ends are absent in annotation
    * ``1`` - one junction end is present in annotation
    * ``2`` - both junction ends are present in annotation
    * ``3`` - both junction ends are present in annotation and correspond to existing intron

* ``splice site`` - donor and acceptor splice site sequences. Usually they are GT and AG but can vary.

Step 4: choose strand
---------------------

Choose correct strand for each junction.

Step 6: filter
--------------

Filter junctions or sites using 3 different criteria:

* total count - for both junctions/sites (default threshold = 1)
* allowed splice sites - all possible splice sites in junctions or GT/AG only (default - all possible sites)
* entropy - entropy for junctions/sites value must be not less than threshold (default = 1.5)

The output is a gzipped tab-delimited file named ``<sample>.J6.gz`` in case of junctions or
``<sample>.(P)S6.gz`` in case of sites.
