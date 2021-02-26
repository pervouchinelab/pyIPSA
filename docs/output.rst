Output
======

pyIPSA stores its output in several folders which store results
of consecutive steps of the pipeline. The structure of the output is the following:

::

    ├── aggregated_junction_stats.tsv
    ├── aggregated_library_stats.tsv
    ├── J1
    |   ├── <sample>.J1.gz
    |   └── <sample>.library_stats.txt
    ├── J2
    |   └── <sample>.J2.gz
    ├── J3
    |   └── <sample>.J3.gz
    ├── J4
    |   ├── <sample>.J4.gz
    |   └── <sample>.junction_stats.txt
    ├── J6
    |   └── <sample>.J6.gz
    ├── R
    |   └── <sample>.R.gz
    ├── S1
    |   └── <sample>.S1.gz
    ├── S2
    |   └── <sample>.S2.gz
    └── S6
    |   └── <sample>.S6.gz

All files are plain-text and tab-separated if they present table data.
Files with junctions and sites are also gzipped.

J1 - Junction and Counts
------------------------

pyIPSA's first step output is 2 files per each sample in J1 folder:
junction counts file and library statistics file. Junction counts file is a
plain-text, tab-separated file without a header. It is named ``<sample>.J1.gz``.
Example of its content:

+--------------+------------+----+----+----+----+
| junction id  |  offset    | F1 | R1 | F2 | R2 |
+==============+============+====+====+====+====+
| chr1_500_700 |    9       | 0  | 5  |  2 |  1 |
+--------------+------------+----+----+----+----+
| chr1_500_700 |    12      | 0  | 0  |  0 |  2 |
+--------------+------------+----+----+----+----+
| chr1_850_950 |    24      | 2  | 0  |  1 |  0 |
+--------------+------------+----+----+----+----+

Each row describes junction and counts from all reads supporting the junction
in given offset. The columns have the following interpretation:

* **junction id** ---
  reference sequence (chromosome), 1-based start position and 1-based end position joined by _

* **offset** ---
  distance from start of the read to the start of the junction

* **F1** ---
  count for all read 1 from the forward strand

* **R1** ---
  count for all read 1 from the reverse strand

* **F2** ---
  count for all read 2 from the forward strand

* **R2** ---
  count for all read 2 from the reverse strand

The second file contains various data about RNASeq library and alignment extracted from BAM.
It is plain-text and named ``<sample>.library_stats.txt``.

..
    TODO: write description of parameters

J2 - Junctions and Aggregated Counts
------------------------------------

The second step's output is one junctions with aggregated counts file
for each sample in J2 folder. It is named ``<sample>.J2.gz`` and has the following content:

..
    TODO: calculate values

+----------------+-------------+-----------------+---------+
| junction id    | total count | staggered count | entropy |
+================+=============+=================+=========+
| chr1_500_700_+ |    10       |        2        |  1.52   |
+----------------+-------------+-----------------+---------+
| chr1_850_950_+ |    3        |        1        |  1.22   |
+----------------+-------------+-----------------+---------+

Each row describes junction and counts from all reads supporting the junction
aggregated by all possible offsets. The columns have the following interpretation:

* **junction id** ---
  reference sequence (chromosome), 1-based start position, 1-based end position and now also strand joined by _

* **total count** ---
  total count of this junction from all reads in all offsets

* **staggered count** ---
  number of offsets that have read with count for this junction

* **entropy** ---
  entropy for offset distribution

J3 - Annotated Junctions
------------------------

The third step's output is annotated junctions file for each sample in J3 folder.
It is named ``<sample>.J3.gz`` and has the following content:

+----------------+-------------+-----------------+---------+-------------------+-------------+
| junction id    | total count | staggered count | entropy | annotation status | splice site |
+================+=============+=================+=========+===================+=============+
| chr1_500_700_+ |    10       |        2        |  1.52   |         3         |    GTAG     |
+----------------+-------------+-----------------+---------+-------------------+-------------+
| chr1_850_950_+ |    3        |        1        |  1.22   |         0         |    GGAG     |
+----------------+-------------+-----------------+---------+-------------------+-------------+

First 4 columns are the same as in J2 file, additional ones are:

* **annotation status** ---
  possible values are:

    * **0** - both junction ends are absent in annotation
    * **1** - one junction end is present in annotation
    * **2** - both junction ends are present in annotation
    * **3** - both junction ends are present in annotation and correspond to existing intron

* **splice site** ---
  donor and acceptor splice site sequences. Usually they are GT and AG but can vary.

J4 - Choose Strand
------------------

The fourth step's output is 2 files for each sample:
annotated junctions file with correct strand chosen for each junctions and
junctions stats file. The first file is named ``<sample>.J4.gz`` and has
the same format as file from J3 step:

+----------------+-------------+-----------------+---------+-------------------+-------------+
| junction id    | total count | staggered count | entropy | annotation status | splice site |
+================+=============+=================+=========+===================+=============+
| chr1_500_700_+ |    10       |        2        |  1.52   |         3         |    GTAG     |
+----------------+-------------+-----------------+---------+-------------------+-------------+
| chr1_850_950_+ |    3        |        1        |  1.22   |         0         |    GGAG     |
+----------------+-------------+-----------------+---------+-------------------+-------------+

But some records from J3 will be missed due to strand choice.
The second file is ``<sample>.junctions_stats.txt``, it just reports
how many read counts support junctions with GTAG or non-GTAG splice sites,
annotated and non-annotated junctions and etc.

J6 - Filter Junctions
---------------------

This step's output is ``<sample>.J6.gz``. It has the same format as files from
J3 and J4 steps:

+----------------+-------------+-----------------+---------+-------------------+-------------+
| junction id    | total count | staggered count | entropy | annotation status | splice site |
+================+=============+=================+=========+===================+=============+
| chr1_500_700_+ |    10       |        2        |  1.52   |         3         |    GTAG     |
+----------------+-------------+-----------------+---------+-------------------+-------------+
| chr1_850_950_+ |    3        |        1        |  1.22   |         0         |    GGAG     |
+----------------+-------------+-----------------+---------+-------------------+-------------+

..
    TODO: update table

The purpose of step is to filter out junctions not passing some of criteria:

* **total count** ---
  must be not less than ``total_count`` value in config file (default is ``1``)

* **entropy** ---
  must be not less than ``entropy`` value in config file (default is ``1.5``)

* **splice site** ---
  must be GTAG only if config file's parameter ``gtag`` set to ``True``
  (default is ``False``)

Gather Junctions
----------------

Some additional files are generated along with J1-J6:

* ``aggregated_library_stats.tsv`` ---
  library parameters (from J1) of all samples present in one table

* ``aggregated_junction_stats.tsv`` ---
  junction stats (from J4) of all samples present in one table

* ``J4/merged_junctions.J4.gz`` ---
  a union of all junctions files from J4 step. Contains all unique junctions
  found in alignments. Computed only if parameter ``pooled`` in config set to ``True``.

S1 (PS1) - Sites and Counts
---------------------------

This step is similar to J1, but now it works with sites, not junctions.
The output file name is ``<sample>.S1.gz``.
``<sample>.PS1.gz`` if junctions were pooled. The format is:

+------------+------------+-------+
|   site id  |   offset   | count |
+============+============+=======+
| chr1_500_+ |    9       |    2  |
+------------+------------+-------+
| chr1_700_+ |    12      |    9  |
+------------+------------+-------+
| chr1_900_+ |    5       |    4  |
+------------+------------+-------+

* **site id** ---
  reference sequence (chromosome), site's 1-based position and strand joined by _

* **offset** ---
  distance from start of the read to the position of site

* **count** ---
  total count of read supporting the site

S2 (PS2) - Aggregate Sites
--------------------------

This step is similar to J2, it aggregates offsets for each junction.
The output file name is ``<sample>.S2.gz`` or
``<sample>.PS2.gz`` if junctions were pooled. The format is:

+------------+-------------+-----------------+---------+
| site id    | total count | staggered count | entropy |
+============+=============+=================+=========+
| chr1_500_+ |      9      |        1        |  0.96   |
+------------+-------------+-----------------+---------+
| chr1_700_+ |      12     |        1        |  1.55   |
+------------+-------------+-----------------+---------+

New columns are:

* **total count** ---
  total count of this site from all reads in all offsets

* **staggered count** ---
  number of offsets that have read with count for this site

* **entropy** ---
  entropy for offset distribution

S6 (PS6) - Filter Sites
-----------------------

The output file name is ``<sample>.S6.gz`` or
``<sample>.PS6.gz`` if junctions were pooled.
The format is the same as S2:

+------------+-------------+-----------------+---------+
| site id    | total count | staggered count | entropy |
+============+=============+=================+=========+
| chr1_500_+ |      9      |        1        |  0.96   |
+------------+-------------+-----------------+---------+
| chr1_700_+ |      12     |        1        |  1.55   |
+------------+-------------+-----------------+---------+


The purpose of this step is to filter out sites not passing some of criteria:

* **total count** ---
  must be not less than ``total_count`` value in config file (default is ``1``)

* **entropy** ---
  must be not less than ``entropy`` value in config file (default is ``1.5``)

R (PR) - Rates
--------------

The output file name is ``<sample>.R.gz`` or ``<sample>.PR.gz``
if junctions were pooled.
The format:

+--------------+-------------+-----------------+-----------+
| site id      | inclusion   | exclusion       | retention |
+==============+=============+=================+===========+
| chr1_500_+_D |      5      |        0        |    6      |
+--------------+-------------+-----------------+-----------+
| chr1_700_+_A |      11     |        0        |    13     |
+--------------+-------------+-----------------+-----------+

* **site id** ---
  reference sequence (chromosome), site's 1-based, and type (D - donor, A - acceptor)
  and strand joined by _

* **inclusion** ---
  number of reads supporting inclusion of this site

* **exclusion** ---
  number of reads supporting exclusion of this site

* **retention** ---
  number of reads support retention of this site
