Running
=======

In order to run **pyIPSA** workflow you have to specify input and output folders in the configuration file.
Open file ``config/config.yaml`` in **pyIPSA** directory and specify your desired input and output folders.
Paths must be absolute or relative to **pyIPSA** directory. Input folder must have at least one alignment file (BAM).

To run **pyIPSA** use the following command while in root directory:

.. code-block:: bash

    $ snakemake --cores <number of cores>

To run in cluster environment using Grid Engine:

.. code-block:: bash

    $ snakemake --cluster qsub --j <number of jobs>

For other running options consult with
`snakemake docs <https://snakemake.readthedocs.io/en/stable/executing/cli.html>`_.

Configuration
-------------

``config/config.yaml`` has many other useful options:

* ``pooled`` - if ``True``, merge junctions from all samples before retrieving sites if True
* ``primary`` - if ``True``, use only primary alignment for multimapped reads
* ``unique`` - if ``True``, do not use multimapped reads
* ``threads`` - number of threads used to read single alignment file
* ``strand_mode`` - to remove in future
* ``min_offset`` - minimal offset when aggregating junctions
* ``min_intron_length`` - minimal allowed length of junction
* ``max_intron_length`` - maximal allowed length of junction
* ``entropy`` - minimal value of entropy used for filtering out junctions or sites
* ``total_count`` - minimal allowed count while filtering out junctions
* ``gtag`` - if ``True``, use only junctions with GT/AG splice sites
* ``genome_filenames`` - stores full names of genomes
* ``genome_urls`` - stores URLs to genome files
* ``annotation_urls`` - stores URLs to annotation files
