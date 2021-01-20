============
Installation
============
pyIPSA requires `snakemake <https://snakemake.github.io/>`_.
You can install it with all dependencies via `conda <https://conda.io/projects/conda/en/latest/index.html>`_ package management system.
To create new conda environment and install snakemake with required packages simply run:

.. code-block:: console

    $ conda create -n ipsa
    $ conda activate ipsa
    $ conda install -c conda-forge -c bioconda snakemake pysam pandas

To install the workflow itself retrieve it from GitHub

.. code-block:: console

    $ git clone https://github.com/Leoberium/pyIPSA.git

Now you can get started.

If you don't want to use conda, see `here <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`_ for alternative ways of snakemake installation.

Dependencies:
    * snakemake
    * pysam
    * pandas
    * numpy

Installation on Arcuda
======================

For Arcuda cluster users the installation is simplier.
Just Load module with Python and install all required packages with these commands:

.. code-block:: console

    $ module load ScriptLang/python/3.8.3
    $ pip3 install --user --upgrade snakemake pandas pysam

Then retrieve the workflow from GitHub:

.. code-block:: console

    $ git clone https://github.com/Leoberium/pyIPSA.git