Installation
============

**pyIPSA** workflow consists of several tools written in Python.
The workflow itself is based on the workflow management system
`snakemake <https://snakemake.github.io/>`_.

Download workflow
-----------------

To download pyIPSA just clone it from `github repository <https://github.com/Leoberium/pyIPSA>`_:

.. code-block:: bash

    $ git clone https://github.com/Leoberium/pyIPSA.git

But it won't work if you don't have all the required packages installed.

Requirements
------------

* Python >= 3.7
* numpy == 1.19.1
* pandas == 1.2.0
* pysam == 0.16.0.1
* snakemake == 5.32.0

Install requirements using ``pip``
----------------------------------

If you have Python environment with ``pip`` just move to directory with pyIPSA and run:

.. code-block:: bash

    $ pip install --upgrade -r requirements.txt

Otherwise you need to install Python and pip firstly.

Installing on Arcuda
----------------------

For Arcuda cluster users the installation is quite similar.
Move to pyIPSA directory, load module with Python and install all required packages:

.. code-block:: bash

    $ module load ScriptLang/python/3.8.3
    $ pip3 install --user --upgrade -r requirements.txt
