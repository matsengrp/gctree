Installation
############

Anaconda install
================
The preferred way to install gctree is with Conda:

.. code-block:: bash

  conda create -n gctree python=3.10
  conda activate gctree
  conda install -c conda-forge gctree

You will likely also want to install PHYLIP (see below)

Pip install
===========

You may install using pip instead:

.. code-block:: bash

  pip install gctree

However, you will then need to separately install PHYLIP


Docker build
============

You may also find the provided Dockerfile useful.
The resulting image runs the test script `tests/test.sh` by default.
This runs the default inference pipeline on the sample data provided with gctree.

.. code-block:: bash

   git clone git@github.com:matsengrp/gctree.git
   docker build gctree -t gctree

To run the test script and verify the image was built successfully:

.. code-block:: bash

   docker run -t gctree

Or use the image interactively with, for example:

.. code-block:: bash

   docker run -i -t gctree bash


PHYLIP Installation
===================

The original use case for gctree is to use genotype abundance information to
rank degenerate maximum parsimony trees. For this, you will need Joe
Felsenstein's PHYLIP package
(https://evolution.genetics.washington.edu/phylip.html).
If you are working in a Conda environment, PHYLIP can be installed with

.. code-block:: bash

  conda install -c bioconda phylip

