Installation
############

gctree
======
The preferred way to install gctree is with Conda:

.. code-block:: bash

  conda create -n gctree python=3.9
  conda activate gctree
  conda install -c conda-forge gctree

This method installs all dependencies, including the PHYLIP package (see below)

You may install using pip instead:

.. code-block:: bash

  pip install gctree

However, you will then need to separately install PHYLIP, as described below:


PHYLIP
======

The original use case for gctree is to use genotype abundance information to
rank degenerate maximum parsimony trees. For this, you will need Joe
Felsenstein's PHYLIP package
(https://evolution.genetics.washington.edu/phylip.html).
If you are working in a Conda environment, PHYLIP can be installed with

.. code-block:: bash

  conda install -c bioconda phylip
