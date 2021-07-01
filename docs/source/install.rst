Installation
############

gctree
======

.. code-block:: bash

  pip install gctree


.. todo::
  setup Conda install


PHYLIP
======

The original use case for gctree is to use genotype abundance information to
rank degenerate maximum parsimony trees. For this, you will need Joe
Felsenstein's PHYLIP package
(https://evolution.genetics.washington.edu/phylip.html).
If you are working in a Conda environment, PHYLIP can be installed with

.. code-block:: bash

  conda install -c bioconda phylip
