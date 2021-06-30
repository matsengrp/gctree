Command line interface

gctree
======

The primary gctree command.

.. argparse::
    :module: gctree.cli
    :func: get_parser
    :prog: gctree
    :nodefault:


Parsimony utilities
===================

Additional command line uilities for processing data for genotype collapse,
inferring parsimony trees with phylip, and prepping parsimony trees for analysis
with gctree.

deduplicate
-----------

.. argparse::
    :module: gctree.deduplicate
    :func: get_parser
    :prog: deduplicate
    :nodefault:

mkconfig
--------

.. argparse::
    :module: gctree.mkconfig
    :func: get_parser
    :prog: mkconfig
    :nodefault:

phylip_parse
------------

.. argparse::
    :module: gctree.phylip_parse
    :func: get_parser
    :prog: phylip_parse
    :nodefault:
