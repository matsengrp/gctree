Quickstart
##########

In this minimal example we show how to use gctree commands to rank parsimony trees from PHYLIP's dnapars program.
See the CLI and API docs for complete documentation of command options, included many not used here.

The data in this example were published in `Tas et al. 2016. Visualizing Antibody Affinity Maturation in Germinal Centers. Science 351 (6277) <http://science.sciencemag.org/content/351/6277/1048>`_ and shown in Fig. 4 (lymph node 2, germinal center 1).

.. image:: _static/gc1.png
  :width: 200


Input data
==========

We start with a fasta alignment file.
The file ``../example/150228_Clone_3-8.fasta`` contains heavy chain V gene sequences from
germinal B cells sorted from a brainbow mouse using multicolor fate mapping.

.. command-output:: tail -30 ../example/150228_Clone_3-8.fasta

In this file the sequence with id ``GL`` is the naive germline sequence, and represents the root of the tree.
It does not refer to an observed sequence, but is included to outgroup root the tree!
Generally this input fasta file should also contain all observed sequences with duplication (some of which may be identical to the root/outgroup sequence).
It is also possible to provide some observed sequences with integers ids to be interpreted as abundances, as in the sequence named ``17`` above.


Deduplication and sequence abundances
=====================================

First we deduplicate the sequences and convert to phylip alignment format, and also determine sequence abundances.
The ``deduplicate`` command writes the phylip alignment of unique sequences to stdout (which we redirect to a file here).
The argument ``--root`` indicates the root id.
The flag ``--id_abundances`` indicates that integer sequence ids should be interepreted as abundances (our ``17`` sequence).
The argument ``--abundance_file`` indicates that sequence abundance information should be written to the specified ``csv`` file.

.. command-output:: deduplicate ../example/150228_Clone_3-8.fasta --root GL --id_abundances --abundance_file abundances.csv > deduplicated.phylip
  :shell:

We now have files ``deduplicated.phylip`` and ``abundances.csv``:

.. command-output:: head deduplicated.phylip
  :shell:

.. command-output:: head abundances.csv
  :shell:


Parsimony trees
===============

We use PHYLIP's ``dnapars`` program (see install docs page) to generate a set of maximum parsimony trees for the deduplicated sequences.
PHYLIP is an interactive command line tool, so we will automatically generate a config file to feed it input, instead of using it interactively.

We generate a config file for ``dnapars`` based on our deduplicated alignment using the ``mkconfig`` command.

.. command-output:: mkconfig deduplicated.phylip dnapars > dnapars.cfg
  :shell:

Run ``dnapars`` using this config file, redirecting stdout to a log file.

.. command-output:: dnapars < dnapars.cfg > dnapars.log
  :shell:

You now have two new files, ``outfile`` and ``outtree``.
Note: if you want to rerun the above ``dnapars`` command, you must delete these two files first!


gctree
======

We're now ready to run ``gctree infer`` to use abundance data (in ``abundances.csv``) to rank the equally parsimonious trees (in ``outfile``).
We can use the optional argument ``--frame`` to indicate the coding frame of the sequence start, so that amino acid substitutions can be annotated on our trees.

.. command-output:: gctree infer outfile abundances.csv --root GL --frame 1
  :shell:
  :ellipsis: 10

A large number of output files with the basename ``gctree.out.*`` are also created.
The SVG image file ``gctree.out.inference.abundance_rank.svg`` shows a distribution of genotype abundances in the original data:

.. image:: gctree.out.inference.abundance_rank.svg
  :width: 600

The SVG image file ``gctree.out.inference.likelihood_rank.svg`` is a rank plot of likelihood over the set of maximum parsimony trees:

.. image:: gctree.out.inference.likelihood_rank.svg
  :width: 600

Then there are files ``gctree.out.inference.[1,2,...].svg`` and ``gctree.out.inference.[1,2,...].nk`` containing an SVG tree image and newick tree file for each of the parsimony trees (numbered according to likelihood rank).
For example here is the top ranked tree ``gctree.out.inference.1.svg``:

.. image:: gctree.out.inference.1.svg
  :width: 1000

You will also see Python pickle files ``gctree.out.inference.[1,2,...].p`` containing a :obj:`gctree.CollapsedTree` object for each tree, which can be loaded and manipulated via the API (e.g. plotted in various ways using :meth:`gctree.CollapsedTree.render`).
