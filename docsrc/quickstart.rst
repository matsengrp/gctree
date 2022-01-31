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
The flag ``--id_abundances`` can be used to indicate that integer sequence ids should be interepreted as abundances.
The argument ``--abundance_file`` indicates that sequence abundance information should be written to the specified ``csv`` file.
The argument ``--idmapfile`` allows us to specify a filename for the output
file containing a mapping of new, unique sequence IDs to original sequence IDs
from the fasta file.

.. command-output:: deduplicate ../example/150228_Clone_3-8.fasta --root GL --abundance_file abundances.csv --idmapfile idmap.txt > deduplicated.phylip
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
  :ellipsis: 0

You now have two new files, ``outfile`` and ``outtree``.
Note: if you want to rerun the above ``dnapars`` command, you must delete these two files first!


gctree
======

We're now ready to run ``gctree infer`` to use abundance data (in ``abundances.csv``) to rank the equally parsimonious trees (in ``outfile``).
We can use the optional argument ``--frame`` to indicate the coding frame of the sequence start, so that amino acid substitutions can be annotated on our trees.

If working in a headless environment, ``gctree infer`` must be run with a tool
like ``xvfb-run`` to provide an X server for rendering the output trees.

.. command-output:: gctree infer outfile abundances.csv --root GL --frame 1 --verbose | tee gctree.inference.log
  :shell:
  :ellipsis: 10

A large number of output files with the basename ``gctree.out.*`` are also created.
The SVG image file ``gctree.out.inference.abundance_rank.svg`` shows a distribution of genotype abundances in the original data:

.. image:: gctree.out.inference.abundance_rank.svg
  :width: 600

Then there are files ``gctree.out.inference.1.svg`` and ``gctree.out.inference.1.nk`` containing an SVG tree image and newick tree file. If more than one parsimony tree is optimal, then up to ten optimal trees will be sampled randomly, and the corresponding files will be numbered **arbitrarily**.
For example here is the top ranked tree ``gctree.out.inference.1.svg``:

.. image:: gctree.out.inference.1.svg
  :width: 1000

You will also see Python pickle files ``gctree.out.inference.[1,2,...].p`` containing a :obj:`gctree.CollapsedTree` object for each optimal tree, which can be loaded and manipulated via the API (e.g. plotted in various ways using :meth:`gctree.CollapsedTree.render`).

Criteria other than branching process likelihoods can be used to break ties
between trees. Providing arguments ``--isotype_mapfile`` and
``--idmapfile`` will allow trees to be ranked by isotype parsimony. Providing
arguments ``--mutability`` and ``--substitution`` allows trees to be ranked
according to a context-sensitive mutation model. By default, trees are ranked
lexicographically, first maximizing likelihood, then minimizing isotype
parsimony and mutabilities, if such information is provided.
Ranking priorities can be adjusted using the argument ``--priority_weights``.

All parsimony trees found by dnapars, as well as branching process parameters
are saved in the file ``gctree.out.inference.parsimony_forest.p``, containing a :class:`gctree.CollapsedForest` object.
This file may be manipulated using ``gctree infer``. For example, to find the optimal tree
according to a linear combination of likelihood, isotype parsimony,
mutabilities, and alleles:

.. command-output:: gctree infer gctree.out.serialized_dag.p --frame 1 --idmap idmap.txt --isotype_mapfile ../example/isotypemap.txt --mutability ../S5F/Mutability.csv --substitution ../S5F/Substitution.csv --priority_weights 2 2 1 0 --outbase newranking --verbose
   :shell:

By default, only the files listed above will be generated, with the optional argument ``--outbase`` specifying how the output files should be named. For a summary of the collection of trees used for ranking, the argument ``--summarize_forest`` is provided. For detailed information about each tree used for ranking, use the argument ``--tree_stats``.

.. image:: newranking.inference.1.svg
   :width: 1000


isotype
=======

If we would like to add observed isotype data to trees output by gctree
inference, we can now do so.
In addition to the outputs from gctree, a file mapping original IDs of observed
sequences to their observed isotypes (like ``example/isotypemap.txt``) is required.

.. command-output:: isotype gctree.inference.log idmap.txt ../example/isotypemap.txt --trees gctree.out.inference.1.p newranking.inference.1.p --out_directory isotyped
  :shell:
  :ellipsis: 10

Trees originally output by gctree are re-rendered with revised labels and node
colors corresponding to observed or inferred isotypes.
For example, here is the top ranked tree above, with isotypes added:

.. image:: isotyped/gctree.out.1.isotype_parsimony.28.svg
  :width: 1000
