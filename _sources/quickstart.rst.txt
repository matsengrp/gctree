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
This input fasta file should also contain all observed sequences with duplication (some of which may be identical to the root/outgroup sequence).

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

You now have two new files, ``outfile`` and ``outtree``.
Note: if you want to rerun the above ``dnapars`` command, you must delete these two files first!


gctree
======

We're now ready to run ``gctree infer`` to use abundance data (in ``abundances.csv``) to rank the equally parsimonious trees (in ``outfile``).
We can use the optional argument ``--frame`` to indicate the coding frame of the sequence start, so that amino acid substitutions can be annotated on our trees.

.. note::
  If working in a headless environment, ``gctree infer`` must be run with a tool
  like XVFB to provide an X server for rendering the output trees.
  Prepend the gctree command with ``xvfb-run -a``.
  Alternatively, we have had success setting the following environment variables instead of using XVFB:
  ``export QT_QPA_PLATFORM=offscreen``,
  ``export XDG_RUNTIME_DIR=/tmp/runtime-runner``.
  You may also want to tell matplotlib to use a headless backend with ``export MPLBACKEND=agg``.

.. command-output:: gctree infer outfile abundances.csv --root GL --frame 1 --verbose
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

All parsimony trees found by dnapars, as well as branching process parameters
are saved in the file ``gctree.out.inference.parsimony_forest.p``, containing a :class:`gctree.CollapsedForest` object.
This file may be manipulated using ``gctree infer``, instead of providing
a dnapars ``outfile``.

.. note::
  Although described below, using context likelihood, mutability parsimony, or isotype parsimony
   as ranking criteria is experimental, and has not yet been shown in a careful
   validation to improve tree inference. Only the default branching process
   likelihood is recommended for tree ranking!

Criteria other than branching process likelihoods can be used to break ties
between trees. Providing arguments ``--isotype_mapfile`` and
``--idmapfile`` will allow trees to be ranked by isotype parsimony. Providing
arguments ``--mutability`` and ``--substitution`` allows trees to be ranked
according to a context-sensitive mutation model. By default, trees are ranked
lexicographically, first maximizing likelihood, then minimizing isotype
parsimony, and finally maximizing a context-based poisson likelihood, if such information is provided.
Ranking priorities can be adjusted using the argument ``--ranking_coeffs``.

For example, to find the optimal tree
according to a linear combination of likelihood, isotype parsimony,
mutabilities, and alleles:

.. command-output:: gctree infer gctree.out.inference.parsimony_forest.p --frame 1 --idmap idmap.txt --isotype_mapfile ../example/isotypemap.txt --mutability ../HS5F_Mutability.csv --substitution ../HS5F_Substitution.csv --ranking_coeffs 1 0.1 0 --outbase newranking --summarize_forest --tree_stats --verbose
   :shell:

The files ``HS5F_Mutability.csv`` and ``HS5F_Substitution.csv`` are a context
sensitive mutation model which can be downloaded from the `Shazam Project <https://bitbucket.org/kleinstein/shazam/src/master/data-raw/`>_.

By default, only the files listed above will be generated, with the optional argument ``--outbase`` specifying how the output files should be named.

.. image:: newranking.inference.1.svg
   :width: 1000

For detailed information about each tree used for ranking, as well as a pairplot like the one below comparing the highest ranked tree to all other ranked trees,use the argument ``--tree_stats``.

.. image:: newranking.tree_stats.pairplot.png
   :width: 1000

Sometimes ranked trees are too numerous, and generating the output of ``--tree_stats`` would require too many resources. For a summary of the collection of trees used for ranking, the argument ``--summarize_forest`` is provided. Most importantly, this option summarizes how much less likely the top ranked tree is, compared to the most likely tree being ranked, for example to validate coefficients passed to ``--ranking_coeffs``.

.. command-output:: cat newranking.forest_summary.log
   :shell:


isotype
=======

If we would like to add observed isotype data to trees output by gctree
inference, we can now do so.
In addition to the outputs from gctree, a file mapping original IDs of observed
sequences to their observed isotypes (like ``example/isotypemap.txt``) is required.

.. command-output:: isotype idmap.txt ../example/isotypemap.txt --trees gctree.out.inference.1.p newranking.inference.1.p --out_directory isotyped
  :shell:
  :ellipsis: 10

Trees originally output by gctree are re-rendered with revised labels and node
colors corresponding to observed or inferred isotypes.
For example, here is the top ranked tree above, with isotypes added:

.. image:: isotyped/gctree.out.inference.1.p.isotype_parsimony.28.svg
  :width: 1000

A note about node names
=======================

The names associated with unobserved nodes (for example, in trees rendered with
`--idlabel`) are arbitrarily chosen, but are guaranteed to correspond
bijectively with unobserved sequences. However, if gctree output contains
multiple trees, **two unobserved nodes which share the same name but occur in
different output trees will not in general possess the same unobserved sequence.**

A note about ambiguous sequence data
====================================

Gctree can handle input sequences which contain ambiguous bases, but handling
of such sequences is experimental. That is, there may be issues with both the
implementation and the underlying methods.

If you use gctree to build trees from ambiguous sequence data, it's probably a good idea
to follow at least these guidelines:

* Make sure that all sites are unambiguous in at least one of the input
  sequences. If a site contains an ambiguous base in all input sequences,
  gctree will choose an arbitrary base for that site. The bases chosen for such
  sites should not be interpreted as meaningful in any way.
* Understand that for any observed sequence disambiguated by gctree, another choice of
  disambiguation for that sequence may exist, which results in a better tree
  with respect to the ranking criteria. The only guarantee is that
  disambiguations of observed sequences are maximally parsimonious given the
  tree topology in which they appear.

Gctree can handle ambiguous input sequences because ``dnapars`` can accept
ambiguous input sequences. Each tree output by ``dnapars`` is then
disambiguated. It is possible that multiple observed sequences may be
disambiguated in identical ways, in which case their corresponding leaf nodes,
and abundances, are merged. The deduplicated sequence ids corresponding to each
node in the final tree output by gctree are retained in the ``original_ids`` node
attribute.

Here's a bit more discussion about how much the disambiguated observed
sequences can be trusted:

**Why not to trust leaf disambiguation:**

As mentioned above, if the same site(s) contain ambiguous bases in all sequences, the disambiguation is completely arbitrary at those sites, but could be mistakenly interpreted as informed.
Also alluded to above, if multiple possible disambiguated leaf sequences exist for a particular dnapars topology, one will be chosen arbitrarily, even though another may be more plausible IRL, or with respect to the ranking criteria that gctree uses.

**Why to trust leaf disambiguation:**

Disambiguation is informed by sequence placement in the tree (by dnapars), which takes into account all that is known about the sequences.
Disambiguation minimizes parsimony score, meaning that
ambiguous sites will be filled in using bases in the most closely related
sequence in which those sites are unambiguous.

Also, different trees, with possibly different disambiguated leaf sequences, are ranked according to all data (isotype, abundance, mutability) provided to gctree, and the disambiguated leaf sequences influence that ranking. This means that ranking may influence the chosen disambiguation of observed sequences to be more plausible.
