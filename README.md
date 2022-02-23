# gctree

Implements phylogenetic inference for data with repeated sequences, as described in:

DeWitt, Mesin, Victora, Minin and Matsen, **Using genotype abundance to improve phylogenetic inference**, [*Molecular Biology and Evolution*](https://matsengrp.github.io/gctree/cite.html).

**Note: full documentation of the gctree package is available at: https://matsengrp.github.io/gctree**

## Installation

### Linux/MacOS

#### Base package install

See the [docs install page](https://matsengrp.github.io/gctree/install.html) to install the gctree Python package and its command line interface.


## Pipeline quick start

### inference
The [*gctree quickstart docs*](https://matsengrp.github.io/gctree/quickstart.html) contain detailed instructions to run the gctree inference pipeline. For quick reference, gctree inference input and output files are described here:

- *input file*: `FASTA` or `PHYLIP` file containing a sequence for each observed individual/cell, and an additional sequence containing the ancestral genotype of all observed sequences (used for outgroup rooting).

- *description of inference output files*: After the inference pipeline has completed, the output directory will contain the following output files:
    - `<input file>.idmap`: text file mapping collapsed sequenced ids to cell ids from the original input file
    - `<input file>.counts`: text file mapping collapsed sequenced ids to their abundances
    - `<input file>.phylip`: phylip alignment file of collapsed sequences for computing parsimony trees
    - `dnapars/`: directory of parsimony tree output from PHYLIP's dnapars
    - `gctree.inference.*.svg`: rendered tree images for each of the parsimony trees
    - `gctree.inference.abundance_rank.pdf`: histogram of genotype abundances
    - `gctree.inference.likelihood_rank.pdf`: rank plot of gctree likelihoods for the parsimony trees
    - `gctree.inference.log`: log file containing parameter fits, numerical likelihood results, and any other program messages
    - `gctree.inference.parsimony_forest.p`: a python pickle file containing the parsimony trees as `CollapsedTree` objects


### simulation
```
gctree simulate -N <integer population size to simulate>
```
