# gctree

Implements phylogenetic inference for data with repeated sequences, as described in:

DeWitt, Mesin, Victora, Minin and Matsen, **Using genotype abundance to improve phylogenetic inference**, [*Molecular Biology and Evolution*](https://matsengrp.github.io/gctree/cite.html).

**Note: full documentation of the gctree package is available at: https://matsengrp.github.io/gctree**

This README provides info on use of SCons pipelines that wrap the base gctree package. **Note: SCons pipelines are deprecated. They do not support more recent functionality, and may be removed in a future release.**

## Installation

### Linux/MacOS

#### Base package install

See the [docs install page](https://matsengrp.github.io/gctree/install.html) to install the gctree Python package and its command line interface. Skip this if you are using the SCons pipelines (see below).


#### SCons pipeline install

`SCons` pipelines can be used to for end-to-end phylogenetic inference from sequence data.
These must be run from the repo directory after cloning.

0. For installing dependencies, [conda](https://conda.io/docs/) environment management is recommended. First install conda or miniconda.
1. Create a python 3 conda environment called gctree from the included environment file:
    ```bash
    conda env create -f environment.yml
    ```
2. Activate the environment:
    ```bash
    conda activate gctree
    ```

## Pipeline quick start

All commands should be issued from within the gctree repo directory.

### inference
- *input file*: `FASTA` or `PHYLIP` file containing a sequence for each observed individual/cell, and an additional sequence containing the ancestral genotype of all observed sequences (used for outgroup rooting).
- *run inference*:
    ```
    scons --inference --outdir=<output directory path> --input=<input FASTA or PHYLIP file> --root_id=<id of ancestral sequence in input file>
    ```
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
scons --simulate  --outdir=<output directory path> --N=<integer population size to simulate>
```

## Pipeline example

### run gctree inference pipeline on the included `FASTA` file

* **Example input data set**

  See the quickstart docs page for a description of these example data.

* **Run inference**  

    From within the `gctree` repository directory:
    ```
    scons --inference --input=example/150228_Clone_3-8.fasta --outdir=_build --id_abundances --root_id=GL --jobs=2
    ```
    This command will produce output in subdirectory `_build/`.


* **Explanation of arguments**    

    `--outdir=_build` specifies that results are to be saved in directory `_build/` (which will be created if it does not exist)

    `--id_abundances` flag means that integer sequence IDs in the input file are interpreted as abundances. The example input `FASTA` includes a sequence with id "17".

    `--root_id=GL` indicates that the root root sequence has id "GL" in the input `FASTA`. This sequence is the germline sequence of the V gene used in the V(D)J rearrangement that defines this clonal family.

    `--jobs=2` indicates that 2 parallel processes should be used

    If running on a remote machine via ssh, it may be necessary to provide the flag `--xvfb` which will allow X rendering of ETE trees without X forwarding.

## Inference pipeline

`scons --inference ...`

### required arguments

`--input=[path]` path to `FASTA` or `PHYLIP` input alignment

`--outdir=[path]` directory for output (created if does not exist)

`--root_id=[string]` ID of root sequence in input file used for outgroup rooting, default 'root'. For BCRs, we assume a known root V(D)J rearrangemnt is an additional sequence in our alignment, regardless of whether it was observed or not. This ancestral sequence must appear as an additional sequence. For applications without a definite root state, an observed sequence can be used to root the tree by duplicating it in the alignment and giving it a new id, which can be passed as this argument.

### optional arguments

`--colorfile=[path]  ` path to a file of plotting colors for cells in the input file. Example, if the input file contains a sequence with ID `cell_1`, this cell could be colored red in the tree image by including the line `cell_1,red` in the color file.

`--bootstrap=[int] ` boostrap resampling, and inference on each, default no bootstrap

`--id_abundances` if this flag is set, parse input IDs that are integers as indicating sequence abundance. Otherwise each line in the input is assumed to indicate an individual (non-deduplicated) sequence. **NOTE:** the example input `FASTA` file `example/150228_Clone_3-8.fasta` requires this option.

## Simulation pipeline

`scons --simulation ...`

### required arguments

`--N=[int]` populaton size to simulate. Note that `N=1` is satisfied before the first time step, so this choice will return the root with no mutation.  

`--outdir=[path]` directory for output (created if does not exist)

### optional arguments

`--root=[string]             ` DNA sequence of root sequence from which to begin simulating, a default is used if omitted

`--mutability=[path]          ` path to S5F mutability file, default 'S5F/mutability'

`--substitution=[path]        ` path to S5F substitution file, default 'S5F/substitution'

`--lambda=[float, float, ...] ` values for Poisson branching parameter for simulation, default 2.0

`--lambda0=[float, float, ...]` values for baseline mutation rate, default 0.25

`--T=[int]                    ` time steps to simulate (alternative to `--N`)

`--nsim=[int]                 ` number of simulation of each set of parameter combination, default 10

`--n=[int]                    ` number of cells to sample from final population, default all

## Optional arguments for both inference and simulation pipelines

`--jobs=[int]  ` number of parallel processes to use

`--srun        ` should cluster jobs be submitted with Slurm's srun?

`--frame=[int] ` codon reading frame, default `None`

`--quick       ` less thorough parsimony tree search (faster, but smaller parsimony forest)

`--idlabel     ` label sequence IDs on tree, and write `FASTA` alignment of distinct sequences. The mapping of the unique names in this `FASTA` file to the cell names in the original input file can be found in the output file with suffix `.idmap`

`--xvfb        ` needed for X rendering in on remote machines. Try setting this option if you get the error:`ETE: cannot connect to X server`

 `--dnaml`    include results for maximum likelihood tree inference using `dnaml` from the PHYLIP package
