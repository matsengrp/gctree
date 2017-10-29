# GCtree

Implements phylogenetic inference for data with repeated sequences, as described in:

DeWitt, Mesin, Victora, Minin and Matsen, *Using genotype abundance to improve phylogenetic inference*, [arXiv:1708.08944](https://arxiv.org/abs/1708.08944).

## DEPENDENCIES
* scons
* Python 2.7, with modules:
  * scipy
  * matplotlib
  * seaborn
  * pandas
  * biopython
  * [ete3](http://etetoolkit.org)
  * [nestly](https://pypi.python.org/pypi/nestly)
* [PHYLIP](http://evolution.genetics.washington.edu/phylip)
* X11 or xvfb-run (for rendering phylogenetic trees using ete3)
* [seqmagick](https://github.com/fhcrc/seqmagick)

### Step-by-step installation instructions

0. For installing dependencies, [conda](https://conda.io/docs/) environment management is recommended. First install conda or miniconda.

1. Create a conda environment (named gctree in this example):
    ```bash
    conda create --name gctree
    ```
2. Activate the environment (note different commands for Linux/MaxOS and Windows):
    - on Linux/MacOS:
    ```bash
    source activate gctree
    ```
    - on Windows:
    ```bash
    activate gctree
    ``` 
3. Install ete3:
    ```bash
    conda install -c etetoolkit ete3 ete3_external_apps
    ```
4. Install python packages:
    ```bash
    conda install biopython matplotlib pandas scipy scons seaborn
    ```
5. Install nestly
    ```bash
    conda install -c conda-forge nestly
    ```
6. Install PHYLIP:
    ```bash
    conda install -c bioconda phylip
    ```
7. Install seqmagick:
    ```bash
    conda install -c cswarth seqmagick
    ```

## SCONS PIPELINES

Two programs are implemented:
- an inference program for experimental data
- a simulation/inference/validation program

All commands should be issued from within the gctree repo directory.

## QUICK START

* **inference:**  `scons --inference --outdir=<output directory path> --fasta=<input fasta file>`
* **simulation:** `scons --simulate  --outdir=<output directory path> --N=<integer population size to simulate>`

### inference output
After the inference pipeline has completed, the output directory will contain the following output files:
- `<input fasta file>.idmap`: text file mapping collapsed sequenced ids to cell ids from the original input file
- `<input fasta file>.counts`: text file mapping collapsed sequenced ids to their abundances
- `<input fasta file>.phylip`: phylip alignment file of collapsed sequences for computing parsimony trees
- `dnapars/`: directory of parsimony tree output from PHYLIP's dnapars
- `gctree.inference.*.svg`: rendered tree images for each of the parsimony trees
- `gctree.inference.abundance_rank.pdf`: histogram of genotype abundances
- `gctree.inference.likelihood_rank.pdf`: rank plot of GCtree likelihoods for the parsimony trees
- `gctree.inference.log`: log file containing parameter fits, numerical likelihood results, and any other program messages
- `gctree.inference.parsimony_forest.p`: a python pickle file containing the parsimony trees


### **example:** run GCtree inference on the included FASTA file
Example input data set `example_input/150228_Clone_3-8.fasta` contains heavy chain V gene sequences from 65 germinal B cells sorted from a brainbow mouse using multicolor fate mapping. These data were published in [Tas et al. 2016. *Visualizing Antibody Affinity Maturation in Germinal Centers.* Science 351 (6277)](http://science.sciencemag.org/content/351/6277/1048)) and shown in Fig. 4 (lymph node 2, germinal center 1).

![](gc1.png)

```
scons --inference --fasta=example_input/150228_Clone_3-8.fasta --outdir=test --converter=tas --naiveID=GL --jobs=2
```
`--outdir=test` specifies that results are to be saved in directory `test/` (which will be created if it does not exist). The `--converter=tas` argument means that integer sequence IDs in the FASTA file are interpreted as abundances. The argument `--naiveID=GL` indicates that the root naive sequence has id "GL". This sequence is the germline sequence of the V gene used in the V(D)J rearrangment that define this clonal family. The argument `--jobs=2` indicates that 2 parallel processes should be used. If running on a remote machine via ssh, it may be necessary to provide the flag `--xvfb` which will allow X rendering of ETE trees without X forwarding.

## **INFERENCE**

`scons --inference ...`

### required arguments

`--fasta=[path]` path to FASTA input alignment
`--outdir=[path]` directory for output (created if does not exist)
`--naiveID=[string]` ID of naive sequence in FASTA file used for outgroup rooting, default 'naive'. For BCRs, we assume a known naive V(D)J rearrangemnt is an additional sequence in our alignment, regardless of whether it was observed or not. This ancestral sequence must appear as an additional sequence. For applications without a definite root state, an observed sequence can be used to root the tree by duplicating it in the alignment and giving it a new id, which can be passed as this argument.

### optional arguments

`colorfile=[path]  ` path to a file of plotting colors for cells in the input FASTA file. Example, if the FASTA contains a sequence with ID `cell_1`, this cell could be colored red in the tree image by including the line `cell_1,red` in the color file.

`--bootstrap=[int] ` boostrap resampling, and inference on each, default no bootstrap

`--converter=[string]` if set to "tas", parse FASTA input IDs that are integers as indicating sequence abundance. Otherwise each line in the FASTA is assumed to indicate an individual (non-deduplicated) sequence. **NOTE:** the example input FASTA file `example_input/150228_Clone_3-8.fasta` requires this option.

## **SIMULATION**

`scons --simulation ...`

### required arguments

`--N=[int]` populaton size to simulate. Note that `N=1` is satisfied before the first time step, so this choice will return the root with no mutation.  

`--outdir=[path]` directory for output (created if does not exist)

### optional arguments

`--naive=[string]             ` DNA sequence of naive sequence from which to begin simulating, a default is used if omitted

`--mutability=[path]          ` path to S5F mutability file, default 'S5F/mutability'

`--substitution=[path]        ` path to S5F substitution file, default 'S5F/substitution'

`--lambda=[float, float, ...] ` values for Poisson branching parameter for simulation, default 2.0

`--lambda0=[float, float, ...]` values for baseline mutation rate, default 0.25

`--T=[int]                    ` time steps to simulate (alternative to `--N`)

`--nsim=[int]                 ` number of simulation of each set of parameter combination, default 10

`--n=[int]                    ` number of cells to sample from final population, default all

## OPTIONAL ARGUMENTS FOR BOTH INFERENCE AND SIMULATION PROGRAMS

`--jobs=[int]  ` number of parallel processes to use

`--srun        ` should cluster jobs be submitted with Slurm's srun?

`--frame=[int] ` codon reading frame, default `None`

`--quick       ` less thorough parsimony tree search (faster, but smaller parsimony forest)

`--idlabel     ` label sequence IDs on tree, and write FASTA alignment of distinct sequences. The mapping of the unique names in this FASTA file to the cell names in the original input FASTA file can be found in the output file with suffix `.idmap`

`--xvfb        ` needed for X rendering in on remote machines

   * Try setting the above option if you get the error:`ETE: cannot connect to X server`

## `gctree.py`
Underlying both pipelines is the `gctree.py` Python library (located in the `bin/` subdirectory) for simulating and compute likelihoods for collapsed trees generated from a binary branching process with mutation and infinite types, as well as forests of such trees. General usage info `gctree.py --help`. There are three subprograms, each of which has usage info:
* `gctree.py infer --help`: takes an `outfile` file made by phylip's `dnapars` as a command line argument, converts each tree therein to a collapsed tree, and ranks by GCtree likelihood.
* `gctree.py simulate --help`: simulate data
* `gctree.py test --help`: performs tests of the likelihood and outputs validation plots.

The under-the-hood functionality of the `gctree.py` library might be useful for some users trying to go beyond the scons pipelines. For example mapping colors to tree image nodes can be achieved with the `--colormap` argument. Colors can be useful for visualizing other cell/genotype properties on the tree.

## FUNCTIONALITY UNDER DEVELOPMENT

### arguments for both inference and simulation programs

`--igphyml`  include results for tree inference with the IgPhyML package

`--dnaml`    include results for maximum likelihood tree inference using `dnaml` from the PHYLIP package

`--nogctree` do not perform gctree inference

### arguments for non-neutral simulation

`--selection`    simulation with affinity selection

`--target_dist`  distance to selection target

`--target_count` number of targets

`--verbose`      verbose printing

`--carry_cap`    carrying capacity of germinal center

`--skip_update`  skip update step

### additional dependencies for development functionality
* IgPhyML (https://github.com/kbhoehn/IgPhyML)
  * Needs to be in $PATH
* perl5, with modules:
  * PDL
  * PDL::LinearAlgebra::Trans

* installing python and perl dependencies *
```
sudo apt-get install python-pip scons
pip install --user ete3 seaborn numpy scipy matplotlib pandas biopython nestly
cpan
> install PDL
> install PDL::LinearAlgebra::Trans
```
