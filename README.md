# GCtree

Implements phylogenetic inference for data with repeated sequences, as described in [[link GCtree paper citation](???)]

## dependencies
* scons
* Python 2, with modules:
  * scipy
  * matplotlib
  * seaborn
  * pandas
  * biopython
  * [ete3](http://etetoolkit.org/download/)
  * [nestly](https://pypi.python.org/pypi/nestly/0.6)
* [PHYLIP](http://evolution.genetics.washington.edu/phylip/getme-new.html)
  * PHYLIP's `dnapars` program is used for generating parsimony trees
* X11 or xvfb-run (for rendering phylogenetic trees using ete3)

## scons pipelines

Two programs are implemented:
- an inference program for experimental data
- a simulation/inference/validation program

All commands should be issued from within the gctree repo directory.

## required arguments for both inference and simulation programs

`--outdir=[path]` directory for output (created if does not exist)


## inference

`scons --inference ...`

### required arguments

`--fasta=[path]` path to FASTA input alignment

### optional arguments

`--naiveID=[string]` ID of naive sequence in FASTA file used for outgroup rooting, default 'naive'

`--bootstrap=[int]` boostrap resampling, and inference on each, default no bootstrap

## simulation

`scons --simulation ...`

### required arguments for simulation program

`--N=[int]` populaton size to simulate

### optional arguments for simulation program

`--naive=[string]             ` DNA sequence of naive sequence from which to begin simulating, a default is used if omitted

`--mutability=[path]          ` path to S5F mutability file, default 'S5F/mutability'

`--substitution=[path]        ` path to S5F substitution file, default 'S5F/substitution'

`--lambda=[float, float, ...] ` values for Poisson branching parameter for simulation, default 2.0

`--lambda0=[float, float, ...]` values for baseline mutation rate, default 0.25

`--T=[int]                    ` time steps to simulate (alternative to `--N`)

`--nsim=[int]                 ` number of simulation of each set of parameter combination, default 10

`--n=[int]                    ` number of cells to sample from final population, default all

## optional arguments for both inference and simulation programs

`--srun        ` should cluster jobs be submitted with Slurm's srun?

`--frame=[int] ` codon reading frame, default `None`

`--quick       ` less thorough parsimony tree search (faster, but smaller parsimony forest)

`--idlabel     ` label sequence IDs on tree, and write FASTA alignment distinct sequences


## `gctree.py`
Underlying both pipelines is the `gctree.py` Python library for simulating and compute likelihoods for collapsed trees generated from a binary branching process with mutation and infinite types, as well as forests of such trees. General usage info `gctree.py --help`. There are three subprograms, each of which has usage info:
* `gctree.py infer --help`: takes an `outfile` file made by phylip's `dnapars` as a command line argument, converts each tree therein to a collapsed tree, and ranks by GCtree likelihood.
* `gctree.py simulate --help`: simulate data
* `gctree.py test --help`: performs tests of the likelihood and outputs validation plots.

# functionality under development

## arguments for both inference and simulation programs

`--igphyml`  include results for tree inference with the IgPhyML package

`--dnaml`    include results for maximum likelihood tree inference using `dnaml` from the PHYLIP package

`--nogctree` do not perform gctree inference

## arguments for non-neutral simulation

`--selection`    simulation with affinity selection

`--target_dist`  distance to selection target

`--target_count` number of targets

`--verbose`      verbose printing

`--carry_cap`    carrying capacity of germinal center

`--skip_update`  skip update step

## additional depencies for development functionality
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
