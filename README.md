# GCtree

Implements phylogenetic inference for data with repeated sequences, as described in:

DeWitt, Mesin, Victora, Minin and Matsen, *Using genotype abundance to improve phylogenetic inference*, [arXiv:1708.08944](https://arxiv.org/abs/1708.08944).

## DEPENDENCIES
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
  * PHYLIP's `dnapars` program is used for generating parsimony trees, so the command-line program `dnapars` should be on your path
* X11 or xvfb-run (for rendering phylogenetic trees using ete3)

**NOTE:** for installing scons, ete, and other python dependencies, [conda](https://conda.io/docs/) is recommended:
```bash
conda install -c etetoolkit ete3 ete3_external_apps
conda install biopython matplotlib nestly pandas scipy scons seaborn
```
Alternatively, an example linux environment spec file is included (`spec-file.txt`), which may be used to create a conda environment.
For example, to create an environment called `gctree`, execute `conda create --name gctree --file spec-file.txt`, then activate the environment with `source activate gctree`.

## SCONS PIPELINES

Two programs are implemented:
- an inference program for experimental data
- a simulation/inference/validation program

All commands should be issued from within the gctree repo directory.

## QUICK START

* **inference:**  `scons --inference --outdir=<output directory path> --fasta=<input fasta file>`
* **simulation:** `scons --simulate  --outdir=<output directory path> --N=<integer population size to simulate>`

### **example:** to run GCtree inference on the included FASTA file on a remote machine
```
scons --inference --fasta=Victora_data/150228_Clone_3-8.fasta --outdir=test --converter=tas --naiveID=GL --xvfb --jobs=10
```
Results are saved in directory `test/`. The `--converter=tas` argument means that integer sequence IDs in the FASTA file are interpreted as abundances. The flag `--xvfb` allows X rendering of ETE trees on remote machines. The argument `--jobs=10` indicates that 10 parallel processes should be used.

## **INFERENCE**

`scons --inference ...`

### required arguments

`--fasta=[path]` path to FASTA input alignment
`--outdir=[path]` directory for output (created if does not exist)

### optional arguments

`--naiveID=[string]` ID of naive sequence in FASTA file used for outgroup rooting, default 'naive'

`colorfile=[path]  ` path to a file of plotting colors for cells in the input FASTA file. Example, if the FASTA contains a sequence with ID `cell_1`, this cell could be colored red in the tree image by including the line `cell_1,red` in the color file.

`--bootstrap=[int] ` boostrap resampling, and inference on each, default no bootstrap

`--converter=[string]` if set to "tas", parse FASTA input IDs that are integers as indicating sequence abundance. Otherwise each line in the FASTA is assumed to indicate an individual (non-deduplicated) sequence. **NOTE:** the included FASTA file `Victora_data/150228_Clone_3-8.fasta` requires this option.

## **SIMULATION**

`scons --simulation ...`

### required arguments

`--N=[int]` populaton size to simulate
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

### additional depencies for development functionality
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
