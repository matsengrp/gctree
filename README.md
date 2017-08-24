# GCtree

Implements phylogenetic inference for data with repeated sequences, as described in [[link GCtree paper citation](???)]

## scons pipelines

Two programs are implemented:
- an inference program for experimental data, initiated with `scons --inference`
- a simulation/inference/validation program, initiated with `scons --simulation`

All commands should be issued from within the gctree repo directory.

**required arguments for both inference and simulation programs**
```
      --outdir=[path]     directory for output (created if does not exist)
```

**required arguments for inference program:**
```
      --fasta=[path]      path to FASTA input alignment
```

**optional arguments for inference program:**
```
      --naiveID=[string]  ID of naive sequence in FASTA file, default 'naive'
      --bootstrap=[int]   boostrap resampling, and inference on each, default no bootstrap
```

**optional arguments for simulation program:**
```   
      --naive=[string]              DNA sequence of naive sequence from which to begin simulating, a default is used if omitted
      --mutability=[path]           path to S5F mutability file, default 'S5F/mutability'
      --substitution=[path]         path to S5F substitution file, default 'S5F/substitution'
      --lambda=[float, float, ...]  values for poisson branching parameter for simulation, default 2.0
      --lambda0=[float, float, ...] values for baseline mutation rate, default 0.25
      --N=[int]                     populaton size to simulate
      --T=[int]                     time steps to simulate (alternatve to --N)
      --nsim=[int]                  number of simulation of each set of parameter combination, default 10
      --n=[int]                     number of cells to sample from final population, default all
```
**NOTE:** one of `--N` or `--T` must be specified.

**optional arguments for both inference and simulation programs**
```
  --srun                      should cluster jobs be submitted with srun?
  --frame=[int]               codon reaading frame, default None
  --quick                     less thourough parsimony tree search (faster, but smaller parsimony forest)
  --idlabel                   label sequence IDs on tree, and write FASTA alignment distinct sequences
  ```

## `gctree.py`
Underlying both pipelines is the `gctree.py` Python library for simulating and compute likelihoods for collapsed trees generated from a binary branching process with mutation and infinite types, as well as forests of such trees. General usage info `gctree.py --help`. There are three subprograms, each of which has usage info:
* `gctree.py infer --help`: takes an `outfile` file made by phylip's `dnapars` as a command line argument, converts each tree therein to a collapsed tree, and ranks by GCtree likelihood.
* `gctree.py simulate --help`: simulate data
* `gctree.py test --help`: performs tests of the likelihood and outputs validation plots.


### Dependencies
* scons
* Python 2, with modules:
  * ete3
  * scipy
  * matplotlib
  * seaborn
  * pandas
  * biopython
* xvfb-run or X11
  * Used for rendering phylogenetic trees using ete3  
* PHYLIP (http://evolution.genetics.washington.edu/phylip/getme-new.html)
  * PHYLIP's `dnapars` program is used for generating parsimony trees
  * Needs to be in $PATH

**For running IgPhyML inference:**
* IgPhyML (https://github.com/kbhoehn/IgPhyML)
  * Needs to be in $PATH
* perl5, with modules:
  * PDL
  * PDL::LinearAlgebra::Trans


### Getting python and perl dependencies
```
sudo apt-get install python-pip scons
pip install --user ete3 seaborn numpy scipy matplotlib pandas biopython argparse
cpan
> install PDL
> install PDL::LinearAlgebra::Trans
```
