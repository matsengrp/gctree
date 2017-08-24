# GCtree

Implements phylogenetic inference for data with repeated sequences, as described in [link GCtree paper citation](???)

## scons pipelines

Two programs are implemented:
- an inference program for experimental data, initiated with `scons --inference`
- a simulation/inference/validation program, initiated with `scons --simulation`

**options for both inference and simulation programs**
```
  --srun                      should cluster jobs be submitted with srun?
  --frame=[int]               codon reaading frame (optional)
  --gctree                    flag for performing gctree inference
  --igphyml                   flag for performing igphyml inference
  --dnaml                     flag for performing dnaml inference
  --outdir=[path]             directory for output (created if does not exist)
  --quick                     less thourough dnapars tree search (faster, but smaller parsimony forest)
  --idlabel                   label sequence IDs on tree, and write FASTA alignment distinct sequences
  ```
NOTE: at least one of `--gctree`, `--igphyml`, and `--dnaml` must be set

**Inference program:**
```
      --fasta=[path]      path to FASTA input alignment
      --naiveID=[string]  ID of naive sequence in FASTA file, default 'naive'
```

**Simulation/validation program:**
```   
      --naive=[string]            DNA sequence of naive sequence from which to begin simulating, a default is used if omitted
      --mutability=[path]         path to S5F mutability file, default 'S5F/mutability'
      --substitution=[path]       path to S5F substitution file, default 'S5F/substitution'
      --lambda=[float]            poisson branching parameter for simulation
      --lambda0=[float]           baseline mutation rate
      --N=[int]     populaton size to simulate
      --T=[int]     time steps to simulate (alternatve to --N
      --nsim=[int]  number of simulation of each set of parameter combination, default 10
      --n=[int]     number of cells to sample from final population
```

## `gctree.py`
library for simulating, and compute likelihoods, for collapsed trees generated from a binary branching process with mutation and infinite types, as well as forests of such trees. General usage info `gctree.py --help`. There are three subprograms, each of which has usage info:
* `gctree.py inference --help`: takes an `outfile` file made by phylip's `dnapars` as a command line argument, converts each tree therein to a collapsed tree, and ranks by likelihood under our model.
* `gctree.py simulation --help`: simulate data under branching model, plus S5F and imperfect sampling.
* `gctree.py validation --help`
* `gctree.py test --help`: performs a test of the likelihood against a by-hand calculation for a simple tree, and simulates a forest of trees and performs MLE, outputting plots validating MLE and gradients.


### Dependencies
* scons
* xvfb-run or X11
  * Used for rendering phylogenetic trees using ete3
* Python 2, with modules:
  * ete3, seaborn, numpy, scipy, matplotlib, pandas, biopython, argparse

**For running GCtree inferrence:**
* PHYLIP (http://evolution.genetics.washington.edu/phylip/getme-new.html)
  * Needs to be in $PATH

**For running IgPhyML inferrence:**
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
