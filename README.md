# gctree

The project is to come up with (relatively simple) models of B cell diversification within germinal centers that can be used as part of inference. When we say relatively simple, we mean analogous to the coalescent rather than something like

[Meyer-Hermann, M., Mohr, E., Pelletier, N., Zhang, Y., Victora, G. D., & Toellner, K.-M. (2012). A theory of germinal center B cell selection, division, and exit. Cell Reports, 2(1), 162–174.] (http://doi.org/10.1016/j.celrep.2012.05.010)

The "bright line" dividing the class of models of interest to us are those for which we can compute likelihoods efficiently.

[Paperpile share] (https://paperpile.com/shared/CF07th)  
[Overleaf document] (https://www.overleaf.com/5906733ngwstr)

There's going to be some simulation, validation, etc, with this project. For that we use [nestly] (http://nestly.readthedocs.io/en/latest/).


## scons pipelines

**Inference pipeline:**
```
scons --inference
  --gctree/--igphyml[inferrence using either gctree, igphyml or both]
  --frame=[DNA reading frame or None, if None stop codons are allowed]
  --outdir=[path to output] --fasta=[path to Tas et al. style fasta]
  --naiveID=[ID of naive sequence in fasta file, default 'naive']
  --srun[optional for SLURM batch job submission]
```

**Simulation/validation pipeline:**
```
scons --simulate 
  --gctree/--igphyml[inferrence using either gctree, igphyml or both]
  --frame=[DNA reading frame or None, if None stop codons are allowed]
  --outdir=[path to output] --naive=[DNA seq of naive sequence from which to start simulating, used default if omitted]
  --mutability=[path to S5F mutability file, default 'S5F/mutability']
  --substitution=[path to S5F substitution file, default 'S5F/substitution']
  --lambda=[poisson branching probability for simulation]
  --lambda0=[baseline mutation rate] --r=[sampling probability]
  --N=[simulation tree size, default None]
  --T=[max generations, default None]
  --n=[number of simulation of each set of parameter combination, default 10]
  --srun[optional for SLURM batch job submission]
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
