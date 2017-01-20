# gctree

The project is to come up with (relatively simple) models of B cell diversification within germinal centers that can be used as part of inference. When we say relatively simple, we mean analogous to the coalescent rather than something like

[Meyer-Hermann, M., Mohr, E., Pelletier, N., Zhang, Y., Victora, G. D., & Toellner, K.-M. (2012). A theory of germinal center B cell selection, division, and exit. Cell Reports, 2(1), 162–174.] (http://doi.org/10.1016/j.celrep.2012.05.010)

The "bright line" dividing the class of models of interest to us are those for which we can compute likelihoods efficiently.

[Paperpile share] (https://paperpile.com/shared/CF07th)  
[Overleaf document] (https://www.overleaf.com/5906733ngwstr)

There's going to be some simulation, validation, etc, with this project. For that we use [nestly] (http://nestly.readthedocs.io/en/latest/).

## scons pipelines

* Inference pipeline: `scons --outdir=[path to output] --fasta=[path to Tas et al. style fasta] --naiveID=[ID of naive sequence in fasta file, default 'naive']`
* Simulation/validation pipeline: `scons --validation --outdir=[path to output] --naive=[DNA seq of naive sequence from which to start simulating, used default if omitted] --mutability=[path to S5F mutability file, default 'S5F/mutability'] --substitution=[path to S5F substitution file, default 'S5F/substitution'] --p=[branching probability for simulation, default 0.49] --lambda0=[baseline mutation rate, default .3] --r=[sampling probability, default 1.] --n=[minimum simulation sequence set size, default 100]`
* gctree.py - library for simulating, and compute likelihoods, for collapsed trees generated from a binary branching process with mutation and infinite types, as well as forests of such trees. General usage info `gctree.py --help`. There are three subprograms, each of which has usage info:
	* `gctree.py inference --help`: takes an `outfile` file made by phylip's `dnapars` as a command line argument, converts each tree therein to a collapsed tree, and ranks by likelihood under our model.
	* `gctree.py simulation --help`: simulate data under branching model, plus S5F and imperfect sampling.
	* `gctree.py validation --help`
	* `gctree.py test --help`: performs a test of the likelihood against a by-hand calculation for a simple tree, and simulates a forest of trees and performs MLE, outputting plots validating MLE and gradients.
