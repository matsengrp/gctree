# gctree

The project is to come up with (relatively simple) models of B cell diversification within germinal centers that can be used as part of inference. When we say relatively simple, we mean analogous to the coalescent rather than something like

[Meyer-Hermann, M., Mohr, E., Pelletier, N., Zhang, Y., Victora, G. D., & Toellner, K.-M. (2012). A theory of germinal center B cell selection, division, and exit. Cell Reports, 2(1), 162–174.] (http://doi.org/10.1016/j.celrep.2012.05.010)

The "bright line" dividing the class of models of interest to us are those for which we can compute likelihoods efficiently.

[Paperpile share] (https://paperpile.com/shared/CF07th)  
[Overleaf document] (https://www.overleaf.com/5906733ngwstr)

There's going to be some simulation, validation, etc, with this project. For that we use [nestly] (http://nestly.readthedocs.io/en/latest/).

## Scripts

* TasParse.py - this reads the Tas data excel file, extracts the sequences from lymph node 2 germinal center 1, aligns with muscle, and prints to phyllip format.
* recurse.py - classes for simulating, and compute likelihoods, for collapsed trees generated from a binary branching process with mutation and infinite types, as well as forests of such trees. The main routine performs some tests using a simulated forest and outputs plots foo.pdf and bar.pdf.
	* foo.pdf shows the 2-norm of the difference between the analytical likelihood gradient and a finite difference approximation over the parameter space. Hopefully it is small everywhere
	* bar.pdf shows the likelihood surface for the randomized forest with the true parameters (red) and MLE parameters (blue)

