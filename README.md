# gctree

The project is to come up with (relatively simple) models of B cell diversification within germinal centers that can be used as part of inference. When we say relatively simple, we mean analogous to the coalescent rather than something like

[Meyer-Hermann, M., Mohr, E., Pelletier, N., Zhang, Y., Victora, G. D., & Toellner, K.-M. (2012). A theory of germinal center B cell selection, division, and exit. Cell Reports, 2(1), 162–174.] (http://doi.org/10.1016/j.celrep.2012.05.010)

The "bright line" dividing the class of models of interest to us are those for which we can compute likelihoods efficiently.

[Paperpile share] (https://paperpile.com/shared/CF07th)  
[Overleaf document] (https://www.overleaf.com/5906733ngwstr)

There's going to be some simulation, validation, etc, with this project. For that we use [nestly] (http://nestly.readthedocs.io/en/latest/).

## Scripts

* TasParse.py - this reads the Tas data excel file, extracts the sequences from lymph node 2 germinal center 1, aligns with muscle, and prints to phyllip format.
* gctree.py - classes for simulating, and compute likelihoods, for collapsed trees generated from a binary branching process with mutation and infinite types, as well as forests of such trees.
	* The main routine takes a phlyip outtree file as a command line argument, converts each tree therein to a collapsed tree, computes there likelihoods, prints the MLE tree parameters, and plots the MLE tree to a file specified with the plot_file argument.
	* In a test mode, the main routine also performs a test of the likelihood against a by-hand calculation for a simple tree, and simulates a forest of trees and performs MLE. It outputs plot foo.pdf.
	* "python gctree.py -h" for usage info
