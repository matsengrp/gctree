#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Infer trees from germinal center data, and validate inference method with simulation
"""
from __future__ import print_function
import os
import sys
import subprocess
import sconsutils
from warnings import warn
from SCons.Script import Environment, AddOption

sconsutils

# Set up SCons environment
environ = os.environ.copy()
env = Environment(ENV=environ)

# Add stuff to PATH
env.PrependENVPath("PATH", "bin")

# Setting up command line arguments/options
AddOption("--inference", action="store_true", help="Run inference")
inference = GetOption("inference")
AddOption(
    "--simulate",
    action="store_true",
    help="validation subprogram, instead of inference",
)
simulate = GetOption("simulate")
AddOption("--srun", action="store_true", help="Should jobs be submitted with srun?")
if GetOption("srun"):
    CommandRunner = env.SRun
else:
    CommandRunner = env.Command
AddOption("--frame", type="int", default=None, help="codon frame")
frame = GetOption("frame")
AddOption("--dnaml", action="store_true", help="use dnaml inference")
dnaml = GetOption("dnaml")
AddOption("--outdir", type="string", help="directory in which to output results")
outdir = GetOption("outdir")
AddOption(
    "--quick", action="store_true", help="less thourough dnapars tree search (faster)"
)
quick = GetOption("quick")
AddOption(
    "--idlabel",
    action="store_true",
    help="label sequence ids on tree, and write associated alignment",
)
idlabel = GetOption("idlabel")
AddOption(
    "--xvfb",
    action="store_true",
    help="use virtual X, for rendering ETE trees on a remote server",
)
xarg = "TMPDIR=/tmp xvfb-run -a " if GetOption("xvfb") else ""
AddOption(
    "--nobuff",
    action="store_true",
    help="use stdbuf to prevent line buffering on linux",
)
AddOption(
    "--disambiguate_with_mutability",
    action="store_true",
    help="use mutability model provided using ``mutability'' and ``substitution'' arguments to attempt to optimally resolve ambiguities in dnapars output trees"
)
disambiguate_with_mutability = GetOption("disambiguate_with_mutability")
AddOption(
    "--mutability",
    type="string",
    metavar="PATH",
    default="S5F/Mutability.csv",
    help="path to S5F mutability data",
)
mutability = GetOption("mutability")
AddOption(
    "--substitution",
    type="string",
    metavar="PATH",
    default="S5F/Substitution.csv",
    help="path to S5F substitution data",
)
substitution = GetOption("substitution")
buffarg = "stdbuf -oL " if GetOption("nobuff") else ""


class InputError(Exception):
    """Exception raised for errors in the input."""

if not simulate and not inference and not GetOption("help"):
    raise InputError(
        'Please provide one of the required arguments. Either "--inference" or "--simulate".'
        'Command line help can then be evoked by "-h" or "--help" and found in the bottom'
        'of the output under "Local Options".'
    )

if simulate:
    AddOption(
        "--root",
        type="string",
        default="ggacctagcctcgtgaaaccttctcagactctgtccctcacctgttctgtcactg"
        "gcgactccatcaccagtggttactggaactggatccggaaattcccagggaataa"
        "acttgagtacatggggtacataagctacagtggtagcacttactacaatccatct"
        "ctcaaaagtcgaatctccatcactcgagacacatccaagaaccagtactacctgc"
        "agttgaattctgtgactactgaggacacagccacatattactgt",
        help="sequence of root from which to simulate",
    )
    root = GetOption("root")
    AddOption(
        "--lambda",
        type="float",
        action="append",
        default=[],
        help="Poisson branching parameter for simulation",
    )
    lambda_list = GetOption("lambda")
    if len(lambda_list) == 0:
        lambda_list = [2.0]
    AddOption(
        "--lambda0",
        type="float",
        action="append",
        default=[],
        help="baseline mutation rate",
    )
    lambda0_list = GetOption("lambda0")
    if len(lambda0_list) == 0:
        lambda0_list = [0.25]
    AddOption("--n", type="int", default=None, help="cells downsampled")
    n = GetOption("n")
    AddOption(
        "--N",
        type="int",
        default=None,
        help="simulation size (number of cells observerved)",
    )
    N = GetOption("N")
    AddOption("--T", type="int", action="append", default=[], help="observation time")
    T = GetOption("T")
    AddOption(
        "--nsim",
        type="int",
        default=10,
        help="number of simulations with each parameter parameter choice",
    )
    nsim = GetOption("nsim")
    AddOption(
        "--experimental",
        type="string",
        action="append",
        default=[],
        help="experimental fastas for comparing summary stats",
    )
    experimental_list = GetOption("experimental")
    AddOption(
        "--root_idexp",
        type="string",
        default="root0",
        help="id of root seq in the experimental data",
    )
    root_idexp = GetOption("root_idexp")

elif inference:
    AddOption(
        "--input",
        dest="input",
        type="string",
        action="append",
        default=[],
        metavar="PATH",
        help="path to input fasta or phylip",
    )
    input_file = GetOption("input")
    if len(input_file) == 1:
        input_file = input_file[0]
        input_file2 = None
    else:
        input_file, input_file2 = input_file
    AddOption(
        "--colorfile",
        dest="colorfile",
        type="string",
        default=None,
        metavar="PATH",
        help="optional two column csv file with colors to associate with each cell",
    )
    colorfile = GetOption("colorfile")
    AddOption(
        "--root_id",
        type="string",
        metavar="seqID",
        default="root",
        help="id of root sequence",
    )
    root_id = GetOption("root_id")
    AddOption(
        "--id_abundances",
        action='store_true',
        help="interpret integer input ids as abundances",
    )
    id_abundances = GetOption("id_abundances")
    AddOption(
        "--bootstrap",
        type="int",
        default=0,
        help="boostrap resampling, and inference on each (default no bootstrap)",
    )
    bootstrap = GetOption("bootstrap")
    AddOption(
        "--extended_parsimony_search",
        action="store_true",
        help="search for more maximum parsimony trees using history DAG",
    )
    extended_parsimony_search = GetOption("extended_parsimony_search")

# First call after all arguments have been parsed
# to enable correct command line help.
if simulate and not GetOption("help"):
    if outdir is None:
        raise InputError("outdir must be specified")
    SConscript(
        "SConscript.simulation",
        exports="env dnaml quick idlabel outdir root mutability substitution lambda_list lambda0_list n frame N T nsim CommandRunner experimental_list root_idexp xarg buffarg",
    )
elif inference and not GetOption("help"):
    if None in [input_file, outdir]:
        raise InputError("input fasta or phylip and outdir must be specified")
    SConscript(
        "SConscript.inference",
        exports="env dnaml quick idlabel frame input_file input_file2 outdir root_id id_abundances CommandRunner bootstrap xarg buffarg colorfile mutability substitution disambiguate_with_mutability extended_parsimony_search",
    )
