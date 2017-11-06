#! /usr/bin/env python
# -*- coding: utf-8 -*-
'''
Infer trees from germinal center data, and validate inference method with simulation
'''
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
env.PrependENVPath('PATH', 'bin')

# Setting up command line arguments/options
AddOption('--inference',
          action='store_true',
          help='Run inference')
inference = GetOption('inference')
AddOption('--simulate',
          action='store_true',
          help='validation subprogram, instead of inference')
simulate = GetOption('simulate')
AddOption('--srun',
          action='store_true',
          help='Should jobs be submitted with srun?')
if GetOption('srun'):
    CommandRunner = env.SRun
else:
    CommandRunner = env.Command
AddOption('--frame',
          type='int',
          default=None,
          help='codon frame')
frame = GetOption('frame')
AddOption('--nogctree',
           action='store_true',
           help='don''t use gctree inference')
gctree = (GetOption('nogctree') != True)
AddOption('--igphyml',
           action='store_true',
           help='use igphyml inference')
igphyml = GetOption('igphyml')
AddOption('--dnaml',
           action='store_true',
           help='use dnaml inference')
dnaml = GetOption('dnaml')
AddOption('--outdir',
          type='string',
          help="directory in which to output results")
outdir = GetOption('outdir')
AddOption('--quick',
           action='store_true',
           help='less thourough dnapars tree search (faster)')
quick = GetOption('quick')
AddOption('--idlabel',
           action='store_true',
           help='label sequence ids on tree, and write associated alignment')
idlabel = GetOption('idlabel')
AddOption('--xvfb',
          action='store_true',
          help='use virtual X, for rendering ETE trees on a remote server')
xarg = 'TMPDIR=/tmp xvfb-run -a ' if GetOption('xvfb') else ''
AddOption('--nobuff',
          action='store_true',
          help='use stdbuf to prevent line buffering on linux')
buffarg = 'stdbuf -oL ' if GetOption('nobuff') else ''

class InputError(Exception):
    """Exception raised for errors in the input."""

if not gctree and not igphyml and not GetOption('help'):
    raise InputError('must set at least one inference method')
if igphyml and frame != 1:
    raise InputError('frame must equal 1 for igphyml')

if not simulate and not inference and not GetOption('help'):
    raise InputError('Please provide one of the required arguments. Either "--inference" or "--simulate".'
                     'Command line help can then be evoked by "-h" or "--help" and found in the bottom'
                     'of the output under "Local Options".')

if simulate:
    AddOption('--naive',
              type='string',
              default='ggacctagcctcgtgaaaccttctcagactctgtccctcacctgttctgtcactg'
                      'gcgactccatcaccagtggttactggaactggatccggaaattcccagggaataa'
                      'acttgagtacatggggtacataagctacagtggtagcacttactacaatccatct'
                      'ctcaaaagtcgaatctccatcactcgagacacatccaagaaccagtactacctgc'
                      'agttgaattctgtgactactgaggacacagccacatattactgt',
              help='sequence of naive from which to simulate')
    naive = GetOption('naive')
    AddOption('--mutability',
              type='string',
              metavar='PATH',
              default='S5F/Mutability.csv',
              help='path to S5F mutability data')
    mutability = GetOption('mutability')
    AddOption('--substitution',
              type='string',
              metavar='PATH',
              default='S5F/Substitution.csv',
              help='path to S5F substitution data')
    substitution = GetOption('substitution')
    AddOption('--lambda',
              type='float',
              action='append',
              default=[],
              help='Poisson branching parameter for simulation')
    lambda_list = GetOption('lambda')
    if len(lambda_list) == 0:
        lambda_list = [2.]
    AddOption('--lambda0',
              type='float',
              action='append',
              default=[],
              help='baseline mutation rate')
    lambda0_list = GetOption('lambda0')
    if len(lambda0_list) == 0:
        lambda0_list = [.25]
    AddOption('--n',
              type='int',
              default=None,
              help='cells downsampled')
    n = GetOption('n')
    AddOption('--N',
              type='int',
              default=None,
              help='simulation size (number of cells observerved)')
    N = GetOption('N')
    AddOption('--T',
              type='int',
              action='append',
              default=[],
              help='observation time')
    T = GetOption('T')
    AddOption('--nsim',
              type='int',
              default=10,
              help='number of simulations with each parameter parameter choice')
    nsim = GetOption('nsim')
    AddOption('--experimental',
              type='string',
              action='append',
              default=[],
              help='experimental fastas for comparing summary stats')
    experimental_list = GetOption('experimental')
    AddOption('--naiveIDexp',
              type='string',
              default='naive0',
              help='id of naive seq in the experimental data')
    naiveIDexp = GetOption('naiveIDexp')
    AddOption('--selection',
              action='store_true',
              help='Simulation with affinity selection.')
    selection = GetOption('selection')
    if selection:
        AddOption('--target_dist',
                  type='int',
                  default=10,
                  help='Distance to selection target.')
        target_dist = GetOption('target_dist')
        AddOption('--target_count',
                  type='int',
                  default=10,
                  help='Number of targets.')
        target_count = GetOption('target_count')
        AddOption('--verbose',
                  action='store_true',
                  help='Verbose printing.')
        verbose = GetOption('verbose')
        AddOption('--carry_cap',
                  type='int',
                  default=1000,
                  help='Number of targets.')
        carry_cap = GetOption('carry_cap')
        AddOption('--skip_update',
                  type='int',
                  default=100,
                  help='Skip update step.')
        skip_update = GetOption('skip_update')
        selection_param = (target_dist, target_count, verbose, carry_cap, skip_update)
    else:
        selection_param = None

elif inference:
    AddOption('--input',
              dest='input',
              type='string',
              action='append',
              default=[],
              metavar='PATH',
              help='path to input fasta or phylip')
    input_file = GetOption('input')
    if len(input_file) == 1:
        input_file = input_file[0]
        input_file2 = None
    else:
        input_file, input_file2 = input_file
    AddOption('--colorfile',
              dest='colorfile',
              type='string',
              default=None,
              metavar='PATH',
              help='optional two column csv file with colors to associate with each cell')
    colorfile = GetOption('colorfile')
    AddOption('--naiveID',
              type='string',
              metavar='seqID',
              default='naive',
              help='id of naive sequence')
    naiveID = GetOption('naiveID')
    AddOption('--converter',
              type='string',
              default=None,
              help='Converter to convert input format e.g. the Victora lab GC fasta format')
    converter = GetOption('converter')
    AddOption('--bootstrap',
              type='int',
              default=0,
              help='boostrap resampling, and inference on each (default no bootstrap)')
    bootstrap = GetOption('bootstrap')

# First call after all arguments have been parsed
# to enable correct command line help.
if simulate and not GetOption('help'):
    if outdir is None:
        raise InputError('outdir must be specified')
    SConscript('SConscript.simulation',
               exports='env gctree igphyml dnaml quick idlabel outdir naive mutability substitution lambda_list lambda0_list n frame N T nsim CommandRunner experimental_list naiveIDexp selection_param xarg buffarg')
elif inference and not GetOption('help'):
    if None in [input_file, outdir]:
        raise InputError('input fasta orp phylip and outdir must be specified')
    SConscript('SConscript.inference', exports='env gctree igphyml dnaml quick idlabel frame input_file input_file2 outdir naiveID converter CommandRunner bootstrap xarg buffarg colorfile')
