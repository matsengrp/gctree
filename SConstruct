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
AddOption('--gctree',
           action='store_true',
           help='use gctree inference')
gctree = GetOption('gctree')
AddOption('--igphyml',
           action='store_true',
           help='use igphyml inference')
igphyml = GetOption('igphyml')
AddOption('--outdir',
          type='string',
          help="directory in which to output results")
outdir = GetOption('outdir')


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
              default=None,
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
              help='experimental fastas for comparing summary stats (CFT)')
    experimental_list = GetOption('experimental')
    AddOption('--naiveIDexp',
              type='string',
              default='naive0',
              help='id of naive seq in the experimental data (CFT)')
    naiveIDexp = GetOption('naiveIDexp')

elif inference:
    AddOption('--fasta',
              dest='fasta',
              type='string',
              metavar='PATH',
              help='path to input fasta')
    fasta = GetOption('fasta')
    AddOption('--naiveID',
              type='string',
              metavar='seqID',
              default='naive',
              help='id of naive sequence')
    naiveID = GetOption('naiveID')
    AddOption('--converter',
              type='string',
              default=None,
              help='Converter to convert input fasta format e.g. the Victora lab GC fasta format')
    converter = GetOption('converter')


# First call after all arguments have been parsed
# to enable correct command line help.
if simulate and not GetOption('help'):
    if outdir is None:
        raise InputError('outdir must be specified')
    SConscript('SConscript.simulation',
               exports='env gctree igphyml outdir naive mutability substitution lambda_list lambda0_list n frame N T nsim CommandRunner experimental_list naiveIDexp')
elif inference and not GetOption('help'):
    if None in [fasta, outdir]:
        raise InputError('input fasta and outdir must be specified')
    SConscript('SConscript.inference', exports='env gctree igphyml frame fasta outdir naiveID converter CommandRunner')
