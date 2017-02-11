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
AddOption('--simulate',
          action='store_true',
          help='validation subprogram, instead of inference')
AddOption('--threads',
          type='int',
          default=1,
          help='Number of cores to run on.')
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
if not gctree and not igphyml:
    raise InputError('must set at least one inference method')

AddOption('--outdir',
          type='string',
          help="directory in which to output results")
outdir = GetOption('outdir')


class InputError(Exception):
    """Exception raised for errors in the input."""

if not GetOption('simulate') and not GetOption('inference'):
    raise InputError('Please provide one of the required arguments. Either "--inference" or "--simulate".'
                     'Command line help can then be evoked by "-h" or "--help" and found in the bottom'
                     'of the output under "Local Options".')

if GetOption('simulate'):
    AddOption('--naive',
              type='string',
              default='ggacctagcctcgtgaaaccttctcagactctgtccctcacctgttctgtcactg'
                      'gcgactccatcaccagtggttactggaactggatccggaaattcccagggaataa'
                      'acttgagtacatggggtacataagctacagtggtagcacttactacaatccatct'
                      'ctcaaaagtcgaatctccatcactcgagacacatccaagaaccagtactacctgc'
                      'agttgaattctgtgactactgaggacacagccacatattactgt',
              help='sequence of naive from which to simulate')
    AddOption('--mutability',
              type='string',
              metavar='PATH',
              default='S5F/Mutability.csv',
              help='path to S5F mutability data')
    AddOption('--substitution',
              type='string',
              metavar='PATH',
              default='S5F/Substitution.csv',
              help='path to S5F substitution data')
    AddOption('--lambda',
              type='float',
              default=None,
              help='Poisson branching parameter for simulation')
    AddOption('--lambda0',
              type='float',
              default=None,
              help='baseline mutation rate')
    AddOption('--r',
              type='float',
              default=None,
              help='sampling probability')
    AddOption('--N',
              type='int',
              default=None,
              help='simulation size (number of cells observerved)')
    AddOption('--T',
              type='int',
              default=None,
              help='observation time')
    AddOption('--n',
              type='int',
              default=10,
              help='number of simulations with each parameter parameter choice')

    naive = GetOption('naive')
    mutability = GetOption('mutability')
    substitution = GetOption('substitution')
    lambda_ = GetOption('lambda')
    lambda0 = GetOption('lambda0')
    r = GetOption('r')
    N = GetOption('N')
    T = GetOption('T')
    n = GetOption('n')
elif GetOption('inference'):
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


# First call after all arguments have been parsed
# to enable correct command line help.
if GetOption('simulate') and not GetOption('help'):
    if outdir is None:
        raise InputError('outdir must be specified')
    SConscript('SConscript.simulation',
               exports='env gctree igphyml outdir naive mutability substitution lambda_ lambda0 r frame N T n')
elif GetOption('inference') and not GetOption('help'):
    if None in [fasta, outdir]:
        raise InputError('input fasta and outdir must be specified')
    SConscript('SConscript.inference', exports='env gctree igphyml frame fasta outdir naiveID')
