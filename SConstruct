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
AddOption('--simulate',
          action='store_true',
          help='validation subprogram, instead of inference')
AddOption('--frame',
          type='int',
          default=None,
          help='codon frame')
frame = GetOption('frame')
AddOption('--outdir',
          type='string',
          help="directory in which to output results")
outdir = GetOption('outdir')

if GetOption('simulate'):
    AddOption('--naive',
              type='string',
              default='ggacctagcctcgtgaaaccttctcagactctgtccctcacctgttctgtcactg'+
                      'gcgactccatcaccagtggttactggaactggatccggaaattcccagggaataa'+
                      'acttgagtacatggggtacataagctacagtggtagcacttactacaatccatct'+
                      'ctcaaaagtcgaatctccatcactcgagacacatccaagaaccagtactacctgc'+
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
              default=.9,
              help='Poisson branching parameter for simulation')
    AddOption('--lambda0',
              type='float',
              default=.1,
              help='baseline mutation rate')
    AddOption('--r',
              type='float',
              default=.5,
              help='sampling probability')
    AddOption('--n',
              type='int',
              default=100,
              help='minimum simulation tree size')
    AddOption('--T',
              type='int',
              default=None,
              help='observation time')

    naive = GetOption('naive')
    mutability = GetOption('mutability')
    substitution = GetOption('substitution')
    lambda_ = GetOption('lambda')
    lambda0 = GetOption('lambda0')
    r = GetOption('r')
    n = GetOption('n')
    T = GetOption('T')
    SConscript('SConscript.simulation',
                exports='env outdir naive mutability substitution lambda_ lambda0 r n frame T')

else:
    AddOption('--fasta',
              dest='fasta',
              type='string',
              metavar='PATH',
              help='path to input fasta')
    AddOption('--naiveID',
              type='string',
              metavar='seqID',
              default='naive',
              help='id of naive sequence')

    fasta = GetOption('fasta')
    naiveID = GetOption('naiveID')
    SConscript('SConscript.inference', exports='env frame fasta outdir naiveID')
