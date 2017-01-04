#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Infer trees from germinal center data
'''
from __future__ import print_function
import os
import re
import sys
import subprocess
import glob
import shutil
import functools
from os import path
from warnings import warn
#from nestly import Nest
#from nestly.scons import SConsWrap
from SCons.Script import Environment
import json

# Setting up command line arguments/options
AddOption('--fasta',
          dest='fasta',
          type='string',
          nargs=1,
          action='store',
          metavar='PATH',
          help='path to input fasta')
AddOption('--outdir',
          dest='outdir',
          type='string',
          nargs=1,
          action='store',
          metavar='DIR',
          default='output',
          help="directory in which to output results")
AddOption('--germline',
          dest='germline',
          type='string',
          nargs=1,
          action='store',
          metavar='GERMLINE ID',
          default='GL',
          help='id of germline sequence')

# Set up SCons environment
environ = os.environ.copy()
env = Environment(ENV=environ)

# Add stuff to PATH
env.PrependENVPath('PATH', 'bin')

# Some environment sanity checks to make sure we have all prerequisits
def cmd_exists(cmd):
    return subprocess.call("type " + cmd, shell=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE) == 0

# only check for dependencies if we aren't in a dry-run.
if not GetOption('no_exec'):
    msg = ""
    if not cmd_exists('dnapars'):
        msg += '''
               Required dependency command,
               `dnapars` not found on PATH
               Consider using,
                $ module use ~matsengrp/modules
                $ module load phylip
                '''

    # if we are missing any prerequisites, print a message and exit
    if len(msg):
        warn(msg)
        sys.exit(1)

# Grab our cli options
fasta = GetOption('fasta')
outdir_base = GetOption('outdir') # we call this outdir_base in order to not conflict with nestly fns outdir arg
germline = GetOption('germline')
print('fasta = {}'.format(fasta))
print('outdir = {}'.format(outdir_base))
print('germline = {}'.format(germline))

phylip = env.Command(path.join(outdir_base, path.basename(fasta)) + '.phylip',
                     fasta,
                     'python bin/TasParse.py $SOURCE ' + germline + ' > $TARGET')

#
#
# # NOTE: below pilfered from cft
#
# # Initialize nestly!
# # ==================
#
# # This lets us create parameter nests and run the given pipeline for each combination of nesting parameters.
# # It also lets us "pop" out of nesting levels in order to aggregate on nested results.
#
# # Our nesting is going to more or less model the nesting of things in our datapath directory.
# # seed > parameter_set > partition
#
# # Here we initialize nestly, and create an scons wrapper for it
# nest = Nest()
# w = SConsWrap(nest, outdir_base, alias_environment=env)
#
#
# # Here we add some aggregate targets which we'll be build upon in order to spit out our final metadata.json file.
#
# w.add_aggregate('metadata', list)
# w.add_aggregate('svgfiles', list)
#
#
#
# # Initialize seed nest
# # --------------------
#
# # The very first nesting is on the seed id.
# # For the sake of testing, we allow for switching between the full set of seeds, as returned by `seeds_fn`, and a small test set (`test_seeds`).
# # Perhaps we'll hook this up to a command line option at some point...
# # In any case, note that w.add lets you either directly pass a collection, or a function of the control dictionary (`c` in the functions below) for determining what nest parameters to add.
#
# test_seeds = ["QB850.424-Vk", "QB850.043-Vk"]
# def seeds_fn(datapath):
#     return os.listdir(path.join(datapath, 'seeds'))
# # You can use this to toggle a small testset vs the entire data set
# seeds = test_seeds
# #seeds = seeds_fn(datapath)
#
# # Initialize our first nest level
# w.add('seed', seeds)
#
#
# # Initialize parameter set nest
# # -----------------------------
#
# # Next we nest on parameter sets; I believe these correspond to sequencing runs.
# # They're the things that look like "Hs-LN2-5RACE-IgG", and are nested within each seed's output directory.
#
# w.add('parameter_set', lambda c: map(lambda x: path.basename(x), glob.glob(path.join(datapath, 'seeds', c['seed'], "*"))))
#
#
# # Some helpers at the seed level
#
# def input_dir(c):
#     "Return the `seed > parameter_set` directory given the closed over datapath and a nestly control dictionary"
#     # This 'seeds' thing here is potentially not a very general assumption; not sure how variable that might be upstream
#     return path.join(datapath, 'seeds', c['seed'], c['parameter_set'])
#
# def path_base_root(full_path):
#     return path.splitext(path.basename(full_path))[0]
#
# def valid_partition(fname):
#     with open(fname, 'r') as partition_handle:
#         partition = csv.DictReader(partition_handle).next()
#         unique_ids = partition['unique_ids'].split(':')
#         return len(unique_ids) >= 2 # Do we want >2 or >=2?
#
#
# def annotations(c):
#     """Return the annotations file for a given control dictionary, sans any partitions which don't have enough sequences
#     for actual analysis."""
#     return map(path_base_root,
#                 # We might eventually want to handle small partitions more manually, so there's data on the other end, but for now, filter is easiest
#                 filter(valid_partition,
#                         glob.glob(path.join(input_dir(c), "*-plus-*.csv"))))
#
# def partitions(c):
#     "Returns the `partition.csv` file path with all the partition information for every partition output by partis."
#     return path.join(input_dir(c), 'partition.csv')
#
# def parameter_dir(c):
#     "The input parameter directory (as in input to partis); see --parameter-dir flag to various partis cli calls."
#     return path.join(datapath, c['parameter_set'])
#
#
#
# # Each c["partition"] value actually points to the annotations for that partition... a little weird but...
# w.add('partition', annotations)
#
#
# def annotation(c):
#     return path.join(input_dir(c), c['partition'] + '.csv')
#
# @w.add_target()
# def process_partis(outdir, c):
#     return env.Command(path.join(outdir, 'metadata.json'),
#             [partitions(c), annotation(c)],
#             'process_partis.py ' +
#                 '--partition ${SOURCES[0]} ' +
#                 '--annotations ${SOURCES[1]} ' +
#                 #'--param_dir ' + parameter_dir(c) + ' '
#                 '--partis_log ' + path.join(input_dir(c), 'partition.log ') +
#                 '--cluster_base cluster ' +
#                 '--output_dir ' + outdir)
#
# # To get the input sequences we need to read in the file output by process_partis to figure out what the input fasta was,
# # and then for convenience copy it over locally
#
# def extract_inseqs(outdir, target, source, env):
#     with open(str(source[0]), "r") as source_handle:
#         metadata = json.load(source_handle)[0]
#         shutil.copy(path.join(outdir, metadata['file']), str(target[0]))
#
# @w.add_target()
# def inseqs(outdir, c):
#     return env.Command(
#         path.join(outdir, "inseqs.fasta"),
#         c['process_partis'],
#         functools.partial(extract_inseqs, outdir))
#
# # Forget why we do this; Chris W. seems to think it may not have been as necessary as original thought.
# @w.add_target()
# def padded_seqs(outdir, c):
#     return env.Command(
#         path.join(outdir, "padded_seqs.fa"),
#         c['inseqs'],
#         "python bin/padseq.py $SOURCE > $TARGET")
#
# # use fasttree to make newick tree from sequences
# @w.add_target()
# def fasttree(outdir, c):
#     return env.SRun(
#         path.join(outdir, "fasttree.nwk"),
#         c['padded_seqs'],
#         "FastTree -nt -quiet $SOURCE > $TARGET")
#
# # calculate list of sequences to be pruned
# @w.add_target()
# def pruned_ids(outdir, c):
#     return env.Command(
#         path.join(outdir, "pruned_ids.txt"),
#         c['fasttree'],
#         "python bin/prune.py --seed " + c['seed'] + " $SOURCE > $TARGET")
#
# # prune out sequences to reduce taxa
# @w.add_target()
# def pruned_seqs(outdir, c):
#     return env.Command(
#         path.join(outdir, "pruned.fa"),
#         [c['pruned_ids'], c['padded_seqs']],
#         "seqmagick convert --include-from-file $SOURCES $TARGET")
#
# # Convert to phylip for dnaml
# @w.add_target()
# def phy(outdir, c):
#     return env.Command(
#         path.join(outdir, "pruned.phy"),
#         c['pruned_seqs'],
#         "seqmagick convert $SOURCE $TARGET")
#
# # Create a config file for dnaml (a persnickety old program with interactive menues...)
# @w.add_target()
# def dnaml_config(outdir, c):
#     return env.Command(
#         path.join(outdir, "dnaml.cfg"),
#         c['phy'],
#         "python bin/mkconfig.py $SOURCE > $TARGET")
#
# # Run dnaml by passing in the "config file" as stdin hoping the menues all stay sane
# # (Aside: There is a program for programatically responding to interactive menus if this gets any hairier)
# @w.add_target()
# def dnaml(outdir, c):
#     "run dnaml (from phylip package) to create tree with inferred sequences at internal nodes"
#     tgt = env.SRun(
#         map(lambda x: path.join(outdir, x), ["outtree", "outfile", "dnaml.log"]),
#         c['dnaml_config'],
#         'cd ' + outdir + ' && dnaml < ${SOURCE.file} > ${TARGETS[2].file}')
#     # Manually depend on phy so that we rerun dnaml if the input sequences change (without this, dnaml will
#     # only get rerun if one of the targets are removed or if the iput dnaml_config file is changed).
#     env.Depends(tgt, c['phy'])
#     return tgt
#
#
# @w.add_target()
# def dnaml_tree(outdir, c):
#     """parse dnaml output into fasta and newick files, and make SVG format tree with ETE package.
#     xvfb-run is needed because of issue https://github.com/etetoolkit/ete/issues/101"""
#     tgt = env.Command(
#             map(lambda x: path.join(outdir, x),
#                 ["dnaml.svg", "dnaml.fa", "dnaml.seedLineage.fa", "dnaml.newick"]),
#             c['dnaml'],
#             # Note: the `-` prefix here tells scons to keep going if this command fails.
#             "- xvfb-run -a bin/dnaml2tree.py --seed " + c['seed'] + " --dnaml ${SOURCES[1]} --outdir ${TARGETS[0].dir} --basename dnaml")
#     # Manually depend on dnaml2tree.py script, since it doesn't fall in the first position within the command
#     # string.
#     env.Depends(tgt, 'bin/dnaml2tree.py')
#     # Do aggregate work
#     c['svgfiles'].append(tgt[0])
#     return tgt
#
# # Might want to switch back to this approach...
# #def extended_metadata(metadata, dnaml_tree_tgt):
#     #"Add dnaml_tree target(s) as metadata to the given metadata dict; used prior to metadata write."
#     #with open(
#     #m = copy.copy(metadata)
#     #m['svg'] = path.relpath(str(dnaml_tree_tgt[0]), outdir_base)
#     #m['fasta'] = path.relpath(str(dnaml_tree_tgt[1]), outdir_base)
#     #m['seedlineage'] = path.relpath(str(dnaml_tree_tgt[2]), outdir_base)
#     #m['newick'] = path.relpath(str(dnaml_tree_tgt[3]), outdir_base)
#     #del m['file']
#     #return m
#
# @w.add_target()
# def cluster_metadata(outdir, c):
#     def dnaml_tgt_relpath(i):
#         return path.relpath(str(c['dnaml_tree'][i]), outdir_base)
#     tgt = env.Command(
#             path.join(outdir, 'extended_metadata.json'),
#             c['process_partis'] + c['dnaml_tree'],
#             # Don't like the order assumptions on dnaml_tgts here...
#             'json_assoc.py $SOURCE $TARGET ' +
#                 # Not 100% on this; Pick optimal attr/key name
#                 #'best_partition ' + c['partition'] + ' ' +
#                 'clustering_step ' + re.compile('run-viterbi-best-plus-(?P<step>.*)').match(c['partition']).group('step') + ' ' +
#                 'svg ' + dnaml_tgt_relpath(0) + ' ' +
#                 'fasta ' + dnaml_tgt_relpath(1) + ' ' +
#                 'seedlineage ' + dnaml_tgt_relpath(2) + ' ' +
#                 'newick ' + dnaml_tgt_relpath(3) + ' ')
#     # Note; we used to delete the 'file' attribute as well; not sure why or if that's necessary
#     c['metadata'].append(tgt)
#     return tgt
#
#
# # Popping out
# # -----------
#
# # Here we pop out to the "seed" level so we can aggregate our metadata
#
# w.pop('seed')
#
#
# # Filtering out bad clusters:
#
# # As mentioned above, dnaml2tree fails for some clusters, so we'd like to filter these clusters out of the final metadata results.
# # We do this based on whether the svg targets of the dnaml2tree command point to actual files or not.
# # Eventually it would be nice to let them pass through with information indicating there was a failure, and handle appropriately in cftweb.
# # We may also want to handle clusters that are too small similarly, but for now we're filtering them out at the very beginning of the pipeline.
#
# # For now, to the stated end, we have this filter predicate
# def tgt_exists(tgt):
#     p = str(tgt)
#     exists = path.exists(p)
#     if not exists:
#         print("Path doesn't exist:", p)
#     return exists
#
# def in_pairs(xs):
#     it = iter(xs)
#     for x in it:
#         yield x, next(it)
#
# def node_metadata(node):
#     with open(str(node), 'r') as handle:
#         return json.load(handle)
#
# def write_metadata(target, source, env):
#     # Here's where we filter out the seeds with clusters that didn't compute through dnaml2tree.py
#     good_metadata = map(lambda x: node_metadata(x[1]), filter(lambda x: tgt_exists(x[0]), in_pairs(source)))
#     with open(str(target[0]), "w") as fh:
#         json.dump(good_metadata, fh, sort_keys=True,
#                        indent=4, separators=(',', ': '))
#
# @w.add_target()
# def metadata(outdir, c):
#     tgt = env.Command(
#         path.join(outdir, "metadata.json"),
#         zip(c['svgfiles'], c['metadata']),
#         write_metadata)
#     env.AlwaysBuild(tgt)
#     return tgt
#
# import textwrap
#
# def print_hints(target, source, env):
#     msg = """\
# 		hint: to run the cft web interface,
# 			$ cd cftweb && python -m cftweb --file {}
# 	""".format(os.path.abspath(str(source[0])))
#     print(textwrap.dedent(msg))
#
# @w.add_target()
# def hints(outdir, c):
#     hints = env.Command(
#         None,
#         c['metadata'],
#         print_hints)
#     env.AlwaysBuild(hints)
#     return hints
