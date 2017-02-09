#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Prepare a configuration file for ancestral state reconstruction under the HLP16 model.
"""

import re
import os
import argparse


class FastaInputError(Exception):
    '''When the fasta file in not reflecting amino acid DNA coding for protein.'''


def which(executable):
    for path in os.environ["PATH"].split(os.pathsep):
        if os.path.exists(os.path.join(path, executable)):
                return os.path.realpath(os.path.join(path, executable))
    return None


def main():

    parser = argparse.ArgumentParser(description='Prepare config file for ASR under the HLP16 model.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--template', help='Config file template', type=str, metavar='FILE', required=True)
    parser.add_argument('--outfile', help='Name of the config file to write.', type=str, metavar='FILENAME', required=True)
    parser.add_argument('--igphyml_exe', help='IgPhyML executable. Will search for the default localt like UNIX `which`.', type=str, required=True)
    parser.add_argument('--model', help='Which model to run? [gy94, hlp16]', type=str, required=True)
    parser.add_argument('--fasta_file', help='To find the length of the amino acid sequence.', type=str, required=True)
    args = parser.parse_args()

    igphyml_path = which(args.igphyml_exe.rstrip('/'))
    assert(igphyml_path is not None)
    IGPHYML_DIR = re.sub(r'/src/\w+$', '', igphyml_path)
    MODEL = args.model
    assert(MODEL in ['gy94', 'hlp16'])
    # Find the length of the translated DNA sequence:
    LEN_AA = None
    with open(args.fasta_file) as fh:
        for l in fh:
            if not l.startswith('>') and l != '':
                this_LEN_AA = len(l.strip()) / 3.0
                if int(this_LEN_AA) != this_LEN_AA or (LEN_AA is not None and this_LEN_AA != LEN_AA):
                    raise FastaInputError('Problem with the input fasta file. Either is the not a multiple of three or it has indels.')
                elif LEN_AA is None:
                    LEN_AA = int(this_LEN_AA)
    assert(LEN_AA is not None)

    # Replace the keyword in the template file:
    with open(args.template) as fh:
        template = fh.read()

    template = template.replace('LEN_AA', str(LEN_AA))
    template = template.replace('IGPHYML_DIR', IGPHYML_DIR)
    template = template.replace('MODEL', MODEL)

    # Write the new config file:
    with open(args.outfile, 'w') as fh_out:
        fh_out.write(template)


if __name__ == "__main__":
    main()
