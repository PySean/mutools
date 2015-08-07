#!/usr/bin/env python
"""
    "reorder_all.py", by Sean Soderman
    Reorders all BAM files in a directory according to the ordering of the
    given FASTA file.
"""

import argparse
import os
import subprocess as sp
import sys

parser = argparse.ArgumentParser(description='BAM file re-orderer')
parser.add_argument('-f', '--fasta', type=str, 
                    help=('The FASTA file to' 
                    'be used for reordering'), required=True)

parser.add_argument('-d', '--directory', type=str, 
                    help='Directory of BAM files to reorder', required=True)

parser.add_argument('-p', '--picard_path', type=str, 
                    help=('The full path of the Picard jarfile. Defaults to a file'
                    'named "picard.jar" in the current working directory'), 
                    default='picard.jar')

args = parser.parse_args()

fasta = args.fasta
directory = args.directory

#The command template to be used for invoking picard on all files in the
#directory. Uses double brace trick for later formatting.
cmd = 'java -jar {path} INPUT={{inbam}} REFERENCE={ref} OUTPUT={{outbam}}'

#Validate command line arguments.
if not os.path.exists(args.picard_path):
    sys.stderr.write(('Please specify an existent path to a picard jar file'
                    'or ensure picard.jar is in your working directory.\n'))
    sys.exit(1)
elif not os.path.exists(fasta):
    sys.stderr.write('Please specify a path to an existent FASTA file\n')
    sys.exit(1)
elif not os.path.exists(directory):
    sys.stderr.write(('Please specify a path to an existent directory of BAM'
                      'files.\n'))
    sys.exit(1)

cmd = cmd.format(path=args.picard_path, ref=fasta)

"""
Walk through the directory tree, invoking ReorderSam on each BAM file with
the supplied FASTA file.
"""
def bamwalk(cmd, directory):
    for dirpath, dirnames, filenames in os.walk(directory):
        for f in filenames:
            if f.endswith('.bam'):
                filepath = os.path.join(dirpath, f)
                command = cmd.format(inbam=path, outbam=path)
                try:
                    sp.check_output(command.split())
                except sp.CalledProcessError as cpe:
                    sys.stderr.write('Issue executing ReorderSam: {}'
                                     .format(cpe))
                    sys.exit(1)
bamwalk(cmd, directory)
