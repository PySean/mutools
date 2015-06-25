#!/usr/bin/env python
"""
    'multimutect.py', by Sean Soderman
    Parallelizer for MuTect.
"""
from synchrom import Synchrom
import argparse
import os #For os.access.
import re 
import subprocess
import sys
try:
    from concurrent.futures import ThreadPoolExecutor 
except ImportError as I:
    sys.stderr.write("Please import the required modules: {}\n".format(I))


def thread_run(t_number, t_object):
    pass
###TODO:Remember to use regular expressions to get at the normal/tumor BAM
#file names.

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='MuTect parallelizer')
    parser.add_argument('-m', '--mupath', type=str, default='mutect.jar',
                        help=('The path to the MuTect jar file.'
                              ' Looks for a file named mutect.jar'
                              ' in the current working directory'
                              ' by default'))
    #Create a group for both the file of bamfiles and cmd line
    #list of bamfiles.
    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument('-b', '--bamlistfile', type=str, 
                        help='File containing tumor:normal pairs for MuTect')
    group.add_argument('-p', '--pairs', type=str, nargs='*',
                        help=('List of arguments specifying tumor:normal'
                            ' filename pairs.'))
                                
    parser.add_argument('-M', '--mutectopts', type=str,
                        help='Extra parameters specific to MuTect')
    #NOTE: Does specifying only the long option also give me a value for
    #the short option? other NOTE: yes.
    parser.add_argument('-f', '--fasta', type=str,
                        help='FASTA formatted reference sequence',
                        required=True)
    #Get the dictionary of command line arguments.
    args = parser.parse_args()
    args_dict = dict(args._get_kwargs())
    #Test the Synchrom object. It works perfectly.
    print args_dict
    s = Synchrom(args)
