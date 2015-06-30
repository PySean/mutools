#!/usr/bin/env python
"""
    'multimutect.py', by Sean Soderman
    Parallelizer for MuTect.
"""
from synchrom import Synchrom
import argparse
import os #For os.path.exists() to check for file existence
import re 
import subprocess
import sys
try:
    from concurrent.futures import ThreadPoolExecutor 
except ImportError as I:
    sys.stderr.write("Please import the required modules: {}\n".format(I))


"""
The big deal function of the entire program.
Has the thread fork a MuTect process and wait for it to complete.
Once this occurs, the thread then acquires a lock on the current
status array and assigns the appropriate exit status of the process
to it. Once a thread returns when all elements of the status
array are DONE, ERROR, or BUSY, the thread creates the data structures
and directory for the next input tumor:normal pair, then goes to
work on the <t_number>'th chromosome.
"""
#TODO: Start fresh.
def thread_run(t_number, synchrom, chromolist):
    pass

"""
Diagnostic function for the program.
Tests whether everything is working as it should, utilising prints.
Takes a Synchrom object and returns None.
"""
def diagnostic(synchrom):
    for item in synchrom.commands:
        dirname = synchrom.output_dirs[-1]
        print(('My command is {}\n'
               ' The directory I am writing to is {}\n'
               #' The filenames will be named after chromosomes: {}\n'
              ).format(item, dirname))


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
                                
    parser.add_argument('-M', '--mutectopts', type=str, default='',
                        help='Extra parameters specific to MuTect')
    parser.add_argument('-f', '--fasta', type=str,
                        help='FASTA formatted reference sequence',
                        required=True)
    parser.add_argument('-i', '--inputdir', type=str, default=os.getcwd(),
                        help=('The name of the directory the input files'
                              ' are located. Default: working directory.'))
    parser.add_argument('-o', '--outputdir', type=str, default='output',
                        help=('The name of the directory the output should go'
                              ' to. Default: a directory called "output"'))
    #Get the dictionary of command line arguments.
    args = parser.parse_args()
    args_dict = dict(args._get_kwargs())
    diagnostic(Synchrom(args))
    #Create the threads, then get the party started!
