#!/usr/bin/env python
"""
    'multimutect.py', by Sean Soderman
    Parallelizer for MuTect.
"""
from multiprocessing import Lock
from synchrom import Synchrom
from chromolist import UNTOUCHED, BUSY, DONE, ERROR, ChromoList
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
def thread_run(t_number, synchrom, chromolist, init_mutex):
    #Preamble code. One thread initializes the objects.
    counter = 0
    sample_number = 0
    chrom_len = 0 #NOTE TODO: Only getting the length for one thread.
    init_mutex.acquire()
    if len(synchrom.cmd_strings) == 0:
        try:
            synchrom.commands.next()
        except StopIteration as S:
            sys.stderr.write('Error: No commands to be executed: {}'.format(S))
            return None
        chromolist.add_chrom_and_array(synchrom.bam_inputs[0])
        init_mutex.release()

    #Initialize length. Should always work as long as the above if executed.
    chrom_len = len(chromolist.chromosomes[sample_number])
    #The main event. Stop when there's no more chromosomes to process.
    while True:
        if (counter + t_number) < chrom_len:
            counter += t_number
        else:
            chromolist.key_array(sample_number, action='lock')
            counter = chromolist.check_ahead(sample_number, counter)
            if counter is None: #Nothing left to fork for this sample
                #Create the new status array if it does not exist.
                if chromolist.status_arrays.get(sample_number) is None:
                    chromolist.add_chrom_and_array(
                                        synchrom.bam_inputs[sample_number])
                #If everything is finished, log statuses + chromosomes.
                elif chromolist.will_log(sample_number):
                    chromolist.log_status_and_chromosomes(
                                        synchrom.output_dirs[sample_number])
                #Increment sample number, reset counter.
                sample_number += 1
                counter = t_number
                chrom_len = len(chromolist.chromosomes[sample_number])
                chromolist.key_array(sample_number, action='unlock')
                
        chromolist.key_array(sample_number, action='lock')
        chrostr, status = chromolist.get_chrostatus(sample_number, counter)
        if status == UNTOUCHED:
            chromolist.set_chrostatus(sample_number, counter, BUSY)
            the_command = synchrom.cmd_strings[sample_number] % (chrostr, 
                                                                chrostr)
            chromolist.key_array(sample_number, action='unlock')
            #exec_command(the_command)

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
