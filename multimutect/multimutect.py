#!/usr/bin/env python
"""
    'multimutect.py', by Sean Soderman
    Parallelizer for MuTect.
"""
#from multiprocessing import Lock, cpu_count
from synchrom import Synchrom
from chromolist import UNTOUCHED, BUSY, DONE, ERROR, ChromoList
import argparse
import multiprocessing
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
#TODO FIXME NOTE: Chromosomes 19, 21, 22, and X were omitted.
def thread_run(t_number, num_t, synchrom, chromolist, init_mutex):
    counter = t_number
    sample_number = 0
    chrom_len = 0
    ##Critical section no. 1##
    #Create initial data structures and directory for 1st sample pair.
    init_mutex.acquire()
    #The memory is indeed, shared. (or is it?)
    if len(synchrom.cmd_strings) == 0:
        try:
            synchrom.commands.next()
        except StopIteration as S:
            sys.stderr.write('Error: No commands to be executed: {}'.format(S))
            return None
        chromolist.add_chrom_and_array(synchrom.bam_inputs[0])
        try:
            os.mkdir(synchrom.output_dirs[0])
        except OSError as O:
            sys.stderr.write('Could not make directory: {}'.format(O))
            return None
    init_mutex.release()
    #Initialize length. Should always work as long as the above if executed.
    chrom_len = len(chromolist.chromosomes[sample_number])
    sys.stderr.write('Thread {} sees this array: {} and these chromosomes: {}\n'
                     .format(t_number, 
                             chromolist.status_arrays[sample_number].get_obj()[:],
                             chromolist.chromosomes[sample_number]))
    #The main event. Stop when there's no more chromosomes to process.
    #TODO: Counter never changes for thread #0!
    while True:
        #if (counter + num_t) < chrom_len:
        #    counter += num_t
        if counter >= chrom_len:
            sys.stderr.write('Do I get here?\n')
            ##Critical section no. 2##
            #Check statuses, create missing data structures.
            chromolist.key_array(sample_number, action='lock')
            counter = chromolist.check_ahead(sample_number, counter)
            if counter is None: #Nothing left to fork for this sample
                #Create the new status array if it does not exist.
                #Also get the next command and make the corresponding
                #directory.
                if chromolist.status_arrays.get(sample_number + 1) is None:
                    try:
                        synchrom.commands.next()
                    except StopIteration as S:
                        return DONE
                    try:
                        os.mkdir(synchrom.output_dirs[sample_number + 1])
                    except OSError as O:
                        sys.stderr.write('Could not make directory: {}'
                                        .format(O))
                    chromolist.add_chrom_and_array(
                                        synchrom.bam_inputs[sample_number + 1])
                #will_log is only true when all threads are done with a sample.
                #This is why I log here.
                elif chromolist.will_log(sample_number):
                    chromolist.log_status_and_chromosomes(
                                        synchrom.output_dirs[sample_number])
                #This thread is moving on, so do its variables.
                sample_number += 1
                counter = t_number
                chrom_len = len(chromolist.chromosomes[sample_number])
            chromolist.key_array(sample_number, action='unlock')
                
        ##Critical section no. 3##
        #Getting/Setting status then forking.
        #sys.stderr.write('I am thread {}\n'.format(t_number))
        #FIXME: This is where all the other threads get stuck.
        chromolist.key_array(sample_number, action='lock')
        chrostr, status = chromolist.get_chrostatus(sample_number, counter)
        exit_status = DONE
        if status == UNTOUCHED:
            chromolist.set_chrostatus(sample_number, counter, BUSY)
            the_command = synchrom.cmd_strings[sample_number] % (chrostr, 
                                                                chrostr)
            chromolist.key_array(sample_number, action='unlock')
            ##TODO: Check what MuTect returns on error.
            #ACTUALLY unnecessary as CalledProcessError exists.
            try:
                subprocess.check_output(the_command.split())
            except subprocess.CalledProcessError as cpe:
                exit_status = ERROR
            chromolist.key_array(sample_number, action='lock')
            chromolist.set_chrostatus(sample_number, counter, exit_status)
            chromolist.key_array(sample_number, action='unlock')
        counter += num_t

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
    args = parser.parse_args()
    if not os.path.exists(args.mupath):
        sys.stderr.write('Error: path to {} does not exist\n'
                         .format(args.mupath))
        sys.exit(1)
    #diagnostic(Synchrom(args))
    #Create the threads and the parent output directory.
    numthreads = multiprocessing.cpu_count() // 2
    os.mkdir(args.outputdir)
    with ThreadPoolExecutor(max_workers=numthreads) as threader:
        synchrom = Synchrom(args)
        chromolist = ChromoList()
        init_mutex = multiprocessing.Lock()
        jobs = {threader.submit(thread_run, i, numthreads, synchrom, 
                                chromolist, init_mutex): i 
                                for i in range(0, numthreads)}
        for future in concurrent.futures.as_completed(jobs.keys()):
            threadnum = jobs[future]
            result = 0
            try:
                result = future.result()
            except Exception as E:
                print('Thread {} terminated abnormally: {}'.format(threadnum,
                                                                   E))
