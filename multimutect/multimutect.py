#!/usr/bin/env python
"""
    'multimutect.py', by Sean Soderman
    Parallelizer for MuTect.
"""
from chromolist import UNTOUCHED, BUSY, DONE, ERROR, ChromoList
from itertools import izip
from synchrom import Synchrom
import argparse
import multiprocessing
import os #For os.path.exists() to check for file existence and os.walk.
import re 
import subprocess
import sys
try:
    from concurrent.futures import ThreadPoolExecutor, as_completed
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
def thread_run(t_number, num_t, synchrom, chromolist, init_mutex):
    counter = t_number
    sample_number = 0
    chrom_len = 0
    ##Critical section no. 1##
    #Create initial data structures and directory for 1st sample pair.
    init_mutex.acquire()
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
    #The main event. Stop when there's no more chromosomes to process.
    while True:
        if counter >= chrom_len:
            counter = counter - t_number
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
                        chromolist.key_array(sample_number, action='unlock')
                        return DONE
                    try:
                        os.mkdir(synchrom.output_dirs[sample_number + 1])
                    except OSError as O:
                        sys.stderr.write('Could not make directory: {}'
                                        .format(O))
                    chromolist.add_chrom_and_array(
                                        synchrom.bam_inputs[sample_number + 1])

                chromolist.key_array(sample_number, action='unlock')
                #This thread is moving on, so do its variables.
                sample_number += 1
                counter = t_number
                chrom_len = len(chromolist.chromosomes[sample_number])
            #This happens either when my counter is none or not.
            else:
                chromolist.key_array(sample_number, action='unlock')
                
        ##Critical section no. 3##
        #Getting/Setting status then forking.
        chromolist.key_array(sample_number, action='lock')
        chrostr, status = chromolist.get_chrostatus(sample_number, counter)
        exit_status = DONE
        if status == UNTOUCHED:
            chromolist.set_chrostatus(sample_number, counter, BUSY)
            the_command = synchrom.cmd_strings[sample_number] % (chrostr, 
                                                            chrostr + '.vcf')
            chromolist.key_array(sample_number, action='unlock')
            try:
                subprocess.check_output(the_command.split())
            except subprocess.CalledProcessError as cpe:
                exit_status = ERROR
            chromolist.key_array(sample_number, action='lock')
            chromolist.set_chrostatus(sample_number, counter, exit_status)
            #Are all threads completely finished with this sample?
            if chromolist.will_log(sample_number):
                chromolist.log_status_and_chromosomes(
                                        synchrom.output_dirs[sample_number],
                                        sample_number)
            chromolist.key_array(sample_number, action='unlock')
        #NOTE: Without this, there was potential for a double lock on the same array.
        else:
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
    #Create a group for both the file of bamfiles and cmd line
    #list of bamfiles.
    file_group = parser.add_mutually_exclusive_group(required=True)

    file_group.add_argument('-b', '--bamlistfile', type=str, 
                        help='File containing tumor:normal pairs for MuTect')

    file_group.add_argument('-p', '--pairs', type=str, nargs='*',
                        help=('List of arguments specifying tumor:normal'
                            ' filename pairs.'))

    parser.add_argument('-f', '--fasta', type=str,
                        help='FASTA formatted reference sequence',
                        required=True)

    #Group for either conf file for MuTect specific args or
    #command line for them.
    cmd_group = parser.add_mutually_exclusive_group(required=False)

    parser.add_argument('-m', '--mupath', type=str, default='mutect.jar',
                        help=('The path to the MuTect jar file.'
                              ' Looks for a file named mutect.jar'
                              ' in the current working directory'
                              ' by default'))

    cmd_group.add_argument('-M', '--mutectopts', type=str, default='',
                        help='Extra parameters specific to MuTect')

    cmd_group.add_argument('-c', '--conf', type=str, default='',
                            help=('File containing extra parameters'
                                  ' specific to MuTect'))

    parser.add_argument('-i', '--inputdir', type=str, default=os.getcwd(),
                        help=('The name of the directory the input files'
                              ' are located. Default: working directory.'))

    parser.add_argument('-o', '--outputdir', type=str, default='output',
                        help=('The name of the directory the output should go'
                              ' to. Default: a directory called "output"'))

    parser.add_argument('--numthreads', type=int, 
                        default=multiprocessing.cpu_count() // 4,
                        help=('The number of threads that will fork mutect'
                              ' processes. Default: The # of cores on your'
                              ' computer / 4, rounded down.'))

    parser.add_argument('--mem', type=int, default=2,
                        help=('The max amount of memory each forked MuTect'
                              ' process can allocate on the Java heap'
                              ' Default: 2'))

    parser.add_argument('--process_whole_bam', action='store_true',
                        help=('Process the entire BAM file at once instead '
                              'of single chromosomes at a time'))
    args = parser.parse_args()
    if not os.path.exists(args.mupath):
        sys.stderr.write('Error: path to {} does not exist\n'
                         .format(args.mupath))
        sys.exit(1)
    #diagnostic(Synchrom(args))
    #Create the threads and the parent output directory.
    numthreads = args.numthreads
    os.mkdir(args.outputdir)
    if args.process_whole_bam != True:
        #Get the the threads going. Output debug info to thread-specific files.
        with ThreadPoolExecutor(max_workers=numthreads) as threader:
            synchrom = Synchrom(args)
            chromolist = ChromoList()
            init_mutex = multiprocessing.Lock()
            jobs = {threader.submit(thread_run, i, numthreads, synchrom, 
                                    chromolist, init_mutex): i 
                                    for i in range(0, numthreads)}
            for future in as_completed(jobs.keys()):
                threadnum = jobs[future]
                result = 0
                try:
                    result = future.result()
                except Exception as E:
                    print('Thread {} terminated abnormally: {}'
                          .format(threadnum,E))

                with open('thread{}.status'.format(threadnum), 'w') as ts:
                    if result == DONE:
                        ts.write('I completed successfully!\n')
                    else:
                        ts.write('I completed unsuccessfully!\n')
    #Simpler forking for whole BAM file processing.
    else:
        #Mini function: execute the command, surround in try except.
        def procfun(dtuple):
            tid, cmd = dtuple
            try:
                cmdlist = cmd.split()
                val = subprocess.check_output(cmdlist)
                print('tid: {}, the cmd is: {}'.format(tid, cmd))
            except subprocess.CalledProcessError as cpe:
                #Log error to a file rather than write to stderr.
                with open('thread{}.err'.format(tid)) as errfile:
                    errfile.write(('I crashed with the command line:{}'
                                   ' {}.{} You may need to use the default '
                                   ' option for individual'
                                   ' chromosome processing option instead, '
                                   ' or a manual command line '
                                   ' accomodating more '
                                   ' memory for the Java heap.'
                                   ).format(os.linesep, cmd, os.linesep))
                return 'Thread {} executed unsuccessfully'.format(tid)
            return 'Thread {} executed successfully'.format(tid)

        #Mini function #2: Generator function for infinite numeric sequence.
        def infinigen():
            i = 0
            while True:
                yield i
                i += 1

        synchrom = Synchrom(args)
        infinity = infinigen()
        pool = multiprocessing.Pool(processes=numthreads)
        pids = 0
        with ThreadPoolExecutor(max_workers=numthreads) as threader:
            results = threader.map(procfun, izip(infinity, synchrom.commands))
            for i in results:
                print i
