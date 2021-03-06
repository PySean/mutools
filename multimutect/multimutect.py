#!/usr/bin/env python
"""
    'multimutect.py', by Sean Soderman
    Parallelizer for MuTect.
"""
from itertools import izip
from synchrom import Synchrom
from time import time
import argparse
import multiprocessing
import os
import re
import subprocess
import sys
try:
    from concurrent.futures import ThreadPoolExecutor
except ImportError as I:
    sys.stderr.write("Please import the required modules: {}\n".format(I))


"""
Diagnostic function for the program.
Tests whether everything is working as it should, utilising prints.
Takes a Synchrom object and returns None.
"""
def diagnostic(synchrom):
    for item in synchrom.commands:
        print('My command is {}'.format(item))


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

    parser.add_argument('--mem', type=int, default=3,
                        help=('The max amount of memory each forked MuTect'
                              ' process can allocate on the Java heap'
                              ' Default: 2'))

    parser.add_argument('--process_whole_bam', action='store_true',
                        help=('Process the entire BAM file at once instead '
                              'of single chromosomes at a time'))
    parser.add_argument('--statistics', type=str,
                        help=('Report statistics on execution time and '
                               ' threads used.'))
    args = parser.parse_args()
    if not os.path.exists(args.mupath):
        sys.stderr.write('Error: path to {} does not exist. cwd: {}\n'
                         .format(args.mupath, os.getcwd()))
        sys.exit(1)
    #Create the threads and the parent output directory.
    numthreads = args.numthreads
    #Mini function: execute the command, surround in try except.
    def procfun(dtuple):
        tid, cmd = dtuple
        try:
            cmdlist = cmd.split()
            val = subprocess.check_output(cmdlist)
            print('tid: {}, the cmd is: {}'.format(tid, cmd))
        except subprocess.CalledProcessError as cpe:
            errfilepath = ''
            if not os.path.exists('errors'):
                try:
                    os.mkdir('errors')
                except OSError as O:
                    sys.stderr.write("Couldn't make directory 'errors':{}{}".
                                     format(O, os.linesep))
            errfilepath = os.path.join('errors', 
                                       'thread{}.err'.format(tid))
            #Log error to a file rather than write to stderr.
            with open(errfilepath, 'w') as errfile:
                errfile.write(('I crashed with the command line:{}'
                               ' {}.{} You may need to use the default '
                               ' option for individual'
                               ' chromosome processing instead, '
                               ' or a command line '
                               ' accomodating more '
                               ' memory for the Java heap.'
                               ' The specific problem was {}\n'
                               ).format(os.linesep, cmd, os.linesep, cpe))
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
    start_time = 0
    end_time = 0
    with ThreadPoolExecutor(max_workers=numthreads) as threader:
        results = threader.map(procfun, izip(infinity, synchrom.commands))
        start_time = time()
        for i in results:
            print i
        end_time = time()
    if args.statistics is not None:
        statfile = args.statistics
        bam_gigs = 0
        cpu_cores = multiprocessing.cpu_count()
        #Gather data for initial run.
        if not os.path.exists(statfile):
            bams = []
            if args.bamlistfile is not None:
                with open(args.bamlistfile, 'r') as blf:
                    bams = [re.split('\s+', b.strip()) for b in blf
                            if re.search('.*bam', b)]
            else:
                bams = [b.split(':') for b in args.pairs]
            #Flatten the list and remove empty strings.
            bams = [i for pair in bams for i in pair if i != '']
            #Prepend input directory name to each bam filename in the list.
            bams = [os.path.join(args.inputdir, b) for b in bams]
            #Attain the size (in bytes) of the processed BAM data.
            bam_gigs = sum([os.stat(b).st_size for b in bams])
        #stats.txt is opened in append mode, as it will take mult.
        #runs to get data for thread performance.
        with open(statfile, 'a') as filestats:
            #Initialize the file if it is of size zero.
            if os.stat(statfile).st_size == 0:
                filestats.write('CPU cores: {}\n'.format(cpu_cores))
                filestats.write('Total BAM data processed: {}\n'
                                .format(bam_gigs))
                filestats.write('Threads\tTime\n')
            thread_and_time = '{}\t{}\n'.format(numthreads, end_time - start_time)
            filestats.write(thread_and_time)
