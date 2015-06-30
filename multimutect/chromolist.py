#!/usr/bin/python
"""
    "chromolist.py", by Sean Soderman
    An object that tightly couples the tuple of chromosomes with the
    list of status arrays.
"""
import sys
from pysam import AlignmentFile
from multiprocessing import Array
"""
Function decorator for the methods get_chrostatus and set_chrostatus.
Catches an IndexError and prints out relevant debug information.
"""
def catch_index_error(func):
    def wrap(self, *args, **kwargs):
        try:
            return func(self, *args, **kwargs)
        except IndexError as I:
            #The first argument is the sample number, second is the
            #chromosome index.
            array_number, cell_and_chromo = args[0:2]
            chromo_str = self.chromosomes[cell_and_chromo]
            sys.stderr.write(('An attempt to access Status array no. {array}'
                              ' status cell/chromosome no. {cell} and' 
                              ' chromosome {chr} was made in function {func}:'
                              '{I}\n'
                             ).format(array = array_number, 
                                      cell = cell_and_chromo,
                                      chr = chromo_str, 
                                      func=func.func_name, I=I))
    return wrap

class ChromoList():
    """
    A class that tightly couples a tuple of chromosome strings
    with a list of arrays containing statuses corresponding 
    to each chromosome.
    """

    """
    Status constants to be used when having a thread assign a status
    to an element within a status array.
    """
    #NOTE: This is a stupid idea. Should be in the main program only.
    UNTOUCHED = 0
    BUSY = 1
    DONE = 2
    ERROR = -1

    """
    A dictionary of tuples used for simple command line construction.
    """
    chromosomes = dict()

    """
    A dictionary of status arrays.
    """
    status_arrays = dict()

    """
    Unfortunately, I must maintain a list of directories to really know
    where I will be outputting logging and chromosome order information
    (for combining the variants later) after a tumor:normal sample pair
    is finished getting processed.
    """
    dirlist = list()

    """
    The latest sample number. _Not_ to be used in the indexing functions,
    only for adding chromosome/array pairs.
    """
    sample_num = int()

    """
    Simply sets the sample number to zero.
    """
    def __init__(self):
        self.sample_num = 0

    """
    Adds both a chromosome and its associated status array
    to their corresponding dictionaries. Increments the
    sample number to reflect this.
    bamname: either the tumor or normal bam file name,
    as both should have the same amount of references.
    """
    def add_chrom_and_array(self, bamname):
        with AlignmentFile(bamname, 'rb') as bamfile:
            self.chromosomes[self.sample_num] = bamfile.references
            self.status_arrays[self.sample_num] = Array('i', 
                                                    len(self.chromosomes))
        self.sample_num += 1

    """
    Atomically retrieves a tuple composed of the chromosome string and status.
    This will allow a thread to know whether another thread is working
    on a chromosome or not, as well as format the string for subprocess
    if the chromosome is untouched.
    sample_number: integer denoting the sample_number'th tumor:normal
    pair.
    chr_ndx: Index specifying chromosome string and status element.
    """
    @catch_index_error
    def get_chrostatus(self, sample_number, chr_ndx):
        status_arr = self.status_arrays[sample_number]
        #status_arr.acquire()
        chr_set = self.chromosomes[sample_number][chr_ndx]
        curr_status = status_arr.get_obj()[chr_ndx]
        #status_arr.release()
        return (chr_set, curr_status)

    """
    Atomically sets the desired status array element to the specified
    status.
    sample_number & chr_ndx: see above
    status: One of the values UNTOUCHED, BUSY, DONE or ERROR.
    """
    @catch_index_error
    def set_chrostatus(self, sample_number, chr_ndx, status):
        status_arr = self.status_arrays[sample_number]
        #status_arr.acquire()
        status_arr.get_obj()[chr_ndx] = status
        #status_arr.release()

    """
    Logs all the information from the status array to a file
    called "status.list" in the corresponding directory of 
    vcf files.

    Should only be used within a critical section, when all elements
    of the status array are ERROR or DONE (In other words, when the current
    sample pair has finished getting processed).
    """
    def log_status_and_chromosomes(self, sample_number):
        dirname = self.dirlist[sample_number]
        status_filename = os.path.join(dirname, 'status.list')
        chromolist_filename = os.path.join(dirname, 'chrs.list')
        #Output diagnostic information from MuTect runs
        #If MuTect's stderr output is relatively small, that will be
        #logged instead of "ERROR", with a pretty header as well.
        with open(status_filename, 'w') as stat:
            status_array = self.status_arrays[sample_number]
            chromosome_list = self.chromosomes[sample_number]
            str_status = {self.DONE: 'DONE', self.ERROR: 'ERROR'}
            chr_number = 0
            for status in status_array.get_obj():
                stat.write(('The MuTect process on chromosome {}'
                            ' completed as {}\n'
                           ).format(chromosome_list[chr_number],
                                    str_status[status]))
        with open(chromolist_filename, 'w') as chromo:
            chromo.write('\n'.join([c for c in chromosome_list]))
