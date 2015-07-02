#!/usr/bin/python
"""
    "chromolist.py", by Sean Soderman
    An object that tightly couples the tuple of chromosomes with the
    list of status arrays.
"""
import os
import sys
from pysam import AlignmentFile
from multiprocessing import Array


"""
Status constants to be used when having a thread assign a status
to an element within a status array.
"""
UNTOUCHED = 0
BUSY = 1
DONE = 2
ERROR = -1

"""
Function decorator for the methods get_chrostatus and set_chrostatus.
Catches an IndexError and prints out relevant debug information.
TODO: Change to KeyError since I switched to using dictionaries.
"""
def catch_key_error(func):
    def wrap(self, *args, **kwargs):
        try:
            return func(self, *args, **kwargs)
        except KeyError as K:
            #The first argument is the sample number, second is the
            #chromosome index.
            array_number = None
            cell_and_chromo = None
            if len(args) >= 2:
                array_number, cell_and_chromo = args[0:2]
            else:
                array_number = args[0]
            chromo_str = self.chromosomes[cell_and_chromo]
            sys.stderr.write(('An attempt to access Status array no. {array}'
                              ' status cell/chromosome no. {cell} and' 
                              ' chromosome {chr} was made in function {func}:'
                              '{K}\n'
                             ).format(array = array_number, 
                                      cell = cell_and_chromo,
                                      chr = chromo_str, 
                                      func=func.func_name, K=K))
    return wrap

class ChromoList():
    """
    A dictionary of tuples used for simple command line construction.
    """
    chromosomes = dict()

    """
    A dictionary of status arrays.
    """
    status_arrays = dict()
    
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
                                                len(bamfile.references))
        self.sample_num += 1

    """
    Retrieves a tuple composed of the chromosome string and status.
    This will allow a thread to know whether another thread is working
    on a chromosome or not, as well as format the string for subprocess
    if the chromosome is untouched.
    sample_number: integer denoting the sample_number'th tumor:normal
    pair.
    chr_ndx: Index specifying chromosome string and status element.
    """
    @catch_key_error
    def get_chrostatus(self, sample_number, chr_ndx):
        status_arr = self.status_arrays[sample_number]
        chr_set = self.chromosomes[sample_number][chr_ndx]
        curr_status = status_arr.get_obj()[chr_ndx]
        return (chr_set, curr_status)

    """
    Sets the desired status array element to the specified
    status.
    sample_number & chr_ndx: see above
    status: One of the values UNTOUCHED, BUSY, DONE or ERROR.
    """
    @catch_key_error
    def set_chrostatus(self, sample_number, chr_ndx, status):
        status_arr = self.status_arrays[sample_number]
        status_arr.get_obj()[chr_ndx] = status

    """
    Checks if all statuses in the specified status array are DONE
    or ERROR.
    """
    #TODO: Check for bad interactions between this function and its
    #decorator. Make sure the log function works as well.
    @catch_key_error
    def will_log(self, sample_number):
        statuses = self.status_arrays[sample_number].get_obj()
        return all([(i == DONE) or (i == ERROR) for i in statuses])

    """
    Locks or unlocks the status array specified by the sample number.
    Returns: None if the array does not exist.
    """
    def key_array(self, sample_number, action='lock'):
        try:
            if action == 'lock':
                self.status_arrays[sample_number].acquire()
            else:
                self.status_arrays[sample_number].release()
        except KeyError:
            sys.stderr.write('Error, status array no. {} does not exist\n'
                             .format(sample_number))

    """
    Checks all elements in the nth (sample_number'th) array
    beginning from chr_ndx + 1 to the end of the array for the
    UNTOUCHED status. 
    Returns the index of this status. Returns None if nothing
    is untouched.
    """
    def check_ahead(self, sample_number, chr_ndx):
        statuses = self.status_arrays[sample_number].get_obj()
        ahead_items = statuses[chr_ndx + 1:]
        for i in range(chr_ndx + 1, len(statuses)):
            if statuses[i] == UNTOUCHED:
                return i
        else:
            return None

    """
    Logs all the information from the status array to a file
    called "status.list" in the corresponding directory of 
    vcf files.

    Should only be used within a critical section, when all elements
    of the status array are ERROR or DONE (In other words, when the current
    sample pair has finished getting processed).
    output_dir: the directory the log and chromosome index will go.
    """
    def log_status_and_chromosomes(self, output_dir, sample_number):
        dirname = output_dir
        status_filename = os.path.join(dirname, 'status.list')
        chromolist_filename = os.path.join(dirname, 'chrs.list')
        #Output diagnostic information from MuTect runs
        #If MuTect's stderr output is relatively small, that will be
        #logged instead of "ERROR", with a pretty header as well.
        with open(status_filename, 'w') as stat:
            status_array = self.status_arrays[sample_number]
            chromosome_list = self.chromosomes[sample_number]
            str_status = {DONE: 'DONE', ERROR: 'ERROR'}
            chr_number = 0
            for status in status_array.get_obj():
                stat.write(('The MuTect process on chromosome {}'
                            ' completed as {}\n'
                           ).format(chromosome_list[chr_number],
                                    str_status[status]))
                chr_number += 1
        with open(chromolist_filename, 'w') as chromo:
            chromo.write('\n'.join([c for c in chromosome_list]))
