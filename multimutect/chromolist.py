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
    UNTOUCHED = 0
    BUSY = 1
    DONE = 2
    ERROR = -1

    #TODO: Since there is no benefit to maintaining two shifting
    #dictionaries other than logging purposes, I am going to
    #add a logging mechanism for my program, basically take the
    #status array and output "There was an error in <sample pair>, 
    #chromosome n" or "Successful run for chromosome n" to a file
    #in the output directory filled with the vcf files.
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
        status_arr.acquire()
        chr_set = self.chromosomes[sample_number][chr_ndx]
        curr_status = status_arr.get_obj()[chr_ndx]
        status_arr.release()
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
        status_arr.acquire()
        status_arr.get_obj()[chr_ndx] = status
        status_arr.release()
