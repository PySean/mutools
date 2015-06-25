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
Catches an IndexError with the 
"""
def catch_index_error(func):
    def wrap(self, *args, **kwargs):
        try:
            return func(self, *args, **kwargs)
        except IndexError as I:
            #The first argument is the index.
            index = args[0]
            array_number = index / len(self.chromosomes)
            cell_and_chromo = index % len(self.chromosomes)
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

    #Tuple utilized for easy command line string construction.
    chromosomes = tuple()

    """
    A list of status arrays.
    """
    status_arrays = list()

    def __init__(self, pair_len):
        self.chromosomes = self._make_chromosomes_()
        self.status_arrays = self._make_arrays_(len(self.chromosomes), 
                                                  pair_len)

    #FIXME: Now is generalized for any kind of genomic data. 
    #filename: The name of a BAM file.
    def _make_chromosomes_(self, filename):
        with AlignmentFile(filename, 'rb') as bamfile:
            return bamfile.references
        """
        biff = ['chrM']
        biff.extend(['chr{}'.format(i) for i in range(1, 23)])
        biff.extend(['chrX', 'chrY'])
        return tuple(biff)
        """

    #Stores several thread-safe Array objects in the list.
    def _make_arrays_(self, chr_len, pair_len):
        return [Array('i', chr_len) for i in range(0, pair_len)]

    """
    Atomically retrieves a tuple composed of the chromosome string and status.
    This will allow a thread to know whether another thread is working
    on a chromosome or not, as well as format the string for subprocess
    if the chromosome is untouched.
    """
    @catch_index_error
    def get_chrostatus(self, ndx):
        status_arr = self.status_arrays[int(ndx / 25)]
        status_arr.acquire()
        chr_set = self.chromosomes[ndx % 25]
        curr_status = status_arr.get_obj()[ndx % 25]
        status_arr.release()
        return (chr_set, curr_status)

    """
    Atomically sets the desired status array element to the specified
    status.
    """
    @catch_index_error
    def set_chrostatus(self, ndx, status):
        status_arr = self.status_arrays[int(ndx / 25)]
        status_arr.acquire()
        status_arr.get_obj()[ndx % 25] = status
        status_arr.release()
