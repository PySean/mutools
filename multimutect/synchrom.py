#!/usr/bin/env python
"""
    "synchrom.py", by Sean Soderman
    Implements a class wrapping a tuple of chromosomes
    from chrM, chr0-22, then chrX-chrY. The class also contains
    a "ChromoList" object that allows a thread to atomically
    get and set job statuses.
"""
from chromolist import ChromoList
import sys


class Synchrom():
    """
    The bulk of the command line.
    Only the chromosome bit is left to the threads to figure out,
    and is denoted by %s. For the -vcf argument, the %s _must_ be 
    constructed by os.path.join for maximum portability.
    """
    cmd_template = 'java -Xmx2g -jar {mupath} --analysis_type MuTect \
--showFullBamList --reference_sequence {fasta} {normal} {tumor} \
--intervals %s -vcf {dirname}%s {mutectopts}'
    """
    List containing a command line template for each thread.#
    The command leaves out the '--intervals <chr>' portion as a {}
    so that the threads can work out which chromosome to process themselves.
    """
    cmd_strings = list()
    """
    A list of thread-safe status arrays. 
    There's one for each tumor/{normal} (or just tumor) pair.
    """
    status_arrays = list()
    """
    The chromolist object for this class.
    Acts as a sort of "state machine" for getting to the next sample.
    """
    chromolist = object()
    """
    Information pertaining to core parameters that musn't be omitted.
    These are:
    fasta: fasta file name
    mutectopts: optional arguments to be passed to MuTect
    mupath: path to MuTect executable.
    """
    fasta, mutectopts, mupath = (str(), str(), str())
    """
    Generator that retrieves commands by reading tumor:normal
    pairs.
    """
    commands = iter()
    def __init__(self, cmd_args):
        self.chromolist = ChromoList()
        self.fasta = cmd_args.fasta
        self.mutectopts = cmd_args.mutectopts
        self.mupath = cmd_args.mupath
        if cmd_args.bamlistfile is not None:
            self.commands = self.get_command(cmd_args.bamlistfile, infile=True)
        else:
            self.commands = self.get_command(cmd_args.pair)

    """
    Retrieves a single command corresponding to the tumor:normal
    pair fetched from the file or command line with _get_pair_.
    
    fasta, mutectopts, and the path to MuTect are never changed
    and needed only in the construction of the command, so I have made 
    them instance variables instead.
    """
    def build_command(self, sample_pair):
        tumor, normal = sample_pair
        normal = normal.strip()
        tumdir, normdir = tumor.split('.')[0] , normal.split('.')[0]
        dirname = tumdir + normdir
        if normal != '':
            normal = '--input_file:normal ' + normal
        #No possibility for tumor to equal '' as the calling function
        #handles this case and returns None as a result.
        tumor = '--input_file:tumor ' + tumor
        #NOTE: It isn't necessary for me to insert the fasta,
        #path, and option data every time as this is currently static.
        #However, if I allow for a mix of species tumor:normal sample
        #pairs, I will need to insert fasta as well.
        cmd = self.cmd_template.format(fasta=self.fasta, normal=normal,
                                           tumor=tumor, dirname=dirname,
                                           mupath=self.mupath,
                                           mutectopts=self.mutectopts)
        return cmd

    """
    Parses a file or cmd line list into tumor:normal pairs. 
    Returns a single command. Must be surrounded in a StopIteration
    try/except. Returns None when a file isn't found or there is no
    tumor sample in the pair.

    This function also adds a chromosome set and status array
    to the chromolist in this object, mapping the sample set to these
    structures.
    """
    def get_command(self, sample_pairs, infile=False):
        err_str = 'Error: argument|line {} has no tumor filename.\n'
        line_number = 0
        if infile == True:
            try:
                sample_pairs = open(sample_pairs, 'r')
            except FileNotFoundError:
                sys.stderr.write(('_parse_pairs_' 
                                  ' Could not open file {}\n'
                                 ).format(sample_pairs))
                return None
        for pair in sample_pairs:
            tumor, normal = pair.split(':')
            if tumor == '':
                sys.stderr.write(err_str.format(line_number))
                return None
            line_number += 1
            self.chromolist.add_chrom_and_array(tumor)
            yield self.build_command(pair)
        if infile == True:
            sample_pairs.close()
