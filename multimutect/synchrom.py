#!/usr/bin/env python
"""
    "synchrom.py", by Sean Soderman
    Implements a class wrapping a tuple of chromosomes
    from chrM, chr0-22, then chrX-chrY. The class also contains
    a "ChromoList" object that allows a thread to atomically
    get and set job statuses.
"""
from chromolist import ChromoList
import os
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
    The chromolist object for this class.
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
    The name of the directory that will contain the output directories.
    Also, the name of the directory containing the input files.
    """
    inputdir = str()
    outputdir = str()

    """
    Generator that retrieves commands by reading tumor:normal
    pairs. Basically "drives" the chromolist object.
    """
    commands = object

    def __init__(self, cmd_args):
        self.chromolist = ChromoList()
        self.fasta = cmd_args.fasta
        self.mutectopts = cmd_args.mutectopts
        self.mupath = cmd_args.mupath
        self.inputdir = cmd_args.inputdir
        self.outputdir = cmd_args.outputdir

        if cmd_args.bamlistfile is not None:
            self.commands = self.get_command(cmd_args.bamlistfile, infile=True)
        else:
            self.commands = self.get_command(cmd_args.pair)
        """
        If an output directory name was specified, create a directory
        with this name. Otherwise, create the default directory called
        "output".
        """
        if self.outputdir is not None:
            try:
                os.mkdir(self.outputdir)
            except OSError as O:
                sys.stderr.write('Error: could not create directory: %s' % O)
                return None
    """
    Parses a file or cmd line list into tumor:normal pairs. 
    Returns a single command. Must be surrounded in a StopIteration
    try/except. Returns None when a file isn't found or there is no
    tumor sample in the pair.
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
                yield None
        for pair in sample_pairs:
            tumor, normal = pair.split(':')
            normal = normal.strip()
            if tumor == '':
                sys.stderr.write(err_str.format(line_number))
                yield None
            line_number += 1
            yield self.build_command( (tumor, normal) )
        #Once all samples have been processed, close the listing file.
        if infile == True:
            sample_pairs.close()
    """
    Retrieves a single command corresponding to the tumor:normal
    pair fetched from the file or command line with get_pair.
    
    fasta, mutectopts, and the path to MuTect are never changed
    and needed only in the construction of the command, so I have made 
    them instance variables instead.

    Side effects:
    
    Adds the chromosomes, status array, and output directory name
    to their respective structures in the ChromoList object.

    Creates a directory (inside output directory) of the form tumor_normal, 
    with the .bam file extensions truncated.
    """
    def build_command(self, sample_pair):
        tumor, normal = sample_pair
        #The directory of vcf files is <tumorbasename>_<normalbasename>,
        #within the parent output directory
        tumdir, normdir = tumor.split('.')[0] , normal.split('.')[0]
        dirname = os.path.join(self.outputdir, (tumdir + '_' + normdir))
        os.mkdir(dirname)
        self.chromolist.add_chrom_and_array(os.path.join(self.inputdir, tumor))
        self.chromolist.dirlist.append(dirname)
        if normal != '':
            normal = '--input_file:normal ' + normal
        #No possibility for tumor to equal '' as the calling function
        #handles this case and returns None as a result.
        tumor = '--input_file:tumor ' + tumor
        #NOTE: It isn't necessary for me to insert the fasta,
        #path, and option data every time as this is currently static.
        #However, if I allow for a mix of species tumor:normal sample
        #pairs, I will need to insert a fasta filename as well.
        cmd = self.cmd_template.format(fasta=self.fasta, normal=normal,
                                           tumor=tumor, dirname=dirname,
                                           mupath=self.mupath,
                                           mutectopts=self.mutectopts)
        return cmd
