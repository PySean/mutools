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
    '''
    cmd_template = 'java -Xmx2g -jar {mupath} --analysis_type MuTect \
--showFullBamList --reference_sequence {fasta} {normal} {tumor} \
--intervals %s -vcf {dirname}%s {mutectopts}'
    '''
    #TODO: Store the path to the directory inside the output directory
    #somewhere (for the logging and stuff).
    cmd_template = 'java -Xmx2g -jar {mupath} --analysis_type MuTect \
--showFullBamList --reference_sequence {fasta} {{normal}} {{tumor}} \
--intervals %s -vcf {{filedir}}%s {mutectopts}'
    """
    List containing a command line template for each thread.#
    The command leaves out the '--intervals <chr>' portion as a %s
    so that the threads can work out which chromosome to process themselves.
    """
    cmd_strings = list()

    """
    List containing the output directories of each sample.
    Required due to logging and chromosome index output
    mechanisms (their output goes in here).
    """
    output_dirs = list()

    """
    The name of the directory that will contain the output directories.
    Also, the name of the directory containing the input files.
    """
    inputdir = str()
    outputdir = str()

    """
    Generator that retrieves commands by reading tumor:normal
    pairs.
    """
    commands = object()

    def __init__(self, cmd_args):
        fasta = cmd_args.fasta
        mutectopts = cmd_args.mutectopts
        mupath = cmd_args.mupath
        #Insert the stuff that won't change into the cmd template.
        self.cmd_template = self.cmd_template.format(fasta=fasta, mutectopts=mutectopts, mupath=mupath)

        self.outputdir = cmd_args.outputdir
        self.inputdir = cmd_args.inputdir

        if cmd_args.bamlistfile is not None:
            self.commands = self.get_command(cmd_args.bamlistfile, infile=True)
        else:
            self.commands = self.get_command(cmd_args.pair)

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
                sys.stderr.write(('get_command' 
                                  ' could not open file {}\n'
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
    them a part of the original command.

    ***REMOVE ALL SIDE EFFECTS*** (implement them elsewhere now)

    Side effects:
    
    Adds the chromosomes, status array, and output directory name
    to their respective structures in the ChromoList object.

    Creates a directory (inside output directory) of the form tumor_normal, 
    with the .bam file extensions truncated.

    Adds the command line to the cmd_strings list.
    """
    def build_command(self, sample_pair):
        tumor, normal = sample_pair
        #The directory of vcf files is <tumorbasename>_<normalbasename>,
        #within the parent output directory
        tumdir, normdir = tumor.split('.')[0] , normal.split('.')[0]
        #Create the output directory name, with a '/' (or '\') at the end.
        #Append it to the list so it can be used later for log/ndx output.
        sample_out = os.path.join(self.outputdir, (tumdir + '_' + normdir), '')
        self.output_dirs.append(sample_out)
        if normal != '':
            normal = os.path.join(self.inputdir, normal)
            normal = '--input_file:normal ' + normal
        #No possibility for tumor to equal '' as the calling function
        #handles this case and returns None as a result.
        tumor = os.path.join(self.inputdir, tumor)
        tumor = '--input_file:tumor ' + tumor
        cmd = self.cmd_template.format(normal=normal, tumor=tumor, 
                                       filedir=sample_out)
        self.cmd_strings.append(cmd)
        return cmd
