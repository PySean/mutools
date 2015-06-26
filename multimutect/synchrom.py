#!/usr/bin/env python
"""
    "synchrom.py", by Sean Soderman
    Implements a class wrapping a tuple of chromosomes
    from chrM, chr0-22, then chrX-chrY. The class also contains
    a "ChromoList" object that allows a thread to atomically
    get and set job statuses.

    TODO: The ChromoList object may also be modified to also store 
    a growing list of command line templates. This is in anticipation
    of the redesign that doesn't use as much memory at once.
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
    Requires a file handler pointing to a file with a list of tumor:{normal}
    pathnames (fully qualified) or a list of files specified on the command
    line as tumor.bam:{normal.bam} (the braces mean optional).
    kwargs contains:
    TODO: Put the following info in a separate markdown file.
    -m
    --mupath: The path to the mutect executable. This is optional, and
    by default the file 'mutect.jar' from the current directory is used.

    -b
    --bamlistfile: If a file was specified (with the -f option), this is an open
    file handler pointing to the file.

    -p
    --pairs: If the tumor:{normal} pairs were specified on the command line
    directly with the -d option, this is a tuple of tumor:{normal} pairs.
    The normal portion is not strictly necessary. This option and "fd"
    are mutually exclusive, and at least one of them is required.

    -r
    --fasta: The name of the reference sequence file.

    -M
    --mutectopts: Extra parameters to MuTect to be concatenated
    to the end of the string.
    """
    def __init__(self, cmd_args):
        self.chromolist = ChromoList()
        tumor_normals = {}
        if cmd_args.bamlistfile is not None:
            tumor_normals = self._parse_pairs_(cmd_args.bamlistfile, 
                                               infile=True)
        elif cmd_args.pairs is not None:
            tumor_normals = self._parse_pairs_(cmd_args.pairs)
        else:
            return None #At least one must be specified.
        #Create the list of command strings from tumor/normal pairs,
        #the specified reference sequence, and unrequired options for MuTect.
        self.cmd_strings = self._build_commands_(tumor_normals,
                                                 cmd_args.fasta,
                                                 cmd_args.mutectopts,
                                                 cmd_args.mupath)

    #TODO: Make this return one tumor:normal pair at a time.
    #Only one line from the file should be read at a time.
    #When a line is read, the underlying ChromoList object
    #should add the chromosome tuple along with its status array.
    """
    Parses a file or cmd line list into tumor:normal pairs. 
    Returns a dictionary of tumor:normal values.
    """
    def _parse_pairs_(self, sample_pairs, infile=False):
        err_str = 'Error: argument {} has no tumor filename.\n'
        samples = {}
        line_number = 0
        if infile == True:
            try:
                sample_pairs = open(sample_pairs, 'r')
            except FileNotFoundError:
                sys.stderr.write(('_parse_pairs_' 
                                  ' Could not open file {}\n'
                                 ).format(sample_pairs))
                return None
        for arg in sample_pairs:
            tumor, normal = arg.split(':')
            if tumor == '':
                sys.stderr.write(err_str.format(line_number))
                return None
            samples[tumor] = normal
            line_number += 1
        if infile == True:
            sample_pairs.close()
        return samples

    #TODO: fix to only return one command...
    """
    Builds up the list of command strings. There is one for each
    tumor, normal pair. The auxiliary strings are above.
    """
    def _build_commands_(self, tumor_normals, fasta, mutectopts, 
                         mupath='mutect.jar'):
        cmds = list(range(0, len(tumor_normals)))
        i = 0
        print(str(tumor_normals)) #DEBUG
        for tumor, normal in tumor_normals.iteritems():
            normal = normal.strip()
            tumdir, normdir = tumor.split('.')[0] , normal.split('.')[0]
            dirname = tumdir + normdir
            normal = '--input_file:normal ' + normal
            tumor = '--input_file:tumor ' + tumor
            #dirname = "".join(tumdir)
            cmds[i] = self.cmd_template.format(fasta=fasta, normal=normal,
                                               tumor=tumor, dirname=dirname,
                                               mupath=mupath,
                                               mutectopts=mutectopts)
            i += 1
        return cmds
