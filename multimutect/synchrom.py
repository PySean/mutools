#!/usr/bin/env python
"""
    "synchrom.py", by Sean Soderman
    Implements a class wrapping a tuple of chromosomes
    from chrM, chr0-22, then chrX-chrY. The class also contains
    a "ChromoList" object that allows a thread to atomically
    get and set job statuses.

    TODO: The ChromoList object may also be modified to also store 
    a growing list of command line templates. This is in anticipation
    of the 
"""
from inheritlist import ChromoList
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
--intervals %s -vcf {dirname}%s {options}'
    """
    The tuple of chromosome segment strings.
    """
    chrm_segments = tuple()
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
    Requires a file handler pointing to a file with a list of tumor:{normal}
    pathnames (fully qualified) or a list of files specified on the command
    line as tumor.bam:{normal.bam} (the braces mean optional).
    kwargs contains:
    TODO: Put the following info in a separate markdown file.
    -m
    --mupath: The path to the mutect executable. This is optional, and
    by default the file 'mutect.jar' from the current directory is used.

    -b
    --bamlist: If a file was specified (with the -f option), this is an open
    file handler pointing to the file.

    -p
    --pairs: If the tumor:{normal} pairs were specified on the command line
    directly with the -d option, this is a tuple of tumor:{normal} pairs.
    The normal portion is not strictly necessary. This option and "fd"
    are mutually exclusive, and at least one of them is required.

    -r
    --fasta: The name of the reference sequence file.

    -o
    --optional: Extra parameters to MuTect to be concatenated
    to the end of the string.
    """
    def __init__(self, **kwargs):
        tumor_normals = {}
        if kwargs.get('fd') is not None:
            tumor_normals = self._parse_pairs_(kwargs['fd'], infile=True)
        elif kwargs.get('pairs') is not None:
            tumor_normals = self._parse_pairs_(kwargs['pairs'])
        else:
            return None #At least one of the two options must be used.
        #Remove pairs for 'fd' and 'pairs' so I can expand kwargs.
        kwargs.pop('fd') if kwargs.get('fd') is not None else None
        kwargs.pop('pairs') if kwargs.get('pairs') is not None else None
        #Create the list of command strings from tumor/normal pairs,
        #the specified reference sequence, and unrequired options for MuTect.
        self.cmd_strings = self._build_commands_(tumor_normals, **kwargs)

    """
    Parses a file or cmd line list into tumor:normal pairs. 
    Returns a dictionary of tumor:normal values.
    """
    def _parse_pairs_(self, sample_pairs, infile=False):
        err_str = 'Error: argument {} has no tumor filename.\n'
        arg_number = 0
        samples = {}
        if infile == True:
            try:
                sample_pairs = open(sample_pairs, 'r')
            except FileNotFoundError:
                sys.stderr.write(('_parse_pairs' 
                                  ' Could not open file {}\n'
                                 ).format(sample_pairs))
        for arg in sample_pairs:
            tumor, normal = arg.split(':')
            if tumor == '':
                sys.stderr.write(err_str.format(line_number))
                return None
            samples[tumor] = normal
        return samples

    """
    Builds up the list of command strings. There is one for each
    tumor, normal pair. The auxiliary strings are above.
    """
    def _build_commands_(self, tumor_normals, fasta, options, 
                            mupath='mutect.jar'):
        cmds = list(range(0, len(tumor_normals)))
        i = 0
        print(str(tumor_normals)) #DEBUG
        for tumor, normal in zip(tumor_normals.keys(), tumor_normals.values()):
            tumdir, normdir = tumor.split('.')[0] , normal.split('.')[0]
            dirname = tumdir + normdir
            normal = '--input_file:normal ' + normal
            tumor = '--input_file:tumor ' + tumor
            #dirname = "".join(tumdir)
            cmds[i] = self.cmd_template.format(fasta=fasta, normal=normal,
                                               tumor=tumor, dirname=dirname,
                                               mupath=mupath, options=options)
            i += 1
        return cmds
