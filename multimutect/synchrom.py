#!/usr/bin/env python
"""
    "synchrom.py", by Sean Soderman
    Implements a class that does pathname and directory
    bookkeeping, as well as building & storing command line templates.
"""
import os
import re
import sys


class Synchrom():
    """
    Standard command template utilized when processing a chromosome at a time
    (using threads).
    """
    cmd_template = ('java -Xmx{mem}g -jar {mupath} --analysis_type MuTect'
    ' --showFullBamList --reference_sequence {fasta} {{normal}} {{tumor}}'
    ' --intervals %s -vcf {{filedir}}%s {mutectopts}')

    """
    Nonstandard command template used when processing entire BAM files at 
    a time. As such, omits the --intervals option.
    """
    ntcmd_template = ('java -Xmx{mem}g -jar {mupath} --analysis_type MuTect'
    ' --showFullBamList --reference_sequence {fasta} {{normal}} {{tumor}}'
    ' -vcf {{filepath}} {mutectopts}')

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

    The input list is required by add_chrom_and_array, as each pair
    will have a (possibly) different chromosome listing.
    """
    bam_inputs = list()
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
        mu_opts = ''
        if cmd_args.conf != '':
            if os.path.exists(cmd_args.conf):
                with open(cmd_args.conf, 'r') as cmd_conf:
                    mu_opts = ''.join([re.sub(os.linesep, ' ', cmd)
                                      for cmd in cmd_conf])
            else:
                sys.stderr.write('Problem: {} does not exist.{}'
                                 .format(cmd_args.conf, os.linesep))
                sys.exit(1)
        else:
            mu_opts = cmd_args.mutectopts
            
        mupath = cmd_args.mupath
        mem = cmd_args.mem
        #Insert the stuff that won't change into the cmd template.
        self.cmd_template = self.cmd_template.format(mem=mem,
                                                     fasta=fasta, 
                                                     mutectopts=mu_opts,
                                                     mupath=mupath)
        self.ntcmd_template = self.ntcmd_template.format(mem=mem,
                                                         fasta=fasta, 
                                                         mutectopts=mu_opts,
                                                         mupath=mupath)
        self.outputdir = cmd_args.outputdir
        self.inputdir = cmd_args.inputdir

        if cmd_args.bamlistfile is not None:
            self.commands = self.get_command(cmd_args.bamlistfile, 
                                             cmd_args.whole,
                                             infile=True)
        else:
            self.commands = self.get_command(cmd_args.pair, cmd_args.whole)

    """
    Parses a file or cmd line list into tumor:normal pairs. 
    Returns a single command. Must be surrounded in a StopIteration
    try/except. Returns None when a file isn't found or there is no
    tumor sample in the pair.
    """
    def get_command(self, sample_pairs, whole, infile=False):
        err_str = 'Error: argument|line {} has no tumor filename.\n'
        line_number = 0
        pairsep = ':'
        if infile == True:
            try:
                sample_pairs = open(sample_pairs, 'r')
                pairsep = '\s+'
            except FileNotFoundError:
                sys.stderr.write(('get_command' 
                                  ' could not open file {}\n'
                                 ).format(sample_pairs))
                yield None
        for pair in sample_pairs:
            #Skip header line
            if not re.search('.*bam', pair):
                continue
            tumor, normal = re.split(pairsep, pair)
            normal = normal.strip()
            if tumor == '':
                sys.stderr.write(err_str.format(line_number))
                yield None
            line_number += 1
            if not whole:
                yield self.build_command( (tumor, normal) )
            else:
                yield self.build_ntcommand( (tumor, normal))
        #Once all samples have been processed, close the listing file.
        if infile == True:
            sample_pairs.close()
    """
    Retrieves a single command corresponding to the tumor:normal
    pair fetched from the file or command line with get_pair.

    fasta, mutectopts, and the path to MuTect are never changed
    and needed only in the construction of the command, so I have made 
    them a part of the original command.

    Adds the output directory to the output_dirs list, input directory to
    the bam_inputs list, and the command line to the cmd_strings list.
    """
    def build_command(self, sample_pair):
        tumor, normal = sample_pair
        filedir = ""
        #The directory of vcf files is <tumorbasename>_<normalbasename>,
        #within the parent output directory
        tumdir, normdir = tumor.split('.')[0] , normal.split('.')[0]
        #Create the output directory name, with a '/' (or '\') at the end.
        #Append it to the list so it can be used later for log/ndx output.
        if normal != '':
            filedir = os.path.join(self.outputdir, (tumdir + '_' + normdir), '')
            normal = os.path.join(self.inputdir, normal)
            normal = '--input_file:normal ' + normal
        else:
            filedir = os.path.join(self.outputdir, tumdir, '')
        self.output_dirs.append(filedir)
        #No possibility for tumor to equal '' as the calling function
        #handles this case and returns None as a result.
        tumor = os.path.join(self.inputdir, tumor)
        self.bam_inputs.append(tumor)
        tumor = '--input_file:tumor ' + tumor
        cmd = self.cmd_template.format(normal=normal, tumor=tumor, 
                                       filedir=filedir)
        self.cmd_strings.append(cmd)
        return cmd

    """
    Builds a command not intended for multithreaded use.
    Omits the allocation steps build_command utilizes.
    """
    def build_ntcommand(self, sample_pair):
        tumor, normal = sample_pair
        filename = ''
        #Construct output filename.
        if normal != '':
            filename = tumor.split('.')[0] + '_' + normal.split('.')[0] \
            + '.vcf'
        else:
            filename = tumor.split('.')[0] + '.vcf'
        #Create output pathname.
        outfile = os.path.join(self.outputdir, filename)
        #Now create the paths to the input tumor/normal files.
        tumor = os.path.join(self.inputdir, tumor)
        tumor = '--input_file:tumor ' + tumor
        if normal != '':
            normal = os.path.join(self.inputdir, normal)
            normal = '--input_file:normal ' + normal
        return self.ntcmd_template.format(normal=normal, tumor=tumor,
                                        filepath=outfile)
