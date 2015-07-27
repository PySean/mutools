#!/usr/bin/env python
"""
    "synchrom.py", by Sean Soderman
    Implements a class that does pathname and directory
    bookkeeping, as well as building & storing command line templates.
"""
import os
import re
import sys
#Same as in mapqto0.py: Check for any package problems concerning pysam.
for i in range(0, len(sys.path)):
    traversed = []
    for dirpath, dirnames, filenames in os.walk(sys.path[i]):
        for dir in dirnames:
            if 'pysam' == dir and not dirpath.endswith('site-packages'):
                if dirpath not in traversed:
                    sys.path.append(sys.path.pop(i))
                    traversed.append(dirpath)
try:
    from pysam import AlignmentFile
except ImportError as I:
    sys.stderr.write('Please install the required modules: {}'
                     .format(I))
    sys.exit(1)

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
                                             cmd_args.process_whole_bam,
                                             infile=True)
        else:
            self.commands = self.get_command(cmd_args.pairs, 
                                            cmd_args.process_whole_bam)
        #Create generator for default case if processing files 
        #by chromosome segments at a time.
        if not cmd_args.process_whole_bam:
            cmd_copy = self.commands
            tumorbam = ''
            line = ''
            #Get the path to the tumor BAM file so protogen can create
            #a list of chromosomes.
            if cmd_args.bamlistfile is not None:
                with open(cmd_args.bamlistfile) as listfile:
                    line = listfile.readline().strip()
                    #Skip any header lines
                    while not re.search('.*bam', line):
                        line = listfile.readline().strip()
                    tumorbam = line.split()[0]
            else:
                tumorbam = cmd_args.pairs[0].split(':')[0]

            tumorbam = os.path.join(self.inputdir, tumorbam)
            self.commands = self.protogen(tumorbam, cmd_copy)

    """
    Generator for default case. Creates n commands for each file, where
    n = the number of chromosomes * the number of files.
    Currently, composes the generator from get_command. May not
    need to do this, but right now it feels natural..
    """
    def protogen(self, tumorbam, cmdgen):
        chrlist = None
        with AlignmentFile(tumorbam, 'rb') as bamfile:
            chrlist = bamfile.references
        cont = True
        cmd = ''
        while cont:
            try:
                cmd = cmdgen.next()
            except StopIteration:
                cont = False
                continue
            for i in range(0, len(chrlist)):
                yield cmd % (chrlist[i], chrlist[i] + '.vcf')

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
            except IOError:
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

    Side effects: Creates an output directory for each BAM file pair.
    Also creates a file 'chrs.list' in the output directory.
    """
    def build_command(self, sample_pair):
        tumor, normal = map(os.path.basename, sample_pair)
        filedir = ""
        #The directory of vcf files is <tumorbasename>_<normalbasename>,
        #within the parent output directory
        tumdir, normdir = tumor.split('.bam')[0], normal.split('.bam')[0]
        #Create the output directory name, with a '/' (or '\') at the end.
        if normal != '':
            filedir = os.path.join(self.outputdir, (tumdir + '_' + normdir), '')
            normal = os.path.join(self.inputdir, normal)
            normal = '--input_file:normal ' + normal
        else:
            filedir = os.path.join(self.outputdir, tumdir, '')
        os.makedirs(filedir.strip('/'))
        #No possibility for tumor to equal '' as the calling function
        #handles this case and returns None as a result.
        tumor = os.path.join(self.inputdir, tumor)
        #Write the chromosome list to the output directory.
        with AlignmentFile(tumor, 'rb') as tumbam:
            with open(os.path.join(filedir, 'chrs.list'), 'w') as chrlist:
                map(chrlist.write, os.linesep.join(tumbam.references))
            
        tumor = '--input_file:tumor ' + tumor
        cmd = self.cmd_template.format(normal=normal, tumor=tumor, 
                                       filedir=filedir)
        return cmd

    """
    Builds a command intended for the processing of an entire
    BAM file.

    Side effect: Creates the output directory specified for the vcf
    outputs.
    """
    def build_ntcommand(self, sample_pair):
        tumor, normal = sample_pair
        filename = ''
        #Construct output filename. 
        #Can't depend on there only being one dot in the filename, so now
        #only splits on the .bam extension!
        if normal != '':
            filename = (tumor.split('.bam')[0] + '_' + 
                       normal.split('.bam')[0]  + '.vcf')
        else:
            filename = tumor.split('.bam')[0] + '.vcf'
        if not os.path.exists(self.outputdir):
            os.mkdir(self.outputdir)
        #Create output pathname.
        outfile = os.path.join(self.outputdir, os.path.basename(filename))
        #Now create the paths to the input tumor/normal files.
        tumor = os.path.join(self.inputdir, tumor)
        tumor = '--input_file:tumor ' + tumor
        if normal != '':
            normal = os.path.join(self.inputdir, normal)
            normal = '--input_file:normal ' + normal
        return self.ntcmd_template.format(normal=normal, tumor=tumor,
                                        filepath=outfile)
