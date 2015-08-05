#!/usr/bin/env python
"""
    "catenate.py", by Sean Soderman
    Concatenates all vcfs under a directory in an order 
    according to the chrs.list file within the directory.
"""
import argparse
import os
import re
import sys
from subprocess import check_output, CalledProcessError
try:
    from arguer import makeparser
except ImportError as I:
    sys.stderr.write('Make sure arguer.py is in my working directory: {}'
                         .format(I))

parser = makeparser(('Concatenates the vcf files within the subdirectories'
                     ' created by multimutect.'))
parser.add_argument('--delete_fragments', action='store_true',
                    help='Delete files utilized for catenation')
parser.add_argument('-l', '--listfile', type=str,
                    help=('A list specifying the samples to'
                          'concatenate. Useful if the BAM files used'
                          'for the analysis were created on a per-'
                          'chromosome basis. Uses chrs.list by default.'),
                    default='chrs.list')
args = parser.parse_args()
"""
Validates the list of chromosome vcf file fragments. Removes paths that
lead to files of size zero, or that do not exist.
"""
def chr_validate(chrlist):
    fil_func = lambda x: os.stat(x).st_size != 0 and os.path.exists(x)
    return filter(fil_func, chrlist)

"""
Concatenates all vcfs under a directory or on a command line.
"""
def vcf_catenate(directory, reference, gatkpath, listfile, delete_fragments):
    cmd = ('java -cp {gatk} org.broadinstitute.gatk.tools.CatVariants'
           ' -assumeSorted -R' 
           ' {ref} -out {{out}} {{vcfs}}').format(gatk=gatkpath, 
                                                ref=reference)
    for dirpath, dirnames, filenames in os.walk(directory):
        for dir in dirnames:
            d_path = os.path.join(dirpath, dir)
            listing = os.path.join(d_path, listfile)
            result_name = dir + '.vcf'
            outpath = os.path.join(dirpath, result_name)
            #Make sure the listing file exists and that this isn't the
            #status directory.
            if os.path.exists(listing):
                with open(listing, 'r') as chrfile:
                    chrlist = [os.path.join(d_path, c.strip()) + '.vcf' 
                               for c in chrfile]
                    #Keep original chrlist for deletions of possibly empty
                    #VCF files.
                    realchrs = chr_validate(chrlist)
                    #The list of vcfs to concatenate, each prepended by
                    #-V.
                    vseries = "".join(['-V ' + c + ' ' for c in realchrs])
                    final_cmd = cmd.format(out=outpath, vcfs=vseries)
                    try:
                        check_output(final_cmd.split())
                    except CalledProcessError as cpe:
                        sys.stderr.write('Problem: {}\n'.format(cpe))
                    #Clean up leftover files after combining them.
                    #Also cleans up .idx files and directories.
                    if delete_fragments == True:
                        map(os.unlink, chrlist)
                        #No index file is created for empty vcfs.
                        map(lambda x: os.unlink(x + '.idx'), realchrs)
                        os.unlink(listing)
                        os.rmdir(d_path)
            else:
                print("I don't exist: {}".format(listing))
"""
Smaller concatenation function for the case of a single directory with
VCFs generated from BAM files each representing single chromosomes.
"""
def minicat(directory, reference, gatkpath, listfile, delete_fragments):
    cmd = ('java -cp {gatk} org.broadinstitute.gatk.tools.CatVariants'
           ' -assumeSorted -R' 
           ' {ref} -out {{out}} {{vcfs}}').format(gatk=gatkpath, 
                                                ref=reference)
    file_list = []
    if os.path.exists(directory) and os.path.exists(listfile):
        with open(listfile, 'r') as lfile:
            #Convert .bam extension to .vcf.
            file_vcfs = [re.sub('\.bam', '.vcf', bam) for bam in lfile]
            file_list = [os.path.join(directory, 
                                      os.path.basename(vcf.strip()))
                        for vcf in file_vcfs]
    else:
        sys.stderr.write(('Error: {} or {} do not exist.'
                         ' Please check your command'
                         ' line').format(directory, listfile))
        sys.exit(1)
    #Create series of VCFs to concatenate
    vseries = "".join(['-V ' + vcf + ' ' for vcf in file_list])
    outfile = directory + '.vcf'
    final_cmd = cmd.format(out=outfile, vcfs=vseries)
    check_output(final_cmd.split())
    
    if delete_fragments == True:
        map(os.unlink, file_list)
        map(lambda x: os.unlink(x + '.idx'), file_list)

cat_func = vcf_catenate
if args.listfile != 'chrs.list':
    cat_func = minicat

#Filter out unused "vcf_files" option (used in combine.py.).
arg_dict = {key: value for (key, value) in vars(args).iteritems()
            if key != 'vcf_files'}
if os.path.exists(args.gatkpath):
    cat_func(**vars(args))
#Try again with the default.
elif os.path.exists('gatk.jar'):
    cat_func(**vars(args))
else:
    sys.stderr.write('Please provide an existent gatk jar filepath.')
    sys.exit(1)
