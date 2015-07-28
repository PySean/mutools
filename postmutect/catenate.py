#!/usr/bin/env python
"""
    "catenate.py", by Sean Soderman
    Concatenates all vcfs under a directory in an order 
    according to the chrs.list file within the directory.
"""
import argparse
import os
import sys
from subprocess import check_output
try:
    from arguer import makeparser
except ImportError as I:
    sys.stderr.write('Make sure arguer.py is in my working directory: {}'
                         .format(I))

parser = makeparser(('Concatenates the vcf files within the subdirectories'
                     ' created by multimutect.'))
parser.add_argument('--delete_fragments', action='store_true',
                    help='Delete files utilized for catenation')
parser.add_argument('-l', '--listing', type=str,
                    help=('A path specifying the list of samples to'
                          'concatenate. Useful if the BAM files used'
                          'for the analysis were created on a per-'
                          'chromosome basis. Uses chrs.list by default.')
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
Concatenates all vcfs under a directory. 
By default, does not remove the vcf files.
At the moment, cleans up after multimutect.py with regards to the 
creation of a status directory as well as moving status files up into
said directory.
"""
def vcf_catenate(parent, reference, gatkpath, listfile, delete):
    cmd = ('java -cp {gatk} org.broadinstitute.gatk.tools.CatVariants'
           ' -assumeSorted -R' 
           ' {ref} -out {{out}} {{vcfs}}').format(gatk=gatkpath, 
                                                ref=reference)
    for dirpath, dirnames, filenames in os.walk(parent):
        for directory in dirnames:
            d_path = os.path.join(dirpath, directory)
            listing = os.path.join(d_path, listfile)
            result_name = directory + '.vcf'
            status_name = directory + '.list'
            outpath = os.path.join(dirpath, result_name)
            #Make sure the listing file exists and that this isn't the
            #status directory.
            if os.path.exists(listing):
                with open(listing, 'r') as chrfile:
                    chrlist = [os.path.join(d_path, c.strip()) + '.vcf' 
                               for c in chrfile]
                    #Keep original chrlist for deletions.
                    realchrs = chr_validate(chrlist)
                    #The list of vcfs to concatenate, each prepended by
                    #-V.
                    vseries = "".join(['-V ' + c + ' ' for c in realchrs])
                    final_cmd = cmd.format(out=outpath, vcfs=vseries)
                    check_output(final_cmd.split())
                    #Clean up leftover files after combining them if the option
                    #is specified.
                    #Also cleans up .idx files and directories.
                    if delete == True:
                        map(os.unlink, chrlist)
                        #No index file is created for empty vcfs.
                        map(lambda x: os.unlink(x + '.idx'), realchrs)
                        os.unlink(listing)
                        os.rmdir(d_path)

if os.path.exists(args.gatkpath):
    vcf_catenate(args.directory, args.reference, args.gatkpath)
elif os.path.exists('gatk.jar'):
    vcf_catenate(args.directory, args.reference, args.gatkpath)
else:
    sys.stderr.write('Please provide an existent gatk jar filepath.')
    sys.exit(1)
