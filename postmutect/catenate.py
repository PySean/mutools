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
def vcf_catenate(parent, reference, gatkpath='gatk.jar'):
    cmd = ('java -cp {gatk} org.broadinstitute.gatk.tools.CatVariants'
           ' -assumeSorted -R' 
           ' {ref} -out {{out}} {{vcfs}}').format(gatk=gatkpath, 
                                                ref=reference)
    for dirpath, dirnames, filenames in os.walk(parent):
        for directory in dirnames:
            d_path = os.path.join(dirpath, directory)
            listing = os.path.join(d_path, 'chrs.list')
            #NOTE: An underscore is no longer at the end of the 
            #directories. However I am going to leave this here just in
            #case.
            result_name = directory.strip('_') + '.vcf'
            status_name = directory.strip('_') + '.list'
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
                    #Clean up leftover files after combining them.
                    #Also cleans up .idx files and directories.
                    map(os.unlink, chrlist)
                    #No index file is created for empty vcfs.
                    map(lambda x: os.unlink(x + '.idx'), realchrs)
                    os.unlink(listing)
            os.rmdir(d_path)

if os.path.exists(args.gatkpath):
    vcf_catenate(args.directory, args.reference, args.gatkpath)
elif os.exists('gatk.jar'):
    vcf_catenate(args.directory, args.reference)
else:
    sys.stderr.write('Please provide an existent gatk jar filepath.')
    sys.exit(1)
