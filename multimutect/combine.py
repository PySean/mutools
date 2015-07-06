#!/usr/bin/env python
"""
    "combine.py", by Sean Soderman
    Combines all vcfs under a directory in an order according to the chrs.list
    file within the directory.
"""
import os
import sys
from subprocess import check_output

#TODO: Either insert this to the end of multimutect, or make it a true
#standalone program with argparse.
if len(sys.argv) < 4:
    sys.stderr.write(('Usage: {} <directory> <reference>' 
                      ' <outputfile> <gatkpath>\n').format(sys.argv[0]))
    sys.exit(1)

"""
Combines all vcfs under a directory. By default, does not remove the vcf files.
"""
def vcf_combine(parent, reference, outfile='out', gatkpath='gatk.jar'):
    cmd = ('java -cp {gatk} org.broadinstitute.gatk.tools.CatVariants'
           ' -assumeSorted -R' 
           ' {ref} -out {out} {{vcfs}}').format(gatk=gatkpath, 
                                                ref=reference, out=outfile)
    for dirpath, dirnames, filenames in os.walk(parent):
        for directory in dirnames:
            d_path = os.path.join(dirpath, directory)
            listing = os.path.join(d_path, 'chrs.list')
            if os.path.exists(listing):
                with open(listing, 'r') as chrfile:
                    chrlist = [os.path.join(d_path, c.strip()) + '.vcf' 
                               for c in chrfile]
                    #The list of vcfs to concatenate, each prepended by
                    #-V.
                    vseries = "".join(['-V ' + c + ' ' for c in chrlist])
                    cmd = cmd.format(vcfs = vseries)
                    print(cmd)
                    check_output(cmd.split())
                    #Clean up leftover files after combining them.
                    print('I would have deleted {}\n'.format(''.join(chrlist)))
                    #map(os.unlink, chrlist)

vcf_combine(sys.argv[1], sys.argv[2], sys.argv[3])
