#!/usr/bin/env python
"""
    "combine.py", by Sean Soderman
    Combines all vcfs under a directory in an order according to the chrs.list
    file within the directory.
"""
import os
import sys
from subprocess import check_output

if len(sys.argv) < 4:
    sys.stderr.write('Usage: {} <directory> <reference> <outputfile>\n'
                     .format(sys.argv[0]))
    sys.exit(1)

"""
Combines all vcfs under a directory. By default, does not remove the vcf files.
"""
def vcf_combine(parent, reference, outfile='out', gatkpath='gatk.jar'):
    cmd = ('java -cp {gatk} -assumeSorted -R' 
          '{ref} -out {out} {{vcfs}}').format(gatk=gatkpath, 
                                                  ref=reference, out=outfile)
    for dirpath, dirnames, filenames in os.walk(parent):
        for directory in dirnames:
            d_path = os.path.join(dirpath, directory)
            listing = os.path.join(d_path, 'chrs.list')
            if os.path.exists(listing):
                with open(listing, 'r') as chrlist:
                    vcat = lambda x: '-V ' + x + ' '
                    #The list of -V commands to execute.
                    vseries = "".join(map(vcat, [c.strip() for c in chrlist]))
                    cmd = cmd.format(vcfs = vseries)
                    check_output(cmd.split())

vcf_combine(sys.argv[1], sys.argv[2])
