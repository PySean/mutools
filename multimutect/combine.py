#!/usr/bin/env python
"""
"combine.py", by Sean Soderman
Combines all vcfs *in* a directory using the GATK tool CombineVariants.
"""
import os
import re
import sys
import subprocess as sp
try:
    from arguer import makeparser
except ImportError as I:
    sys.stderr.write('Make sure arguer.py is in my working directory: {}'
                     .format(I))


#TODO: Add option for deleting all vcfs in the directory.
parser = makeparser('Combines all vcfs in a directory with CombineVariants')
parser.add_argument('-o', '--outfile', type=str,
                    help='The resulting combined output file',
                    default='outfile.vcf')
args = parser.parse_args()

"""
Contains a little too much to be contained in a lambda expression.
Returns an argument string with the hyphens separated by 
"""
def vformat(filename): 
    base_filename = os.path.basename(filename).split('.')[0]
    no_hyphen_filename = re.sub('-', '_', base_filename)
    return '-V:{} '.format(no_hyphen_filename)

def vcf_combine(dirname, reference, outfile, gatkpath):
    cmd = ('java -jar {gatk} -T CombineVariants -R {ref}' 
           ' -nt 4 {{vcfs}} -o {outfile} --genotypemergeoption UNSORTED')
    cmd = cmd.format(gatk=gatkpath, ref=reference, outfile=outfile)
    
    vcfs = filter(lambda x: x.endswith('.vcf'), os.listdir(dirname))
    path_vcfs = [os.path.join(dirname, s) for s in vcfs]
    #Use the basenames of the vcf files for more descriptive 'set='
    #areas in the output VCF file.
    vstring = ' '.join(map(lambda x: vformat(x) + x, path_vcfs))
    try:
        sp.check_output(cmd.format(vcfs=vstring).split())
    except sp.CalledProcessError as cpe:
        sys.stderr.write('CombineVariants ran into a problem: {}'.format(cpe))
        sys.exit(1)

vcf_combine(args.directory, args.reference, args.outfile, args.gatkpath)
