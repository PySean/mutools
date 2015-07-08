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


parser = makeparser('Combines all vcfs in a directory with CombineVariants')
parser.add_argument('-o', '--outfile', type=str,
                    help='The resulting combined output file',
                    default='outfile.vcf')

parser.add_argument('-D', '--delete', action='store_true',
                    help=('If this is supplied, delete all vcfs in the'
                    ' directory'))
                    
parser.add_argument('-l', '--listing', type=str,
                    help=('If this is supplied, combine VCF files in the order'
                    'given in the listing file.'))

args = parser.parse_args()

"""
Contains a little too much to be contained in a lambda expression.
Returns an argument string with the hyphens separated by underscores.
"""
def vformat(filename): 
    base_filename = os.path.basename(filename).split('.')[0]
    no_hyphen_filename = re.sub('-', '_', base_filename)
    return '-V:{} '.format(no_hyphen_filename)

"""
Combines all vcfs in a directory with CombineVariants.
"""
def vcf_combine(directory, reference, outfile, gatkpath, delete, listing):
    cmd = ('java -jar {gatk} -T CombineVariants -R {ref}' 
           ' -nt 4 {{vcfs}} -o {outfile} --genotypemergeoption UNSORTED')
    cmd = cmd.format(gatk=gatkpath, ref=reference, outfile=outfile)
    
    vcfs = []
    if listing is not None:
        with open(listing, 'r') as bamlist:
            vcfs = [x.strip() for x in bamlist if re.search('.bam\s*$', x)]
            vcfs = map(lambda x: re.sub('.bam', '.vcf', x), vcfs)
    else:
        vcfs = filter(lambda x: x.endswith('.vcf'), os.listdir(directory))

    path_vcfs = [os.path.join(directory, s) for s in vcfs]
    #Use the basenames of the vcf files for more descriptive 'set='
    #areas in the output VCF file.
    vstring = ' '.join(map(lambda x: vformat(x) + x, path_vcfs))
    try:
        sp.check_output(cmd.format(vcfs=vstring).split())
    except sp.CalledProcessError as cpe:
        sys.stderr.write('CombineVariants ran into a problem: {}\n'
                         .format(cpe))
        sys.exit(1)
    if delete == True:
        map(lambda x: os.unlink(x), path_vcfs)

vcf_combine(**vars(args))
