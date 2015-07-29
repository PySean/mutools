#!/usr/bin/env python
"""
"combine.py", by Sean Soderman
Combines all vcfs *in* a directory using the GATK tool CombineVariants.

To accomodate the case of only a few VCF files needing combining,
the user may also specify VCFs directly on the command line.
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
parser.add_argument('-D', '--delete_input_vcfs', action='store_true',
                    help=('If this is supplied, delete all vcfs in the'
                    ' directory'))
                    
parser.add_argument('-l', '--listing', type=str,
                    help=('If this is supplied, combine VCF files in the order'
                    ' given in the listing file.'))
parser.add_argument('-w', '--without_nonecol', action='store_true',
                    help=('If this option is specified, the column containing'
                          ' "none" will be omitted.'))

args = parser.parse_args()

"""
Contains a little too much to be contained in a lambda expression.
Returns an argument string with the hyphens replaced by underscores.
"""
def vformat(filename): 
    base_filename = os.path.basename(filename).split('.')[0]
    no_hyphen_filename = re.sub('-', '_', base_filename)
    return '-V:{} '.format(no_hyphen_filename)

"""
Combines all vcfs in a directory with CombineVariants.
"""
def vcf_combine(directory, vcf_files, reference, outfile, gatkpath, 
                delete_input_vcfs, listing, without_nonecol):
    #Use the cwd if VCF files are specified on the command line.
    if directory is None:
        directory = '.'
    cmd = ('java -jar {gatk} -T CombineVariants -R {ref}' 
           ' -nt 4 {{vcfs}} -o {outfile}'
           ' -dt NONE --genotypemergeoption UNSORTED')
    cmd = cmd.format(gatk=gatkpath, ref=reference, outfile=outfile)
    
    vcfs = []
    if listing is not None:
        with open(listing, 'r') as bamlist:
            bamonly = [x for x in bamlist if re.search('\.bam\s*$', x)]
            #Make respective tumors/normals lists.
            tumors, normals = zip(*[x.split() for x in bamonly])
            #Substitute .bam in the tumor element for an underscore
            #and append the corresponding normal if y isn't blank, else
            #just use the tumor filename (.bam -> .vcf in either case)
            vcfs = map(lambda tumor, normal: (re.sub('.bam', '_', tumor) + 
                       re.sub('.bam', '.vcf', normal)
                       if not normal.isspace() 
                       else re.sub('.bam', '.vcf', tumor)), tumors, normals)
            vcfs = [x.strip() for x in vcfs]
    elif directory != '.':
        vcfs = filter(lambda x: x.endswith('.vcf'), os.listdir(directory))
    else: #VCF files were specified on the cmd line.
        vcfs = filter(lambda x: os.path.exists(x), os.listdir(directory))

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
    if delete_input_vcfs:
        map(lambda x: os.unlink(x), path_vcfs)
    if without_nonecol is True:
        rm_nonecol(outfile, 'newfile.vcf')
"""
Deletes the column called "none" in the resulting combined VCF.
This results from tumor only samples being provided.
Overwrites the vcf file with the new vcf file without the none column.
"""
def rm_nonecol(infile, noneless):
    with open(infile, 'r') as orig, open(noneless, 'w') as new:
        none_pos = -1
        new_line = []
        vcf_match = re.compile('\s*#*[0-9_a-zA-Z/.:;,=-]+')
        for line in orig:
            if re.match('##', line):
                new.write(line)
            else:
                new_line = re.findall(vcf_match, line)

            if re.match('#CHROM', line) and not re.search('none', line):
                sys.stderr.write(("I am sorry, but you don't"
                                  " even need to use this switch, as your file"
                                  " has no 'none' column in it.\n"))
                sys.exit(1)
            elif re.search('#CHROM.*none', line):
                none_pos = map(lambda x: x.strip(), new_line).index('none')
            if none_pos != -1:
                new_str = ''.join([x for x in new_line 
                                  if new_line.index(x) != none_pos])
                #Portable newline catenation if the newline was removed.
                if not new_str.endswith(os.linesep):
                    new_str += os.linesep
                new.write(new_str)
    #Replace the input file with the file omitting the none column.
    os.rename(noneless, infile)
vcf_combine(**vars(args))
