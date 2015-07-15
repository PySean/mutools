#!/usr/bin/env python
"""
    "addgroups.py", by Sean Soderman
    Adds read group information to all BAM files in a directory.
"""
import argparse
import os
import subprocess
import sys


parser = argparse.ArgumentParser(description=('Tool for adding read group'
                                             ' information to all files within'
                                             ' a directory'))
parser.add_argument('-d', '--inputdir', type=str,
                    help=('A directory containing the BAM files to add read'
                         'groups to.'), required=True, metavar='input_dir')

parser.add_argument('--info', type=str, 
                    help=('Arguments for '
                          ' AddOrReplaceReadGroups to use.'
                          ' Inserts dummy data by default.'
                          ' The arguments may be in any order you desire,'
                          ' but must be in the form rg_field=value.'),
                    nargs=5, 
                    default=['rgid=group1', 'rglb=lib1', 'rgpl=illumina',
                             'rgpu=unit1', 'rgsm=sample1'],
                    metavar=('rgid', 'rglb', 'rgpl', 'rgpu', 'rgsm'))

parser.add_argument('-p', '--picard_path', type=str,
                    help=('The complete path to the picard jar file.'
                          ' Default: a file called picard.jar in '
                          ' the current working directory'),
                    default='picard.jar', metavar='path_to_picard')
args = parser.parse_args()

picard_cmd = ('java -jar {picard_path} INPUT={{bam}} OUTPUT={{bam}}'
              ' RGID={rgid} RGLB={rglb} RGPL={rgpl} RGPU={rgpu}'
              ' RGSM={rgsm}')
#Transform RG info from command line into dictionary after popping
#it from the argument dict, then add it back to arg dict.
argdict = vars(args)
rgs = {opt: arg for opt, arg in [t.split('=') for t in argdict.pop('info')]}

argdict.update(rgs)

picard_cmd = picard_cmd.format(**argdict)

#I use a generator just in case this is a massive directory...
bams = (os.path.join(args.inputdir, x) for x in 
        os.listdir(args.inputdir) if x.endswith('bam'))

for bam in bams:
    #print(picard_cmd.format(bam=bam)) #Just for testing..
    subprocess.check_output(picard_cmd.format(bam=bam).split())
