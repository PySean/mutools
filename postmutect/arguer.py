#!/usr/bin/env python
"""
"arguer.py", by Sean Soderman
Wrap the common command line options in a single function.
This does *not* parse the arguments, this is so that the calling program may
add more arguments if it desires.
"""
import argparse

def makeparser(description):
    parser = argparse.ArgumentParser(description=description)
    input_group = parser.add_mutually_exclusive_group(required=True)

    input_group.add_argument('-d', '--directory', type=str,
                        help='The input directory containing vcf files',)
    parser.add_argument('-r', '--reference', type=str,
                        help='The reference genome for the BAM files',
                        required=True)

    parser.add_argument('-g', '--gatkpath', type=str,
                        help=('The path to the gatk jar file. This looks for'
                        ' gatk.jar in the current working directory by'
                        ' default'), default='gatk.jar'
                       )
    return parser
