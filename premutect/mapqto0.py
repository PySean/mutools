#!/usr/bin/env python
"""
'mapqto0.py', by Sean Soderman
Sets the MAPQ value of all unmapped reads to zero for each bam file in the
specified directory.
"""
import sys
import os
#Some preamble code:
#Check for packages that have pysam as a subdirectory and push them to
#the end of sys.path to avoid any possibility of the wrong pysam version
#being used.
for i in range(0, len(sys.path)):
    traversed = []
    for dirpath, dirnames, filenames in os.walk(sys.path[i]):
        for dir in dirnames:
            if 'pysam' == dir and not dirpath.endswith('site-packages'):
                if dirpath not in traversed:
                    sys.path.append(sys.path.pop(i))
                    traversed.append(dirpath)
try:
    from pysam import AlignmentFile
except ImportError as I:
    sys.stderr.write('Please install the necessary packages: {}'
                     .format(I))
    sys.exit(1)
    

if len(sys.argv) < 2:
    usestr = "Usage: {} <bamdir> \n"
    sys.stderr.write(usestr.format(sys.argv[0]))
    sys.exit(1)

def umappedq2zero(bamdir):
    """
    Reads in a BAM file, setting the MAPQ value for an alignment segment
    to zero if it is unmapped.
    Opens up both infile and outfile and outputs these modified
    reads to outfile.
    """
    if not os.path.exists(bamdir):
        sys.stderr.write('Sorry, but the specified directory does not exist.')
        sys.exit(1)
    for bam in [bam for bam in os.listdir(bamdir) if bam.endswith('.bam')]:
        inbam = AlignmentFile(bam, 'rb')
        #Template is specified to maintain the same header information.
        outbam = AlignmentFile('temp.bam', 'wb', template = inbam)
        #Construct reads iterator using fetch.
        reads = inbam.fetch(until_eof = True)
        for read in reads:
            if read.is_unmapped == True:
                read.mapping_quality = 0
            outbam.write(read) #Don't omit any reads!
        #Overwrite the original with the new file with MAPQs set to zero.
        os.rename('temp.bam', os.path.join(bamdir, bam))

if __name__ == '__main__':
    umappedq2zero(sys.argv[1])
