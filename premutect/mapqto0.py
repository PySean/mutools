#!/usr/bin/python2.7
"""
'mapqto0.py', by Sean Soderman
Sets the MAPQ value of all unmapped reads to zero for each bam file in the
specified directory.
"""
import pysam
import sys

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
    if not os.exists(bamdir):
        sys.stderr.write('Sorry, but the specified directory does not exist.')
        sys.exit(1)
    for bam in [bam for bam in os.listdir(bamdir) if bam.endswith('.bam')]:
        inbam = pysam.AlignmentFile(bam, 'rb')
        #Template is specified to maintain the same header information.
        outbam = pysam.AlignmentFile('temp.bam', 'wb', template = inbam)
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
