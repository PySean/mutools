#mutools

##Description
A suite of tools for working with MuTect.

##Dependencies
* [pysam](https://pypi.python.org/pypi/pysam)
* [argparse](https://pypi.python.org/pypi/argparse/1.3.0)
   \*if your python version is < 2.7
* [futures](https://pypi.python.org/pypi/futures/3.0.3)
* *multimutect requires that python 2 be used, as pysam is written in it.*

##Required Files
* [MuTect](https://github.com/broadinstitute/mutect)
* [GATK](https://www.broadinstitute.org/gatk/download/)
* [Picard](http://broadinstitute.github.io/picard/)

##Typical Workflow
You'd most likely make use of the tools in premutect to prepare
BAM files first. First you'd use addgroups.py to add (or replace) read groups 
within the BAM files of the directory you specified. Then you would sort them using
reorder\_all.py. You must have a FASTA file sorted in 
coordinate order for this to work, however.

Next, you would build bam indices using reindex.py. And finally,
if MuTect complains about unmapped reads (or if you are simply unsure about them)
you can use mapqto0.py to set MAPQ to 0 for all unmapped reads.

This should be sufficient to prepare your BAM files for processing by MuTect.

The postprocessing steps are different depending on how you go about running
multimutect. If you didn't use the --process\_whole\_bam option, you will need
to run catenate.py on your directory before combine.py, otherwise combine.py
is enough. If you've just processed a directory filled with chromosomal regions
(i.e., chromosome1.bam, chromosome2.bam...)
with the --process\_whole\_bam option, however, you will only need to run catenate.py
with the --listfile option. Either way, appropriate processing will 
produce a single VCF file with columns corresponding to each sample used 
in MuTect's analysis.

Unless you use an option provided by combine.py, you will also have a leftover
directory of VCF files. I leave this choice up to the user, as multiple VCF
files may be desired and the combine phase may be omitted entirely.

Further documentation for each tool is contained within README files inside their
respective directories.
