multimutect
===========

This program runs multiple instances of the MuTect mutation analysis tool on
different chromosomal regions of your input BAM files (by default).

The BAM files used in the analysis can be specified by either a list of tumor:normal
pairs on the command line, like this:

tumor1.bam:normal1.bam tumor2.bam:normal2.bam ... tumorN.bam:normalN.bam

Or in a text file where the tumor/normal pairs are separated by whitespace instead:

tumor1.bam  normal1.bam 
tumor2.bam  normal2.bam

Of course, the names do not have to be tumorX.bam and normalX.bam.

Summary
-------
`
./multimutect.py [--help] [-m path_to_mutect]
                 -b list_of_bams | -p tumor1.bam:[normal1.bam] tumor2.bam...
                 [-M mutect_options] -f fasta_file [-i input_directory]
                 [-o output_directory] [--process_whole_bam]
`

Usage Notes
-----------
Either the -b or -p options *must* be used. The -f option is always required.
When -M is specified, you must give the MuTect specific options 
in single or double quotes, otherwise there is no way multimutect 
can tell whether the options are specific to it or MuTect. I 
recommend single quotes as they suppress  all command line expansion.

The -i option specifies the path to the directory in which all BAM files
reside. By default, this is the directory in which you run multimutect.

The -o option mirrors the -i one, as it specifies an output directory.
It creates a directory called "output" as the output directory, by default.

--process\_whole\_bam has multimutect process entire BAM files with each
thread, rather than chromosome regions. This option is useful for when
you have smaller bam files.

Example Usage
-------
TODO

###Memory Usage Notes (NOTE: Perhaps I should allow the usr to spec numthreads?)
This program will only work (usefully) if you have > 4 cores on your machine.
The reason behind this is that I divide the number of cores on the machine
by four. Each thread from multimutect will run MuTect processes which fork
a considerable number of children processes. This will cause a lot of swapping
to happen if I created a thread for each of your machine's cores.

That being said you can manually change the number of threads multimutect
will create, however you will be doing this at your own risk.
