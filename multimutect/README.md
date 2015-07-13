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
...
tumorN.bam  normalN.bam

Of course, the names do not have to be tumorX.bam and normalX.bam.

Also, for either use case a normal bam file in the pair is optional, as MuTect
works fine with tumor only data.

Summary
-------
```
./multimutect.py [--help] [-m path_to_mutect]
                 -b list_of_bams | -p tumor1.bam:[normal1.bam] tumor2.bam...
                 [-M mutect_options] -f fasta_file [-i input_directory]
                 [-o output_directory] [--process_whole_bam]
```

-  `-h`

   `--help`: Displays basic usage information.
 
- `-m`

  `--mupath`: Specifies the full pathname of the MuTect jar file.
              *Default*: a file called "mutect.jar" in the current 
              working directory.

- `-b`

  `--bamlistfile`: *Required*. Specifies a file listing the BAM files 
                   multimutect will run MuTect on. The format of the 
                   file is described above.

- `-p`

  `--pairs`: *Required*. A colon separated list of tumor:normal pairs as 
             described above. **NOTE**: Only this option or -b 
             must be specified for the list of tumor:normal pairs.

- `-M`

  `--mutectopts`: Extra parameters specific to MuTect. Must be surrounded
                  by single or double quotes.

- `-f`

  `--fasta`: *Required*. FASTA formatted reference sequence. 

- `-i`

  `--inputdir`: The name of the directory the input BAM files are located.
                *Default*: The current working directory.

- `-o`

  `--outputdir`: The name of the directory multimutect will create & 
                 output VCF files to.
                 *Default*: A directory called "output" is created.
- `-w`

  `--process_whole_bam`: Each thread will process an entire BAM file at once
                         instead of multiple chromosomes at a time. This is
                         a good idea for smaller BAM files.

-  `--mem`: The maximum amount of heap memory the Java interpreter should
            use. It is a good idea to increase this when the no downsampling
            MuTect option is specified. Keep in mind that each thread will
            be forking processes that will each have a maximum of this number
            to allocate on the heap. 
            *Default*: 2.

-  `--numthreads`: The number of threads multimutect will create.
                   *Default*: Number of cores on your machine / 4.
                   **Use this option at your own risk!**.


Example Usage
-------

Here, two pairs are specified for processing, one without a normal sample. 
Here the option -dt NONE is given to MuTect to tell it not to use downsampling.

```
./multimutect.py -m /opt/bin/mutect-1.1.7.jar -p tumor.bam:normal.bam tumor2.bam:
-M "-dt NONE" -f hg19.fa -i myBAMs -o BAMresults --process_whole_bam
```

###Memory Usage Notes
This program will only work (usefully) if you have > 4 cores on your machine.
The reason behind this is that I divide the number of cores on the machine
by four. Each thread from multimutect will run MuTect processes which fork
a considerable number of children processes. This will most likely cause
a lot of swapping to happen if I created a thread for each of your machine's cores.

That being said you can manually change the number of threads multimutect
will create, however you will be doing this at your own risk.
