premutect
=========
These are some Python scripts that aid in the preprocessing of BAM files for
MuTect.

##addgroups
This script allows one to easily add read groups to all BAM files within
the specified directory.

##Summary
```
./addgroups.py -d input_dir 
              [-h] [--info rgid rglb rgpl rgpu rgsm]
              [-p path_to_picard]
```
- ``
  
  ``:

- `-d`

  `--inputdir`: *Required*: Specifies the directory containing the input BAM
                files to be modified.

- `--info`: Arguments that map exactly to the required arguments for
            AddOrReplaceReadGroups, save for the input and output
            files. They must be provided in the same format as well,
            but in lowercase.

            rglb=lib1 rgid=group1 rgpl=illumina rgsm=sample1  rgpu=unit1

            Take note that these parameters can be specified *in any order*.
- `-p`
  
  `--picard_path`: The complete path to the picard jar file, defaults to a file
                   called "picard.jar" in the working directory.

###Example Usage
```
./addgroups.py --inputdir bam_files rgid=grouper rglb=libby rgpu=unit7 rgsm=sample11
               --picard_path /opt/picard.jar
```
Here, I've opted to provide my own read-group information. The abbreviated
rg.. commands correspond to the following items:

* rgid: Read Group ID.
* rglb: Read Group Library.
* rgpl: Read Group platform.
* rgpu: Read Group platform unit.
* rgsm: Read Group sample name.

##reorder_all
This script reorders all BAM files in the specified directory according to
the ordering of chromosomes in the supplied FASTA file.

##Summary
```
./reorder_all.py -f fasta -d directory
              [-h] [-p path_to_picard]
```

- `-f`
  
  `--fasta`: *Required*: The FASTA file that will be used for reordering.

- `-d`
  
  `--directory`: *Required*: The directory of BAM files to reorder.

- `-p`
  
  `--picardpath`: The full pathname of Picard. Defaults to a file named picard.jar
   in the working directory.

####reindex and mapqto0
These two scripts are rather short. Both take the directory of
bam files as their only positional arguments (invoked as ./reindex.py wham.bam)
