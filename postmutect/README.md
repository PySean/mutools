#postmutect
Here are two scripts that allow for the necessary postprocessing steps of
multimutect's output.

If you used the default options with multimutect (that is, allowing for 
chromosome-by-chromosome processing), you need only use combine.py on the
resulting directory. Otherwise, you will need to use catenate.py on it first
to glue all the chromosome-specific vcf files together.

##catenate
Concatenates all vcf "pieces" generated from a BAM file together.

##Summary
```
catenate.py -d input_dir -r fasta
              [-h] [-g gatkpath] [--delete_fragments] 
              [-l sample_list_file]
```
- `-d`

  `--directory`: The input directory containing the chromosome segment VCFS
  for each BAM sample pair.


- `-r`

  `--reference`: The reference genome for the BAM files, in FASTA format.

- `-g`

  `--gatkpath`: The path to the gatk jar file. This is a file named 'gatk.jar'
  in the current working directory, by default.

- `--delete_fragments`: Deletes all files within the directory used
  for catenation.
 
- `-l`

  `--listfile`: A list specifying the the samples to concatenate. 
   Useful for when you have VCF files generated from BAM files that correspond
   to each chromosome region.
   *Default*: A file called "chrs.list" contained within each sample's VCF
   output directory.

##combine
The command you'd use after running multimutect with the --process\_whole\_bam
option, or the one you would use after catenate if you didn't use that option.

##Summary
```
combine.py (-d input_dir | -v vcf_files) -r fasta
              [-h] [-g gatkpath] [--delete_input_vcfs] 
              [-l sample_list_file] [-w]
```

- `-d`

  `--directory`: The directory of vcf files to combine.

- `-r`

   `--reference`: The reference genome used for the BAM files.

- `-g`

   `--gatkpath`: The path to the gatk jar file.
   *Default*: Looks for a file called gatk.jar in the current working directory.

- `-v`

   `--vcf_files`: A list of vcfs on the command line to combine. There can be 
   as many of these as you wish.

- `-o`

  `--outfile`: The desired name of the combined output file.
  *Default:* outfile.vcf.

- `-D`

  `--delete_input_vcfs`: Delete all input vcfs in the input directory, after
  they have been used. *Use with caution*.

- `-l`

  `--listing`: If this is supplied, combine VCF files in the order given
  in the file containing tumor normal pairs.

- `-w`

  `--without_nonecol`: If this is specified, the column containing 'none'
  will be omitted. *This is probably what you want, most of the time.*
