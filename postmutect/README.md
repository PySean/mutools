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
./catenate.py -d input_dir -r fasta
              [-h] [-g gatkpath]
```
- `-d`

  `--directory`: The output directory containing the chromosome segment VCFS
  for each BAM sample pair.


- `-r`

  `--reference`: The reference genome for the BAM files, in FASTA format.

- `-g`

  `--gatkpath`: The path to the gatk jar file. This is a file named 'gatk.jar'
  in the current working directory, by default.

##combine
The command you'd use after running multimutect with the --process\_whole\_bam
option, or the one you would use after catenate if you didn't use that option.

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