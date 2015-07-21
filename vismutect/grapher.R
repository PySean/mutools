#!/usr/bin/Rscript
#
# "grapher.R", by Sean Soderman
# Program that will read in + visualize data from VCF files, pertaining to
# allele frequency. It will count certain mutations from certain contexts
# specified as options on the command line).
#

#library('ggplot2')

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
   print(paste('Usage:', './grapher.R', '<vcf>', '<refseq>', '<fai>'))
   stop(status=1)
}

#Function that reads a file connection to a fai file
#Returns: A list wrapping a list indexed by contig names
#that are each associated with a cumulative byte offset, as well
#as the line length in bytes for this fasta file.
fai_fields <- function(fai_file) {
   contigs <- list()
   col_names <- c('contig', 'sz', 'pos', 'line_bp', 'line_bytes')
   fai_table <- read.table(fai_file, sep='\t', 
                           col.names=col_names)

   #Since it's read in as a factor, convert the contig col to char vectors.
   fai_table$contig <- as.character(fai_table$contig)
   accum <- 0
   #Utilized for calculating the progressive file offset from the
   #contig column. 2 is added to account for the newline and '>'.
   getaccum <- function(x) {
      if (is.character(x)) {
         accum <<- accum + length(strsplit(x, split='')[[1]]) + 2
         accum
      }
      else{
         NA
      }
   }
   contigs[fai_table$contig] <- sapply(fai_table$contig, getaccum)
   close(fai_file)
   return(list(contigs=contigs, 
               line_nts=fai_table$line_nts[1], 
               line_bytes=fai_table$line_bytes[1]))
}

#Utility function for making regex matches.
#Returns: the character vector composed of desired matches.
ezmatch <- function(regex, string) {
   match_obj <- regexec(regex, string)
   matches <- regmatches(string, match_obj)
   if (length(matches) > 0) {
      return(matches[[1]])
   }else {
      return(NULL)
   }
}

#Function that (currently) retrieves chromosomes, along with these fields:
#CHROM
#POS
#REF
#ALT
#Allele frequency (NOTE: I don't know if there is a way to make VarScan
#VCF files generate this field. At the moment this will only do for MuTect)
vcf_info <- function(vcf_line) {
   fields <- list(chrom='', pos=0, ref='', alt='', freq='')
   mstr <- '(\\S*)\\s+(\\S*)\\s+\\S*\\s+(\\S*)\\s+(\\S*).*AF=([0-9.,]*);'
   matches <- ezmatch(mstr, vcf_line) 
   fields$chrom <- matches[2]
   fields$pos <- as.numeric(matches[3])
   fields$ref <- matches[4]
   fields$alt <- matches[5]
   frequencies <- sapply(strsplit(matches[6], ','), as.numeric)
   fields$freq <- as.vector(frequencies)
   return (fields)
}

#Seeks through the FASTA file for desired nucleotide triple, consisting
#of the leading nucleotide, the ref nt a SNP/indel was compared with,
#and the trailing nucleotide.
#Returns: A string that looks like I_J, where I is the chromosome in front
#and J is the chromosome behind the reference chromosome.
get_context <- function(genome_pos, offset, line_bytes, line_nts, ref_file) {
   #Calculate offset to the chromosome behind the SNP.
   chr_off <- floor(genome_pos / line_nts) + (genome_pos %% line_nts) +
              offset - 2 #Get single context nucleotide before SNP/indel.
   #Seek to the position. If there's a newline, step back twice and try again.
   #Make sure to skip newlines.
   seek(ref_file, where=chr_off, origin='start')
   char <- rawToChar(readBin(ref_file, what=raw(), n=1))
   if (char == '\n') 
      seek(ref_file, where=-2, origin='current')
   counter <- 0
   context_string <- character(3)
   while (counter < 3) {
      char <- rawToChar(readBin(ref_file, what=raw(), n=6))
      if (char != '\n') {
         context_string[counter] <- ifelse(counter == 1, char, '_')
         counter <- counter + 1
      }
   }
   return(context_string)
}

#Somewhat conspicuous name, this function reads through a VCF file so it
#can then seek to the location of the indel and increment "bins" based on
#what the change is.
#Returns: A data structure containing the counts and frequencies
#of mutations. 
gather_data <- function(fai_data, vcf_name, ref_name) {
   contigs <- fai_data$contigs
   line_bytes <- fai_data$line_bytes
   line_nts <- fai_data$line_nts
   vcf_file <- file(vcf_name, 'r')
   ref_file <- file(ref_name, 'rb', raw=TRUE)
   line <- character()
   repeat { #Skip header lines.
      line <- readLines(vcf_file, n=1)
      if (ezmatch('^#*', line) == '')
         break
   }
   #Read in pertinent information.
   repeat {
      if (length(line) == 0)
         break #TODO: Return approp. data structure.
      info <- vcf_info(line)
      offset <- contigs[[info$chrom]]
      context <- get_context(info$pos, offset, line_bytes, line_nts,  ref_file)
      line <- readLines(vcf_file, n=1)
   }
   close(vcf_file)
   close(ref_file)
}

vcf <- args[1]
refseq <- args[2]
fai_file <- args[3]
gather_data(fai_fields(file(fai_file, 'r')), vcf, refseq)
