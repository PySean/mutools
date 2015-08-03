#!/usr/bin/Rscript
#
#"grapher.R", by Sean Soderman
#Program that will read in + visualize data from VCF files, pertaining to
#allele frequency. It will count mutations from any context.
#The desired allele frequency threshold can be specified on the command line
#as a range, ex: 0 <= AF <= .5.
#By default, the alt allele with the maximum allelic frequency is 
#used. Otherwise, SNPS falling out of the specified range are ignored.

library('ggplot2')
library('reshape2')

args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
   print(paste('Usage:', './grapher.R', 'vcf', 'refseq', 'fai', '[frequency]'))
   stop(status=1)
}

#Function that reads a file connection to a FASTA index file.
#Returns: A descriptive data frame representative of the data in the file.
fai_fields <- function(fai_file) {
   contigs <- list()
   col_names <- c('contig', 'sz', 'pos', 'line_nts', 'line_bytes')
   fai_table <- read.table(fai_file, sep='\t', 
                           col.names=col_names)
   #Since it's read in as a factor, convert the contig col to char vectors.
   fai_table$contig <- as.character(fai_table$contig)
   close(fai_file)
   return(fai_table)
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
   fields <- list(chrom='', pos=0, ref='', alt='', freqs=0)
   mstr <- '(\\S*)\\s+(\\S*)\\s+\\S*\\s+(\\S*)\\s+(\\S*).*AF=([0-9.,]*);'
   matches <- ezmatch(mstr, vcf_line) 
   fields$chrom <- matches[2]
   fields$pos <- as.numeric(matches[3])
   fields$ref <- matches[4]
   fields$alt <- matches[5]
   fields$freqs <- as.numeric(strsplit(matches[6], ',')[[1]])
   return (fields)
}

#Creates the "bin" data structure necessary for graphing.
#Returns: A list of data frames. 
#The list is indexed by the "change" string,
#matrices are indexed by the preceding nucleotide (row) and the succeeding
#nucleotide (col)
make_bins <- function() {
   nucleotides <- c('A', 'C', 'G', 'T')
   matlist <- list()
   matlab <- NA
   snps <- character(12) #4 choose 3 possibilities for DNA snps.
   contexts <- character(16) #4 choose 4 possilities for nt contexts of len=1.

   #Initialize & create the matrix that will get copied copiously.
   #Create the factors (SNP possibilities)
   snp_ndx <- 1
   for (i in 1:length(nucleotides)) {
      ref <- nucleotides[i]
      for (j in nucleotides[nucleotides != ref]){
         alt <- j
         snps[snp_ndx] <- paste0(ref, '>', alt)
         snp_ndx <- snp_ndx + 1
      }
   }
   #Now for the column names.
   context_ndx <- 1
   for (i in 1:length(nucleotides)) {
      preceding_nt <- nucleotides[i]
      for (j in 1:length(nucleotides)) {
         succeeding_nt <- nucleotides[j]
         contexts[context_ndx] <- paste0(preceding_nt, '_', succeeding_nt)
         context_ndx <- context_ndx + 1
      }
   }
   #Create blank data frame.
   #Since d.frames don't have "nrow" and ncol, I have to create a matrix 1st.
   matlab <- matrix(data=rep(rep(0, 12), 16), nrow=12, ncol=16)
   matlab <- as.data.frame(matlab)
   #Add the names, then bind the snp factors in as the leftmost column.
   names(matlab) <- contexts
   matlab <- cbind(snps=as.factor(snps), matlab)
   return(matlab)
}

#Seeks through the FASTA file for desired nucleotide triple, consisting
#of the leading nucleotide, the ref nt a SNP/indel was compared with,
#and the trailing nucleotide.
#Returns: A character vector of the SNP, with context.
get_context <- function(genome_pos, offset, line_bytes, line_nts, ref_file) {
   #Calculate offset to the chromosome behind the SNP.
   genome_rem <-  genome_pos %% line_nts
   chr_off <- (floor((genome_pos - 1) / line_nts) * line_bytes) + 
              ifelse(genome_rem, genome_rem, line_nts) +
              offset - 2 #Get single context nucleotide before SNP.

   seek(ref_file, where=chr_off, origin='start')
   counter <- 1
   context_string <- character(3)
   #Check for initial newline and seek back twice if one is encountered.
   char <- rawToChar(readBin(ref_file, what=raw(), n=1))
   if (char == '\n') {
      seek(ref_file, where=-2, origin='current')
   } else {
      seek(ref_file, where=-1, origin='current')
   }
   while (counter <= 3) {
      char <- rawToChar(readBin(ref_file, what=raw(), n=1))
      if (char != '\n') {
         #In the rare case an ambiguous nt is picked up, swap it with another.
         if (char == 'N' || char == 'n') {
            char <- sample(c('A', 'G', 'C', 'T'), 1)
         }
         context_string[counter] <- char
         counter <- counter + 1
      } 
   }
   return(context_string)
}

#Checks the frequency of the alt allele(s), returning the allele
#with the highest frequency that falls within the frequency range.
check_freq <- function(inequality, alts, freq) {
      if ( (length(inequality) == 0) || inequality[1] == 'freq') {
         print(paste('Please specify the range for allelic frequencies',
                     'appropriately.'))
         stop(status=1)
      }
      #Evaluate expression for each alt allele.
      #If more than two satisfy it, return the greater of the two.
      expr_left <- parse(text=paste(inequality[2], 'freq'))
      expr_right <- parse(text=paste('freq', inequality[3]))
      passing_freqs <- freq[eval(expr_left) & eval(expr_right)]
      #If no allele passes the frequency test, return an empty vector.
      if (length(passing_freqs) > 0)
         max_freq <- max(passing_freqs)
      else
         return(numeric(0))
      #If two alleles are at the max and have the same frequency,
      #return one of them at random.
      return(sample(x=alts[max_freq == freq], size=1))
}
#Somewhat conspicuous name, this function reads through a VCF file so it
#can then seek to the location of the SNP and increment "bins" based on
#what the change is.
#Returns: A data structure containing the counts and frequencies
#of mutations. 
gather_data <- function(fai_data, vcf_name, ref_name, range) {
   snp_bins <- make_bins()
   inequality <- ezmatch('(.*)freq(.*)', range)
   #line_bytes <- fai_data$line_bytes
   #line_nts <- fai_data$line_nts
   vcf_file <- file(vcf_name, 'r')
   ref_file <- file(ref_name, 'rb', raw=TRUE)
   line <- character()
   repeat { #Skip header lines.
      line <- readLines(vcf_file, n=1)
      if (ezmatch('^#*', line) == '')
         break
   }
   #Read in pertinent information, perform binning.
   repeat {
      if (length(line) == 0)
         break 
      info <- vcf_info(line)
      offset <- fai_data[fai_data$contig == info$chrom,]$pos
      line_nts <- fai_data[fai_data$contig == info$chrom,]$line_nts
      line_bytes <- fai_data[fai_data$contig == info$chrom,]$line_bytes
      context <- get_context(info$pos, offset, line_bytes, line_nts,  ref_file)
      context <- paste0(context[1], '_', context[3])
      #Some SNPS have multiple possibilities.. account for them.
      alts <- strsplit(info$alt, split=',')[[1]]
      #By default pick the alt allele with the highest
      #allelic frequency. If the command line option for AF was used however,
      #use the highest allelic frequency used in that range.
      alt <- ifelse(range != '', check_freq(inequality, alts, info$freqs),
                    max(alts))
      if (length(alt) == 0) #Frequency interval unsatisfied, skip line.
         next
      snp <- paste0(info$ref, '>', alt)
      snp_bins[snp_bins$snps == snp,][[context]] <- 
              snp_bins[snp_bins$snps == snp,][[context]] + 1
      line <- readLines(vcf_file, n=1)
   }
   close(vcf_file)
   close(ref_file)
   return(snp_bins)
}


vcf <- args[1]
refseq <- args[2]
fai_file <- args[3]
freq_range <- args[4]

histdata <- gather_data(fai_fields(file(fai_file, 'r')), vcf, refseq,
                        ifelse(is.na(freq_range), '', freq_range))
#Transform data into "long" format suitable for graphing.
histmelt <- melt(histdata, id.vars='snps', variable.name='context',
                 value.name='count')
legend <- ifelse(!is.na(freq_range), sub(pattern='freq', 
                                     x=freq_range, 
                                     replacement=' AF '), '')
#Graph the data!
if (legend != '') {
qplot(x=context, y=count, data=histmelt, facets=snps ~ . ,
       geom='bar', fill=snps, stat='identity') + 
       theme(axis.text.y=element_text(size=9), legend.position='none') + 
       ggtitle(legend)
} else {
qplot(x=context, y=count, data=histmelt, facets=snps ~ . ,
       geom='bar', fill=snps, stat='identity') + 
       theme(axis.text.y=element_text(size=9), legend.position='none') +
       ggtitle('SNPs and their contexts')
}
