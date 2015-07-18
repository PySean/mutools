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
   print(paste0('Usage: ', './grapher.R', '<vcf>', '<refseq>', '<fai>'))
   stop(status=1)
}

#Function that greps for the appropriate fields in the fa.fai file.
#The fields are (in this order):
#contig, size(in bp), location (in bytes), line length(bp), line length(bytes)
fai_fields <- function(fai_line) {
   fields <- list(contig='', size=0, location=0, lenbp=0, lenbyte=0)
   mstr <- '(\\w+)\\s+(\\w+)\\s+(\\w+)\\s+(\\w+)\\s+(\\w+)\\s*'
   mobj <- regexec(mstr, fai_line)
   #Unfortunately regmatches wraps the char vector in a list..
   matches <- regmatches(fai_line, mobj)[[1]]
   #Return a list containing each item, named appropriately.
   #A bit clunky but I can't used named captures with regexec.
   fields$contig <- matches[2]
   fields$size <- as.numeric(matches[3])
   fields$location <- as.numeric(matches[4])
   fields$lenbp <- as.numeric(matches[5])
   fields$lenbyte <- as.numeric(matches[6])
   return (fields)
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
   mstr <- ' (\\S*)\\s+(\\S*)\\s+\\S*\\s+(\\S*)\\s+(\\S*).*AF=([0-9.,]*);'
   mobj <- regexec(mstr, vcf_line)
   matches <- regmatches(vcf_line, mobj)[[1]]
   fields$chrom <- matches[2]
   fields$pos <- as.numeric(matches[3])
   fields$ref <- matches[4]
   fields$alt <- matches[5]
   frequencies <- sapply(strsplit(matches[6], ','), as.numeric)
   fields$freq <- as.vector(frequencies)
   return (fields)
}
