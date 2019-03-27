

library(rtracklayer)
library(GenomicFiles)
library(Rsamtools)
library(Biostrings)

args <- c('pipeline/my_gencode.vM12.annotation.gff3','pipeline/tRNAs.gtf')
args <- commandArgs(trailingOnly = TRUE) 
annot <- args[1] %>% import
genome <- FaFile(args[2])
trnas <- args[3]


#takes a length, a region, a genome, optional scoring later?



simulate_reads <- function(regions, readlength, genome){



	#reads

	readfile <- 
}





#call it on our regions