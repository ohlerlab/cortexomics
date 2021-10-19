
library(R.utils)
library(GenomicRanges)
library(Biostrings)
library(tidyverse)
args = R.utils::commandArgs(trailingOnly=TRUE,asValues=TRUE,defaults =list(
  pepfile='test_peps.fa',
  protfile='test.prots.fa',
  output='test_pepmatch.tsv'
))
attach(args)

#make sure inputs look right
stopifnot(file.exists(pepfile))
stopifnot(str_detect(pepfile,'.(fa|fasta)$'))
stopifnot(file.exists(protfile))
stopifnot(str_detect(protfile,'.(fa|fasta)$'))

#read in as AA string objects (amino acid sequences)
peptides = readAAStringSet(pepfile)
protseqs = readAAStringSet(protfile)

peptide=peptides[[1]]
#loop over each peptide
pepmatch = mclapply(peptides,mc.cores=detectCores(),function(peptide){
    #get a bioconductor object containin the location of each match
    #in each peptide
    match = vmatchPattern(peptide,fixed=TRUE,protseqs)
    #convert this to a data frame
    match%>%as("GRanges")%>%as.data.frame%>%
      select(protein=seqnames,start,end)
})
#combine all the dataframes
pepmatch <- pepmatch%>%bind_rows(.id='peptide')

#output as a file
dir.create(dirname(output),showWarn=F,rec=TRUE)
pepmatch%>%
  select(peptide,protein,start,end)%>%
  write_tsv(output)

