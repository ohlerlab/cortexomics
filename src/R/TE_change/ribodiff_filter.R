#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
message('...done')

# defaultargs <- c(
#   ribodiffdone='ribodiff/.done',
#   scaledata='exprdata/cent_scaled_exprdata.txt',
#   lfoldchanges='exprdata/limma_fold_changes.txt',
# )
# args <- coalesce(

args <- commandArgs(trailingOnly=TRUE)
ribodifffolder <- args[1]%>%dirname


ALPHA=0.05

idfile <- args[2]

exprtablefiles <- args[-c(1,2)]
stopifnot(length(exprtablefiles)>0)

outputfolder<-exprtablefiles[1]%>%dirname

ribodifffiltgenes <- Sys.glob(file.path(ribodifffolder,'riboseqres_*.txt'))%>%
  map(.%>%{suppressMessages(read_tsv(.))}%>%filter(padj<ALPHA)%>%.[[1]])%>%unlist%>%unique

ids <- read.table(idfile)%>%distinct
ribodifffiltgenes <- ids[[2]][ids[[1]] %in% ribodifffiltgenes]
assert_that('Satb2' %in% ribodifffiltgenes)


dir.create(outputfolder,showWarn=FALSE,rec=TRUE)

filepostpend<-function(file,pend)paste0(tools::file_path_sans_ext(file),pend,'.',tools::file_ext(file))


for (exprtablefile in exprtablefiles){
  exprtable<-read_tsv(exprtablefile)
  gcol=FALSE
  if('gene' %in% colnames(exprtable))   gcol = 'gene'
  if('gene_name' %in% colnames(exprtable))     gcol = 'gene_name'
  if(gcol==FALSE) next
  passesfilt <- exprtable[[gcol]]%in%ribodifffiltgenes
  assert_that(!all(passesfilt==FALSE))
  message(paste0('Got ',sum(passesfilt),' genes passing filter for ',exprtablefile))
  exprtablefilt <- exprtable[passesfilt,]
  write_tsv(exprtablefilt,filepostpend(exprtablefile,'_ribodfilt'))

}

