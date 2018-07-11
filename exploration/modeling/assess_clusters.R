#!/usr/bin/env Rscript
message('loading libraries')
library(warn.conflicts = FALSE,quietly=TRUE,tidyverse)
library(warn.conflicts = FALSE,quietly=TRUE,magrittr)
library(warn.conflicts = FALSE,quietly=TRUE,stringr)
library(warn.conflicts = FALSE,quietly=TRUE,data.table)
library(warn.conflicts = FALSE,quietly=TRUE,assertthat)
library(warn.conflicts = FALSE,quietly=TRUE,limma)
library(warn.conflicts = FALSE,quietly=TRUE,DESeq2)
message('...done')

defaultargs <- c(
  clusterfolder='clusters/limma_fold_changes/tsne/',
  outputfolder='clsuter_assessement/limma_fold_changes/tsne/'
)
args <- coalesce(
  commandArgs(trailingOnly=TRUE)[1:length(defaultargs)]%>%setNames(names(defaultargs)),
  defaultargs
)
for(i in names(args)) assign(i,args[i])

args[1]%>%
	list.files(full=T,pattern='clusters')