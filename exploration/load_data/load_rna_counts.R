source = function(path,...){
  require(stringr)
  require(magrittr)
  path%<>%str_replace('/home/dharnett','~')
  message(path)
  base::source(path,...)
}
root='~/projects/cortexomics'

library(rtracklayer)
library(stringr)
library(magrittr)
library(assertthat)
library(tidyverse)
library(data.table)
select =dplyr::select

#this is join, won't multiply columns if repeated, but instead overwrites them
#and also ensures that exactly one match per row is present. 
left_join_ov <- function(LHS,RHS,...,.allow_dups=FALSE,.allow_missing=FALSE){
  require(assertthat)
  require(rlang)
  #arguments as strings
  joincols<-map_chr(quos(...),rlang::quo_text)
  for(j in joincols)  assert_that(has_name(LHS , j ))
  for(j in joincols)  assert_that(has_name(RHS , j ))
  #columns to keep in x - delete those in RHS
  xcols<-colnames(LHS)%>%setdiff(colnames(RHS))%>%c(joincols,.)
  LHS <- select(LHS,!!!xcols)
  #now join
  if(!.allow_missing) stopifnot(nrow(anti_join(LHS,RHS,by=joincols))==0)
  out = left_join(LHS,RHS,by=joincols)
  if(!.allow_dups) stopifnot(nrow(out)==nrow(LHS))
  out
}

# rids <- read_tsv('~/projects/cortexomics/ext_data/riboprotids.tsv')
# ridssplit<-rids%>%filter(!`Mito-RP`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist

file.path(root,"data/kallisto/")%>%normalizePath
#first load kallisto counts
kallistodir = file.path(root,"data/kallisto/")%T>%{stopifnot(file.exists(.))}

kallistocounts <- kallistodir %>%
  list.files(full=TRUE, patt="abundance.tsv",recurs=TRUE)%>%
  # .[[1]]%>%
  setNames(.,basename(dirname(dirname(.))))%>%
  map(read_tsv,col_types=cols())%>%
  map(select,target_id,tpm)%>%
  bind_rows(.id = 'dataset')%T>%
  {assert_that(is.data.frame(.)&(nrow(.)>1e3)&(ncol(.)==3))}
kallistocounts%<>%select(signal = tpm,transcript_id=target_id,dataset)

#we need to ensure each entry of our kallisto table has a unique gene ID
stopifnot(kallistocounts$transcript_id%>%str_count('ENSMUST\\d+')%>%equals(1)%>%all)

#get annotation
annotation_gr <- rtracklayer::import(file.path(root,'data','my_gencode.vM12.annotation.gff3'))

# annotation_gr <- rtracklayer::import('/fast/projects/cubit/0.12.0/static_data/annotation/GENCODE/M12/GRCm38/gencode.vM12.annotation.gff3')

#join info on transcript names and gene names etc to counts
kallistocounts<- kallistocounts %>%left_join_ov(transcript_id,RHS=annotation_gr%>%
                                                  mcols%>%
                                                  as.data.frame%>%
                                                  select(transcript_id,gene_id=matches('gene(ID|_id)'),gene_name,transcript_name=Name)
)

#
kallistocounts<-kallistocounts%>%separate(dataset,c('time','assay','rep'),remove=FALSE)
write_tsv(kallistocounts,file.path(root,'exploration/tables/kallistocounts.tsv'))

kallistocounts<-read_tsv(file.path(root,'exploration/tables/kallistocounts.tsv'),col_types=cols(
  transcript_id = col_character(),
  signal = col_double(),
  dataset = col_character(),
  time = col_character(),
  assay = col_character(),
  rep = col_integer(),
  gene_id = col_character(),
  gene_name = col_character(),
  transcript_name = col_character()
))