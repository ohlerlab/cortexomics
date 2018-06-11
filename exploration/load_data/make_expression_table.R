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


######Now load ms data

delete_japanese<-function(dblstring) str_extract(dblstring,'(NA)|(NaN)|(\\d+\\.?\\d*)')%>%as.numeric


ms_data=
  read_tsv('/fast/groups/cubi/projects/2017-10-10_cortexomics/gdrive/cortexomics_ms_total/325_new_2_all_log2_LFQ_n7464.txt',comment = '#')%>%
  {colnames(.) = str_replace_all(colnames(.),' ','_');.}%>%
  select(gene_name=Gene_names,everything())%>%
  mutate_at(vars(matches('^(E\\d\\d|P0)')),funs(delete_japanese))
#label the LFQ columns
colnames(ms_data) <- str_replace(colnames(ms_data),'^(E\\d\\dp?\\d?|P0)_n','LFQ_\\1_total_rep')

for(lfqcol in colnames(ms_data)%>%str_subset('LFQ')) ms_data[[lfqcol]] %<>% {2^.}
# ms_data %<>% mutate_at(vars(matches('LFQ')),funs(~ 2 ^.))

# read_tsv('/fast/groups/cubi/projects/2017-10-10_cortexomics/gdrive/cortexomics_ms_cyt+80+Poly/proteinGroups.txt',comment = '#')
# read_tsv('/fast/groups/cubi/projects/2017-10-10_cortexomics/gdrive/cortexomics_ms_cyt+80+Poly/MS_Ribo+Total_Proteome_raw data.txt',comment = '#')%>%
#     colnames%>%str_subset('E13')
# 
# 
# cbind(ms_data%>%filter(gene_name=='Pa2g4')%>%select(matches('LFQ'))%>%t,
#       ms_data%>%filter(gene_name=='Pa2g4')%>%select(matches('iBAQ_'))%>%t)%>%set_colnames(c('LFQ','iBAQ'))
# 
# 
# cbind(ms_data%>%filter(gene_name=='Pa2g4')%>%select(matches('LFQ'))%>%t,
#       ms_data%>%filter(gene_name=='Pa2g4')%>%select(matches('iBAQ_'))%>%t)%>%set_colnames(c('LFQ','iBAQ'))%>%{cor(.[,1],.[,2],method='s')}
# 
# cbind(ms_data%>%filter(gene_name=='Pa2g4')%>%select(matches('LFQ'))%>%t%>%{2^.},
#       ms_data%>%filter(gene_name=='Pa2g4')%>%select(matches('iBAQ_'))%>%t)%>%set_colnames(c('LFQ','iBAQ'))%>%
#   as.data.frame%>%
#   gather(coltype,signal)%>%qplot(data=.,color=coltype,y=signal)


ms_data_spec=
  read_tsv('/fast/groups/cubi/projects/2017-10-10_cortexomics/gdrive/cortexomics_ms_cyt+80+Poly/proteinGroups.txt',comment = '#')%>%
  {colnames(.) = str_replace_all(colnames(.),' ','_');.}%>%
  select(gene_name=Gene_names,everything())%>%
  mutate_at(vars(matches('^(E\\d\\d|P0)')),funs(delete_japanese))

colnames(ms_data_spec) <- str_replace(colnames(ms_data_spec),'^(E\\d\\d|P0)','LFQ_intensity_\\1')
colnames(ms_data_spec) <- str_replace(colnames(ms_data_spec),'input_rep','cyto_rep')
colnames(ms_data_spec) <- str_replace(colnames(ms_data_spec),'^LFQ_intensity','LFQ')


intersect(colnames(ms_data),colnames(ms_data_spec))
setdiff(colnames(ms_data),colnames(ms_data_spec))
setdiff(colnames(ms_data_spec),colnames(ms_data))

setdiff(ms_data$gene_name,ms_data_spec$gene_name)
setdiff(ms_data_spec$gene_name,ms_data$gene_name)

stopifnot(!anyDuplicated(ms_data$`Protein_IDs`))
stopifnot(!anyDuplicated(ms_data_spec$`Protein_IDs`))
# 
# #join, use majority protein IDs
# ms_data_all <- full_join(ms_data,ms_data_spec,by='Majority_protein_IDs')
# ms_data_all%<>%select(Protein_IDs=Majority_protein_IDs,everything())

ms_data_all <- full_join(ms_data,ms_data_spec,by='Protein_IDs')
ms_data_all%<>%select(Protein_IDs=Protein_IDs,everything())


ms_data$gene_name%>%str_subset(regex('foxp2|CAGH44|SPCH1|TNRC10|2810043D05Rik|D0Kist7',ignore_case=TRUE))
ms_data$gene_name%>%str_subset(regex('fezf2|fezl|Zfp312',ignore_case=TRUE))
ms_data$gene_name%>%str_subset(regex('cux2|cutl2|cux-2',ignore_case=TRUE))


#get data in tall format, one signal value per row
ms_tall <- 
  ms_data_all %>%
  # select(gene_name,matches('iBAQ_'))%>%
  filter(!is.na(Protein_IDs))%>%
  select(Protein_IDs,matches('_(E|P)[0-9p\\.]+_[^_]+_rep\\d'))%>%
  select(-matches('Identification_type'))%>%
  gather(dataset,signal,-Protein_IDs)%>%
  identity


#parse the dataset column in into multiple fields and join back to the main df
ms_tall<-
  ms_tall%>%
  mutate(dataset = str_replace(dataset,'(.*)_\\[%\\]$','Perc_\\1')) %>%
  {left_join(by='dataset',.,
             distinct(.,dataset)%>%
               tidyr::extract(dataset,
                              into =  c('sigtype','time','fraction','replicate'),
                              regex = '(.*)_([EP][0-9p\\.]+)_([^_]+)_(rep\\d$)',
                              remove=FALSE
               )
  )
  }


#check out the na columns
ms_tall%>%filter(is.na(as.numeric(signal)))%>%distinct(sigtype,signal)
#mutate the columns that are 
ms_tall%<>%mutate(signal = as.numeric(signal))



#convert 0 values to NAs
intenssigs<-c('LFQ','iBAQ','Intensity')
ms_tall <- ms_tall %>%
  mutate(signal = ifelse((signal==0)&(sigtype%in%intenssigs) ,NA,signal))%>%
  group_by(time,replicate)%>%
  # mutate(signal = (1e6*signal) / sum (na.omit(signal)) )
  identity
# 

# ms_tall%>%filter(fraction=='total',time=='E16',gene_name==unames[5])


#make outlying low values into NAs
# ms_tall%>%group_by(fraction,Protein_IDs)%>%mutate()
# ms_tall$nearzero <- log10(ms_tall$signal) < -4 & (ms_tall$sigtype=='LFQ')
# stopifnot( between ( mean(na.omit(ms_tall%>%filter(sigtype=='LFQ')%>%.$nearzero)) , 0 , 0.04))
# ms_tall%<>%mutate(signal=ifelse(nearzero ,NA,signal))

#clean this up with gene names
ms_tall%<>%left_join(ms_data_all%>%select(Protein_IDs,gene_name.x),by='Protein_IDs')
colnames(ms_tall)%<>%str_replace('gene_name.x','gene_name')

ms_tall%>%ungroup%>%filter(gene_name=='Pa2g4',sigtype=='LFQ')%>%qplot(data=.,x=time,y=signal,ylab='LFQ',main='MS for Pa2g4 in different fractions')+facet_grid(~fraction)

#i need to characterize genes by their presence of data
ms_tall<-ms_tall%>%
  group_by(Protein_IDs,fraction)%>%
  mutate(fraction_missing = sum(!is.na(signal))<=1)%>%
  group_by(Protein_IDs,fraction,time)%>%
  mutate(frac_time_missing = sum(!is.na(signal))<=1)

#how many genes are just missing most or all of their proteomics data?
ms_tall%>%group_by(Protein_IDs,)%>%filter(!fraction_missing)%>%distinct(Protein_IDs,gene_name,fraction)%>%
  group_by(Protein_IDs,gene_name)%>%tally%>%.$n%>%table

#get protein groups not measured
nonmeasured_pids <- ms_tall%>%group_by(Protein_IDs)%>%filter(all(fraction_missing))%>%distinct(Protein_IDs)%>%.$Protein_IDs

#now filter out those missing genes
ms_tall %<>% filter(!Protein_IDs %in% nonmeasured_pids)

totalmslfq<-ms_tall%>%
  filter(fraction=='total')%>%select(-fraction)%>%
  filter(sigtype=='LFQ')%>%select(-sigtype)%>%
  mutate(replicate==str_replace(replicate,'rep',''))%>%
  identity 

########Normalized iBAQs

message('normalizing iBAQ')
mscolinfo = ms_tall%>%ungroup%>%distinct(dataset,sigtype,time,fraction,replicate)


rowcol = "Protein_IDs"
sigcol = "signal"
colcols = c('time','fraction','replicate')

#get matrix of all data for size factors
ibaqcols <-mscolinfo%>%filter(sigtype=='iBAQ')%>%.$dataset

library(DESeq2)
ibaq_sizefactors = ms_data_all%>%select(Protein_IDs,one_of(ibaqcols))%>%
  {set_rownames(as.matrix(.[,-1]),.[[1]])}%>%
  DESeq2::estimateSizeFactorsForMatrix(.)
ibaq_sizefactors = ibaq_sizefactors%>%stack%>%set_colnames(c('sizefactor','tmp'))%>%
  separate(tmp,into=c('time','fraction','replicate'))

ibnormms_tall <- ms_tall%>%filter(sigtype=='iBAQ')%>%left_join(ibaq_sizefactors)%>%mutate(sigtype='norm_iBAQ',signal=signal/sizefactor)%>%select(-sizefactor)
#now add in the 
ms_tall = bind_rows(ms_tall%>%filter(!sigtype=='norm_iBAQ'),ibnormms_tall)

#output files for all of the ms signal types.
sigtypes = ms_tall$sigtype%>%unique
stopifnot(sigtypes%>%unique%>%str_detect('norm_iBAQ')%>%any)
for(sigtype in sigtypes){
  mstallfile <- file.path(root,'tables','ms',sigtype,'ms_tall.tsv')
  dir.create(dirname(mstallfile),showWarn=FALSE)
  write_tsv(ms_tall%>%filter(sigtype==sigtype),mstallfile)
}
