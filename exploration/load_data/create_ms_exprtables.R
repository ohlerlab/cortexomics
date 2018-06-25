#!/usr/bin/env Rscript
suppressMessages(library(limma))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
defaultargs <- c(
	totalmsfile='/fast/groups/cubi/projects/2017-10-10_cortexomics/gdrive/cortexomics_ms_total/325_new_2_all_log2_LFQ_n7464.txt',
	specmsfile='/fast/groups/cubi/projects/2017-10-10_cortexomics/gdrive/cortexomics_ms_cyt+80+Poly/proteinGroups.txt',
	outfolder=file.path('/fast/groups/cubi/projects/2017-10-10_cortexomics','pipeline','ms_tables')
)

args <- coalesce(
	commandArgs(trailingOnly=TRUE)[1:length(defaultargs)]%>%setNames(names(defaultargs)),
	defaultargs
)
for(i in names(args)) assign(i,args[i])
#function to delete japanese fromt he tables
delete_japanese<-function(dblstring) str_extract(dblstring,'(NA)|(NaN)|(\\d+\\.?\\d*)')%>%as.numeric

#process koshi's table
ms_data=
  read_tsv(totalmsfile,comment = '#')%>%
  {colnames(.) = str_replace_all(colnames(.),' ','_');.}%>%
  select(gene_name=Gene_names,everything())%>%
  mutate_at(vars(matches('^(E\\d\\d|P0)')),funs(delete_japanese))
#label the LFQ columns
colnames(ms_data) <- str_replace(colnames(ms_data),'^(E\\d\\dp?\\d?|P0)_n','LFQ_\\1_total_rep')

#convert them back from their log2 values
for(lfqcol in colnames(ms_data)%>%str_subset('LFQ')) ms_data[[lfqcol]] %<>% {2^.}

#now read in the specific fracitons in the MS
ms_data_spec=
  read_tsv(specmsfile,comment = '#')%>%
  {colnames(.) = str_replace_all(colnames(.),' ','_');.}%>%
  select(gene_name=Gene_names,everything())%>%
  mutate_at(vars(matches('^(E\\d\\d|P0)')),funs(delete_japanese))

#format the column names
colnames(ms_data_spec) <- str_replace(colnames(ms_data_spec),'^(E\\d\\d|P0)','LFQ_intensity_\\1')
colnames(ms_data_spec) <- str_replace(colnames(ms_data_spec),'input_rep','cyto_rep')
colnames(ms_data_spec) <- str_replace(colnames(ms_data_spec),'^LFQ_intensity','LFQ')


assert_that(!anyDuplicated(ms_data$`Protein_IDs`))
assert_that(!anyDuplicated(ms_data_spec$`Protein_IDs`))

#join, use majority protein IDs
ms_data_all <- full_join(ms_data,ms_data_spec,by='Majority_protein_IDs')
ms_data_all%<>%select(Protein_IDs=Majority_protein_IDs,everything())
# ms_data_all <- full_join(ms_data,ms_data_spec,by='Protein_IDs')
# ms_data_all%<>%select(Protein_IDs=Protein_IDs,everything())


commoncols<-intersect(colnames(ms_data),colnames(ms_data_spec))
totalcols<-setdiff(colnames(ms_data),colnames(ms_data_spec))
speccols<-setdiff(colnames(ms_data_spec),colnames(ms_data))

totalungenes<-setdiff(ms_data$gene_name,ms_data_spec$gene_name)
specungenes<-setdiff(ms_data_spec$gene_name,ms_data$gene_name)

message(str_interp('found ${length(totalungenes)} genes that were unique to the total MS data, e.g.\n ${sample(totalungenes,10)}'))
message(str_interp('found ${length(specungenes)} genes that were unique to the total MS data, e.g.\n ${sample(specungenes,10)}'))

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


#signal column should be numeric
ms_tall%<>%mutate(signal = as.numeric(signal))

intenssigs<-c('LFQ','iBAQ','Intensity')
ms_tall <- ms_tall %>%
  mutate(signal = ifelse((signal==0)&(sigtype%in%intenssigs) ,NA,signal))%>%
  group_by(time,replicate)%>%
  # mutate(signal = (1e6*signal) / sum (na.omit(signal)) )
  identity


#make outlying low values into NAs
ms_tall$nearzero <- log10(ms_tall$signal) < -4 & (ms_tall$sigtype=='LFQ')
stopifnot( between ( mean(na.omit(ms_tall%>%filter(sigtype=='LFQ')%>%.$nearzero)) , 0 , 0.04))
ms_tall%<>%mutate(signal=ifelse(nearzero ,NA,signal))
#clean this up with gene names
ms_tall%<>%left_join(ms_data_all%>%select(Protein_IDs,gene_name.x),by='Protein_IDs')
colnames(ms_tall)%<>%str_replace('gene_name.x','gene_name')

#CATAGORIZE GENES by amount of missing data
ms_tall<-ms_tall%>%
  group_by(Protein_IDs,fraction)%>%
  mutate(fraction_missing = sum(!is.na(signal))<=1)%>%
  group_by(Protein_IDs,fraction,time)%>%
  mutate(frac_time_missing = sum(!is.na(signal))<=1)


# ms_tall%>%ungroup%>%filter(gene_name=='Pa2g4',sigtype=='LFQ')%>%qplot(data=.,x=time,y=signal,ylab='LFQ',main='MS for Pa2g4 in different fractions')+facet_grid(~fraction)

#how many genes are just missing most or all of their proteomics data?
ms_tall%>%group_by(Protein_IDs,)%>%filter(!fraction_missing)%>%distinct(Protein_IDs,gene_name,fraction)%>%
  group_by(Protein_IDs,gene_name)%>%tally%>%.$n%>%table

#get protein groups not measured
nonmeasured_pids <- ms_tall%>%group_by(Protein_IDs)%>%filter(all(fraction_missing))%>%distinct(Protein_IDs)%>%.$Protein_IDs

message(str_interp('filtered out ${length(nonmeasured_pids)} protein groups\
 that were unique to the total MS data, e.g.\n ${sample(specungenes,10)}'))


#now filter out those missing genes
ms_tall %<>% filter(!Protein_IDs %in% nonmeasured_pids)


####################################

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
ms_tall$replicate%<>%str_replace('rep','')
stopifnot(sigtypes%>%unique%>%str_detect('norm_iBAQ')%>%any)
for(sigtypei in sigtypes){
	for(fractioni in unique(ms_tall$fraction)){
	  mstallfile <- file.path(outfolder,str_interp('ms_${sigtypei}_${fractioni}_ms_tall.tsv'))
	  dir.create(dirname(mstallfile),showWarn=FALSE)
	  write_tsv(ms_tall%>%filter(sigtype==sigtypei,fraction==fractioni),mstallfile%T>%message)
	}
}

