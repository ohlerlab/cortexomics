#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(limma))
suppressMessages(library(DESeq2))
message('...done')
filter<-dplyr::filter
select<-dplyr::select
slice<-dplyr::slice

LOWCOUNTLIM <- 10

lsf.str(dplyr)
as.character(lsf.str('package:dplyr')[i])


args <- c(
	countfile='feature_counts/all_feature_counts',
	msfile=file.path('ms_tables/ms_LFQ_total_ms_tall.tsv'),
	transformdexprfile=file.path('exprdata/transformed_data.txt'),
  transformd_scale_cent_exprfile=file.path('exprdata/cent_scaled_exprdata.txt'),
  designmatrixfile=file.path('exprdata/designmatrix.txt'),
  normcountstable='exprdata/allcounts_snorm.tsv'
)
args <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
for(i in names(args)) assign(i,args[i])

# message(str_interp('filtered out ${length(nonmeasured_pids)} protein groups\
#  that were unique to the total MS data, e.g.\n ${sample(specungenes,10)}'))

countstable <- data.table::fread(countfile)

#carry out individual processing of the data sources
ids <- fread('ids.txt')%>%set_colnames(c('feature_id','gene_name'))%>%distinct
countstable <- left_join(countstable,ids)
multidgenes <- countstable$gene_name%>%table%>%keep(~ . > 1)
countstable <- dplyr::filter(countstable,!gene_name %in% multidgenes) 
countstable %<>% select(gene_name,everything())

    
#convert to matrix
countsmatrix <- countstable%>% { set_rownames(as.matrix(.[-(1:2)]),.[[1]]) }

#filter out stuff with very low counts
lowmediancounts <- countsmatrix %>% apply(1,median) %>%`<`(LOWCOUNTLIM)
countsmatrix <- countsmatrix[!lowmediancounts,]

# pdf(file.path('../plots','tmp.pdf')%T>%{normalizePath(.)%>%message})
# countsmatrix%>%apply(1,median)%>%add(0.5)%>%log10%>%hist(breaks=20)
# dev.off()

#and then transform the counts
countsmatrix<-DESeq2::vst(countsmatrix)

countsmatrix_snorm <- countsmatrix %>% {sweep(.,2,STATS = DESeq2::estimateSizeFactorsForMatrix(.),FUN='/')}

countsmatrix_snorm %>%{cbind(gene_name=rownames(.),as_data_frame(.))} %>% write_tsv(normcountstable)

mstable=data.table::fread(msfile)
#some formatting differences
mstable$dataset%<>%str_replace('p5','5')
mstable$dataset%<>%str_replace('_rep','_')
mstable$dataset%<>%str_replace('^[^_]+_','')#no need to annotate what signal type it is
mstable$dataset%<>%str_replace('total','MS')#no need to annotate what signal type it is



message('Taking genes for which every timepoint has at least some information')
#for each gene, take the protein with most signal
mstable_comp <- mstable
mstable_comp %<>% group_by(Protein_IDs)%>%filter(!any(frac_time_missing))


n_filtered_protids <- n_distinct(mstable$Protein_IDs) - n_distinct(mstable_comp$Protein_IDs)
n_filtered_gidids <- n_distinct(mstable$gene_name) - n_distinct(mstable_comp$gene_name)

message(str_interp('filtered out ${n_filtered_protids} protein groups\
for ${n_filtered_gidids} genes which leaves us with ${n_distinct(mstable_comp$gene_name)} gene ids'))


mstable_gene <- 
  mstable_comp%>%semi_join(., 
    group_by(.,gene_name,Protein_IDs)%>%
    summarize(msig=median(signal,na.rm=TRUE))%>%
    ungroup%>%
    arrange(desc(msig))%>%
    distinct(gene_name,.keep_all=TRUE)
  )


msmatrix<-mstable_gene%>%
  ungroup%>%
  select(gene_name,dataset,signal)%>%
  filter(!is.na(gene_name))%>%
  spread(dataset,signal)%>%
  { set_rownames(as.matrix(.[,-1]),.[[1]]) }%>%
  log2


assert_that(! msmatrix%>%colnames%>%str_detect('E17p5')%>%any)
assert_that(! msmatrix%>%colnames%>%str_detect('rep\\d+')%>%any)

#we first need to produce some plots demonstraitng that the cdata are homoskedastic


#print a pot showing homoskedasticity to the 
ms_meanvar_plotname <- basename(msfile)%>%paste0(.,'.pdf')
pdf(file.path('../plots','mean_variance_plots',ms_meanvar_plotname)%T>%message)
vsn::meanSdPlot(msmatrix)
dev.off()

count_meanvar_plotname <- basename(countfile)%>%paste0(.,'.pdf')
pdf(file.path('../plots','mean_variance_plots',count_meanvar_plotname)%>%normalizePath%T>%message)
vsn::meanSdPlot(countsmatrix)
dev.off()

#join the data sources together
exprmatrix <- inner_join(
  countsmatrix %>% as_data_frame %>% cbind(gene_name = row.names(countsmatrix), .),
  msmatrix %>% as_data_frame %>% cbind(gene_name = row.names(msmatrix), .)
)%>%  { set_rownames(as.matrix(.[,-1]),.[[1]]) }

#normalize the various 'libraries' with DESeq norm factors
exprmatrix %<>% sweep(.,2,STATS = estimateSizeFactorsForMatrix(.),FUN='/')
#centered and scaled data
cent_scaled_exprmatrix<-
  exprmatrix%>% 
  sweep(.,1,STATS = rowMeans(.),FUN='-')%>%
  sweep(.,1,STATS = rowSds(.),FUN='/')


#and export
dir.create('exprdata')
exprmatrix %>% as_data_frame %>% cbind(gene_name = row.names(exprmatrix))%>%write_tsv(transformdexprfile)

cent_scaled_exprmatrix %>% as_data_frame %>% cbind(gene_name = row.names(exprmatrix))%>%write_tsv(transformd_scale_cent_exprfile)

designmatrix<-exprmatrix%>%colnames%>%as_data_frame%>%set_colnames('dataset')%>%separate(dataset,c('time','assay','rep'),remove=FALSE)

designmatrix%>%write_tsv(designmatrixfile)
