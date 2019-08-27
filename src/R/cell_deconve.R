#Cell deconve

library(data.table)
library(tidyverse)
library(data.table)
library(magrittr)
library(data.table)
library('DeconRNASeq')

signatures <- '/fast/groups/ag_ohler/work/dharnet_m/cortexomics/ext_data/scharma_etal_proteomic/'

scharmaxls <- '/fast/groups/ag_ohler/work/dharnet_m/cortexomics//ext_data/scharma_etal_proteomic/nn.4160-S15.xlsx'
scharmafile <- '/fast/groups/ag_ohler/work/dharnet_m/cortexomics/pipeline/../ext_data/scharma_etal_proteomic/nn.4160-S15.munged.txt'

exprdata <- here('pipeline/exprdata/transformed_data.txt')%>%fread
gnamestbl <- here('pipeline/ids.txt')%>%fread
exprdata %<>% inner_join(gnamestbl)


if(!file.exists(scharmafile)){
  cells<-tidyxl::xlsx_cells(scharmaxls)
  #munge it
  library(unpivotr)
  scharma_munged <- cells%>%
    select(-sheet,-address,-comment,-height,-width)%>%
    behead('NNW','cellsource')%>%
    behead('NNW','sigtype')%>%
    behead('N','colname')%>%
    mutate(colname = paste0(sigtype,':',colname))%>%
    select(data_type,numeric,character,row,colname)%>%
    spatter(colname)
  write_tsv(scharma_munged,scharmafile)
  scharma <- scharma_munged
}else{
  scharma<-read_tsv(scharmafile)
}
#filter only very differential genes
scharma <- scharma[scharma$`Log2 Fold Expression:>ten fold expressed in atleast one cell type`%in%'+',]
colnames(scharma)%<>%str_replace('NA:Gene names','gene_name')
scharma <- scharma%>%inner_join(gnamestbl)%>%distinct(gene_id,.keep_all = T)

LFQsigs <- scharma%>%select(matches('Log 2 LFQ Intensity.*Isolated'))%>%{colnames(.)%<>%str_replace('Log 2 LFQ Intensity:Isolated','');.}%>%
    as.matrix%>%
    {rownames(.)<-scharma$gene_id;.} %>%
    as.data.frame

LFQsigs_cultured <- scharma%>%select(gene_id,matches('Log 2 LFQ Intensity'))%>%select(-matches('Isolated'))%>%{colnames(.)%<>%str_replace('Log 2 LFQ Intensity:','');.}%>%
  gather(celltype,signal,-gene_id)%>%
  mutate(signal = ifelse(is.nan(signal),0,signal))%>%
  spread(celltype,signal)%>%
  {set_rownames(as.matrix(select(.,-gene_id)),.$gene_id)} %>%
  as.data.frame

tegenes <- Sys.glob('/fast/groups/ag_ohler/work/dharnet_m/cortexomics/pipeline/ribodiff//riboseqres*.txt')%>%map(fread)%>%bind_rows%>%filter(padj<0.05)%>%.$geneIDs
ntegenes <- Sys.glob('/fast/groups/ag_ohler/work/dharnet_m/cortexomics/pipeline/ribodiff//riboseqres*.txt')%>%
  map(.%>%fread%>%{colnames(.)[7]<-'l2fc';.})%>%
  bind_rows%>%
  filter(padj>0.05)%>%
  filter(between(l2fc,-0.025,0.025))%>%
  .$geneIDs

#LFQsigs%>%rownames%>%is_in(tegenes)%>%table

#exprdata<-exprdata%>%select(-gene_id,-gene_name)%>%as.data.frame%>%set_rownames(exprdata$gene_id)
#impute data by mean and put back in df form
exprdata <- 
    exprdata%>%
    select(-gene_name)%>%
    gather(dataset,signal,-gene_id)%>%separate(dataset,into=c('time','assay','rep'))%>%
    group_by(time,assay)%>%
    filter(!sum(is.na(signal))==3)%>%
    mutate(signal = ifelse(is.na(signal),mean(signal,na.rm=T),signal))%>%
    unite(dataset,time,assay,rep)%>%
    spread(dataset,signal)%>%
    {set_rownames(as.data.frame(select(.,-gene_id)),.$gene_id)}

exprdata %>%head

#exprdata %>%select('')
stopifnot(mean(rownames(exprdata)%in%rownames(LFQsigs))>0.9)
stopifnot(mean(rownames(LFQsigs)%in%rownames(exprdata))>0.5)

deconresults <- DeconRNASeq(exprdata,LFQsigs)
deconresults_cult <- DeconRNASeq(exprdata,LFQsigs_cultured)




extractdecon <- .%>%set_rownames(colnames(exprdata))%>%as.data.frame%>%rownames_to_column('sample')%>%
  gather(celltype,proportion,-sample)%>%
  separate(sample,into=c('time','assay','rep'))%>%
  mutate(time = factor(time,unique(time)))

plotdecon <- . %>%{ggplot(data=.,aes(color=celltype,y=proportion,x=as.numeric(time)))+
    geom_point()+geom_smooth()+scale_x_continuous(name='Time Point',labels=levels(plotdf$time))+
    facet_grid(~assay)+theme_bw()+
    ggtitle('')}

plotdf<-deconresults$out.all%>%extractdecon 
plotdecon(plotdf)

plotdf<-deconresults_cult$out.all%>%extractdecon 
plotdecon(plotdf)

deconresults_cult <- DeconRNASeq(exprdata[rownames(exprdata)%in%tegenes,],LFQsigs_cultured)
plotdf_te<-deconresults_cult$out.all%>%extractdecon 

deconresults_cult <- DeconRNASeq(exprdata[rownames(exprdata)%in%ntegenes,],LFQsigs_cultured)
plotdf_nte<-deconresults_cult$out.all%>%extractdecon 

pdf(here('plots/deconv/deconv.pdf')%T>%{dir.create(dirname(.));message(normalizePath(.))})
ggpubr::ggarrange(plotdecon(plotdf_te)+ggtitle('TE Changing Genes'),plotdecon(plotdf_nte)+ggtitle('non TE Changing Genes'),nrow=2)
dev.off()

exprdata[gnamestbl%>%filter(gene_name=='Aif1')%>%pluck('gene_id'),]



## Please refer our demo
source("DeconRNASeq.R")
### multi_tissue: expression profiles for 10 mixing samples from multiple tissues
data(multi_tissue)

datasets <- x.data[,2:11]
signatures <- x.signature.filtered.optimal[,2:6]
#proportions <- fraction

#DeconRNASeq(datasets, signatures, proportions, checksig=FALSE, known.prop = TRUE, use.scale = TRUE)
#
datasets%>%rownames
DeconRNASeq(datasets,signatures)

gnamestbl%>%filter(gene_name=='iba1')

exprdata[]


exprdata

#' Initial tests with count data vs LFQs form Scharma, on log scale, are weird - mostly neurons, but this prop goes down in favor of oligodendrocytes, mainly

