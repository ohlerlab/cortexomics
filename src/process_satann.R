library(svglite)
library(readr)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(stringr)
library(data.table)
library(assertthat)
library(parallel)

filter<-dplyr::filter
slice<-dplyr::slice
filter<-dplyr::filter
filter<-dplyr::filter

source('../src/R/Rprofile.R')

#initiaal form when testing
#let's look at the scores for our
satannfiles <- Sys.glob('SaTAnn/*/*Final_ORFs*')%T>%{stopifnot(length(.)>0)}
satannorfs <- 
	# Sys.glob('SaTAnn/*/*Final_ORFs*')%>%
	satannfiles%>%
	setNames(.,basename(dirname(.)))%>%
	mclapply(load_objs)

#now filter out the weird GRanges columns, for now, and aggregate into one data table 
all_orfs <- satannorfs%>%map(.%>%.$ORFs_tx)
all_orfs_gen <- satannorfs%>%map(.%>%.$ORFs_gen)%>%GRangesList%>%unique

for (i in seq_along(all_orfs))all_orfs[[i]]$Protein%<>%as.character

orfs_dt <- all_orfs%>%map(.%>%{
	issimplelist <- mcols(.)%>%vapply(is,TRUE,'atomic')
	.[,issimplelist]%>%GR2DT
})%>%bind_rows(.id='sample')



##Outputing merged files for uORFs etc.

#get the canonical annotation
annocds <- read_compressed_gfile('my_gencode.vM12.annotation.gtf','CDS')
exons <- read_compressed_gfile('my_gencode.vM12.annotation.gtf','exon')
exons %<>% split(.,.$transcript_id)

#genomic seqinfo object
gseqinfo<-fread('cat *.chromsizes')%>%set_colnames(c('seqnames','seqlengths'))%>%do.call(Seqinfo,.)

#transcript seqinfo object
trseqinfo<-exons%>%width%>%sum%>%enframe%>%set_colnames(c('seqnames','seqlengths'))%>%do.call(Seqinfo,.)

#the genomic locations of all ORFs found by the pipeline
orfgendt<-all_orfs_gen%>%as.list%>%map(.%>%{.$ORF_id_tr<-names(.);GR2DT(.)})%>%bind_rows%>%unique

#filter the table of ORFs to uORFs, then join them 
uorfdt <- orfs_dt %>% 
	filter(ORF_category_Tx=='uORF')  %>% 
	select(-seqnames,-start,-end,-width,-strand)%>%
	left_join(orfgendt)

#now export to a file 
uorfdt_cdsfilt <- uorfdt%>%
	DT2GR(gseqinfo)%>%
	subsetByOverlaps(annocds,invert=TRUE)

#export cdsfilt
uorfdt_cdsfilt %>% export('SaTAnn/uORFs.gtf')



# orfs_dt%<>%group_by(sample)%>%mutate(fdr =  p.adjust(pval, method='fdr'))
# orfs_dt %<>% dplyr::filter(fdr<0.05)%>%ungroup
orfs_dt%>%nrow

#TODO - eveyrthing in linc

n_genes_translated <- orfs_dt %>% group_by(sample)%>%dplyr::summarise(n_genes=n_distinct(gene_id))

readthroughs <- satannorfs%>%map_df(.%>%.$ORFs_readthroughs%>%length)%>%stack%>%set_colnames(c('n_readthroughs','sample'))
readthroughs$sample%<>%as.character
readthroughs$n_readthroughs%<>%as.numeric

genetypehits <- satannorfs%>%lapply(
	.%>%.$ORFs_tx%>%.[T,c('gene_biotype','gene_id')]%>%mcols%>%as.data.frame%>%distinct%>%group_by(gene_biotype)%>%tally%>%
	mutate(gene_biotype = ifelse(gene_biotype%>%str_detect('pseudogene'),'pseudogene',gene_biotype))%>%
	group_by(gene_biotype)%>%summarise(n=sum(n))%>%
	dplyr::filter(gene_biotype %in% c('antisense','lincRNA','protein_coding','pseudogene')))%>%
	bind_rows(.id='sample')%>%
	spread(gene_biotype,n)

#but does our uorf 


satann_summary_table <- genetypehits%>%left_join(readthroughs)

satann_summary_table%>%write_tsv('satann_summary.tsv'%T>%{message(normalizePath(.))})

orfs_dt%>%distinct(gene_id,ORF_category_Tx,.keep_all=T)%>%.$ORF_category_Tx%>%table%>%stack

satannorfs[[1]]$ORFs_tx[,c('gene_biotype','gene_id')]%>%mcols%>%as.data.frame