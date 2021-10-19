sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
proDAfitms<-readRDS('data/proDAfitms.rds')
sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
sel_ms_mat<-readRDS('data/sel_ms_mat.rds')

countpred_df<-readRDS('data/countpred_df.rds')
countcontr_df<-readRDS('data/countcontr_df.rds')
tx_countdata<-readRDS('data/tx_countdata.rds')

ms_metadf <- readRDS('data/ms_metadf.rds')

library(tidyverse)

library(txtplot)

countpred_df %>% 
	filter(contrast%>%str_detect('E175_ribo'))%>%
	select(contrast,gene_id,logFC)%>%
	separate(contrast,c('time','assay'))%>%
	left_join(sel_prodpreds%>%select(time,gene_id,diff),by=c('gene_id','time'))%>%
	filter(logFC > quantile(logFC,.1))%>%
	filter(diff > quantile(diff,.1,na.rm=T))%>%
	mutate(SDR = diff - logFC)%>%
	mutate(SDR = log10(2^SDR))%>%
	.$SDR%>%{.=.-mean(.,na.rm=T);.}%>%na.omit%>%sd


allratios = countpred_df %>% 
	filter(contrast%>%str_detect('E175_ribo'))%>%
	inner_join(ms_metadf%>%filter(!is_rpl|is_rps)%>%select(gene_id))%>%
	select(contrast,gene_id,logFC)%>%
	separate(contrast,c('time','assay'))%>%
	left_join(sel_prodpreds%>%select(time,gene_id,diff),by=c('gene_id','time'))%>%
	# filter(logFC > quantile(logFC,.1))%>%
	# filter(diff > quantile(diff,.1,na.rm=T))%>%
	mutate(SDR = diff - logFC)%>%
	mutate(SDR = log10(2^SDR))%>%
	.$SDR

rpratios = countpred_df %>% 
	filter(contrast%>%str_detect('E175_ribo'))%>%
	inner_join(ms_metadf%>%filter(is_rpl|is_rps)%>%select(gene_id))%>%
	select(contrast,gene_id,logFC)%>%
	separate(contrast,c('time','assay'))%>%
	left_join(sel_prodpreds%>%select(time,gene_id,diff),by=c('gene_id','time'))%>%
	# filter(logFC > quantile(logFC,.1))%>%
	# filter(diff > quantile(diff,.1,na.rm=T))%>%
	mutate(SDR = diff - logFC)%>%
	mutate(SDR = log10(2^SDR))%>%
	.$SDR

t.test(allratios,rpratios)

	txtdensity

countpred_df %>% 
	filter(contrast%>%str_detect('E175_total'))%>%
	select(contrast,gene_id,logFC)%>%
	separate(contrast,c('time','assay'))%>%
	left_join(sel_prodpreds%>%select(time,gene_id,diff),by=c('gene_id','time'))%>%
	filter(logFC > quantile(logFC,.1))%>%
	filter(diff > quantile(diff,.1,na.rm=T))%>%
	mutate(SDR = diff - logFC)%>%
	mutate(SDR = log10(2^SDR))%>%
	.$SDR%>%{.=.-mean(.,na.rm=T);.}%>%na.omit%>%sd

	txtdensity

