################################################################################
########Various diagnostic plots attempting to understand the MuSiC results
################################################################################
  

source('src/Figures/Figure5/1_MuSiC.R')

#let's plot music proportions
mbasis <- music_basis(sscountexprdata,clusters='cellType',sample='sampleID')
bulkabund_long <- exprs(bestcountexprdata)%>%as_tibble(rownames='gene_name')%>%gather(dset,bulk_abundance,-gene_name)
sc_abund_long <- mbasis[[1]]%>%as_tibble(rownames='gene_name')%>%gather(ftset,sc_abundance,-gene_name)

pred_expr_df <- music_reslist[[1]][[1]] %>% as.data.frame%>%rownames_to_column('dset')%>%gather(ftset,music_prop,-dset)%>%
  left_join(sc_abund_long,by='ftset')%>%
  group_by(dset,gene_name)%>%
  summarise(pred_abundance = sum(sc_abundance * music_prop))%>%
  left_join(bulkabund_long%>%select(dset,gene_name,bulk_abundance))

stem_prop_df  <- music_reslist[[1]][[1]] %>% as.data.frame%>%rownames_to_column('dset')%>%gather(ftset,music_prop,-dset)%>%
  left_join(sc_abund_long,by='ftset')%>%
  mutate(isstem=str_detect(ftset,'1H'))%>%
  group_by(dset,gene_name,isstem)%>%
  summarise(pred_abundance = sum(sc_abundance * music_prop))%>%
  group_by(dset,gene_name)%>%summarise(prop_stem = pred_abundance[isstem]/sum(pred_abundance))
 

pred_expr_df%<>%mutate(techangegene = gene_name%in%techangegenes)

#now plot
plotfile<- here(paste0('plots/','musicplots','.pdf'))
pdf(plotfile)
pred_expr_df%>%filter(is.finite(log2(pred_abundance)),is.finite(log2(bulk_abundance)))%>%
  arrange(desc(techangegene))%>%
  filter(bulk_abundance>32)%>%
  filter(dset%>%str_detect('E13|E175'),dset%>%str_detect('_1'))%>%
  ggplot(data=.,aes(x=log2(pred_abundance),y=log(bulk_abundance)))+
  scale_x_continuous('predicted abundance (MuSic)')+
  scale_y_continuous('Actual Bulk abundance')+
  facet_wrap(dset~.)+
  geom_point(aes(color=techangegene),size=I(.2))+
  theme_bw()
dev.off()
normalizePath(plotfile)
stopifnot(exists('stem_prop_df'))

if(!file.exists(here('data/sel_prodpreds.rds'))){
  source(here('src/Figures/Figure0/2_proDA.R'))
}
sel_prodpreds<-readRDS(here('data/sel_prodpreds.rds'))
if(exists('countpred_df')){
  load('data/1_integrate_countdata.R')
}

prediction_df = bind_rows(
  sel_prodpreds%>%
    mutate(assay='MS')%>%
    select(gene_id,time,assay,estimate,CI.L,CI.R),
  countpred_df%>%
    separate(contrast,c('time','assay'))%>%
    select(gene_id,time,assay,estimate=logFC,CI.L,CI.R)
)

stem_prop_df%<>%mutate(techangegene = gene_name%in%techangegenes)
stem_prop_df %<>% separate(dset,c('time','assay','rep'),remove=F)%>%
  inner_join(prediction_df%>%filter(assay=='TE')%>%
  mutate(gene_name=gid2gnm[[gene_id]]),by=c('time','gene_name'))


#now plot
plotfile<- here(paste0('plots/','stemprop_vs_te','.pdf'))
pdf(plotfile)
print(
stem_prop_df %>%
  filter(dset%>%str_detect('E16_total_1'),dset%>%str_detect('_1'))%>%
  ggplot(data=.,aes(x=(prop_stem),y=estimate))+
  scale_x_continuous('predicted proportion of RNA in stem cells')+
  scale_y_continuous('TE')+
  facet_wrap(dset~.)+
  # geom_smooth(method='lm')+
  # geom_point(aes(color=techangegene),size=I(.2))+
  geom_point(aes(color=str_detect(gene_name,'^Rp[sl]')),size=I(.2))+
  # geom_point(aes(color=str_detect(gene_name,'^Rp[sl]')),size=I(.2))+
  theme_bw()
  )
dev.off()
normalizePath(plotfile)




stem_prop_df%>%
  filter(rep==1,assay.x=='total')%>%
  group_by(gene_name)%>%
  mutate(estimate = estimate - estimate[time=='E13'])%>%
  {cor.test(.$estimate,.$prop_stem,method='pearson')}


################################################################################
########Try to understand
################################################################################
library(ggrepel)

library(LSD)
source('Applications/LSD/R/LSD.heatscatter.R')
#if it was two cell types we could just plot like, 
comptp1='24H'
comptp2='96H'
telleystemcountmat = cbind(
  exprs(sscountexprdata)[,sscountexprdata@phenoData$ftdiff==comptp1]%>%rowSums,
  exprs(sscountexprdata)[,sscountexprdata@phenoData$ftdiff==comptp2]%>%rowSums
)
telleystemcountmat <- telleystemcountmat[telleystemcountmat%>%rowMins%>%`>`(32),]
telleystemcountmat = telleystemcountmat%>%sweep(.,2,FUN='/',STAT= colSums(.))
stemcountweights = log2(telleystemcountmat) %>% {.[,2]-.[,1]}
#

#now plot
plotfile<- here(paste0('plots/','telleyweightsvste','.pdf'))
pdf(h=14,w=14,plotfile)
# print(
stemcountweights%>%enframe('gene_name','diff_weight')%>%
# wteltcoregsdiff%>%enframe('gene_name','diff_weight')%>%
  inner_join(prediction_df%>%filter(assay=='TE')%>%mutate(gene_name=gid2gnm[[gene_id]]),by=c('gene_name')) %>%
  # filter(gene_name %in% telleycoregenesdiff)%>%
  mutate(lbl = ifelse(diff_weight>3,gene_name,''))%>%
  left_join(techangedf%>%transmute(gene_name,dTE = (up==1|(down==1))))%>%
  filter(!is.na(dTE))%>%{
  . = filter(.,time=='E13');
  heatscatter((.$diff_weight),.$estimate,ggplot=TRUE)+
  geom_smooth(method='lm')+
  # geom_label_repel()+
  scale_x_continuous('Ratio of Expression in Neurons to Stem Cells as per Telley')+
  # scale_x_continuous('Expression weight Neurons to Stem Cells as per Telley')+
  scale_y_continuous('TE')+
  facet_wrap(time+assay~.)+
  geom_point(size=I(1.0))+
  theme_bw()+
  ggtitle(paste0(comptp1,' vs. ',comptp2))
  }
dev.off()
normalizePath(plotfile)


techangedf<-read_tsv('tables/xtailTEchange.tsv') 
#now plot

plotfile<- here(paste0('plots/','telleyweightsvsdte','.pdf'))
pdf(plotfile)
print(
stemcountweights%>%enframe('gene_name','diff_weight')%>%
  inner_join(countcontr_df%>%filter(assay=='TE',!is.na(time))%>%mutate(gene_name=gid2gnm[[gene_id]]),by=c('gene_name')) %>%
  ggplot(data=.,aes(x=(diff_weight),y=logFC))+
  scale_x_continuous('Ratio of Expression in Neurons to Stem Cells as per Telley')+
  scale_y_continuous('dTE')+
  facet_wrap(time+assay~.)+
  geom_point(size=I(1.0))+
  theme_bw()
  )
dev.off()
normalizePath(plotfile)


################################################################################
########now weight birthdate
################################################################################

telleybdcountmat = cbind(
  exprs(sscountexprdata)[,sscountexprdata@phenoData$fttime=='E15']%>%rowSums,
  exprs(sscountexprdata)[,sscountexprdata@phenoData$fttime=='E12']%>%rowSums
)
telleybdcountmat <- telleybdcountmat[telleybdcountmat%>%rowMins%>%`>`(32),]
bdcountweights = log2(telleybdcountmat) %>% {.[,2]-.[,1]}
#now plot
plotfile<- here(paste0('plots/','telleyweightsvste','.pdf'))
pdf(plotfile)
print(
bdcountweights%>%
  enframe('gene_name','BD_weight')%>%
  mutate(techangegene = gene_name%in%techangegenes)%>%
  inner_join(prediction_df%>%filter(assay=='TE')%>%
  mutate(gene_name=gid2gnm[[gene_id]]),by=c('gene_name')) %>%
  ggplot(data=.,aes(x=(BD_weight),y=estimate))+
  scale_x_continuous('Ratio of Expression in E15 to 12 Cells as per Telley')+
  scale_y_continuous('TE')+
  facet_wrap(time+assay~.)+
  geom_point(aes(color=techangegene),size=I(.2))+
  theme_bw()
  )
dev.off()
normalizePath(plotfile)

#now plot
plotfile<- here(paste0source("", chdir = TRUE)('plots/','bdweight_vs_te','.pdf'))
pdf(plotfile)
print(
wteltcoregs%>%
enframe('gene_name','bd_weight')%>%
  inner_join(prediction_df%>%filter(assay=='TE')%>%
    mutate(gene_name=gid2gnm[[gene_id]]),by=c('gene_name'))%>%
  ggplot(data=.,aes(x=(bd_weight),y=estimate))+
  scale_x_continuous('Birthdate_weight')+
  scale_y_continuous('TE')+
  facet_wrap(time~.)+
  geom_point(aes(color=str_detect(gene_name,'^Rp[sl]')),size=I(.2))+
  theme_bw()
  )
dev.off()
normalizePath(plotfile)

