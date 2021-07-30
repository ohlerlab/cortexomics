base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
  base::source("src/Figures/Figure0/0_load_annotation.R")
}
rename<-dplyr::rename
first<-dplyr::first
last<-dplyr::last

  gid2gnm<-load_hashmap('gid2gnm.hmp')
  gnm2gid<-load_hashmap('gnm2gid.hmp')
  gnm2gid<-load_hashmap('gnm2gid.hmp')


sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
sel_ms_mat<-readRDS('data/sel_ms_mat.rds')
countpred_df<-readRDS('data/countpred_df.rds')
tx_countdata<-readRDS('data/tx_countdata.rds')
ms_metadf<-readRDS('data/ms_metadf.rds')
allvoom <- readRDS(here('data/allvoom.rds'))
prediction_df = bind_rows(
  sel_prodpreds%>%
    mutate(assay='MS')%>%
    select(gene_id,time,assay,estimate,CI.L,CI.R),
  countpred_df%>%
    separate(contrast,c('time','assay'))%>%
    select(gene_id,time,assay,estimate=logFC,CI.L,CI.R)
)
#
exprdf = bind_rows( 
  allvoom$E%>%
    as.data.frame%>%
    rownames_to_column('gene_id')%>%
    gather(dataset,signal,-gene_id)%>%
    separate(dataset,into=c('time','assay','replicate')),
  sel_ms_mat[,]%>%
    as.data.frame%>%
    rownames_to_column('gene_id')%>%
    filter(!is.na(gene_id))%>%
    gather(dataset,signal,-gene_id)%>%
    separate(dataset,into=c('time','assay','replicate'))
)
exprdf%>%head
exprdf$gene_name = gid2gnm[[exprdf$gene_id]]
prediction_df$gene_name = gid2gnm[[prediction_df$gene_id]]


#' Thi sis a title with inline R code `r foo`

#' First we load the list of protein IDs, handpicked by Matt, using only the small
#' or large subunits - no mitochondrial riboproteins  
#define ambigous protein groups as those which have elements that appear in more than one protein group
# allpgroups <- mstall$Protein_IDs%>%unique
# multids<-allpgroups%>%unique%>%str_split_fast(';')%>%unlist%>%table%>%keep(~ . > 1)%>%names
# all_ambig_pgroups<-allpgroups%>%sep_element_in(multids)
# library(data.table)

rpexprdf = exprdf%>%inner_join(ms_metadf%>%filter(is_rpl|is_rps))
rpexprdf$cat = case_when(
  rpexprdf$gene_id %in% (ms_metadf%>%filter(is_rpl)%>%.$gene_id) ~ 'RPL',
  rpexprdf$gene_id %in% (ms_metadf%>%filter(is_rps)%>%.$gene_id) ~ 'RPS',
  TRUE ~ 'other'
)
rpexprdf$cat%>%table
rpexprdf%<>%bind_rows(rpexprdf%>%filter(assay=='ribo')%>%mutate(assay='TE',signal = signal - (rpexprdf%>%filter(assay=='total')%>%.$signal)))

#possibly useful formodeling?
ggdf%>%filter(assay=='ribo')%>%group_by(time)%>%summarise_at(vars(signal),median)%>%.$signal%>%
  saveRDS('data/rp_ribo_offsets.rds')

ggdf = rpexprdf%>%
  arrange(time,assay=='MS',assay=='TE',assay=='ribo')%>%
  mutate(assay=as_factor(assay))%>%
  group_by(gene_name,assay,time,cat)%>%
  summarise(signal=mean(signal,na.rm=TRUE))%>%
    group_by(gene_name,assay,cat)%>%
    mutate(signal = signal - median(signal[time=='E13']))

medggdf = ggdf%>%group_by(assay,cat,time)%>%summarise(signal = median(signal))

#now plot
plotfile<- here('plots/figures/figure1/RP_traj.pdf')
pdf(plotfile,w=14,h=7)
ggdf%>%
    # filter(signal < -3)
  ggplot(.,aes(y=signal,x=as.numeric(as.factor(time)),group=gene_name))+
  geom_line(alpha=I(.8),color=I("grey"))+
  # scale_color_discrete(name='colorname',colorvals)+
  scale_x_continuous(paste0('Time'),labels = rpexprdf$time%>%unique)+
  scale_y_continuous(paste0('log2(CPM/iBAQ) relative to Mean'),limits=c(-1,1))+
  ggtitle(paste0('Ribosomal Proteins - Expression Trajectory'))+
  facet_grid(cat~assay,scales='free')+
  geom_line(data=medggdf,aes(group=cat),color=I('black'))+
  theme_bw()
dev.off()
normalizePath(plotfile)

