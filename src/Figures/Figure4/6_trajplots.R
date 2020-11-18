{
base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
  base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
}
pdf <- grDevices::pdf
ms_metadf<-readRDS('data/ms_metadf.rds')

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
allvoom <- readRDS(here('data/allvoom.rds'))
prediction_df = bind_rows(
  sel_prodpreds%>%
    mutate(assay='MS')%>%
    select(gene_id,time,assay,estimate,CI.L,CI.R),
  countpred_df%>%
    separate(contrast,c('time','assay'))%>%
    select(gene_id,time,assay,estimate=logFC,CI.L,CI.R)
)

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
}



gname='Satb2'
#now plot
# dir.create(plotfile%>%dirname,rec=TRUE)
#
assaynames = c('total'='RNA-seq','ribo'='Ribo-seq','TE'='TE','MS'='Mass-Spec')
stagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
tpnames = names(stagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
make_trajplot = function(gname){
  stopifnot(exists('exprdf'))
  stopifnot(exists('prediction_df'))
  prediction_df%<>%mutate(passay=assaynames[assay])
  exprdf%<>%mutate(passay=assaynames[assay])
ggpred <- prediction_df%>%filter(gene_name==gname)%>%filter(!assay=='TE')%>%
  arrange(passay=='MS',passay=='Ribo-seq')%>%
  ungroup%>%
  mutate(passay=as_factor(passay))%>%
  group_by(gene_name,assay)
t0sig = ggpred%>%ungroup%>%filter(time==time[1])%>%distinct%>%select(assay,gene_name,t0=estimate)
ggpred = ggpred%>%left_join(t0sig)%>%mutate_at(vars(estimate,CI.L,CI.R),~ .-t0)
ggexpr = exprdf%>%
  filter(gene_name==gname)%>%
  arrange(passay=='Mass-Spec',passay=='Ribo-seq')%>%
  ungroup%>%
  mutate(passay=as_factor(passay))%>%
  group_by(gene_name,assay)%>%
  left_join(t0sig)%>%
  mutate(signal = signal - t0)
#
ggexpr%>%
    # filter(signal < -3)
  ggplot(.,aes(y=signal,x=as.numeric(as.factor(time)),group=gene_name))+
  geom_point()+
  geom_line(data=ggpred,aes(y=estimate))+
  geom_ribbon(data=ggpred,aes(y=estimate,ymax=CI.R,ymin=CI.L),,fill='darkgreen',alpha=I(0.5))+
  # geom_line(alpha=I(.8),color=I("black"))+
  # scale_color_discrete(name='colorname',colorvals)+
  scale_x_continuous(paste0('Time'),labels = exprdf$time%>%unique%>%tpnames[.])+
  scale_y_continuous(paste0('log2(CPM/iBAQ) relative to T0'))+
  ggtitle(str_interp('${gname} - Expression Trajectory'))+
  facet_wrap(~passay,ncol=3)+
  # geom_line(data=medggdf,aes(group=cat),color=I('black'))+
  theme_bw()
}
pdf(plotfile,w=3*4,h=4*4)
plotfile<- here('plots/figures/figure4/trajectory.pdf')
ggarrange(plotlist=map(c('Satb2','Nes','Flna','Tle4'),make_trajplot),nrow=4)
dev.off()
normalizePath(plotfile)
