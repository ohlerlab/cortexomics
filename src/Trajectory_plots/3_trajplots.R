{
base::source(here::here('src/Rprofile.R'))
if(!exists("cdsgrl")) {
  base::source("src/Preprocess/0_load_annotation.R")
}
pdf <- grDevices::pdf
ms_metadf<-readRDS('data/ms_metadf.rds')

assaynames = c('total'='RNA-seq','ribo'='Ribo-seq','MS'='Mass-Spec','TE'='TE')
stagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
tpnames = names(stagecols)%>%setNames(c('E13','E145','E16','E175','P0'))


rename<-dplyr::rename
first<-dplyr::first
last<-dplyr::last


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
exprdf$gene_name = gid2gnmv[exprdf$gene_id]
prediction_df$gene_name = gid2gnmv[prediction_df$gene_id]
}

teexprdf = left_join(
  exprdf%>%filter(assay=='ribo'),
  exprdf%>%filter(assay=='total'),
  suffix=c('_ribo','_total'),
  by=c("gene_id", "time", "replicate", "gene_name")
)%>%mutate(signal = signal_ribo - signal_total,assay='TE')%>%
  select(one_of(c("gene_id", "time", "replicate", "gene_name",'assay','signal')))
exprdf = bind_rows(exprdf,teexprdf)

gname='Satb2'

{
#now plot
# dir.create(plotfile%>%dirname,rec=TRUE)
#
getp <- function(sym){
  e=parent.frame(2)
  get(sym,env=e)
}
# myfun = function(x=getp('x')){
#   x+1
# }
# x=3
# myfun()
#see becker_modelling for this data
bmodelvals <- read_tsv(here('data/bmodelvals.tsv'))
assays2plot=c('total','MS','ribo')
myymin=-3
myymax=3
show_model=TRUE
modelcols <- c(
  msdev = '#c64730',
  production = '#fda27c',
  linear = '#ba5e47',
  degredation = '#1d4c9b',
  stationary = '#6289cc'
)
make_trajplot = function(gname,
  myymin=NULL,
  myymax=NULL,
  exprdf=getp('exprdf'),
  prediction_df=getp('prediction_df'),
  tpnames=getp('tpnames'),
  assaynames=getp('assaynames'),
  assays2plot=c('total','MS','ribo'),
bmodelvals=getp('bmodelvals'),
show_model=getp('show_model'),
modelcols=getp('modelcols')){
  #
  stopifnot(exists('exprdf'))
  stopifnot(exists('prediction_df'))
  if(!gname %in% exprdf$gene_name) return(NULL)
    #
ggpred <- prediction_df%>%
  filter(gene_name%in%gname)%>%
  mutate(passay=assaynames[assay])%>%
  filter(assay%in%assays2plot)%>%
  arrange(passay=='TE',passay=='MS',passay=='Ribo-seq')%>%
  ungroup%>%
  mutate(passay=as_factor(passay))%>%
  group_by(gene_name,assay)
t0sig = ggpred%>%ungroup%>%filter(time==time[1])%>%distinct%>%select(assay,gene_name,t0=estimate)
ggpred = ggpred%>%left_join(t0sig)%>%mutate_at(vars(estimate,CI.L,CI.R),~ .-t0)
ggexpr = exprdf%>%
  filter(assay%in%assays2plot)%>%
  filter(gene_name%in%gname)%>%
  mutate(passay=assaynames[assay])%>%
  arrange(passay=='TE',passay=='MS',passay=='Ribo-seq')%>%
  ungroup%>%
  mutate(passay=as_factor(passay))%>%
  group_by(gene_name,assay)%>%
  left_join(t0sig)%>%
  mutate(signal = signal - t0)
  ggexpr$passay%<>%factor(levels=assaynames)
#
drange = c(ggexpr$signal,ggpred$CI.R,ggpred$CI.L)%>%range
drange = drange%>%divide_by(0.5)%>%floor%>%{.[2]=.[2]+1;.}%>%multiply_by(0.5)
drange[1]%<>%min(.,-1)
drange[2]%<>%max(.,1)
if(!is.null(myymin))drange[1] = myymin
if(!is.null(myymax))drange[2] = myymax
#
gmodeldat <- bmodelvals%>%
  filter(gene%in%gname)%>%
  mutate(gene_name=gene)%>%
  mutate(passay=assaynames['MS'])
  gmodeldat$passay%<>%factor(levels=assaynames)
  gmodeldat$estimate = gmodeldat$estimate - (t0sig%>%filter(assay=='MS')%>%.$t0)
  best_model = gmodeldat$model[1]
  best_model = paste0('Kinetic Class - ',best_model)
#
  labldf=gmodeldat%>%group_by(gene_name)%>%
    dplyr::slice(1)%>%mutate(labl=paste0('Kinetic Class - ',model))%>%
    mutate(modelcol=modelcols[model])
outplot <- ggexpr%>%
  ggplot(.,aes(y=signal,x=as.numeric(as.factor(time)),group=gene_name))+
  geom_point()+
  geom_line(data=ggpred,aes(y=estimate,group=gene_name))+
  geom_ribbon(data=ggpred,aes(y=estimate,ymax=CI.R,ymin=CI.L,group=gene_name),fill='darkgrey',alpha=I(0.5))+
  # geom_line(alpha=I(.8),color=I("black"))+
  scale_x_continuous(paste0('Time'),labels = exprdf$time%>%unique%>%tpnames[.])+
  scale_y_continuous(paste0('log2(Counts/iBAQ)\n vs E12.5'))+
  coord_cartesian(ylim=drange)+
  ggtitle(str_interp('${gname} - Expression Trajectory'))+
  facet_wrap(~passay,ncol=4)+
  theme_bw()
  if(show_model){
    outplot = outplot + 
        geom_line(data=gmodeldat,aes(y=estimate,group=gene_name,color=model),linetype='dashed',size=I(2))+
        geom_text(data=labldf,aes(label=labl,x=-Inf,y=Inf,color=model),hjust=0,vjust=1)+
        scale_color_discrete(modelcols,guide=F)
    outplot
  }
  else{
    outplot
  }
}
}
#
{
datafile <- here('src/Shiny/data/make_trajplots_arglist.rds')
trajplot_arglist <- list(
  make_trajplot=make_trajplot ,
  exprdf=exprdf,
  prediction_df=prediction_df,
  tpnames=tpnames,
  assaynames=assaynames[c('total','ribo','MS','TE')],
  bmodelvals=bmodelvals,modelcols=modelcols,show_model=TRUE,
  assays2plot=c('total','ribo','MS','TE'))
saveRDS(trajplot_arglist,datafile)
saveRDS(trajplot_arglist,'data/make_trajplots_arglist.rds')
}
#
pdf<-grDevices::pdf
{
  'plots/Trajectory_plots/trajectory.pdf'%>%dirname%>%dir.create
plotfile<- here('plots/Trajectory_plots/trajectory.pdf')
pdf(plotfile,w=4*3,h=4*1)
print(ggarrange(plotlist=map(c('Nes'),make_trajplot,show_model=TRUE),nrow=1))
dev.off()
normalizePath(plotfile)
}



pdf<-grDevices::pdf
{
  'plots/Trajectory_plots/trajectory.pdf'%>%dirname%>%dir.create
plotfile<- here('plots/Trajectory_plots/trajectory.pdf')
pdf(plotfile,w=4*3,h=4*5)
print(ggarrange(plotlist=map(c('Satb2','Nes','Flna','Tle4','Bcl11b','Pum2'),make_trajplot),nrow=5))
dev.off()
normalizePath(plotfile)

plotfile<- here('plots/Trajectory_plots/trajectory_te.pdf')
pdf(plotfile,w=4,h=4*5)
print(ggarrange(plotlist=map(c('Satb2','Nes','Flna','Tle4','Bcl11b','Pum2'),make_trajplot,assays2plot='TE'),nrow=5))
dev.off()
normalizePath(plotfile)
}

{
plotfile<- here('plots/Trajectory_plots/pumsatoverlay.pdf')
pdf(plotfile,w=4,h=4*2)
print(ggarrange(plotlist=list(make_trajplot(c('Satb2','Pum2'),assays2plot='MS'))),nrow=2)
dev.off()
normalizePath(plotfile)
}


make_trajplots<-function(gnames,make_trajplot,...){
  ngenes <- length(gnames)
  ggarrange(plotlist=map(gnames,make_trajplot,...),nrow=ngenes)
}
make_trajplots <- do.call(what=partial,args=c(list(make_trajplots),trajplot_arglist))
{
igenes <- c('Satb2','Nes')
plotfile<- here('plots/trajapptest.pdf')
pdf(plotfile,w=4*4,h=4*length(igenes))
print(make_trajplots(igenes))
dev.off()
normalizePath(plotfile)
}

