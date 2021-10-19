################################################################################
################################################################################
# devtools::install_git('https://github.com/renozao/xbioc')
# devtools::install_git('https://xuranw.github.io/MuSiC/index.html')
base::source(here::here('src/Rprofile.R'))
if(!exists("cdsgrl")) {
  base::source("src/Figures/load_annotation.R")
}

# load
library(xbioc)
library(MuSiC)
library(naniar)
base::source('src/Rprofile.R')

allxtail = Sys.glob('pipeline/xtail/*')%>%map_df(.id='time',fread)%>%group_by(gene_name)
allxtail$gene_id = gnm2gid[[ allxtail$gene_name]]
techangedf <- allxtail%>%group_by(gene_id,gene_name)%>%
  mutate(sig = (adj_p_value < 0.05)& (abs(log2fc)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (log2fc > 0))),
    down = as.numeric(any(sig & (log2fc < 0)))
  )
techangegenes = techangedf%>%filter(up==1|down==1)%>%.$gene_name
dtegenes = techangegenes

if(!exists('ss_emat')) ss_emat <- projmemoise(fread)(Sys.glob(here('ext_data/GSE11*tsv.gz')))
ss_emat <- ss_emat[ss_emat[[1]]%in%allxtail$gene_name,]
if(!exists('tx_countdata'))tx_countdata <- readRDS(here('data/tx_countdata.rds'))

{

	featuredata = tibble(gene_name=ss_emat[[1]],gene_id=gnm2gid[[ss_emat[[1]]]])%>%set_rownames(.,.$gene_name)
	sscoldata <- ss_emat[,-1]%>%colnames%>%str_split('[\\.\\_]')%>%simplify2array%>%t%>%
	      .[,1:2]%>%set_colnames(c('fttime','ftdiff'))%>%as.data.frame
	sscoldata%<>%mutate(cellType = paste0(fttime,'_',ftdiff))%>%as.data.frame
	rownames(sscoldata) = ss_emat[,-1]%>%colnames
	sscountexprdata <- ExpressionSet(
	  ss_emat[,-1]%>%as.matrix,
	  AnnotatedDataFrame(sscoldata%>%as.data.frame),
	  AnnotatedDataFrame(featuredata)
	)
	sscountexprdata%>%saveRDS(here('data/telleyexprset.rds'))
	}

#
telleycoregenes = 'ext_data/telley_weights_comb.xlsx'%>%readxl::read_excel(.)%>%.[[1]]%>%tail(-1)%>%intersect(rownames(sscountexprdata))
wteltcoregs = 'ext_data/telley_weights_comb.xlsx'%>%readxl::read_excel(.)%>%tail(-1)%>%{setNames(as.numeric(.[[2]]),.[[1]])}%>%.[telleycoregenes]
#
telleycoregenesdiff = 'ext_data/telley_weights_comb.xlsx'%>%readxl::read_excel(.)%>%.[[3]]%>%tail(-1)%>%intersect(rownames(sscountexprdata))
wteltcoregsdiff = 'ext_data/telley_weights_comb.xlsx'%>%readxl::read_excel(.)%>%tail(-1)%>%{setNames(as.numeric(.[[4]]),.[[3]])}%>%.[telleycoregenesdiff]
telleycoregenesboth = 'ext_data/telley_weights_comb.xlsx'%>%readxl::read_excel(.)%>%{c(.[[1]],.[[3]])}%>%tail(-1)%>%intersect(rownames(sscountexprdata))


clusters.type = pData(sscountexprdata)$cellType%>%unique%>%setNames(.,.)
cl.type = as.character(sscountexprdata$cellType)
for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}

#get our best genes for the expression data
bestcountexprdata = tx_countdata$counts%>%set_rownames(.,gid2gnm[[rownames(.)]])
bestcountexprdata = ExpressionSet(bestcountexprdata,featureData=data.frame(gene_name=rownames(bestcountexprdata))%>%set_rownames(.$gene_name)%>%AnnotatedDataFrame)

pData(sscountexprdata)$sampleID = rep(c(1,2),nrow(pData(sscountexprdata)))[1:nrow(pData(sscountexprdata))]

#define our gene sets
genessets <- list(
  allgenes=rownames(sscountexprdata),
  dte_genes = rownames(sscountexprdata)%>%intersect(dtegenes),
  non_dte_genes = rownames(sscountexprdata)%>%setdiff(dtegenes),
  telleycoregenes = telleycoregenes,
  telleycoregenesdiff = telleycoregenesdiff,
  telleycoregenesboth = telleycoregenesboth
)

#Now run music on our data
exprs <- Biobase::exprs
music_reslist = genessets%>%lapply(function(geneset){
  music_prop(bulk.eset = bestcountexprdata, sscountexprdata[geneset,],clusters = 'cellType',clusters.type=clusters.type,sample='sampleID')
})
dir.create(here('plots/Figure5'),showWarn=F)
#now plot all sets
music_reslist %>% seq_along%>%lapply(function(i){
  musicres <- music_reslist[[i]]
  tesetname = names(music_reslist)[[i]]
  #
  gtdf1 <- musicres[[1]] %>% as.data.frame%>%rownames_to_column('dset')%>%gather(ftset,music_prop,-dset)%>%
    separate(ftset,c('ft_time','ft_diff'))%>%
    group_by(dset,ft_time,ft_diff)%>%
    separate(dset,c('time','assay','rep'))%>%
    spread(assay,music_prop)%>%
    mutate(est_gte = ribo/total)
  #
  telleyestgteplots<-lapply(list(quo(ft_diff),quo(ft_time)),function(telleyaxis){
  #   
  gtdf <- gtdf1%>%
    group_by(!!telleyaxis,time,rep)%>%
    summarise(ribo=sum(ribo),total=sum(total))%>%
    mutate(est_gte = ribo/total)%>%
    mutate(log2_est_gte = log2(est_gte))%>%
    ungroup%>%
    mutate(ntime = as.numeric(as.factor(time)))
  gtdf$log2_est_gte%<>%replace(.,!is.finite(.),NA)
  #
  mval = gtdf$log2_est_gte%>%min
  #
  ldf <- gtdf%>%ungroup%>%group_by(ntime,!!telleyaxis)%>%
    mutate(log2_est_gte = replace(log2_est_gte,!is.finite(log2_est_gte),mval*1.1))  %>%
    summarise(log2_est_gte=mean((log2_est_gte)))
  #
  #
  gtdf%>%  ggplot(aes(x=ntime,color=!!telleyaxis,y=log2_est_gte))+
    geom_point(size=I(4))+
    geom_line(data=ldf)+
    geom_miss_point(shape=I(4),size=I(3))+
    scale_x_continuous(name=paste0('Bulk Time Point'))+
    scale_y_continuous(name=paste0('Estimated Relative Global TE (Ribo/RNA)'))+
    theme_bw()
  })
  #
  telleygtepanel = telleyestgteplots%>%ggarrange(nrow=2,plotlist=.)
  plotfile = paste0('plots/Figure5/','telley_MuSiC_est',tesetname,'.pdf')
  pdf(width=10,h=8,file=plotfile)
  print(telleygtepanel)
  dev.off()
  plotfile%>%normalizePath(mustWork=TRUE)%>%message
}) 


for(geneset in names(music_reslist)){
  cosdistres <- music_reslist[[geneset]][[1]]%>%as.data.frame%>%rownames_to_column('dset')%>%gather(ftset,music_prop,-dset)
  cosdistres<-cosdistres%>%separate(dset,c('time','assay','rep'))%>%group_by(time,assay,ftset)%>%summarise(music_prop=mean(music_prop))
  cosdistres%>%
    separate(ftset,c('ft_time','ft_diff'))%>%
    group_by(time,assay,ft_diff)%>%summarise(music_prop=sum(music_prop))%>%
    group_by(time,assay,ft_diff)%>%group_slice(2) 
  #now plot
  cspftdiff=cosdistres%>%
    separate(ftset,c('ft_time','ft_diff'))%>%
    group_by(time,assay,ft_diff)%>%summarise(music_prop=sum(music_prop))%>%
    # ungroup%>%mutate(music_prop=music_prop/min(music_prop))%>%
    # ggplot(.,aes(x=as.numeric(as.factor(time)),y=ft_diff,size=music_prop,color=assay))+
    ggplot(.,aes(x=time,y=music_prop,fill=ft_diff))+
    stat_identity(geom='bar',position='stack')+
    # geom_point()+
    facet_grid(assay~.)+
    # scale_color_discrete(name='',colorvals)+
    scale_x_discrete(paste0('Time Point (bulk)'))+
    scale_y_continuous(paste0('MuSic Estimated Proportion'))+
    scale_fill_discrete(paste0('Telley time since Flashtag\n(96H = Neuron)'))+
    ggtitle(paste0('MuSiC results using Telley FT groups'))+
    theme_bw()
  cspftt=cosdistres%>%
    separate(ftset,c('ft_time','ft_diff'))%>%
    group_by(time,assay,ft_time)%>%summarise(music_prop=sum(music_prop))%>%
    # ungroup%>%mutate(music_prop=music_prop/min(music_prop))%>%
    # ggplot(.,aes(x=as.numeric(as.factor(time)),y=ft_diff,size=music_prop,color=assay))+
    ggplot(.,aes(x=time,y=music_prop,fill=ft_time))+
    stat_identity(geom='bar',position='stack')+
    # geom_point()+
    facet_grid(assay~.)+
    # scale_color_discrete(name='',colorvals)+
    scale_x_discrete(paste0('Time Point (bulk)'))+
    scale_y_continuous(paste0('MuSic Estimated Proportion'))+
    scale_fill_discrete(paste0('Telley Time of Collection'))+
    ggtitle(paste0('MuSiC results using Telley FT groups'))+
    theme_bw()
  pdf(width=10,h=4,file=paste0('plots/','telley_MuSiC_bulk_both_',geneset,'.pdf'))
  # ggarrange(ncol=2,plotlist=list(cspftdiff))
  print(ggarrange(ncol=2,plotlist=list(cspftdiff,cspftt)))
  dev.off()
  paste0('plots/Figure5/','telley_MuSiC_bulk_both_',geneset,'.pdf')%>%normalizePath%>%message
}

cosdistres%>%write_tsv('tables/musicdata.tsv')


pred_expr_df <- music_reslist[[1]][[1]] %>% 
  as.data.frame%>%
  rownames_to_column('dset')%>%
  gather(ftset,music_prop,-dset)%>%
  left_join(sc_abund_long,by='ftset')%>%
  group_by(dset,gene_name)%>%
  summarise(pred_abundance = sum(sc_abundance * music_prop))%>%
  left_join(bulkabund_long%>%select(dset,gene_name,bulk_abundance))

#now plot
dir.create('plots/Figures/FigureS6',showWarn=F,rec=TRUE)
plotfile<- here(paste0('plots/Figures/FigureS6','musicplots','.pdf'))
pdf(plotfile)
# lm(data=pred_expr_df%>%filter(is.finite(log2(sc_abundance)),is.finite(log2(bulk_abundance))),log2(bulk_abundance)~log2(sc_abundance))%>%plot
pred_expr_df%>%filter(is.finite(log2(pred_abundance)),is.finite(log2(bulk_abundance)))%>%
  arrange(desc(techangegene))%>%
  filter(bulk_abundance>32)%>%
  filter(dset%>%str_detect('E13|E175'),dset%>%str_detect('_1'))%>%
  ggplot(data=.,aes(x=log2(pred_abundance),y=log(bulk_abundance)))+
  scale_x_continuous('predicted abundance (MuSic)')+
  scale_y_continuous('Actual Bulk abundance')+
  facet_wrap(dset~.)+
  # geom_smooth(method='lm')+
  geom_point(aes(color=techangegene),size=I(.2))+
  # geom_point(aes(color=str_detect(gene_name,'^Rp[sl]')),size=I(.2))+
  theme_bw()
dev.off()
normalizePath(plotfile)



################################################################################
########dTE vs 
################################################################################
xtailte = readxl::read_xlsx('tables/S2.xlsx',sheet='xtail')
limmafc = readxl::read_xlsx('tables/S2.xlsx','limma_count_lfc',col_types='text')
absTEs = limmafc%>%filter(assay=='TE')%>%select(time,gene_id,logFC)%>%spread(time,logFC)
colnames(absTEs)%<>%str_replace('<NA>','E12.5')
for(tp2 in colnames(absTEs)%>%setdiff(c('E12.5','gene_id'))) absTEs[[tp2]] = as.numeric(absTEs[[tp2]])+as.numeric(absTEs[['E12.5']] )
absTEs%<>%gather(time,TE,-gene_id)
absTEs%<>%inner_join(xtailte%>%distinct(gene_name,gene_id))
absTEs$TE%<>%as.numeric

#now plot
plotfile<- here(paste0('plots/','telleyweight_vs_TE','.pdf'))
pdf(plotfile,w=16,h=4)
ggdf = wteltcoregsdiff%>%enframe('gene_name','telley_weight_diff')%>%
  inner_join(absTEs%>%
  select(time,gene_name,TE))
sigtests = ggdf%>%group_by(time)%>%
  split(.,.$time)%>%
  # lapply(.%>%{split(.$TE,.$telley_weight_diff>0)}%>%
  # {wilcox.test(.[[1]],.[[2]])})%>%
  lapply(.%>%{cor.test(.$TE,.$telley_weight_diff)})%>%
  map_df(.id='time',tidy)
sigtests$label = sigtests$p.value%>%round(3)%>%as.character%>%{paste0('Wilcox Test\np=',.)}
ggdf %>% ggplot(.,aes(x=telley_weight_diff>0,y=TE))+
  geom_boxplot()+
  facet_wrap(~time,ncol=4)+
  scale_x_discrete(paste0('Neuron Associated'))+
  scale_y_continuous(paste0('TE'))+
  coord_cartesian(ylim=c(-2,2))+
  ggtitle(paste0('Change in TE vs Telley Weights'))+
  geom_text(show.legend=F,data=sigtests,hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=label))+
  theme_bw()
dev.off()
message(normalizePath(plotfile))

#
#
#now plot
plotfile<- here(paste0('plots/','telleyweight_Time_vs_TE','.pdf'))
pdf(plotfile,w=16,h=4)
ggdf = wteltcoregs%>%enframe('gene_name','telley_weight_time')%>%
  inner_join(absTEs%>%
  select(time,gene_name,TE))
sigtests = ggdf%>%group_by(time)%>%
  split(.,.$time)%>%
  lapply(.%>%{split(.$TE,.$telley_weight_time>0)}%>%
  {wilcox.test(.[[1]],.[[2]])})%>%
  map_df(.id='time',tidy)
sigtests$label = sigtests$p.value%>%round(3)%>%as.character%>%{paste0('Wilcox Test\nP=',.)}
ggdf%>%
  ggplot(.,aes(x=telley_weight_time>0,y=TE))+
  # geom_point()+
  geom_boxplot()+
  # geom_smooth(method='lm')+
  facet_wrap(~time,ncol=4)+
  scale_x_discrete(paste0('Late Stage Associated'))+
  scale_y_continuous(paste0('TE'))+
  coord_cartesian(ylim=c(-2,2))+
  ggtitle(paste0('Change in TE vs Telley Weights'))+
  geom_text(show.legend=F,data=sigtests,hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=label))+
  theme_bw()
dev.off()
message(normalizePath(plotfile))







#now plot
plotfile<- here(paste0('plots/','telleyweight_vs_dTE','.pdf'))
pdf(plotfile,w=16,h=4)
ggdf = wteltcoregsdiff%>%enframe('gene_name','telley_weight_diff')%>%
  inner_join(limmafc%>%
  select(time,gene_name,log2fc))
sigtests = ggdf%>%group_by(time)%>%
  split(.,.$time)%>%
  lapply(.%>%{split(.$log2fc,.$telley_weight_diff>0)}%>%
  {wilcox.test(.[[1]],.[[2]])})%>%
  map_df(.id='time',tidy)
sigtests$label = sigtests$p.value%>%round(3)%>%as.character%>%{paste0('Wilcox Test\np=',.)}
ggdf %>% ggplot(.,aes(x=telley_weight_diff>0,y=log2fc))+
  geom_boxplot()+
  facet_wrap(~time,ncol=4)+
  scale_x_discrete(paste0('Neuron Associated'))+
  scale_y_continuous(paste0('log2(delta TE)'))+
  coord_cartesian(ylim=c(-2,2))+
  ggtitle(paste0('Change in TE vs Telley Weights'))+
  geom_text(show.legend=F,data=sigtests,hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=label))+
  theme_bw()
dev.off()
message(normalizePath(plotfile))
#
#
#now plot
plotfile<- here(paste0('plots/','telleyweight_Time_vs_dTE','.pdf'))
pdf(plotfile,w=16,h=4)
ggdf = wteltcoregs%>%enframe('gene_name','telley_weight_time')%>%
  inner_join(limmafc%>%
  select(time,gene_name,log2fc))
sigtests = ggdf%>%group_by(time)%>%
  split(.,.$time)%>%
  lapply(.%>%{split(.$log2fc,.$telley_weight_time>0)}%>%
  {wilcox.test(.[[1]],.[[2]])})%>%
  map_df(.id='time',tidy)
sigtests$label = sigtests$p.value%>%round(3)%>%as.character%>%{paste0('Wilcox Test\nP=',.)}
ggdf%>%
  ggplot(.,aes(x=telley_weight_time>0,y=log2fc))+
  # geom_point()+
  geom_boxplot()+
  # geom_smooth(method='lm')+
  facet_wrap(~time,ncol=4)+
  scale_x_discrete(paste0('Late Stage Associated'))+
  scale_y_continuous(paste0('log2(delta TE)'))+
  coord_cartesian(ylim=c(-2,2))+
  ggtitle(paste0('Change in TE vs Telley Weights'))+
  geom_text(show.legend=F,data=sigtests,hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=label))+
  theme_bw()
dev.off()
message(normalizePath(plotfile))









