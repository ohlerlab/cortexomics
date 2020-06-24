
# # install devtools if necessary
# install.packages('devtools')

# # install the MuSiC package
# devtools::install_github('xuranw/MuSiC')

# load
library(xbioc)
library(MuSiC)
source('src/R/Rprofile.R')
metainfo<-'data/metainfo.tsv'%>%fread

if(!exists('ss_emat')) ss_emat <- projmemoise(fread)(Sys.glob(here('ext_data/GSE11*')))

if(!exists('countexprdata')) load('data/integrate_exprdata2.Rdata')

{
	featuredata = tibble(gene_name=ss_emat[[1]])%>%left_join(metainfo%>%filter(isbest)%>%distinct(gene_name,protein_id,transcript_id,uprotein_id))%>%
	  set_rownames(.$gene_name)
	sscoldata <- ss_emat[,-1]%>%colnames%>%str_split('[\\.\\_]')%>%simplify2array%>%t%>%
	      .[,1:2]%>%set_colnames(c('fttime','ftdiff'))%>%as.data.frame
	sscoldata%<>%mutate(cellType = paste0(fttime,'_',ftdiff))%>%as.data.frame
	rownames(sscoldata) = ss_emat[,-1]%>%colnames
	sscountexprdata <- ExpressionSet(
	  ss_emat[,-1]%>%as.matrix,
	  AnnotatedDataFrame(sscoldata%>%as.data.frame),
	  AnnotatedDataFrame(featuredata)
	)
	# pData(Mousesub.eset)$clusterType = factor()
	sscountexprdata%>%saveRDS(here('pipeline/exprdata/telleyexprset.rds'))
	}


countexprdata%>%colnames



clusters.type = pData(sscountexprdata)$cellType%>%unique%>%setNames(.,.)
cl.type = as.character(sscountexprdata$cellType)
for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}

#get our best genes for the expression data
bestcountexprdata = countexprdata[rownames(countexprdata)%in%(metainfo%>%filter(isbest)%>%.$protein_id)]
rownames(bestcountexprdata) = tibble(protein_id=rownames(bestcountexprdata))%>%left_join(metainfo%>%filter(isbest)%>%distinct(gene_name,protein_id,transcript_id,uprotein_id))%>%
	  .$gene_name

pData(sscountexprdata)$sampleID = rep(c(1,2),nrow(pData(sscountexprdata)))[1:nrow(pData(sscountexprdata))]

dtegenes = metainfo%>%filter(dTE)%>%.$gene_name

genessets <- list(
  allgenes=rownames(sscountexprdata),
  dte_genes = rownames(sscountexprdata)%>%intersect(dtegenes),
  non_dte_genes = rownames(sscountexprdata)%>%setdiff(dtegenes)
)

music_reslist = genessets%>%lapply(function(geneset){
  music_prop(bulk.eset = bestcountexprdata, sc.eset = sscountexprdata[geneset,],clusters = 'cellType',clusters.type=clusters.type,sample='sampleID')
})


#now plot all sets
music_reslist %>% seq_along%>%lapply(function(i){
  musicres <- music_reslist[[i]]
  tesetname = names(music_reslist)[[i]]


library(naniar)

gtdf1 <- musicres[[1]] %>% as.data.frame%>%rownames_to_column('dset')%>%gather(ftset,music_prop,-dset)%>%
  separate(ftset,c('ft_time','ft_diff'))%>%
  group_by(dset,ft_time,ft_diff)%>%
  separate(dset,c('time','assay','rep'))%>%
  spread(assay,music_prop)%>%
  mutate(est_gte = ribo/total)


telleyestgteplots<-lapply(list(quo(ft_diff),quo(ft_time)),function(telleyaxis){
   
gtdf <- gtdf1%>%
  group_by(!!telleyaxis,time,rep)%>%
  summarise(ribo=sum(ribo),total=sum(total))%>%
  mutate(est_gte = ribo/total)%>%
  mutate(log2_est_gte = log2(est_gte))%>%
  ungroup%>%
  mutate(ntime = as.numeric(as.factor(time)))
gtdf$log2_est_gte%<>%replace(.,!is.finite(.),NA)

mval = ldf$log2_est_gte%>%min
ldf <- gtdf%>%ungroup%>%group_by(ntime,!!telleyaxis)%>%
  mutate(log2_est_gte = replace(log2_est_gte,!is.finite(log2_est_gte),mval*1.1))  %>%
  summarise(log2_est_gte=mean((log2_est_gte)))

gtdf%>%  ggplot(aes(x=ntime,color=!!telleyaxis,y=log2_est_gte))+
  geom_point(size=I(4))+
  geom_line(data=ldf)+
  geom_miss_point(shape=I(4),size=I(3))+
  scale_x_continuous(name=paste0('Bulk Time Point'))+
  scale_y_continuous(name=paste0('Estimated Relative Global TE (Ribo/RNA)'))+
  theme_bw()

})
telleygtepanel = telleyestgteplots%>%ggarrange(nrow=2,plotlist=.)
plotfile = paste0('plots/','telley_MuSiC_est',tesetname,'.pdf')
pdf(width=10,h=8,file=plotfile)
print(telleygtepanel)
dev.off()
plotfile%>%normalizePath(mustWork=TRUE)%>%message
})  


#now plot
cspftdiff=cosdistres%>%
  #group_by(dset,ft_diff)%>%summarise(music_prop=sum(music_prop))%>%
  separate(dset,c('time','assay'))%>%
  # ungroup%>%mutate(music_prop=music_prop/min(music_prop))%>%
  # ggplot(.,aes(x=as.numeric(as.factor(time)),y=ft_diff,size=music_prop,color=assay))+
  ggplot(.,aes(x=as.numeric(as.factor(time)),y=music_prop,fill=ft_diff))+
  stat_identity(geom='bar',position='stack')+
  # geom_point()+
  facet_grid(assay~.)+
  # scale_color_discrete(name='',colorvals)+
  scale_x_continuous(paste0('Time Point (bulk)'))+
  scale_y_discrete(paste0('Telley time since Flashtag\n(96H = Neuron)'))+
  ggtitle(paste0('MuSiC results using Telley FT groups'))+
  theme_bw()
#now plot
cspftt=cosdistres%>%
  #group_by(dset,ft_time)%>%summarise(music_prop=sum(music_prop))%>%
  separate(dset,c('time','assay'))%>%
  # ggplot(.,aes(x=as.numeric(as.factor(time)),y=ft_time,size=music_prop,color=assay))+
  # geom_point()+
  # ungroup%>%mutate(music_prop=music_prop/min(music_prop))%>%
  ggplot(.,aes(x=as.numeric(as.factor(time)),y=music_prop,fill=ft_time))+
  stat_identity(geom='bar',position='stack')+
  # geom_point()
  facet_grid(assay~.)+
  # scale_color_discrete(name='',colorvals)+
  scale_x_continuous(paste0('Time Point (bulk)'))+
  scale_y_discrete(paste0('Telley Collection Stage'))+
  ggtitle(paste0('MuSiC results using Telley FT groups'))+
  theme_bw()
pdf(width=10,h=4,file=paste0('plots/','telley_MuSiC_bulk_both','.pdf'))
ggarrange(ncol=2,plotlist=list(cspftdiff,cspftt))
dev.off()
paste0('plots/','telley_MuSiC_bulk_both','.pdf')%>%normalizePath



