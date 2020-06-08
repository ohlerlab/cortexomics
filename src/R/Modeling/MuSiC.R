/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/R/Modeling/telleycellsim.Rmd

# # install devtools if necessary
# install.packages('devtools')

# # install the MuSiC package
# devtools::install_github('xuranw/MuSiC')

# load
library(xbioc)
library(MuSiC)
metainfo<-'data/metainfo.tsv'%>%fread

load('data/integrate_exprdata2.Rdata')

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


ss_emat <- projmemoise(fread)(Sys.glob(here('ext_data/GSE11*')))

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

music_res = music_prop(bulk.eset = bestcountexprdata, sc.eset = sscountexprdata,clusters = 'cellType',clusters.type=clusters.type,sample='sampleID')

cosdistres <- music_res[[1]] %>% as.data.frame%>%rownames_to_column('dset')%>%gather(ftset,music_prop,-dset)%>%
  separate(ftset,c('ft_time','ft_diff'))%>%
  mutate(dset=str_replace(dset,'_\\d$',''))%>%
  group_by(dset,ft_time,ft_diff)%>%
  summarise(music_prop=mean(music_prop))%>%
  identity
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



Est.mouse.bulk = music_prop.cluster(bulk.eset = Mouse.bulk.eset, sc.eset = Mousesub.eset, group.markers = IEmarkers, clusters = 'cellType', group = 'clusterType', samples = 'sampleID', clusters.type = clusters.type)

load('https://xuranw.github.io/MuSiC/data/IEmarkers.RData')

