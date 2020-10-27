# load
library(xbioc)
# library(MuSiC)
base::source('src/R/Rprofile.R')

if(!exists('ss_emat')) ss_emat <- projmemoise(fread)(Sys.glob(here('ext_data/GSE11*.gz')))

rename = dplyr::rename
sscoldata <- ss_emat[,-1]%>%colnames%>%str_split('[\\.\\_]')%>%simplify2array%>%t%>%
	      .[,1:2]%>%
	      set_colnames(c('fttime','ftdiff'))%>%as.data.frame%>%
			mutate(cellType = paste0(fttime,'_',ftdiff))%>%
			set_rownames(ss_emat[,-1]%>%colnames)%>%
			mutate(cell = ss_emat[,-1]%>%colnames)

ube3names = ss_emat%>%.[[1]]%>%str_subset('Ube3')

libsizes = ss_emat%>%gather(cell,count,-V1)%>%
	left_join(sscoldata)%>%
	rename('gene_name':=V1)%>%
	group_by(fttime,ftdiff)%>%
	summarise(libsize=sum(count))

ube3dat = ss_emat%>%.[.[[1]]%in%ube3names]%>%gather(cell,count,-V1)%>%
	left_join(sscoldata)%>%
	rename('gene_name':=V1)%>%
	group_by(gene_name,fttime,ftdiff)%>%
	summarise(count=sum(count))

ube3dat%<>%left_join(libsizes)%>%mutate(ncount = count/libsize)

#now plot
plotfile<- here(paste0('plots/','ube3counts','.tiff'))
grDevices::tiff(plotfile)
ube3dat%>%
	group_by(gene_name)%>%mutate(relative_normalized_read_count = ncount/ncount[(fttime=='E12')&(ftdiff=='1H')])%>%
	ggplot(.,aes(x=gene_name,y=relative_normalized_read_count,fill=gene_name))+
	stat_identity(geom='bar')+
	scale_fill_discrete(name='Gene')+
	scale_x_discrete(paste0('Gene'))+
	scale_y_continuous(paste0('Reads Per Million (relative to E12/1H)'))+
	facet_grid(ftdiff~fttime)+
	ggtitle(paste0('Ube3 Family - read counts in Telley Data'))+
	theme_bw()+
	theme(axis.text.x=element_text(vjust=0.5,angle=45))
dev.off()
normalizePath(plotfile)