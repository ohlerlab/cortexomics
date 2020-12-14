################################################################################
########
################################################################################


################################################################################
################################################################################
base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
}

	
igene = 'Satb2'
igenetrs = ids_nrgname %>% 
	safe_filter(gene_name=='Satb2')%>%
	.$transcript_id%>%unique

dsets = iso_tx_countdata$abundance%>%colnames
colgrps = dsets%>%str_extract('.*?(?=_\\d+)')
colgrp = colgrps[[1]]

ucolgrps = unique(colgrps)%>%setNames(.,.)
relativeusedf = map_df(.id='dset',ucolgrps,function(colgrp){
	sab = iso_tx_countdata$abundance[,colgrps==colgrp]%>%rowSums
	sab%>%
		enframe('transcript_id','abundance')%>%
		left_join(ids_nrgname%>%select(transcript_id,gene_name))%>%
		group_by(gene_name)%>%mutate(abundance=abundance/sum(abundance))
})


ign_usedf = relativeusedf%>%
	safe_filter(gene_name==igene)%>%
	separate(dset,c('time','assay'))


#now plot
plotfile<- here(paste0('plots/','tmp','.pdf'))
pdf(plotfile)
ign_usedf%>%
	ggplot(data=.,aes(x=time,fill=transcript_id,y=abundance))+
	geom_bar(stat='identity')+
	# scale_fill_discrete(name='colorname',colorvals)+
	facet_grid(assay~.)+
	scale_x_discrete(paste0('xname'))+
	scale_y_continuous(paste0('yname'))+
	ggtitle(paste0('title'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))


iso_tx_countdata$abundance