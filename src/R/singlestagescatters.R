pdf = grDevices::pdf
# src/R/Load_data/integrate_exprdata2.R
# src/R/Modeling/proDA.R



#get a set of genes with MS and high expression
#For now just use the best ones from integrate

#get the tracks for these
#Within our time point, compare RNAseq, Riboseq, Ribo/RNAseq (TE), 'entropy (maybe weighted)'
#Geometric mean of the CDS, with without clipping, and maybe some other shit too.


#let's simply try with the existing data
# load('data/integrate_exprdata2.Rdata')
base::source('src/R/Rprofile.R')
library(GenomicAlignments)
library(GenomicFeatures)
msid2upid = metainfo%>%
	distinct(ms_id,uprotein_id)%>%
	{safe_hashmap(.[['ms_id']],.[['uprotein_id']])}
msid2pid = metainfo%>%
	distinct(ms_id,uprotein_id)%>%
	{safe_hashmap(.[['ms_id']],.[['uprotein_id']])}
gname2pid = metainfo%>%
	distinct(gene_name,protein_id)%>%
	{safe_hashmap(.[['gene_name']],.[['protein_id']])}
upid2msid = metainfo%>%
	distinct(uprotein_id,ms_id)%>%
	{safe_hashmap(.[['uprotein_id']],.[['ms_id']])}
upid2pid = metainfo%>%
	distinct(uprotein_id,protein_id)%>%
	{safe_hashmap(.[['uprotein_id']],.[['protein_id']])}
trid2pid = metainfo%>%
	distinct(transcript_id,protein_id)%>%
	{safe_hashmap(.[['transcript_id']],.[['protein_id']])}
pid2gname = metainfo%>%
	distinct(protein_id,gene_name)%>%
	{safe_hashmap(.[['protein_id']],.[['gene_name']])}
pid2msid = metainfo%>%
	distinct(protein_id,ms_id)%>%
	{safe_hashmap(.[['protein_id']],.[['ms_id']])}
# distinct(data.frame(protein_id=gtf_gr$protein_id,gene_id=gtf_gr$gene_id))%>%filter(!is.na(protein_id))%>%
		# {safe_hashmap(.[['protein_id']],.[['gene_name']])}
# distinct(data.frame(protein_id=gtf_gr$protein_id,gene_id=gtf_gr$gene_id))%>%filter(!is.na(protein_id))

bestmsmat=matchedms_mat[upid2msid[[best_uprotein_ids]],]
#
bestprotids<-bestmsmat%>%rownames%>%str_replace('_\\d+','')
#
#make combined matrix
best_protids <- best_uprotein_ids%>%str_replace('_\\d+$','')
best_msids = upid2msid[[best_uprotein_ids]]
#
scatterdatamat <- cbind(allcountmat[best_protids,]%>%add(.1)%>%log2,bestmsmat)%>%set_rownames(rownames(bestmsmat))
scatterdatamat%<>%set_rownames(best_uprotein_ids)
#
for(tp in tps){
	for(rep in 1:2){
		scatterdatamat%<>%cbind(scatterdatamat[,paste0(tp,'_total_',rep)] - scatterdatamat[,paste0(tp,'_ribo_',rep)])
		colnames(scatterdatamat)[ncol(scatterdatamat)] = paste0(tp,'_TE_',rep)
	}
}


	
#
countpredmat = prediction_df_countonly%>%mutate(effect=paste0(time,'_',assay))%>%
  # filter(time!='E13')%>%
  filter(uprotein_id%in%best_protids)%>%
  select(uprotein_id,effect,logFC)%>%
  spread(effect,logFC)%>%
  {set_rownames(as.matrix(.[,-1]),.[[1]])}
stopifnot(exists('proda_ms_pred'))
#
scatterdatamat%<>%
	cbind(countpredmat[best_protids,],proda_ms_pred[best_msids,])
##
for(tp in tps){
	scatterdatamat%<>%cbind(scatterdatamat[,paste0(tp,'_total')] - scatterdatamat[,paste0(tp,'_ribo')])
	colnames(scatterdatamat)[ncol(scatterdatamat)] = paste0(tp,'_TE')
}
#
scatterdatamat%<>%{proDA::median_normalization(.)}


#add wave and spec individual stages
lspecmat = bamstatdf%>%select(sample,lspec,protein_id)%>%mutate(sample=str_replace(sample,'ribo','spec'))%>%spread(sample,lspec)%>%as.data.frame%>%{set_rownames(as.matrix(.[,-1]),.[,1])}
scatterdatamat%<>%cbind(lspecmat[upid2pid[[rownames(.)]],])
#add wave and spec averages
lwavemat = bamstatdf%>%select(sample,lwave,protein_id)%>%mutate(sample=str_replace(sample,'ribo','wave'))%>%spread(sample,lwave)%>%as.data.frame%>%{set_rownames(as.matrix(.[,-1]),.[,1])}
scatterdatamat%<>%cbind(lwavemat[upid2pid[[rownames(.)]],])

colnames(lspecmat)

lwavemat = lwavedf%>%mutate(stage=paste0(stage,'_','wave'))%>%spread(stage,value)%>%as.data.frame%>%{set_rownames(as.matrix(.[,-1]),.[,1])}
lspecmat = lspecdf%>%mutate(stage=paste0(stage,'_','spec'))%>%spread(stage,value)%>%as.data.frame%>%{set_rownames(as.matrix(.[,-1]),.[,1])}

scatterdatamat%<>%cbind(lwavemat[upid2pid[[rownames(.)]],],lspecmat[upid2pid[[rownames(.)]],])

#now plot
compcolpairs = fread('
one	two
E13_total_1 E13_MS_1
E13_ribo_1 E13_MS_1
E13_TE_1 E13_MS_1
P0_total_1 P0_MS_1
P0_ribo_1 P0_MS_1
P0_TE_1 P0_MS_1
E13_total E13_MS
E13_ribo E13_MS
P0_total P0_MS
P0_ribo P0_MS
E13_TE E13_MS
P0_TE P0_MS
E13_wave E13_MS
P0_wave P0_MS
E13_spec E13_MS
P0_spec P0_MS
'
)

#now plot
compcolpairs = fread('
one	two
E13_total_1 E13_MS_1
E13_ribo_1 E13_MS_1
E13_TE_1 E13_MS_1
E13_wave_1 E13_MS_1
E13_spec_1 E13_MS_1
E13_total E13_MS
E13_ribo E13_MS
E13_TE E13_MS
E13_wave E13_MS
E13_spec E13_MS
'
)
setdiff(compcolpairs[[1]],colnames(scatterdatamat))
setdiff(compcolpairs[[2]],colnames(scatterdatamat))
stopifnot(all(compcolpairs[[1]]%in%colnames(scatterdatamat)))
stopifnot(all(compcolpairs[[2]]%in%colnames(scatterdatamat)))
compcolpairs=compcolpairs%>%t%>%as.data.frame%>%as.list


library(lsd)
base::source('https://raw.githubusercontent.com/stineb/lsd/master/r/lsd.heatscatter.r')

plotlist = list()
for(compcolpair in compcolpairs){
	compcol = compcolpair[1]
	compcol2 = compcolpair[2]
#
ggdf <- scatterdatamat[,c(compcol,compcol2)]%>%as_tibble
#
plotfile<- here(paste0('plots/singlestagescatter/',compcol,compcol2,'.pdf'))
scattertitle = cor.test(ggdf[[1]],ggdf[[2]])%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n ${.$conf.low} - ${.$conf.high}')}
ecompcol = rlang::parse_expr(compcol)
ecompcol2 = rlang::parse_expr(compcol2)
ggdf = ggdf%>%{.[(.[2]>10),]}
pdf('tmp.pdf')
scatterplot=
	# ggplot(.,aes(!!ecompcol,!!ecompcol2))+
	# geom_point()+
	heatscatter(ggdf[[1]],ggdf[[2]],ggplot=TRUE)+
	# geom_smooth(formula= y ~ 1+x)+
	# scale_color_discrete(name='colorname',colorvals)+
	scale_x_continuous(paste0(compcol))+
	scale_y_continuous(paste0(compcol2))+
	ggtitle(paste0(scattertitle))+
	theme_bw()
dev.off()
plotlist%<>%append(list(scatterplot)%>%setNames(paste0(compcol,'_',compcol2)))
grDevices::pdf(plotfile)
print(scatterplot)
dev.off()
normalizePath(plotfile)%>%message
#
}

pdf('plots/singlestagescatter/grid.pdf'%>%normalizePath%T>%message)
print(ggarrange(plotlist=plotlist))
dev.off()
stop()
scatterdatamat%>%colnames
scatterdatamat[,'E145_MS']


load.image('riboreprocess_comps.RData')


#Now some other plots
#spec vs wave
#spec vs density
#ribowave vs density

#f




