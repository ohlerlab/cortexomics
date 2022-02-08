################################################################################
################################################################################
base::source(here::here('src/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("src/Figures/load_annotation.R")
}
gnm2gid = ids_nrgname%>%distinct(gene_id,gene_name)%>%
	{safe_hashmap(.[[1]],.[[2]])}
gid2gnm = ids_nrgname%>%distinct(gene_id,gene_name)%>%
	{safe_hashmap(.[[2]],.[[1]])}
#We are going off of- https://peerj.com/articles/270/#fig-4

################################################################################
########Within stage variance decomp
################################################################################


allxtail = Sys.glob('pipeline/xtail/*')%>%map_df(.id='time',fread)%>%group_by(gene_name)
allxtail$gene_id = gnm2gid[[ allxtail$gene_name]]
techangedf <- allxtail%>%group_by(gene_id,gene_name)%>%
  mutate(sig = (adj_p_value < 0.05)& (abs(log2fc)>log2(2)))%>%
  summarise(
    up = as.numeric(any(sig & (log2fc > 0))),
    down = as.numeric(any(sig & (log2fc < 0)))
  )
teupgenes = techangedf%>%filter(up==1)%>%.$gene_name
tedowngenes = techangedf%>%filter(down==1)%>%.$gene_name

sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
proDAfitms<-readRDS('data/proDAfitms.rds')
sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
sel_ms_mat<-readRDS('data/sel_ms_mat.rds')

countpred_df<-readRDS('data/countpred_df.rds')
countcontr_df<-readRDS('data/countcontr_df.rds')
tx_countdata<-readRDS('data/tx_countdata.rds')

countpred_df%<>%separate(contrast,c('time','assay'))
countpred_df%<>%mutate(se = (CI.R-CI.L)/3.92)

#okay so pretty much linear relationship between Mass spec and riboseq
#but NOT between mass spec and RNAseq!

tx_countdata$abundance%>%
	as.data.frame%>%rownames_to_column('gene_id')%>%
	gather(dataset,TPM,-gene_id)%>%
	separate(dataset,c('time','assay','rep'))%>%
	group_by(time,assay)%>%
	filter(time=='E13',assay=='ribo')%>%
	left_join(sel_prodpreds,by=c('gene_id','time'))%>%
	select(gene_id,diff,TPM)%>%
	filter(!is.na(diff),!is.na(TPM))%>%
	filter(TPM > quantile(TPM,.1))%>%
	# {txtplot(log2(.$TPM),.$diff)}
	lm(data=.,diff~log2(TPM))%>%confint

# countpred_df%>%
# 	separate(contrast,c('time','assay'))%>%
# 	filter(time=='E13',assay=='ribo')%>%
# 	left_join(sel_prodpreds,by=c('gene_id','time'))%>%
# 	select(gene_id,diff,logFC)%>%
# 	filter(!is.na(diff),!is.na(logFC))%>%
# 	filter(logFC > quantile(logFC,.1))%>%
# 	# {txtplot(log2(.$logFC),.$diff)}
# 	lm(data=.,diff~logFC)%>%confint

stopifnot('se' %in% colnames(countpred_df))

tps = unique(countpred_df$time) 

withinstage_varexpldf = map_df(.id='time',tps%>%setNames(.,.),function(stage){
map_df(.id='assay',c('total','ribo')%>%setNames(.,.),function(cassay){
	varmodeldf  = countpred_df%>%
		# separate(contrast,c('time','assay'))%>%
		filter(time==stage,assay==cassay)%>%
		left_join(sel_prodpreds,by=c('gene_id','time'))%>%
		filter(!is.na(diff),!is.na(logFC))%>%
		# filter(logFC > quantile(logFC,.1))%>%
		identity

	varmodeldf%>%select(gene_id,diff,logFC)%>%
		filter(logFC > quantile(logFC,.1))
		# {txtplot(log2(.$logFC),.$diff)}

	varmodel = varmodeldf%>%	lm(data=.,diff~logFC)
	b_all = varmodel$coef[2]
	b_r = 1
	allvar = sum(varmodel$residuals^2)
	var_p = sum(varmodeldf$se.y^2)
	var_r = sum(varmodeldf$se.x^2)
	var_pdt = allvar - (b_all/b_r)^2 * var_r - var_p
	var_mp = sum(varmodeldf%>%lm(data=.,diff~1)%>%.$residuals%>%.^2)
	var_explained = (var_mp - var_p - var_pdt) / (var_mp - var_p)
	names(var_explained)=NULL
	c(var_explained = var_explained,count_variance = var_r / (var_mp - var_p))
})})


#or use ms
contrdf<-readRDS('data/contrdf.rds')
contrdf$gene_name <- gid2gnm[[contrdf$gene_id]]
mschangedf = contrdf%>%
  group_by(gene_name)%>%
  mutate(sig = (adj_pval < 0.05)& (abs(diff)>log2(2)))%>%
  summarise(
    up = as.numeric(any(sig & (diff > 0))),
    down = as.numeric(any(sig & (diff < 0)))
 ) 
msupgenes = mschangedf%>%filter(up==1)%>%.$gene_name
msdowngenes = mschangedf%>%filter(down==1)%>%.$gene_name

genesets = list(all = sel_prodpreds$gene_name,
	TE_up = teupgenes,
	TE_down=tedowngenes,
	MS_up = msupgenes,
	MS_down = msdowngenes
)


# kineticres = 
cassay='ribo_kinetic' 
cassay='ribo' 
geneset=genesets[[1]]

betweenstagevardf =
map_df(.id='geneset',genesets,function(geneset){
# map_df(.id='assay',c('total','ribo','ribo_kinetic')%>%setNames(.,.),function(cassay){
map_df(.id='assay',c('total','ribo')%>%setNames(.,.),function(cassay){
		counttraj  = countpred_df%>%
			# separate(contrast,c('time','assay'))%>%
			filter(assay==cassay%>%str_replace('_kinetic',''))
		counttraj$gene_name = gid2gnm[[counttraj$gene_id]]
		countses = counttraj%>%select(gene_id,time,se)
	if(!cassay=='ribo_kinetic'){
	}else{

			kinetictraj = bmodelopts[gns_by_mod$production]%>%map(~.[['riboseq']][['production']]$par$prot)%>%
				discard(is.null)%>%
				map_df(.id='gene_name',.%>%set_colnames(tps)%>%as.data.frame)%>%
				pivot_longer(-gene_name,names_to='time',values_to='logFC')
			kinetictraj$logFC%<>%log
			alldrawses = bmodelopts%>%map('riboseq')%>%map('production')%>%map('cov')%>%discard(is.null)%>%map(diag)%>%map(~.[str_subset(names(.),'prot\\[')])%>%
					map_df(.id='gene_name',enframe,'time','se')
			alldrawses$time %<>% as.factor%>%as.numeric%>%tps[.]
			kinetictraj$gene_id = gnm2gid[[kinetictraj$gene_name]]
			# kinetictraj%<>%left_join(countses)	
			kinetictraj%<>%left_join(alldrawses)	
			kinetictraj$se = kinetictraj$se*(kinetictraj$logFC^(-2))
		counttraj <- kinetictraj%>%group_by(gene_id)%>%filter(!any(is.na(se)))	

	}

	counttraj$gene_name<-NULL
	varmodeldf <- counttraj%>%left_join(sel_prodpreds,by=c('gene_id','time'))%>%
		filter(!is.na(diff),!is.na(logFC))%>%
		# filter(logFC > quantile(logFC,.1))%>%
		identity
	varmodeldf%<>%group_by(gene_id)%>%mutate(logFC= logFC-logFC[1])
	varmodeldf%<>%group_by(gene_id)%>%mutate(diff= diff-diff[1])
	varmodeldf%<>%filter(gene_name%in%geneset)
	varmodeldf%>%{quicktest(.$diff,.$logFC)}

	varmodel = varmodeldf%>%lm(data=.,diff~logFC)
	b_all = varmodel$coef[2]
	b_r = 1
	allvar = sum(varmodel$residuals^2)
	var_p = sum(varmodeldf$se.y^2)
	var_r = sum(varmodeldf$se.x^2)
	var_pdt = allvar - (b_all/b_r)^2 * var_r - var_p
	var_mp = sum(varmodeldf%>%lm(data=.,diff~1)%>%.$residuals%>%.^2)
	var_explained = (var_mp - var_p - var_pdt) / (var_mp - var_p)
	names(var_explained)=NULL
	c(var_explained = var_explained,count_variance = var_r / (var_mp - var_p))

})
})

#now plot
plotfile<- here(paste0('plots/','figures/figure3/within_tp_var_explained','.pdf'))
cairo_pdf(plotfile)
withinstage_varexpldf%>%
	ggplot(.,aes(y=var_explained,x=assay,fill=assay))+
	facet_grid(.~time)+
	stat_identity(geom='bar')+
	scale_color_manual(name='assay',c('total'='blue','ribo'='green','ribo_traj'='red'))+
	scale_x_discrete(paste0('Assay'))+
	scale_y_continuous(paste0('% Variance Explained'))+
	ggtitle(paste0('Between - Timepoint Variance Explained'))+
	theme_bw()
dev.off()
normalizePath(plotfile)


#now plot
plotfile<- here(paste0('plots/','figures/figure3/tp_var_explained','.pdf'))
cairo_pdf(plotfile)
betweenstagevardf%>%
	ggplot(.,aes(y=var_explained,x=assay,fill=assay))+
	facet_grid(.~geneset)+
	stat_identity(geom='bar')+
	scale_color_manual(name='assay',c('total'='blue','ribo'='green'))+
	scale_x_discrete(paste0('Assay'))+
	scale_y_continuous(paste0('% Variance Explained'))+
	ggtitle(paste0('Between - Timepoint Variance Explained'))+
	theme_bw()
dev.off()
normalizePath(plotfile)

#not sure how to deal with our many variances as opposed to the single epsilon
#in the ref... I guess we just 


