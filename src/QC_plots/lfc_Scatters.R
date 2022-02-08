################################################################################
################################################################################

base::source(here::here('src/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("src/Figures/load_annotation.R")
}
if(!exists('tx_countdata')) {
	base::source("src/Figures/1_integrate_countdata.R")
	# load('data/1_integrate_countdata.R')
}

gnm2gid <- ids_nrgname%>%distinct(gene_id,gene_name)%>%
	{setNames(.$gene_id,.$gene_name)}

#
protgnmtrids<-readRDS('data/protgnmtrids.rds')

proteinmsdata<-readRDS('data/proteinmsdata.rds')
proteinmsdata$ms_id = proteinmsdata[[1]]
proteinmsdatamatch <- proteinmsdata
proteinmsdatamatch$cdsidx = match(proteinmsdatamatch$gene_name,fmcols(cdsgrl,gene_name))
proteinmsdatamatch <- proteinmsdatamatch%>%filter(!is.na(cdsidx))
ms_transcripts <- protgnmtrids$tr_id
#select iBAQ and input cols
sigcols <- proteinmsdatamatch%>%colnames%>%str_subset('iBAQ')%>%str_subset('input')
proteinmsdatamatch%<>%select(ms_id,gene_name,g_id,one_of(sigcols))
ms_mat <- proteinmsdatamatch[,-c(1:3)]%>%as.matrix
#
rownames(ms_mat) <- proteinmsdatamatch$ms_id
colnames(ms_mat) <- colnames(ms_mat)%>%
	str_replace('iBAQ\\.','')%>%
	str_replace('input','MS')%>%
	str_replace('rep','')%>%
	str_replace('p(\\d)','\\1')
ms_mat = proDA::median_normalization(ms_mat)

################################################################################
########Import it all to get tr length scaled counts
################################################################################
{
library(tximport)
library(tidyverse)

rnasalmonfiles = Sys.glob(here('pipeline/salmon/data/*total*/quant.sf'))
dpoutfiles = Sys.glob(here('pipeline/deepshapeprime/fakesalmonfiles/*ribo*/*'))

allquantfiles = c(rnasalmonfiles,dpoutfiles)
names(allquantfiles) <- allquantfiles%>%dirname%>%basename

dptrs = dpoutfiles[[1]]%>%fread%>%.$Name
salmontrs = allquantfiles[[1]]%>%fread%>%.$Name%>%str_extract('ENSMUST\\w+')
# inclusiontable(dptrs,salmontrs)
trs = intersect(dptrs,salmontrs)


####KEY STEP
ms_trs = intersect(trs,ms_transcripts)

ms_tx_countdata = tximport(files=allquantfiles,
	ignoreTxVersion=TRUE,
	tx2gene=mcols(cds)[,c('transcript_id','gene_id')]%>%as.data.frame%>%distinct(transcript_id,gene_id),
	type='salmon',
	countsFromAbundance='scaledTPM',
	importer=function(file){
		read_tsv(file,col_types=cols())%>%
			mutate(Name=str_extract(Name,'ENSMUST\\w+'))%>%
			filter(Name%in%ms_trs)%>%arrange(match(Name,trs))
})

randomround = function(x)floor(x)+rbinom(length(x),1,x%%1)
ms_tx_countdata$counts%<>%as.data.frame%>%
	mutate_all(randomround)%>%
	as.matrix%>%
	set_rownames(rownames(ms_tx_countdata$counts))

}

{

sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
proDAfitms<-readRDS('data/proDAfitms.rds')
sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
contrdf<-readRDS('data/contrdf.rds')
sel_ms_mat<-readRDS('data/sel_ms_mat.rds')



countcontr_df <- readRDS(here('data/countcontr_df.rds'))


}

# stop()

stage_count='E13'
stage_ms='E13'


tps = sel_ms_mat%>%colnames%>%str_extract('[^_]+')%>%unique
tps %<>%setNames(.,.)
cassays = c('total','ribo')%>%setNames(.,.)

allxtail = Sys.glob('pipeline/xtail/*')%>%map_df(.id='time',fread)%>%group_by(gene_name)
allxtail$gene_id = gnm2gid[ allxtail$gene_name]

txnchangedf <- countcontr_df%>%filter(assay=='all')%>%group_by(gene_id)%>%
  mutate(sig = (adj.P.Val < 0.05)& (abs(logFC)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (logFC > 0))),
    down = as.numeric(any(sig & (logFC < 0)))
  )

ribochangedf <- countcontr_df%>%filter(assay=='ribo')%>%group_by(gene_id)%>%
  mutate(sig = (adj.P.Val < 0.05)& (abs(logFC)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (logFC > 0))),
    down = as.numeric(any(sig & (logFC < 0)))
  )

techangedf <- allxtail%>%group_by(gene_id,gene_name)%>%
  mutate(sig = (adj_p_value < 0.05)& (abs(log2fc)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (log2fc > 0))),
    down = as.numeric(any(sig & (log2fc < 0)))
  )

#or use ms
mschangedf = contrdf%>%
  group_by(gene_id)%>%
  mutate(sig = (adj_pval < 0.05)& (abs(diff)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (diff > 0))),
    down = as.numeric(any(sig & (diff < 0)))
  )

mschangedf%>%group_by(up,down)%>%tally
techangedf%>%group_by(up,down)%>%tally
ribochangedf%>%group_by(up,down)%>%tally
techangedf%>%group_by(up,down)%>%tally

genesets = list(
	all = techangedf$gene_id,
	txn_up = txnchangedf%>%filter(up==1)%>%.$gene_id,
	txn_down = txnchangedf%>%filter(down==1)%>%.$gene_id,
	te_up = techangedf%>%filter(up==1)%>%.$gene_id,
	te_down = techangedf%>%filter(down==1)%>%.$gene_id,
	ms_up = mschangedf%>%filter(up==1)%>%.$gene_id,
	ms_down = mschangedf%>%filter(down==1)%>%.$gene_id
)

genesets%>%saveRDS('data/genesets.rds')


avtpms = ms_tx_countdata$abundance[mssharedgenes,]%>%
	as.data.frame%>%
	rownames_to_column('gene_id')%>%
	gather(contrast,TPM,-gene_id)%>%
	separate(contrast,into=c('time','assay','rep'))%>%
	group_by(time,assay,gene_id)%>%summarise(TPM=mean(TPM))%>%
	unite(contrast,c('time','assay'))
ms_metadf <- readRDS('data/ms_metadf.rds')

itime = 'E13'

msscatterlist = lapply(cassays,function(countassay){
	lapply(genesets,function(geneset){
	lapply(tps,function(stage_count){
	lapply(tps,function(stage_ms){
	countname = paste0('log10(TPM) ',countassay,'_',stage_count)
	MSname = paste0('MS_',stage_ms)
	misscol = paste0(stage_ms,'_anymissing')
	istagepresent <- ms_metadf%>%filter(!(!!sym(misscol)))

	scatdf = avtpms%>%
		filter(!str_detect(contrast,'_TE'))%>%
		separate(contrast,c('time','assay')	)%>%
		filter(time==stage_count)%>%
		filter(assay==countassay)%>%
		filter(gene_id%in%geneset)%>%
		left_join(sel_prodpreds%>%filter(time==stage_ms)%>%select(gene_id,diff))%>%
		filter(gene_id %in% mssharedgenes)%>%
		inner_join(istagepresent)%>%
		filter(TPM>1)%>%
		mutate(TPM = log10(TPM))%>%
		select(gene_id,countvals=TPM,msvals=diff)
	#
	testtidy = tidy(cor.test(scatdf$countvals,scatdf$msvals))
	acor_label=paste0('rho = ',round(testtidy$conf.low,3),' - ',round(testtidy$conf.high,3))
	acor_label=paste0('rho = ',round(testtidy$estimate,3),'\n','pval = ',round(testtidy$p.value,3))
	# pdf('tmp.pdf')
	# gplot = heatscatter((scatdf$countvals),scatdf$msvals,ggplot=TRUE)+
	# 	scale_x_continuous(countname)+
	# 	scale_y_continuous(MSname)+
	# 	ggtitle(acor_label)
	# dev.off()
	# list(gplot,testtidy)
	list(NA,testtidy)
})
})
})
})


{

countmscors = msscatterlist%>%map_df(.id='assay',.%>%map_df(.id='geneset',.%>%map_df(.id='stage_count',.%>%map_df(.id='stage_ms',~.[[2]]$estimate))))

#now plot
plotfile<- here(paste0('plots/','QC_plots/cortilesribo','.pdf'))
pdf(plotfile,w=10,h=15)
print(countmscors%>%
	ggplot(data=.,aes(fill=cor,x=stage_count,y=stage_ms))+
	geom_tile()+
	geom_text(aes(label = round(cor,2)))+
	facet_grid(geneset~assay)+
	scale_x_discrete(paste0('Stage Ribo-seq'))+
	scale_y_discrete(paste0('Stage MS'))+
	# scale_fill_gradient2(low='#8904B1',mid='white',high='#FF8000',limits = range(countmscors$cor))+
	scale_fill_gradient2(low='purple',mid='white',high='orange',midpoint = mean(range(countmscors$cor)),limits = range(countmscors$cor))+
	ggtitle(paste0('Stage-stage correlations MS vs Ribo-seq'))+
	theme_bw()
)
dev.off()
message(normalizePath(plotfile))

techangedf = read_tsv(here('tables/xtailTEchange.tsv'))
#now plot
plotfile<- here(paste0('plots/','QC_plots/cortilesribo','.pdf'))
pdf(plotfile,w=10,h=15)
print(countmscors%>%
	ggplot(data=.,aes(fill=cor,x=stage_count,y=stage_ms))+
	geom_tile()+
	geom_text(aes(label = round(cor,2)))+
	facet_grid(geneset~assay)+
	scale_x_discrete(paste0('Stage Ribo-seq'))+
	scale_y_discrete(paste0('Stage MS'))+
	# scale_fill_gradient2(low='#8904B1',mid='white',high='#FF8000',limits = range(countmscors$cor))+
	scale_fill_gradient2(low='purple',mid='white',high='orange',midpoint = mean(range(countmscors$cor)),limits = range(countmscors$cor))+
	ggtitle(paste0('Stage-stage correlations MS vs Ribo-seq'))+
	theme_bw()
)
dev.off()
message(normalizePath(plotfile))

}

{
################################################################################
########Time correlations
################################################################################
countpred_df<-readRDS('data/countpred_df.rds')
tx_countdata<-readRDS('data/tx_countdata.rds')


cests = countpred_df%>%filter(str_detect(contrast,'total'))%>%distinct(gene_id,contrast,logFC)%>%spread(contrast,logFC)
rcests = countpred_df%>%filter(str_detect(contrast,'ribo'))%>%distinct(gene_id,contrast,logFC)%>%spread(contrast,logFC)
pests = sel_prodpreds%>%distinct(gene_id,time,diff)%>%spread(time,diff)

rnacors = pests%>%left_join(cests)%>%.[,-1]%>%apply(1,function(x)cor(x[1:5],x[6:10]))
ribocors = pests%>%left_join(rcests)%>%.[,-1]%>%apply(1,function(x)cor(x[1:5],x[6:10]))

corggdf = bind_rows(.id='assay',rnaseq=data.frame(cor=rnacors),ribo=data.frame(cor=ribocors))
#now plot
plotfile<- here(paste0('plots/','timecordensity','.pdf'))
pdf(plotfile)
print(corggdf%>%
	ggplot(.,aes(x=cor,color=assay))+
	geom_density()+
	# facet_grid(assay~.)+
	# scale_color_discrete(name='assay',colorvals)+
	scale_x_continuous(paste0('Temporal Correlation'))+
	ggtitle(paste0('Correlations for Riboseq vs RNAseq'))+
	theme_bw()
)
dev.off()
normalizePath(plotfile)


gnm2gid <- ids_nrgname%>%distinct(gene_id,gene_name)%>%
	{setNames(.$gene_id,.$gene_name)}
cassaynames = c('RNAseq','Riboseq','TE')%>%setNames(c('all','ribo','TE'))

stepcountcontrdf <- readRDS(here('data/stepcountcontrdf.rds'))
stepstepprotcontrdf<-readRDS('data/stepprotcontrdf.rds')

xtail_stepwise <- Sys.glob('pipeline/xtail/xtail*v*')%>%
	setNames(.,str_extract(.,'(?<=_)[^_]+'))%>%
	map_df(.id='time',read_tsv)%>%
	mutate(gene_id = gnm2gid[gene_name])%>%
	mutate(assay='TE')%>%
	select(gene_id,time,logFC=log2fc,assay,adj.P.Val=adj_p_value)

testepcountcontrdf <- stepcountcontrdf%>%
	bind_rows(xtail_stepwise)
stepcountcontrdf%>%colnames%>%dput
cassay='ribo'

stop()

cassay='TE'
v_mrna=FALSE
v_msR=TRUE

for(cassay in c('all','ribo','TE')){
	for(v_mrna in c(TRUE,FALSE)){
		for(v_msR in c(TRUE,FALSE)){
			if((v_mrna) & cassay!='TE') next

	library(broom)
	yaxisfoldchanges = stepstepprotcontrdf%>%
		mutate(psig = adj_pval<0.05)%>%
		select(gene_id,time,prot_logFC=diff,psig)
	yaxassay='MS'
	if((cassay=='TE')&(v_mrna) ){
		yaxisfoldchanges = testepcountcontrdf%>%
			filter(assay=='all')%>%
			mutate(psig = adj.P.Val<0.05)%>%
		  select(gene_id,time,prot_logFC=logFC,psig)
		 yaxassay='RNA-seq'
	}
	if(v_msR){
		if(cassay!='TE')next
		if(v_mrna)next
		denom = testepcountcontrdf%>%
			filter(assay=='all')%>%
			mutate(rsig = adj.P.Val<0.05)%>%
		  select(gene_id,time,rna_logFC=logFC)
		yaxisfoldchanges = yaxisfoldchanges%>%inner_join(denom,by=c('gene_id','time'))%>%mutate(RM_diff=prot_logFC-rna_logFC)
		yaxisfoldchanges <- yaxisfoldchanges%>%mutate(prot_logFC = RM_diff)
		yaxisfoldchanges <- yaxisfoldchanges%>%mutate(psig = FALSE)
		yaxassay = 'MS_RNA_Ratio'
	}
	#
	steplcdf = testepcountcontrdf%>%
		mutate(csig = adj.P.Val<0.05)%>%
		select(gene_id,time,logFC,assay,csig)%>%
		filter(assay==cassay)%>%
		inner_join(yaxisfoldchanges,by=c('gene_id','time'))%>%
		group_by(gene_id)
	#
	csigname = paste0(cassaynames[cassay],' Sig')
	steplcdf%<>%mutate(sigset = case_when(
		csig & psig ~ 'Both Sig',
		psig ~ str_interp('${yaxassay} Sig'),
		csig ~ csigname,
		TRUE ~ 'None'
	))
	#
	lfccorlabels = steplcdf%>%group_by(time)%>%filter(is.finite(logFC),is.finite(prot_logFC))%>%summarise(tidy(cor.test(logFC,prot_logFC)))%>%
		mutate(label = paste0('rho = ',round(conf.low,3),'-',round(conf.high,3),'\np = ',format(p.value,sci=T,digits=2)%>%str_replace('.*0.0.*',' < 1.0e-300')
	))
	sisetcolors = c('None'='#D8D8D8',csigname='#000000',
		'MS Sig'='#04B404','Both Sig'='#e81887')
	if(cassay=='TE')sisetcolors['MS Sig']='#241e20'
	if(cassay=='TE')sisetcolors['csigname']='#159244'
	names(sisetcolors)%<>%str_replace('MS Sig',str_interp('${yaxassay} Sig'))
	# BOTH sig (); Riboseq or RNAseq sig (); MS sig (); not sig ()
	names(sisetcolors)%<>%str_replace('csigname',csigname)
	#steplcdf%>%arrange(-match(sigset,names(sisetcolors)))%>%.$sigset%>%setdiff(names(sisetcolors))
	#now plot
	p<- steplcdf%>%
		arrange(match(sigset,names(sisetcolors)))%>%
		# ungroup%>%distinct(sigset)
		mutate(sigset = sigset%>%as_factor)%>%
		ggplot(aes(x=logFC,y=prot_logFC,color=sigset))+
		facet_grid(.~time)+
		geom_text(data=lfccorlabels,aes(label=label,color=NULL),y=I(-Inf),x=I(-Inf),
				vjust=0,hjust=0,legend=F,show_guide=F)+
		geom_point(size=I(0.2),aes(alpha=sigset!='None'))+
		scale_color_manual(name='colorname',values=sisetcolors)+
		scale_x_continuous(paste0(cassay,' Fold Change'),limits=c(-5,5))+
		scale_y_continuous(paste0(yaxassay,' Fold Change'),limits=c(-5,5))+
		ggtitle(paste0('Fold Change Comparison - ',yaxassay,' vs ',cassaynames[cassay]))+
		scale_alpha_discrete(guide=F)+
		theme_bw()
	if((cassay=='TE')&(v_mrna)) p = p+coord_flip()	
	#
	plotfile<- here(paste0('plots/','QC_plots/lfc_cors_',cassay,ifelse(v_mrna,'vmrna',''),ifelse(v_msR,'v_msR',''),'.pdf'))
	pdf(plotfile,w=16,h=4)
	print(p)
	dev.off()
	message(normalizePath(plotfile))
	}
}
}

}