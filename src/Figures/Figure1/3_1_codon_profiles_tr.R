{
library(abind)
library(tidyverse)
roundup <- function(x,n) n*ceiling(x/n)
rounddown <- function(x,n) n*floor(x/n)
number_ticks <- function(limits,n=3){
	out = c(seq(min(rounddown(limits,n)),- n,n),seq(0,max(roundup(limits,n)),by=n))
	out = out[between(out,limits[1],limits[2])]
	out
}

STARTCDSSIZE = 60
STOPCDSSIZE = 60
FPEXT = 36
STARTWINDSIZE = STARTCDSSIZE+FPEXT
TPUTREXT = 36
STOPWINDSIZE = STOPCDSSIZE+TPUTREXT

MINCDSSIZE = STARTCDSSIZE+STOPCDSSIZE+1

TOTBINS = (STARTWINDSIZE)+(STOPWINDSIZE)+1
}
# mystagecols = stagecols%>%setNames(names(stageconv))

################################################################################
########Do with deepshape etc instead
################################################################################
{
library(GenomicFeatures)
base::source(here::here('src/R/Rprofile.R'))
if(!exists('iso_tx_countdata')) load('data/1_integrate_countdata.R')
if(!exists('cdsgrl')) base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
# base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
trid2gid = cds%>%mcols%>%as.data.frame%>%select(transcript_id,gene_id)%>%{safe_hashmap(.[[1]],.[[2]])}


displaystagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
stageconv = names(displaystagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
stagecols <- displaystagecols%>%setNames(names(stageconv))

CODOONSOFINTEREST<-DNAStringSet(c('TCA','TCG','TCC','TCT'))
codons2scan <- c(CODOONSOFINTEREST,reverseComplement(CODOONSOFINTEREST))
# codons2scan <- c(CODOONSOFINTEREST,(CODOONSOFINTEREST))
ctrlcodons<-c('TTT','GGG','AAA','CCC')
codons2scan <- c(CODOONSOFINTEREST,ctrlcodons)
#codons2scan <- CODOONSOFINTEREST

FLANKCODS<-15

allbamtbls = Sys.glob('pipeline/deepshapebamdata/*.bam.reformat')
names(allbamtbls) = allbamtbls%>%str_extract('[^/]*?(?=.bam.ref)')
mainbamtbls <- allbamtbls%>%.[str_detect(.,'ribo_\\d+')]
bamtbl=mainbamtbls[1]
mainsamps = names(mainbamtbls)

# ribocovtrs = iso_tx_countdata$abundance[,mainsamps]%>%
# 	as.data.frame%>%
# 	rownames_to_column('transcript_id')%>%
# 	left_join(tx2genemap%>%set_colnames(c('transcript_id','gene_id')))%>%
# 	pivot_longer(-one_of('gene_id','transcript_id'))%>%
# 	separate(name,c('time','assay','rep'))%>%
# 	group_by(gene_id,transcript_id)%>%summarise(value=median(value))%>%
# 	group_by(gene_id)%>%
# 	filter(gene_id %in% highcountgenes)%>%
# 	slice(which.max(value))%>%.$transcript_id

# names(fpcovlist[[1]][['29']])%>%saveRDS('data/ribocovtrs.rds')
ribocovtrs <- readRDS(here('data/ribocovtrs.rds'))
ribocovtrs%>%n_distinct

cds2use <- cdsgrl[ribocovtrs]

cdsseq <- cds2use%>%{GenomicFeatures::extractTranscriptSeqs(.,x=fafile)}
allcodons=getGeneticCode()

highcountcovtrs = ribocovtrs[trid2gid[[ribocovtrs]]%in%highcountgenes]
trlens = exonsgrl%>%width%>%sum
trlensgr = trlens[ribocovtrs]%>%enframe('seqnames','end')%>%mutate(start=1)%>%GRanges

offsets <- read_tsv('ext_data/offsets_manual.tsv')
}

trseqinfo = Seqinfo(seqnames=ribocovtrs,seqlengths=trlens[ribocovtrs])

cdsstarts = cdsgrl[highcountcovtrs]%>%sort_grl_st%>%resize_grl(1)%>%unlist%>%
	pmapToTranscripts(exonsgrl[names(.)]%>%sort_grl_st)%>%
	{setNames(start(.),as.character(seqnames(.)))}
cdsends = cdsgrl[highcountcovtrs]%>%sort_grl_st%>%resize_grl(1,'end')%>%unlist%>%
	pmapToTranscripts(exonsgrl[names(.)]%>%sort_grl_st)%>%
	{setNames(start(.),as.character(seqnames(.)))}
trcds = GRanges(names(cdsstarts),IRanges(cdsstarts,cdsends))

if(!file.exists(here('data/fpcovlist.rds'))){
	fpcovlist = allbamtbls%>%mclapply(mc.cores=4,function(bamtbl){
		bamtbl%>%
			str_interp('grep -e $')%>%
			fread(select=c(2,6,7))%>%
			set_colnames(c('transcript_id','start','readlen'))%>%
			mutate(transcript_id=trimids(transcript_id))%>%
			filter(transcript_id%in%highcountcovtrs)%>%
			filter(between(readlen,25,31))%>%
			{GRanges(.$transcript_id,IRanges(.$start,w=1),readlen=.$readlen)}%>%
			{seqlevels(.) = highcountcovtrs ;.}%>%
			{seqlengths(.) = trlens[highcountcovtrs] ;.}%>%
			split(.,.$readlen)%>%
			lapply(coverage)
	})
	saveRDS(fpcovlist,here('data/fpcovlist.rds'))
}else{
	fpcovlist<-readRDS(here('data/fpcovlist.rds'))
	stopifnot(names(fpcovlist)==names(allbamtbls))
	stopifnot(all(ribocovtrs%in%names(fpcovlist[[1]][['29']])))
	inclusiontable(ribocovtrs,names(fpcovlist[[1]][['29']]))
}


if(!file.exists(here('data/psitecovlist.rds'))){
	psitecovlist = allbamtbls%>%mclapply(mc.cores=1,function(bamtbl){
		bamtbl%>%
			str_interp('grep -e $')%>%
			fread(select=c(2,6,7))%>%
			set_colnames(c('transcript_id','start','length'))%>%
			mutate(transcript_id=trimids(transcript_id))%>%
			filter(transcript_id%in%ribocovtrs)%>%
			inner_join(offsets%>%select(offset,length))%>%	
			mutate(start = start+offset)%>%
			{GRanges(.$transcript_id,IRanges(.$start,w=1))}%>%
			{seqlevels(.) = ribocovtrs ;.}%>%
			{seqlengths(.) = trlens[ribocovtrs] ;.}%>%
			coverage
	})
	saveRDS(psitecovlist,here('data/psitecovlist.rds'))
}else{
	psitecovlist<-readRDS(here('data/psitecovlist.rds'))
	stopifnot(names(psitecovlist)==names(allbamtbls))
	stopifnot(all(ribocovtrs%in%names(psitecovlist[[1]])))
}

if(!file.exists(here('data/psitecovnorm.rds'))){
	psitecovnorm = lapply(psitecovlist,function(psitecov){psitecov/(psitecov%>%sum%>%pmax(.,1))})
	saveRDS(psitecovnorm,here('data/psitecovnorm.rds'))
}else{
	psitecovnorm<-readRDS(here('data/psitecovnorm.rds'))
}

get_utr_exts <- function(cov,trspacecds,FPEXT,TPUTREXT){
	trlens = cov%>%runLength%>%sum
	seqnms = names(trlens)
	fputrlens = start(trspacecds[seqnms])-1
	fputrext = pmax(0,FPEXT - fputrlens)%>%setNames(seqnms)

	tputrlens = trlens - end(trspacecds[seqnms])
	tputrext = pmax(0,TPUTREXT - tputrlens)%>%setNames(seqnms)

	list(fputrext,tputrext)
}


ext_cov <- function(cov,fputrext,tputrext){
	gr = as(cov,"GRanges")
	seqnms = names(fputrext)
	stopifnot(identical(seqnms,names(tputrext)))
	stopifnot(identical(seqnms,names(cov)))
	seqlengths(gr)[seqnms]%<>%add(fputrext+tputrext)
	gr %<>% shift(fputrext)
	coverage(gr,weight='score')
}
{
#filter for only long enough ones, if we're doing windows
trspacecds = pmapToTranscripts(cdsgrl[ribocovtrs],exonsgrl[ribocovtrs])
trspacecds%<>%unlist

ltrspacecds = trspacecds
longcdstrs = names(ltrspacecds)[ltrspacecds%>%width%>%`>`(MINCDSSIZE)]
ltrspacecds = ltrspacecds[longcdstrs]
#see how many we eliminate for length reasons
message(psitecovlist[[1]]%>%length)
message(length(ltrspacecds))
#spand the cds as necessary

cov = psitecovlist[[1]][longcdstrs]

c(fputrext,tputrext) %<-% (psitecovlist[[1]][longcdstrs]%>%get_utr_exts(ltrspacecds,FPEXT,TPUTREXT))
stopifnot(all(names(fputrext)==longcdstrs))
stopifnot(all(names(tputrext)==longcdstrs))

if(!file.exists(here('data/epsitecovlist.rds'))){
	epsitecovlist  <- psitecovlist%>% 
		lapply(.%>%.[longcdstrs])%>%
		lapply(ext_cov,fputrext,tputrext)
	saveRDS(epsitecovlist,here('data/epsitecovlist.rds'))
}else{
	epsitecovlist<-readRDS(here('data/epsitecovlist.rds'))

}
stopifnot(all(names(epsitecovlist[[1]])==longcdstrs))
# stopifnot(seqnames(epsitecovlist[[1]])==longcdstrs)
#also lift cds to this new space
eltrspacecds <- ltrspacecds[longcdstrs]
seqlengths(eltrspacecds)[longcdstrs] %<>% add(fputrext+tputrext)
eltrspacecds %<>% shift(fputrext)
eltrspacecds%<>%resize(width(.)+FPEXT,'end')%>%resize(width(.)+TPUTREXT,'start')
#also the cds in this space

# highcovtrs = psitecovlist[[1]]%>%sum%>%sort%>%tail(500)%>%names%>%intersect(names(ltrspacecds))
# highcovtrs = names(ltrspacecds)

cov = epsitecovlist[[1]]
}
# file.remove(here('data/bindatamats.rds'))
if(!file.exists(here('data/bindatamats.rds'))){
	bindatamats = lapply(epsitecovlist,function(cov){
		midwind = eltrspacecds%>%
			resize(width(.)-(STOPWINDSIZE),'end')%>%
			resize(width(.)-(STARTWINDSIZE),'start')
		midmat = cov[midwind]%>%
			sum%>%#compress these variable length Rles down to 1 value per gene
			matrix#make a 1 column matrix
		# midmat = midmat%>%mean
		stwind = eltrspacecds%>%resize(STARTWINDSIZE,ignore.strand=TRUE)
		startmat = cov[stwind]%>%as.matrix
	#	%>%colMeans(na.rm=T)
		endwind = eltrspacecds%>%resize(STOPWINDSIZE,fix='end',ignore.strand=TRUE)
		endmat = cov[endwind]%>%as.matrix
		#%>%colMeans(na.rm=T)
		#sum codons
		# startmat%<>%matrix(nrow=3)%>%colSums
		# endmat%<>%matrix(nrow=3)%>%colSums
		message('.')
		#output as a data frame
		out = cbind(startmat,midmat,endmat)
		# c(startmat,midmat,endmat)%>%enframe('start','signal')%>%mutate(section=c(rep('AUG',150/3),'middle',rep('end',150/3)))
	})
	saveRDS(bindatamats,here('data/bindatamats.rds'))
}else{
	bindatamats<-readRDS(here('data/bindatamats.rds'))
}

names(bindatamats)
#functions for norm
rownorm <-function(x) x %>%sweep(.,MARGIN=1,F='/',STAT=rowSums(.)%>%{pmax(.,min(.[.>0])/10)})

#matrix <- matrix(1,6,24)
codmerge <- function(matrix,size=3) {
	matrix%>%t%>%rowsum(.,ceiling((1:nrow(.))/size),na.rm=T)%>%t
}

#name and check
{
matgenes = longcdstrs
bindatamats%<>%lapply(function(mat)set_rownames(mat,longcdstrs))
stopifnot(bindatamats[[1]]%>%dim%>%identical(c(length(longcdstrs),as.integer(TOTBINS))))
#why are the start and end mats so similiar?
startinds = 1:(STARTWINDSIZE)
endinds = (1+STARTWINDSIZE):TOTBINS
}


# bindatamats[[1]][,1:(STARTWINDSIZE*3)]%>%rownorm%>%apply(2,is.na)%>%colMeans%>%txtplot
# bindatamats[[1]][1,]

#function that outputs the new dimesnions of each tr given the required extension, the
#length of each ones existing fputr, and the length of each ones existing tputr
#function that gets the utr lengths from the cds and the exonsgrl

#function that shifts a grange in the actal trspace to a grange in the extended tr space
#codmerge to get codons not locations
get_metasignaldf<-function(bindatamats,genelist=TRUE){
	metasignaldf = bindatamats%>%
		map_df(.id='sample',.%>%.[genelist,]%>%
			rownorm%>%
			# {./sum(.)}%>%
			colMeans(na.rm=T)%>%
			enframe('start','signal'))%>%
		group_by(sample)%>%
		mutate(section=c(rep('AUG',STARTWINDSIZE),'middle',rep('end',STOPWINDSIZE)))
	#assign stage
	metasignaldf$stage <- metasignaldf$sample%>%str_extract('[^_]+')
	# metasignaldf%<>%mutate(start = as.numeric(start) - ifelse(section=='AUG',STARTWINDSIZE+1,TOTBINS- STOPWINDSTART ))
	metasignaldf%<>%mutate(start = as.numeric(start) - ifelse(section=='AUG',1+FPEXT,1+FPEXT+STARTCDSSIZE+STOPCDSSIZE-2 ))
	metasignaldf
}
get_metaplot <- function(metasignaldf,ylims=NULL){
	if(!'fraction'%in%colnames(metasignaldf)) metasignaldf$fraction='total'
	if(!'sample'%in%colnames(metasignaldf)) metasignaldf$sample=metasignaldf$stage
	
	metasignaldf%>%
	filter(section!='middle')%>%
	split(.,.$section)%>%map( .%>%
		# slice_by(sample,8)%>%
		# filter(position%%3 == 0)%>%
		{
			isfirst = .$section[1]=='AUG'
			qplot(data=.,group=sample,color=stage,x=start,y=signal,geom='blank')+
			geom_line()+
			scale_x_continuous(name='position (bp)',
				limits=if(.$section[1]=='AUG') c(-FPEXT,STARTCDSSIZE) else c(-STOPCDSSIZE,TPUTREXT) ,
				# limits=if(isfirst) c(0,(STARTWINDSIZE)-1) else c(1-(STOPWINDSIZE),0) ,
				minor_breaks=number_ticks,breaks=partial(number_ticks,n=12)
			)+
			facet_grid( fraction ~ section,scale='free_x')+
			scale_color_manual(values=stagecols)+
			scale_y_continuous(name='Mean Psite Count / CDS Total',limits=ylims)+
			theme_bw()+
			# theme_minimal()+
			theme(panel.grid = element_line(color=I('grey')))
		}
	)%>%
	ggarrange(plotlist=.,ncol=2,common.legend=T)
}
metasignaldf_stgrp <- get_metasignaldf(bindatamats[mainsamps],longcdstrs) %>% 
	group_by(stage,section,start)%>%
	summarise(signal=mean(signal))
{
library(rlang)
plotfile<-'plots/figures/figure1/fig1c_myribowaltz_allsec_stageov.pdf'%T>%pdf(h=6,w=12)
rwplot <- metasignaldf_stgrp%>%get_metaplot(c(0,0.012))
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message
}



nonmainsamps = names(bindatamats)%>%setdiff(mainsamps)

metasignaldf_stgrp <- get_metasignaldf(bindatamats[nonmainsamps],longcdstrs) %>% 
	mutate(fraction = str_extract(sample,'80S|Poly'))%>%
	mutate(stage = sample%>%str_extract('(?<=80S|Poly).*?(?=_)'))%>%
	group_by(stage,section,start,fraction)%>%
	summarise(signal=mean(signal))
{
library(rlang)
plotfile<-'plots/figures/figure1/fig1c_myribowaltz_frac_metaplots.pdf'%T>%pdf(h=6,w=12)
rwplot <- metasignaldf_stgrp%>%get_metaplot(c(0,0.006))
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message
}

#Let's define a subset of the data - those with the apparently biggest shift over time.
#Let's also define a 
inf.omit = function(x) x[is.finite(x)]


# if(F){

# {

# earlymats = (bindatamats[[1]]+bindatamats[[2]])%>%rownorm
# latemats = (bindatamats[[9]]+bindatamats[[10]])%>%rownorm
# earlystdens = earlymats[,startinds[1:15]]%>%rowSums
# latestdens = latemats[,startinds[1:15]]%>%rowSums
# has_finitedens = is.finite(log2(earlystdens) - log2(latestdens))
# startchange = (latestdens+1e-12) / (earlystdens+1e-12)
# # has_finitedens = is.finite(startchange)
# moststart_incr = startchange[has_finitedens]%>%sort%>%tail(1000)%>%names
# moststart_decr = startchange[has_finitedens]%>%sort%>%head(1000)%>%names
# #
# counts = (rowSums((bindatamats[[1]]+bindatamats[[2]]))+rowSums(bindatamats[[9]]+bindatamats[[10]]))%>%log2
# isaboveq <- function(x,p) x > quantile(x,p,na.rm=T)
# isbelowq <- function(x,p) x < quantile(x,p,na.rm=T)

# plotstchange = startchange%>%log2

# zeroclass = case_when(
# 	plotstchange>20 ~ 'Zero to Smth',
# 	plotstchange< -20~ 'Smth to Zero',
# 	TRUE ~ "Smth to Smth"
# )




# }

# {
# #now plot
# plotfile<- here(paste0('plots/','startchange_vs_signal','.pdf'))
# pdf(plotfile)
# print(
# 	qplot(data=data.frame(counts,plotstchange,zeroclass),x=counts,y=plotstchange,color=zeroclass)+
# 	facet_grid(zeroclass~.,scales='free_y')
# )
# # qplot(x=counts,y=plotstchange,color=(plotstchange>4) & (plotstchange/counts)<0.6)
# dev.off()
# message(normalizePath(plotfile))
# }

# {
# metasignaldf_stgrp <- get_metasignaldf(bindatamats,matgenes[counts>quantile(counts,.9,na.rm=T)]) %>% 
# 	group_by(stage,section,start)%>%
# 	summarise(signal=mean(signal)) 
# library(rlang)
# plotfile<-'plots/figures/figure1/metaplottest.pdf'%T>%pdf(h=6,w=12)
# rwplot <- metasignaldf_stgrp%>%get_metaplot
# print(rwplot)
# dev.off()
# normalizePath(plotfile)%>%message
# }

# metasignaldf_stgrp <- get_metasignaldf(bindatamats,moststart_incr) %>% 
# 	group_by(stage,section,start)%>%
# 	summarise(signal=mean(signal)) 
# {
# library(rlang)
# plotfile<-'plots/figures/figure1/metaplottestmoststart_incr.pdf'%T>%pdf(h=6,w=12)
# rwplot <- metasignaldf_stgrp%>%get_metaplot
# print(rwplot)
# dev.off()
# normalizePath(plotfile)%>%message
# }


# metasignaldf_stgrp <- get_metasignaldf(bindatamats,moststart_decr) %>% 
# 	group_by(stage,section,start)%>%
# 	summarise(signal=mean(signal)) 
# {
# library(rlang)
# plotfile<-'plots/figures/figure1/metaplottestmoststart_decr.pdf'%T>%pdf(h=6,w=12)
# rwplot <- metasignaldf_stgrp%>%get_metaplot
# print(rwplot)
# dev.off()
# normalizePath(plotfile)%>%message
# }

# bindatamats[[1]]%>%.[,1: STARTWINDSIZE]%>% codmerge%>%codmerge(4)

# highmincovtrs = psitecovlist%>%map(.%>%mean)%>%purrr::reduce(pmin)%>%sort%>%tail(100)%>%names%>%intersect(names(ltrspacecds))


# }



# ################################################################################
# ########Try pos specific analysis of the frac data
# ################################################################################
# {

# #compress our bins into 4 codon units
# data_binned <- bindatamats[nonmainsamps]%>%lapply(function(mat){cbind(
# 	mat%>%.[,1:(FPEXT+STOPCDSSIZE)]%>%codmerge%>%codmerge(4)%>%rowSums,
# 	mat%>%.[,FPEXT+STOPCDSSIZE+1,drop=FALSE],
# 	mat%>%.[,(FPEXT+STOPCDSSIZE+2):TOTBINS,drop=FALSE]%>%codmerge%>%codmerge(4)%>%rowSums
# )})

# binneddata = data_binned%>%imap(~.x%>%as.data.frame%>%
# 	set_colnames(c('start','mid','stop'))%>%
# 	rownames_to_column('tr_id')%>%
# 	gather(bin,val,-tr_id)%>%
# 	mutate(bin = bin%>%
# 		as.character%>%
# 		str_replace('V9','main')%>%
# 		str_replace('(.*)\\.1','e_\\1')%>%
# 		str_replace('^(\\d)','s_\\1'))%>%
# 	{colnames(.)[3]=.y;.}
# 		)
	
# binneddata = binneddata%>%purrr::reduce(.f=left_join,by=c('tr_id','bin'))

# binneddata$bin%>%unique%>%print

# ribofrac_design = data.frame(sample=nonmainsamps)%>%	
# 	mutate(fraction = str_extract(sample,'80S|Poly'))%>%
# 	mutate(time = sample%>%str_extract('(?<=80S|Poly).*?(?=_)'))%>%
# 	mutate(rep = sample%>%str_extract('[0-9]+$'))%>%
# 	arrange(fraction=='80S',time=='P0',time=='E16')%>%
# 	set_rownames(.$sample)
# ribofrac_design$fraction%<>%as_factor
# ribofrac_design$timefrac = paste0(ribofrac_design$time,ribofrac_design$fraction)

# library(DEXSeq)

# dxd  <- DEXSeqDataSet(
# 	binneddata%>%select(-tr_id,-bin)%>%mutate_at(vars(everything()),list(~replace_na(.,0)))%>%as.matrix,
# 	featureID = paste0(binneddata$tr_id,binneddata$bin	),
# 	groupID = paste0(binneddata$tr_id),
# 	design = ~ sample+timefrac*exon,
# 	sampleData =ribofrac_design
# )

# sizeFactors(dxd) <- dxd%>%counts%>%estimateSizeFactorsForMatrix
# dxd%<>%estimateDispersions
# dxd = testForDEU(dxd)
# dxd = estimateExonFoldChanges( dxd, fitExpToVar="timefrac")

# dxr1 = DEXSeqResults( dxd )%>%as.data.frame

# library(txtplot)
# dxr1$padj%>%na.omit%>%txtdensity

# }

# dxr1[rownames(dxr1)%>%str_detect('start'),]$padj%>%na.omit%>%txtdensity

# dxr1$gene_name = gid2gnm[[trid2gid[[dxr1$groupID]]]]

# dxr1[dxr1$gene_name%>%str_detect('Satb2'),]







