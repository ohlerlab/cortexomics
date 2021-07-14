################################################################################
########Metagene plots and positions specific analysis
################################################################################

{	

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

################################################################################
########Do with deepshape etc instead
################################################################################
{
library(GenomicFeatures)
base::source(here::here('src/R/Rprofile.R'))
if(!exists('cdsgrl')) base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
if(!exists('iso_tx_countdata')) load('data/1_integrate_countdata.R')
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
stopifnot(length(allbamtbls)==22)
names(allbamtbls) = allbamtbls%>%str_extract('[^/]*?(?=.bam.ref)')
mainbamtbls <- allbamtbls%>%.[str_detect(.,'ribo_\\d+')]
bamtbl=mainbamtbls[1]
mainsamps = names(mainbamtbls)

if(!file.exists(here('data/ribocovtrs.rds'))){
	ribocovtrs = iso_tx_countdata$abundance[,mainsamps]%>%
		as.data.frame%>%
		rownames_to_column('transcript_id')%>%
		left_join(tx2genemap%>%set_colnames(c('transcript_id','gene_id')))%>%
		pivot_longer(-one_of('gene_id','transcript_id'))%>%
		separate(name,c('time','assay','rep'))%>%
		group_by(gene_id,transcript_id)%>%summarise(value=median(value))%>%
		group_by(gene_id)%>%
		filter(gene_id %in% highcountgenes)%>%
		slice(which.max(value))%>%.$transcript_id
	saveRDS(ribocovtrs,here('data/ribocovtrs.rds'))
}else{
	ribocovtrs<-readRDS(here('data/ribocovtrs.rds'))
}


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
}

################################################################################
########Collect transript level info for top transcripts
################################################################################
	
if(!file.exists(here('data/fpcovlist.rds'))){
	fpcovlist = allbamtbls%>%mclapply(mc.cores=8,function(bamtbl){
		bamtbl%>%
			str_interp('grep -e $')%>%
			fread(select=c(2,6,7))%>%
			set_colnames(c('transcript_id','start','readlen'))%>%
			mutate(transcript_id=trimids(transcript_id))%>%
			filter(transcript_id%in%highcountcovtrs)%>%
			filter(between(readlen,19,31))%>%
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
	trlens <- seqlengths(trspacecds)
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
	gr %<>% shift(fputrext[as.character(seqnames(gr))])
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
	
	bindatamats <- get_cds_bin_counts(epsitecovlist,trsp,STOPWINDSIZE,STARTWINDSIZE)
	saveRDS(bindatamats,here('data/bindatamats.rds'))
}else{
	bindatamats<-readRDS(here('data/bindatamats.rds'))
}

names(bindatamats)
#functions for norm
rownorm <-function(x) x %>%sweep(.,MARGIN=1,F='/',STAT=rowSums(.)%>%{pmax(.,min(.[.>0])/10I)})

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