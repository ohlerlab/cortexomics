################################################################################
########This figure produces locus plots with mRNA,Riboseq and MS data
################################################################################
	
#function that takes all the granges in a region, and compresses them down so introns
#are all a fixed distance
library(Gviz)
library(rtracklayer)
library(tidyverse)
library(magrittr)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)

peptidemsfile = here('ext_data/MS_Data_New/Ages_Brain_PEP_summ.txt')

mymem<-projmemoise(function(){rnorm(1)})

if(!exists("fafileob")) {
  base::source("src/Figures/Figure0/0_load_annotation.R")
}


tpcols <- c(E13='#214098',E145='#2AA9DF',E16='#F17E22',E175='#D14E28',P0='#ED3124')


options(ucscChromosomeNames=FALSE) 


#this is the function that shifts a bunch of ranges to minimize the space between them, down to a ceiling maxwidth
maxwidth=1
get_comp_shift<-function(rangestoshift,maxwidth=1){
	rangegaps <- rangestoshift%>%{strand(.)<-'*';.}%>%gaps%>%.[-1]
	rangegaps%>%width
	#maxwidth <- 4
	toshift<-rep(0,length(rangestoshift))
	i=3
	for(i in seq_along(rangegaps)){
		adj <- width(rangegaps[i])-maxwidth
		if(adj>0){
			rightofgap <- start(rangestoshift)>end(rangegaps[i])
			toshift[rightofgap] <- toshift[rightofgap] + adj
		}
	}
	toshift %<>% {. + min(start(rangestoshift))-1}
	#rangestoshift<-GenomicRanges::shift(rangestoshift,-rangestoshift$toshift)
	toshift
}




########

offsets<-read_tsv(here('ext_data/offsets_manual.tsv'))
totbams = Sys.glob('pipeline/star/data/*total*/*[12].bam')
ribobams = Sys.glob('pipeline/star/data/*ribo*/*[12].bam')
totbams
stopifnot(length(totbams)>0)
stopifnot(length(ribobams)>0)
ribobams%<>%setNames(.,basename(dirname(.)))
totbams%<>%setNames(.,basename(dirname(.)))

if(!exists('peptidemsdata')){
	peptidemsdata = fread(peptidemsfile)%>%
	rename('gene_name'=Gene.names)

	peptidemsdata%>%colnames
	normpepmat = peptidemsdata%>%select(matches('Intensity.*input'))%>%as.matrix%>%
		proDA::median_normalization(.)

	normpepmat = colnames(normpepmat)%>%str_extract('(?<=\\.).*?(?=_input)')%>%
		str_replace('p','')%>%
		{split(seq_along(.),.)}%>%
		lapply(.,function(colgrp){
			normpepmat[,colgrp]%>%rowMedians(na.rm=T)
		})%>%
		simplify2array
	rownames(normpepmat) <- peptidemsdata$Sequence

	extrseq = function(gr)GenomicFeatures::extractTranscriptSeqs(gr,x=fafileob)


	libstg = c(ribobams,totbams)%>%dirname%>%basename%>%str_split_fixed('_',3)%>%.[,1]
	libstg%<>%setNames(c(ribobams,totbams)%>%dirname%>%basename)
	allvoom<-readRDS('data/allvoom.rds')
	libsizes = allvoom@.Data[[1]]%>%{setNames(.$lib.size,rownames(.))}
	rbamcov = allvoom@.Data[[1]]$lib.size
	offsets%<>%mutate(comp='nucl')

	cdsgrl = cdsgrl[str_order_grl(cdsgrl)]
	proteinsequences = GenomicFeatures::extractTranscriptSeqs(cdsgrl,x=fafileob)
	aaproteinsequences = translate(proteinsequences)
	trlist = ids_nrgname%>%distinct(gene_name,transcript_id)%>%{split(.[[2]],.[[1]])}
	iso_tx_countdata <- readRDS(here('data/iso_tx_countdata.rds'))

}

igene='Flna'

COMPRESS=T
# for(COMPRESS in c(TRUE,FALSE)){
for(COMPRESS in c(TRUE)){
for(igene in c('Satb2','Flna','Nes','Bcl11b','Tle4'))
# for(igene in c('Kat6a','Arhgef12'))
# for(igene in c('Slc31a2','Dlg2','Dock9','Msra','Sbno2'))
# for(igene in c('Tuba1a','Tuba1b',))
# for(igene in c('Satb2','Bcl11b'))
{
{

types = c('UTR','CDS')
if(!COMPRESS) types = c(types,'gene')


igeneanno <- fread(str_interp('grep -ie ${igene} ${gtf}'))%>%
	as.matrix%>%apply(1,paste0,collapse='\t')%>%paste0(collapse='\n')%>%import(text=.,format='gtf')%>%
	subset(gene_name==igene)%>%
	subset(type%in%types)
stopifnot(length(igeneanno)>0)
if(COMPRESS){
	igeneanno$toshift <- get_comp_shift(igeneanno,maxwidth=10)
}else{
	igeneanno$toshift <- 0
}
igeneanno$transcript_id%<>%trimids


exontrack =
  igeneanno%>%
  shift(.,-.$toshift)%>%
  # {.$transcript_id %<>% str_replace('\\.[0-9]+$','');.}%>%
  subset(type%in%setdiff(types,'gene'))%>%
  .[,c('type','transcript_id')]%>%
  {.$feature<-as.character(.$type);.}%>%
  {.$transcript<-as.character(.$transcript_id);.}%>%
  Gviz::GeneRegionTrack(.,thinBoxFeature=c("UTR"))
  
library(Gviz)
}


get_cov_track <- function(igeneanno,bams=ribobams,offsets_=offsets){
	st = as.character(strand(igeneanno[1]))
	ichr = as.character(seqnames(igeneanno))[[1]]
	riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=200,which=igeneanno)
  	psites <- lapply(bams,function(bam)readGAlignments(bam,param=riboparam))

	psites%<>%lapply(coverage)
	normcov = map2(psites,libsizes[names(psites)],~ (.x / .y)*1e6)
	tpnormcov = normcov%>%{split(.,libstg[names(.)])} %>%lapply(function(trackpair) purrr::reduce(.f=`+`,trackpair))
	tpnormcov = tpnormcov%>%map(GRanges)%>%GRangesList%>%unlist%>%subset(score!=0)
	bp1tpnormcov=tpnormcov@ranges%>%as('IntegerList')%>%setNames(.,seq_along(.))%>%stack%>%as.data.frame
	bp1tpnormcov$stage = names(tpnormcov)[as.integer(bp1tpnormcov$name)]
	bp1tpnormcov$score = tpnormcov$score[as.integer(bp1tpnormcov$name)]
	bp1tpnormcov = bp1tpnormcov%>%select(start=value,stage,score)%>%
		spread(stage,score)%>%
		mutate(seqnames=ichr,end=start)%>%
		GRanges(.)
	

	ov = findOverlaps(bp1tpnormcov,igeneanno)
	ov = ov[!duplicated(queryHits(ov)),]
	bp1tpnormcov = GenomicRanges::shift(bp1tpnormcov[ov@from],-igeneanno$toshift[ov@to])
	strand(bp1tpnormcov) = st
	# bp1tpnormcov
	stopifnot(all(tps%in%colnames(mcols(bp1tpnormcov))))
	# for(tp in tps) mcols(bp1tpnormcov)[[tp]]%<>%{log2(.+(.5*min(.,na.rm=T)))}
	DataTrack(bp1tpnormcov,type='hist',col = tpcols,groups=tps)
}

 
peptidemsdata_=peptidemsdata;cdsgrl_=cdsgrl;aaproteinsequences_ = aaproteinsequences
get_gene_pepgr = function(igene,peptidemsdata_=peptidemsdata,cdsgrl_=cdsgrl,aaproteinsequences_ = aaproteinsequences){
	gpep_df = peptidemsdata_%>%filter(gene_name==igene)
	if(nrow(gpep_df)==0) return(NULL)
	peptides=gpep_df$Sequence
	ig_seqs = aaproteinsequences_[trlist[[igene]]]
	#
	peptidegr = peptides%>%setNames(.,.)%>%map_df(.id='peptide',function(peptide){
		vmatchPattern(peptide,ig_seqs)%>%
		as.data.frame%>%
		transmute(start,end,width,seqnames=names(ig_seqs)[group])
	})%>%GRanges

	gpeptidegr = peptidegr%>%{
		x=.
		x = resize(x,1,'start',ignore.strand=TRUE)
		end(x)  = ((start(x)-1)*3)+1
		start(x)  = end(x)
		# end(x)  = start(x)
		out = mapFromTranscripts(x,cdsgrl_)
		out$transcript_id = names(cdsgrl_)[out$transcriptsHits]
		out$peptide = x$peptide[out$xHits]
		out
		}
	# gpep	
	gpeptidegr$xHits=NULL
	gpeptidegr$transcriptsHits=NULL
	mcols(gpeptidegr)%<>%cbind(normpepmat[match(gpeptidegr$peptide,rownames(normpepmat)),])
	gpeptidegr
}
# get_gene_pepgr('Satb2')

get_peptide_track<-function(igeneanno){
	igeneanno%<>%subset(type!='gene')
	igene = igeneanno$gene_name[1]
	pepgr = get_gene_pepgr(igene)
	if(is.null(pepgr))return(DataTrack(NULL))
	pepgr = pepgr[!duplicated(pepgr$peptide),]
	# igeneanno = annotation%>%subset(gene_name=igene)
	ov = findOverlaps(pepgr,igeneanno)
	ov = ov[pepgr$transcript_id[queryHits(ov)]== igeneanno$transcript_id[subjectHits(ov)],]
	ov = ov[!duplicated(queryHits(ov)),]
	pepgr = GenomicRanges::shift(pepgr[ov@from],-igeneanno$toshift[ov@to])
	# start(pepgr)%>%txtplot
	pepgr$transcript_id=NULL
	pepgr = unique(pepgr)
	pepgr = sort(pepgr)
	# pepgr[60:70]
	# pepgr = pepgr[70:80]
	DataTrack(pepgr[,tps],cols = tpcols,groups=tps)	
	# pepgr
}

motifpattern='TGTANATA'

getmotiftrack = function(motifpattern,igeneanno){
	igeneanno%<>%subset(type!='gene')
	spligeneanno = igeneanno%>%split(.,.$transcript_id)
	spligeneanno%<>%sort_grl_st
	trs = names(spligeneanno)
	igeneannoseq = spligeneanno%>%extrseq

	motifgr = vmatchPattern(motifpattern,fixed=F,igeneannoseq)%>%
		as.data.frame%>%
		transmute(start,end,width,seqnames=names(igeneannoseq)[group])
	if(nrow(motifgr)==0) return(AnnotationTrack(NULL))
	motifgr = motifgr%>% GRanges%>%
		mapFromTranscripts(spligeneanno)
		
	ov = findOverlaps(motifgr,igeneanno)
	ov = ov[trs[motifgr$transcriptsHits][queryHits(ov)]== igeneanno$transcript_id[subjectHits(ov)],]
	ov = ov[!duplicated(queryHits(ov)),]
	motifgr = GenomicRanges::shift(motifgr[ov@from],-igeneanno$toshift[ov@to])
	motifgr = unique(motifgr)
	AnnotationTrack(motifgr)
}


peaktrack = function(gr,igeneanno){		
	ov = findOverlaps(gr,igeneanno,select='first')
	if(is.na(ov)) return(AnnotationTrack(NULL))
	gr = GenomicRanges::shift(gr,-igeneanno$toshift[ov])
	gr = unique(gr)
	AnnotationTrack(gr)
}
getstart_track = function(igeneanno){
	spligeneanno = igeneanno%>%subset(type=='CDS')%>%split(.,.$transcript_id)%>%.[str_order_grl(.)]
	spligeneanno = spligeneanno%>%resize_grl(1,'start')%>%unlist
	spligeneanno %<>% GenomicRanges::shift(-spligeneanno$toshift)
	AnnotationTrack(spligeneanno,background.panel = "grey")
}

tps = names(tpcols)

zhangetal_pum2clip = GRanges('chr1:56794013-56794101:-')
pum1clip = GRanges('chr1:56796631-56796712')
dir.create('plots/Figure7')

nametrack = function(x,tnm) x%>%{.@name=tnm;.}

gviztracks <- list(
	get_cov_track(igeneanno,totbams)%>%nametrack('RNA-Seq (RPM)'),
	get_cov_track(igeneanno,ribobams)%>%nametrack('Ribo-Seq (RPM)'),
	get_peptide_track(igeneanno)%>%nametrack('log2(Normalized Intensity)'),
 	exontrack%>%nametrack('Transcripts'),
 	getstart_track(igeneanno)%>%nametrack('AUG'),
 	getmotiftrack('TGTANATA',igeneanno)%>%nametrack('Pum2 Motifs'),
 	peaktrack(zhangetal_pum2clip,igeneanno)%>%nametrack('Zhang et al Clip'),
 	peaktrack(pum1clip,igeneanno)%>%nametrack('PUM1 Clip')
)

dir.create('data/Shiny_track_data/')
if(COMPRESS) saveRDS(gviztracks,str_interp('data/Shiny_track_data/${igene}.rds'))

{
#now plot
options(ucscChromosomeNames=FALSE)
compressstr = if(COMPRESS) '' else '_NONcompressed_'
plotfile<- here(paste0('plots/Figure7/locusplot_',compressstr,igene,'.pdf'))
pdf(plotfile,w=24,h=12)
plotTracks(gviztracks,col = tpcols)
dev.off()
message(normalizePath(plotfile))
}
}
}

gviztracklist = Sys.glob('data/Shiny_track_data/*.rds')%>%map(readRDS)


