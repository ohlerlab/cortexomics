################################################################################
########Apply the model
################################################################################


# suppressMessages({library(svglite)})
suppressMessages({library(readr)})
suppressMessages({library(Biostrings)})
suppressMessages({library(Rsamtools)})
#suppressMessages({library(psd)})
suppressMessages({library(txtplot)})
suppressMessages({library(rtracklayer)})
suppressMessages({library(stringr)})
suppressMessages({library(data.table)})
suppressMessages({library(assertthat)})
suppressMessages({library(parallel)})
suppressMessages({library(dplyr)})
suppressMessages({library(purrr)})
suppressMessages({library(here)})
suppressMessages({library(magrittr)})
suppressMessages({library(stringr)})
suppressMessages({library(tidyverse)})
suppressMessages({library(GenomicFiles)})
# suppressMessages({library(bamsignals)})
suppressMessages({library(zeallot)})
suppressMessages({library(stringr)})
suppressMessages({library(here)})
suppressMessages({library(assertthat)})
suppressMessages({library(multitaper)})
suppressMessages({library(GenomicAlignments)})
suppressMessages({library(GenomicFeatures)})
suppressMessages({library(here)})

library(parallel)

reduce <- GenomicRanges::reduce

MAPQTHRESH <- 200
USEPHASE <- FALSE
USERIBOSEQC <- FALSE
source(here('src/R/Rprofile.R'))



for(fname in lsf.str('package:GenomicRanges')) assign(fname,get(fname,'package:GenomicRanges'))
for(fname in lsf.str('package:dplyr')) assign(fname,get(fname,'package:dplyr'))


argv <- c(
	bam = here('pipeline/star/data/P0_ribo_1/P0_ribo_1.bam')%T>%{stopifnot(file.exists(.))},
	gtf = here('pipeline/my_gencode.vM12.annotation.gtf'),
	REF = here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'),
	shiftmodel = 'pipeline/seqshift_reads/data/P0_ribo_1/seqshiftmodel.rds',
	outfolder = 'riboseq_quant/data/P0_ribo_1/'
)

sampleparams <- fread(here('pipeline/sample_parameter.csv'))

argv[] <- commandArgs(trailing=TRUE)

for (nm in names(argv)) assign(nm,argv[[nm]])


message(outfolder)
#get exons
 if(!exists('gtf_gr')) gtf_gr<-import(con=gtf,format='gtf')
# exons <- gtf_gr%>%subset(type=='exon')

# if(!is('exons','GRanges')) exons <- gtf_gr%>%subset(type=='exon')
# if(!is('cds','GRanges')) cds <- gtf_gr%>%subset(type=='CDS')
# if(!exists('startcodsa')) startcods <- gtf_gr%>%subset(type=='start_codon')
exons <- gtf_gr%>%subset(type=='exon')
cds <- gtf_gr%>%subset(type=='CDS')
	

if(!file.exists(shiftmodel)){
	message('no psite model found')
	psite_model <- NULL
}else{
	psite_model <- readRDS(shiftmodel)
}


# exons2use[2:3]%>%revElements(.,any(strand(.)=='-'))

# cds	%>%split(.,.$protein_id)%>%head(2)%>%lapply(head,2)%>%GRangesList%>%unlist%>%{.[isfpmost(.)]%<>%clip_start(5);.}


#NOte - startflank is in basepairs, STARTCLIP in codons

# exons2use[2:3]%>%revElements(.,any(strand(.)=='-'))


# tmpgrl <- exons2use[2:3]
# any(strand(exons2use)=='+')%>%which%>%head(1)


# grl<-tmpgrl
# tmpgrl<-tmpgrl%>%unlist%>%.[,NULL]%>%split(.,names(.))
resize_grl<-function(grl,nbp,fixend='end'){
	stopifnot(identical(grl,sort(grl)))
	if(is.null(names(grl))) names(grl)<-seq_along(grl)
	#if we fix the start we want to select the end
	if(fixend=='start'){
		strandtorev <- '-' 
		endinds <-  grl%>%{as(ifelse(any(strand(.)%in%strandtorev),1,elementNROWS(.)),'IntegerList')}

	}else if (fixend=='end'){
		#and if the end then we want the start
		strandtorev <- '+'
		endinds <-  grl%>%{as(ifelse(any(strand(.)%in%strandtorev),1,elementNROWS(.)),'IntegerList')}

	}else{stop('fixend needs to be start or end')}
	cuminds <- cumsum(elementNROWS(grl))

	ulgrl <- unlist(grl,use.names=TRUE)
	endinds <- unlist(lag(cuminds,1,0)+endinds)
	#but this isn't safe since the gene might end up out of bounds
	#add the difference between it and the limit, rounded down to the nearest 3, up to 15
	stopifnot(all(start(ulgrl[endinds]) > nbp))
	ulgrl[endinds] %<>% resize(width(.)+nbp,fixend)
	grl <- split(setNames(ulgrl,NULL),names(ulgrl))
	stopifnot(is(grl,'GRangesList'))
	grl
}

strandedorder <- function(grl) {
	order <- order(tmpgrl)
	order[any(strand(tmpgrl)%in%'-')]%<>%revElements
	order
}


################################################################################
library(Biobase)
library(data.table)

countfile='pipeline/exprdata/countexprset.rds'
countexprdata <- readRDS(countfile)
count_tbl <- assayData(countexprdata)[['exprs']][countexprdata@featureData$is_gid_highest,]%>%
  as.data.table%>%mutate(feature_id=fData(countexprdata)$gene_name[fData(countexprdata)$is_gid_highest])
highexprgeneinds <- assayData(countexprdata)[['exprs']][countexprdata@featureData$is_gid_highest,]%>%rowMedians%>%order%>%tail(30e3)
highprotein_ids<-countexprdata@featureData$protein_id[countexprdata@featureData$is_gid_highest][highexprgeneinds]



cdstouse <- cds%>%subset(protein_id%in%highprotein_ids)

firstcds <- cdstouse%>%subset(strand=='+')%>%split(.,.$protein_id)%>%.[as(rep(1,length(.)),'IntegerList')]%>%unlist
lastcds <- cdstouse%>%subset(strand=='+')%>%split(.,.$protein_id)%>%revElements%>%.[as(rep(1,length(.)),'IntegerList')]%>%unlist

startwindows <- firstcds%>%resize(50,'start')%>%resize(50+20,'end')
endwindows <- lastcds%>%resize(50,'end')%>%resize(50+20,'start')
endwindows%>%length

bams <- here('pipeline/star/data/*/*.bam')%>%Sys.glob%>%str_subset(neg=TRUE,'transcript')%>%str_subset('_ribo_|total')
ribobams<-bams%>%str_subset('ribo')

allpsitedfs <- mclapply(mc.cores=20,bams,function(bam){

	ovcounts <- summarizeOverlaps(c(startwindows,endwindows),BamFile(bam))

	highreads <- assay(ovcounts)[1:length(startwindows)] %>%as.vector%>%`>`(10)

	highwinds <- c(startwindows,endwindows)[rep(highreads,2)]

	riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=MAPQTHRESH,which=highwinds)
	reads <- readGAlignments(bam,param=riboparam)
	mcols(reads)$cdsshift <- get_cds_offsets(reads,psite_model$offsets,psite_model$compartments)

	#reads <- reads%>%subset(width==27)
	reads <- reads%>%subset(width %in% psite_model$offsets$length)
	mcols(reads)$length <- width(reads)
	reads%<>%subset(!is.na(cdsshift))
	psites <- apply_psite_offset(reads,c('cdsshift'))%>%as("GRanges")
	mcols(psites)$length <- mcols(reads)$length

	#
	startwindpsites <- psites%>%mapToTranscripts(highwinds%>%head(sum(highreads)))
	startwindpsites$length <- psites$length[startwindpsites$xHits]
	startwindpsites<-startwindpsites%>%as.data.frame%>%group_by(start,length)%>%tally%>%ungroup%>%mutate(section='AUG',start=start-21)
	endwindpsites <- psites%>%mapToTranscripts(highwinds%>%tail(sum(highreads)))
	endwindpsites$length <- psites$length[endwindpsites$xHits]
	endwindpsites<-endwindpsites%>%as.data.frame%>%group_by(start,length)%>%tally%>%ungroup%>%mutate(section='STOP',start=start-51) 
	allpsites<-rbind(startwindpsites,endwindpsites)
	allpsites%<>%mutate(section=factor(section,unique(section)),bam=basename(bam))
	message('done')
	message(bam)
	allpsites
})
allpsitedfs%<>%setNames(bams)
allpsitedfs%<>%.[ribobams]
allpsitedfs_bind <- allpsitedfs%>%bind_rows
allpsitedfs_bind$bam%>%unique
allpsitedfs_bind %<>% mutate(stage=str_extract(bam,'[^_]+'))
stagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
stageconv = names(stagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
allpsitedfs_bind$stage %<>% stageconv[.]

# save.image('data/myribowaltzplots.Rdata')
# load('data/myribowaltzplots.Rdata')



{

bamcols <- allpsitedfs_bind$bam%>%unique%>%str_extract('[^_]+')%>%stageconv[.]%>%stagecols[.]%>%setNames(allpsitedfs_bind$bam%>%unique)
roundup <- function(x,n) n*ceiling(x/n)
rounddown <- function(x,n) n*floor(x/n)
number_ticks <- function(limits,n=3){
	out = c(seq(min(rounddown(limits,n)),- n,n),seq(0,max(roundup(limits,n)),by=n))
	out = out[between(out,limits[1],limits[2])]
	out
}
totalcounts<-allpsitedfs_bind%>%group_by(stage,section,start)%>%tally%>%.$n%>%sum

plotfile<-'plots/figures/figure1/fig1c_myribowaltz_allsec_stageov.pdf'%T>%pdf(h=4,w=8)
rwplot <- qplot(data=allpsitedfs_bind%>%group_by(stage,section,start)%>%tally%>%mutate(n=n/totalcounts),color=stage,x=start,y=n,geom='line')+
	scale_x_continuous(name='position',limits=c(-24,24),minor_breaks=number_ticks,breaks=partial(number_ticks,n=12))+
	# scale_x_continuous(name='position')+
	# scale_color_discrete(guide=TRUE)+
	facet_grid(  ~ section)+
	scale_color_manual(values=stagecols)+
	scale_y_continuous(name='normalized bin count')+
	theme_bw()
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message


plotfile<-'plots/figures/figure1/fig1c_myribowaltz_allsec_stagesep.pdf'%T>%pdf(h=4,w=8)
rwplot <- qplot(data=allpsitedfs_bind%>%group_by(stage,bam,section,start)%>%tally%>%mutate(n=n/totalcounts),group=bam,color=stage,x=start,y=n,geom='line')+
	scale_x_continuous(name='position',limits=c(-24,24),minor_breaks=number_ticks,breaks=partial(number_ticks,n=12))+
	# scale_color_discrete(guide=TRUE)+
	facet_grid( stage ~ section)+
	scale_color_manual(values=stagecols)+
	scale_y_continuous(name='normalized bin count')+
	theme_bw()
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message

allpsitedfs_bind$bam%>%unique%>%is_in(names(bamcols))

plotfile<-'plots/figures/figure1/fig1c_myribowaltz_aug_stageov.pdf'%T>%pdf(h=4,w=8)
qplot(data=allpsitedfs_bind%>%filter(section=='AUG',as.numeric(start) > -10,start<10)%>%group_by(bam,section,start)%>%tally%>%mutate(n=n/totalcounts),color=bam,x=start,y=n,geom='line')+
	scale_x_continuous(name='position',limits=c(-9,9),minor_breaks=number_ticks,breaks=partial(number_ticks,n=3))+
	# scale_color_discrete(guide=TRUE)+
	facet_grid( . ~ section)+
	scale_color_manual(values=bamcols)+
	scale_y_continuous(name='normalized bin count')+
	theme_bw()
dev.off()
normalizePath(plotfile)%>%message

plotfile<-'plots/figures/figure1/fig1c_myribowaltz_stop_stageov.pdf'%T>%pdf(h=4,w=8)
qplot(data=allpsitedfs_bind%>%filter(section=='STOP',as.numeric(start) > -10,start<10)%>%group_by(bam,section,start)%>%tally%>%mutate(n=n/totalcounts),color=bam,x=start,y=n,geom='line')+
	scale_x_continuous(name='position',limits=c(-9,9),minor_breaks=number_ticks,breaks=partial(number_ticks,n=3))+
	# scale_color_discrete(guide=TRUE)+
	facet_grid( . ~ section)+
	scale_color_manual(values=bamcols)+
	scale_y_continuous(name='normalized bin count')+
	theme_bw()
dev.off()
normalizePath(plotfile)%>%message
}

# # cds2use < cds%>%split(.,.$protein_id)%>%sample(50)
# cds2use <- cds%>%split(.,.$protein_id)
# p2trdf<-with(cds,data_frame(protein_id,transcript_id))%>%distinct
# cdstrs<-data_frame(protein_id=names(cds2use))%>%safe_left_join(p2trdf)%>%.$transcript_id
# exons2use <- exons%>%split(.,.$transcript_id)%>%.[cdstrs]
# names(exons2use)<-names(cds2use)


# exons2use <- sort(exons2use)
# exons2use <- exons2use%>%resize_grl(9,'end')
# exons2use <- exons2use%>%resize_grl(9,'start')



# cdstousetrna<-cds%>%split(.,.$protein_id)%>%.[highprotein_ids]%>%unlist

# codons <- getSeq(cdstousetrna,x=FaFile(REF))%>%split(cdstousetrna$protein_id)%>%lapply(.%>%unlist%>%xscat)%>%DNAStringSet%>%vmatchPattern('TCT',.)%>%unlist%>%subset((start %% 3) == 1)

# codonsgr <- GRanges(names(codons),codons)%>%mapFromTranscripts(cdstousetrna%>%split(.,.$protein_id))

# codonsgr_exp <- codonsgr%>%resize(3+9+9,'center')


# signalovercodons<-ribobams%>%mclapply(function(bam){
# 	riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=MAPQTHRESH,which=codonsgr_exp)
# 	reads <- readGAlignments(bam,param=riboparam)
# 	mcols(reads)$cdsshift <- get_cds_offsets(reads,psite_model$offsets,psite_model$compartments)
# 	#reads <- reads%>%subset(width==27)
# 	reads <- reads%>%subset(width %in% psite_model$offsets$length)
# 	mcols(reads)$length <- width(reads)
# 	reads%<>%subset(!is.na(cdsshift))
# 	psites <- apply_psite_offset(reads,c('cdsshift'))%>%as("GRanges")
# 	mcols(psites)$length <- mcols(reads)$length
# 	startwindpsites <- psites%>%mapToTranscripts(codonsgr_exp)
# 	start(startwindpsites)
# })


################################################################################
########Check if the different stop codon classes look different
################################################################################
	
stopcodontypes <- lastcds%>%resize(3,'end')%>%strandshift(3)%>%getSeq(FaFile(REF),.)%>%as.character%>%setNames(lastcds$protein_id)
endwindowtypes<-endwindows%>%split(.,stopcodontypes[.$protein_id])%>%.[c('TAA','TAG','TGA')]
endwindowtypes%>%lapply(length)
stoptypes<-names(endwindowtypes)

typesepallpsitedfs <- mclapply(mc.cores=20,ribobams[],function(bam){
	stoptype=names(endwindowtypes)[1]
	lapply(names(endwindowtypes),function(stoptype){
		endwindowsi<-endwindowtypes[[stoptype]]

		ovcounts <- summarizeOverlaps(endwindowsi,BamFile(bam))

		highreads <- assay(ovcounts)[1:length(endwindowsi)] %>%as.vector%>%`>`(10)

		highwinds <- c(endwindowsi)[highreads]

		riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=MAPQTHRESH,which=highwinds)
		reads <- readGAlignments(bam,param=riboparam)
		mcols(reads)$cdsshift <- get_cds_offsets(reads,psite_model$offsets,psite_model$compartments)

		#reads <- reads%>%subset(width==27)
		reads <- reads%>%subset(width %in% psite_model$offsets$length)
		mcols(reads)$length <- width(reads)
		reads%<>%subset(!is.na(cdsshift))
		psites <- apply_psite_offset(reads,c('cdsshift'))%>%as("GRanges")
		mcols(psites)$length <- mcols(reads)$length
		#
		endwindpsites <- psites%>%mapToTranscripts(highwinds%>%tail(sum(highreads)))
		endwindpsites$length <- psites$length[endwindpsites$xHits]
		endwindpsites<-endwindpsites%>%as.data.frame%>%group_by(start,length)%>%tally%>%ungroup%>%mutate(section='STOP',start=start-51) 
		allpsites<-rbind(endwindpsites)
		allpsites%<>%mutate(section=factor(section,unique(section)),bam=basename(bam))
		message('done')
		message(bam)
		allpsites
	})	
})
typesepallpsitedfs%<>%setNames(ribobams[])
# typesepallpsitedfs%<>%.[ribobams]
typesepallpsitedfs_bind <- typesepallpsitedfs%>%lapply(bind_rows,.id='stoptype')%>%bind_rows
typesepallpsitedfs_bind %<>% mutate(stage=str_extract(bam,'[^_]+'))
if(all(typesepallpsitedfs_bind$stoptype%in%1:100)) typesepallpsitedfs_bind$stoptype = stoptypes[as.numeric(typesepallpsitedfs_bind$stoptype)]
typesepallpsitedfs_bind$section='STOP'
totalcounts<-typesepallpsitedfs_bind%>%group_by(stage,section,start)%>%tally%>%.$n%>%sum
typesepallpsitedfs_bind$bam%<>%str_replace('.bam$','')

plotfile<-'plots/figures/figure2/fig2_myribowaltz_stop_stageov_stoptypesep.pdf'%T>%pdf(h=4,w=8)
qplot(data=typesepallpsitedfs_bind%>%filter(section=='STOP',as.numeric(start) > -10,start<10)%>%group_by(stoptype,bam,section,start)%>%tally%>%mutate(n=n/totalcounts),color=bam,x=start,y=n,geom='line')+
	scale_x_continuous(name='position',limits=c(-9,9),minor_breaks=number_ticks,breaks=partial(number_ticks,n=3))+
	# scale_color_discrete(guide=TRUE)+
	facet_grid( stoptype ~ section)+
	scale_color_manual(values=bamcols)+
	scale_y_continuous(name='normalized bin count')+
	theme_bw()
dev.off()
normalizePath(plotfile)%>%message

