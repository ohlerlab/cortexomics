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
suppressMessages({library(ggpubr)})
library(parallel)

stagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
stageconv = names(stagecols)%>%setNames(c('E13','E145','E16','E175','P0'))

reduce <- GenomicRanges::reduce
#
MAPQTHRESH <- 200
USEPHASE <- FALSE
USERIBOSEQC <- FALSE
source(here('src/R/Rprofile.R'))



roundup <- function(x,n) n*ceiling(x/n)
rounddown <- function(x,n) n*floor(x/n)
number_ticks <- function(limits,n=3){
	out = c(seq(min(rounddown(limits,n)),- n,n),seq(0,max(roundup(limits,n)),by=n))
	out = out[between(out,limits[1],limits[2])]
	out
}


displaystagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
stageconv = names(stagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
stagecols = displaystagecols%>%setNames(names(stageconv))

for(fname in lsf.str('package:GenomicRanges')) assign(fname,get(fname,'package:GenomicRanges'))
for(fname in lsf.str('package:dplyr')) assign(fname,get(fname,'package:dplyr'))


argv <- c(
	bam = here('pipeline/star/data/P0_ribo_1/P0_ribo_1.bam')%T>%{stopifnot(file.exists(.))},
	gtf = here('pipeline/my_gencode.vM12.annotation.gtf'),
	REF = here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'),
	shiftmodel = 'pipeline/seqshift_reads/data/P0_ribo_1/seqshiftmodel.rds',
	outfolder = 'riboseq_quant/data/P0_ribo_1/'
)
#
sampleparams <- fread(here('pipeline/sample_parameter.csv'))
#
argv[] <- commandArgs(trailing=TRUE)
#
for (nm in names(argv)) assign(nm,argv[[nm]])
if(!file.exists(shiftmodel)){
	message('no psite model found')
	psite_model <- NULL
}else{
	psite_model <- readRDS(shiftmodel)
}




message(outfolder)
#get exons
 if(!exists('gtf_gr')) gtf_gr<-import(con=gtf,format='gtf')
# exons <- gtf_gr%>%subset(type=='exon')

# if(!is('exons','GRanges')) exons <- gtf_gr%>%subset(type=='exon')
# if(!is('cds','GRanges')) cds <- gtf_gr%>%subset(type=='CDS')
# if(!exists('startcodsa')) startcods <- gtf_gr%>%subset(type=='start_codon')
exons <- gtf_gr%>%subset(type=='exon')
cds <- gtf_gr%>%subset(type=='CDS')
	

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


cdsids2use<-
cds2use <- cds%>%subset(protein_id%in%cdsids2use)%>%split(.,.$protein_id)
firstcds <- cds2use%>%.[as(rep(1,length(.)),'IntegerList')]%>%unlist
lastcds <- cds2use%>%revElements%>%.[as(rep(1,length(.)),'IntegerList')]%>%unlist
#

startwindows <- firstcds%>%resize(100,'start')%>%resize(100+30,'end')
endwindows <- lastcds%>%resize(100,'end')%>%resize(100+30,'start')
#
bams <- here('pipeline/star/data/*/*.bam')%>%Sys.glob%>%str_subset(neg=TRUE,'transcript')%>%str_subset('_ribo_|total')
ribobams<-bams%>%str_subset('ribo')%>%setNames(.,.)

stopifnot(startwindows%>%strand%>%`==`('-')%>%any)
# #
# allpsitedfs <- mclapply(mc.cores=1,ribobams[1],function(bam){
# 	#
# 	ovcounts <- summarizeOverlaps(cds2use,BamFile(bam))
# 	#
# 	assay(ovcounts)%>%dim
# 	length(cds2use)
# 	highreads <- assay(ovcounts)[1:length(cds2use)] %>%as.vector%>%`>`(32)
# 	ovcounts <- assay(ovcounts)[,1]%>%setNames(cds2use)
# 	highreads%>%is.na%>%table
# 	c(startwindows,endwindows)%>%length
# 	#
# 	highwinds <- c(startwindows,endwindows)[rep(highreads,2)]
# 	#now scan the start and stop
# 	riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=MAPQTHRESH,which=highwinds)
# 	reads <- readGAlignments(bam,param=riboparam)
# 	mcols(reads)$cdsshift <- get_cds_offsets(reads,psite_model$offsets,psite_model$compartments)
# 	#reads <- reads%>%subset(width==27)
# 	reads <- reads%>%subset(width %in% psite_model$offsets$length)
# 	mcols(reads)$length <- width(reads)
# 	reads%<>%subset(!is.na(cdsshift))
# 	psites <- apply_psite_offset(reads,c('cdsshift'))%>%as("GRanges")
# 	mcols(psites)$length <- mcols(reads)$length
# 	#
# 	# highwinds%>%head(sum(highreads))%>%length
# 	startwindpsites <- psites%>%mapToTranscripts(highwinds%>%head(sum(highreads)))	
# 	startwindpsites<-startwindpsites%>%coverage%>%.[sum(.)>0]%>%{ovcounts}%>%Reduce(f='+')
# 	startwindpsites<-startwindpsites%>%as.data.frame%>%rownames_to_column('start')%>%mutate(section='AUG',start=as.numeric(start)-21)
# 	endwindpsites <- psites%>%mapToTranscripts(highwinds%>%tail(sum(highreads)))
# 	endwindpsites<-endwindpsites%>%coverage%>%.[sum(.)>0]%>%{./sum(.)}%>%Reduce(f='+')
# 	endwindpsites<-endwindpsites%>%as.data.frame%>%rownames_to_column('start')%>%mutate(section='STOP',start=as.numeric(start)-51)

# 	allpsites<-rbind(startwindpsites,endwindpsites)

# 	allpsites%<>%mutate(section=factor(section,unique(section)),bam=basename(bam))
# 	message('done')
# 	message(length(cds2use))
# 	allpsites
# })

#25 to 10
# offsets = fread(header=TRUE"
# offset compartment length
# 	11 nucl            26
# 	11 nucl            27
# 	12 nucl            28
# 	11 nucl            29
# 	10 nucl            25
# 	12 nucl            30
# 	12 nucl            31
# 	"
# )%T>%write_tsv(here('myoffsets.tsv'))
offsets = fread('ext_data/offsets_manual.tsv')
#Now run through the bam file and get a profile for each codon

bams2scan<-bams%>%
	# .[c(1:2,length(.)-1,length(.))] %>%
	identity

widths <- (25:31)%>%setNames(.,paste0('rl',.))
getreadfps <- function(bam,mapqthresh=200){
	  riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=cds2use%>%unlist)
	  reads <- readGAlignments(bam,param=riboparam)
	  reads%<>%as('GRanges')
	  reads$length <- width(reads)
	  reads%<>%subset(length %in% widths)
	  reads %<>% resize(1,'start')
	  reads
}
fpsitelist <- mclapply(bams2scan,mymemoise(getreadfps))
fpsitelist%<>%setNames(bams2scan%>%basename%>%str_replace('\\.\\w+$',''))
stopifnot('length'%in%colnames(mcols(fpsitelist[[1]])))

cds2useext <- ext_grl(sort(cds2use),45,fixend='start')%>%ext_grl(45,fixend='end')

{
	reads=fpsitelist[[1]]
startstopprofiles<-mclapply(mc.cores=1,fpsitelist,(function(reads){
	# reads%<>%resize(.,.$length)
	#
	mcols(reads)$cdsshift <- offsets$offset[match(reads$length,offsets$length)]
	# reads <- reads%>%subset(width==27)
	# reads <- reads%>%subset(width %in% psite_model$offsets$length)
	# mcols(reads)$length <- width(reads)
	reads%<>%subset(!is.na(cdsshift))
	psites <- strandshift(reads,reads$cdsshift)
	#
	cdscounts <- countOverlaps(cds2useext,psites)
	zerocds <- names(cdscounts)[cdscounts < 32]
	#
	#
	cdscounts_densitys <- cdscounts%>%divide_by(sum(width(cds2useext)))%>%setNames(names(cds2useext))
	#
	stopifnot(all(names(cdscounts_densitys)==names(cds2useext)))
	#
	startwindows->sectionwind
	mclapply(mc.cores=1,list('AUG'=startwindows,'STOP'=endwindows),function(sectionwind){
		#
		message(length(sectionwind))
		#
		length_i=widths[4]
		lapply(widths%>%.[.%in%reads$length],function(length_i){
			#
			message(length_i)
			pospsites<-psites%>%subset(length==length_i)%>%subset(strand=='+')
			negpsites<-psites%>%subset(length==length_i)%>%subset(strand=='-')

			#
			poswinds<-sectionwind%>%subset(strand=='+')#%>%resize(.,width(.-((FLANKCODS-3)*3)),'center')
			negwinds<-sectionwind%>%subset(strand=='-')#%>%resize(.,width(.-((FLANKCODS-3)*3)),'center')
			#
			profiles <- 
			c(
				pospsites%>%coverage%>%.[poswinds],
				negpsites%>%coverage%>%.[negwinds],
				NULL
			)
			#
			nposwinds<-poswinds%>%length
			profiles[-(1:nposwinds)]%<>%revElements
			# profilemat<-profiles%>%simplify2array%>%t
			windnames<-c(names(poswinds),names(negwinds))
			# sapply(i:n_col,function(i) profiles[as(rep(i,np),'IntegerList')]%>%sum%>%add(1)%>%log)
			#now rather than covert to a matrix, we can do this to 
			profiles <- profiles[!windnames%in%zerocds]
			normdensities <- cdscounts_densitys[windnames[!windnames%in%zerocds]]
			np<-length(profiles)
			if(np==0){return(data.frame(position=1:width(poswinds[1]),signal=0))}
			n_col<-length(profiles[[1]])
			logsumpvect<-sapply(1:n_col,function(i){
				# i=1
				profiles[as(rep(i,np),'IntegerList')]%>%sum%>%divide_by(normdensities)%>%mean(na.rm=TRUE)
				# profiles[as(rep(i,np),'IntegerList')]%>%sum%>%sum
			})
			message('.')
			# if(any(!is.finite(logsumpvect)))
			# invisible(logsumpvect)
			enframe(logsumpvect,'position','signal')
		})%>%bind_rows(.id='readlen')
	})%>%bind_rows(.id='section')
}))%>%bind_rows(.id='sample')

startstopprofiles[[2]]
startstopprofiles%<>%mutate(start = as.numeric(position) - ifelse(section=='AUG',31,101))
}

#
{
allreadlens2plot<-c('rl29','rl27','rl26','rl28','rl25','rl30','rl31')
allpsitedfs_bind <- startstopprofiles
allpsitedfs_bind$stage %<>% stageconv[.]
allpsitedfs_bind%<>%filter(readlen%in%allreadlens2plot)
# allpsitedfs_bind%<>%mutate(start = start + ifelse(readlen=='rl28',3,0))
grpreadlens=allreadlens2plot%>%setdiff(c(''))
allpsitedfs_bind%<>%mutate(readlen2grp = ifelse(readlen %in%grpreadlens,'all',readlen))
allpsitedfs_bind%<>%select(-readlen)
# allpsitedfs_bind%<>%group_by(sample,section,start)%>%mutate(signal=sum(signal))
allpsitedfs_bind%<>%group_by(sample,section,start,readlen2grp)%>%mutate(signal=sum(signal))
# allpsitedfs_bind%<>%group_by(sample,section,start,readlen2grp)%>%mutate(signal=sum(signal))
allpsitedfs_bind%<>%group_by(sample,section,readlen2grp)%>%mutate(signal=signal / median(signal))
allpsitedfs_bind %<>% ungroup%>%mutate(stage=str_extract(sample,'[^_]+'))
#
allpsitedfs_bind_stgrp <- allpsitedfs_bind %>% group_by(stage,section,readlen2grp,start)%>%summarise(signal=mean(signal)) 
samplecols <- allpsitedfs_bind$sample%>%unique%>%str_extract('[^_]+')%>%stageconv[.]%>%stagecols[.]%>%setNames(allpsitedfs_bind$sample%>%unique)
}
{
	library(rlang)
	# allpsitedfs_bind%>%slice_by(sample,1)
	plotfile<-'plots/figures/figure1/rltest.pdf'%T>%pdf(h=6,w=12)
	rwplot <- allpsitedfs_bind_stgrp%>%
		split(.,.$section)%>%map( .%>%
			# slice_by(sample,8)%>%
			# filter(position%%3 == 0)%>%
			{
				qplot(data=.,color=stage,x=start,y=signal,group=readlen2grp,geom='blank')+
				geom_line(aes(linetype=readlen2grp))+
				scale_x_continuous(name='position',
					# limits=if(.$section[1]=='AUG') c(-30,100) else c(-100,30) ,
					limits=if(.$section[1]=='AUG') c(-30,100) else c(-100,30) ,
					minor_breaks=number_ticks,breaks=partial(number_ticks,n=12))+
				# scale_x_continuous(name='position',minor_breaks=number_ticks,breaks=partial(number_ticks,n=12))+
				# scale_x_continuous(name='position')+
				# scale_color_discrete(guide=TRUE)+
				# facet_grid( sample+readlen ~ section)+
				facet_grid( readlen2grp+stage ~ section,scale='free_x')+
				scale_color_manual(values=stagecols)+
				scale_y_continuous(name='normalized bin count')+
				theme_bw()
			}
		)%>%
		{.[[1]]$show.legend=F;.}%>%
		ggarrange(plotlist=.,ncol=2)
	print(rwplot)
	dev.off()
	normalizePath(plotfile)%>%message
}

{
	library(rlang)
	# allpsitedfs_bind%>%slice_by(sample,1)
	plotfile<-'plots/figures/figure1/fig1c_myribowaltz_allsec_stageov.pdf'%T>%pdf(h=6,w=12)
	rwplot <- allpsitedfs_bind_stgrp%>%
		split(.,.$section)%>%map( .%>%
			# slice_by(sample,8)%>%
			# filter(position%%3 == 0)%>%
			{
				isfirst = .$section[1]=='AUG'
				qplot(data=.,color=stage,x=start,y=signal,geom='blank')+
				geom_line()+
				scale_x_continuous(name='position',
					# limits=if(.$section[1]=='AUG') c(-30,100) else c(-100,30) ,
					limits=if(isfirst) c(-30,100) else c(-100,30) ,
					minor_breaks=number_ticks,breaks=partial(number_ticks,n=12))+
				# scale_x_continuous(name='position',minor_breaks=number_ticks,breaks=partial(number_ticks,n=12))+
				# scale_x_continuous(name='position')+
				# scale_color_discrete(guide=TRUE)+
				# facet_grid( sample+readlen ~ section)+
				facet_grid( ~ section,scale='free_x')+
				scale_color_manual(values=stagecols)+
				scale_y_continuous(name='Mean Psite Count / CDS Total',limits=c(0,20))+
				theme_bw()
			}
		)%>%
		ggarrange(plotlist=.,ncol=2,common.legend=T)
	print(rwplot)
	dev.off()
	normalizePath(plotfile)%>%message
	library(rlang)
	# allpsitedfs_bind%>%slice_by(sample,1)
	plotfile<-'plots/figures/figure1/fig1c_myribowaltz_allsec_stagesep.pdf'%T>%pdf(h=6,w=12)
	rwplot <- allpsitedfs_bind_stgrp%>%
		split(.,.$section)%>%map( .%>%
			# slice_by(sample,8)%>%
			# filter(position%%3 == 0)%>%
			{
				isfirst = .$section[1]=='AUG'
				qplot(data=.,color=stage,x=start,y=signal,geom='blank')+
				geom_line()+
				scale_x_continuous(name='position',
					# limits=if(.$section[1]=='AUG') c(-30,100) else c(-100,30) ,
					limits=if(isfirst) c(-30,100) else c(-100,30) ,
					minor_breaks=number_ticks,breaks=partial(number_ticks,n=12))+
				# scale_x_continuous(name='position',minor_breaks=number_ticks,breaks=partial(number_ticks,n=12))+
				# scale_x_continuous(name='position')+
				# scale_color_discrete(guide=TRUE)+
				# facet_grid( sample+readlen ~ section)+
				facet_grid( stage ~ section,scale='free_x')+
				scale_color_manual(values=stagecols)+
				scale_y_continuous(name='Mean Psite Count / CDS Total',limits=c(0,20))+
				theme_bw()
			}
		)%>%
		ggarrange(plotlist=.,ncol=2,common.legend=T)
	print(rwplot)
	dev.off()
	normalizePath(plotfile)%>%message

	plotfile<-'plots/figures/figure1/fig1c_myribowaltz_aug_stageov.pdf'%T>%pdf(h=4,w=8)
	qplot(data=allpsitedfs_bind%>%filter(section=='AUG',as.numeric(start) > -10,start<10),color=sample,x=start,y=signal,geom='line')+
		scale_x_continuous(name='position',limits=c(-9,9),minor_breaks=number_ticks,breaks=partial(number_ticks,n=3))+
		# scale_color_discrete(guide=TRUE)+
		facet_grid( . ~ section)+
		scale_color_manual(values=samplecols)+
		scale_y_continuous(name='Mean Psite Count / CDS Total',limits=c(0,20))+
		theme_bw()
	dev.off()
	normalizePath(plotfile)%>%message

	plotfile<-'plots/figures/figure1/fig1c_myribowaltz_stop_stageov.pdf'%T>%pdf(h=4,w=8)
	qplot(data=allpsitedfs_bind%>%filter(section=='STOP',as.numeric(start) > -10,start<10),color=sample,x=start,y=signal,geom='line')+
		scale_x_continuous(name='position',limits=c(-9,9),minor_breaks=number_ticks,breaks=partial(number_ticks,n=3))+
		# scale_color_discrete(guide=TRUE)+
		facet_grid( . ~ section)+
		scale_color_manual(values=samplecols)+
		scale_y_continuous(name='Mean Psite Count / CDS Total',limits=c(0,20))+
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

