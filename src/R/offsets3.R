#' ---
#' title: "Sequence specific offset determination"
#' author: "Dermot Harnett"
#' ---


#+ setup, include=FALSE, echo=FALSE, eval=T
knitr::opts_chunk$set(root.dir = here::here(),eval=FALSE,cache=FALSE,echo=FALSE,warning = FALSE,message = FALSE,include=FALSE)
isknitr<-isTRUE(getOption('knitr.in.progress'))


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
# suppressMessages({library(riboWaltz)})
suppressMessages({library(purrr)})
suppressMessages({library(here)})
suppressMessages({library(magrittr)})
suppressMessages({library(stringr)})
suppressMessages({library(tidyverse)})
suppressMessages({library(GenomicAlignments)})
suppressMessages({library(GenomicFeatures)})
suppressMessages({library(GenomicFiles)})
suppressMessages({library(bamsignals)})
suppressMessages({library(zeallot)})
suppressMessages({library(here)})

reduce <- GenomicRanges::reduce

MAPQTHRESH <- 50
USEPHASE <- FALSE
USERIBOSEQC <- FALSE
source(here('src/R/Rprofile.R'))



for(fname in lsf.str('package:dplyr')) assign(fname,get(fname,'package:dplyr'))

source(here('src/R/Rprofile.R'))

argv <- c(
	bam = here('pipeline/star/data/E13_ribo_2/E13_ribo_2.bam')%T>%{stopifnot(file.exists(.))},
	gtf = here('pipeline/my_gencode.vM12.annotation.gtf'),
	REF = here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'),
	outfolder = 'seqshift_reads/data/E13_ribo_2/'
)

sample_params<-here::here('src/sample_parameter.csv')%>%fread%>%filter(assay=='ribo')%>%.$sample_id


# bams <- paste0('star/data/',sample_params,'/',sample_params,'.bam')%>%map(. %T>%{stopifnot(file.exists(.))})

argv[] <- commandArgs(trailing=TRUE)

for (nm in names(argv)) assign(nm,argv[[nm]])




#don't memoise if running on
if(!interactive()) mymemoise<-identity

# purely(function(){
STOPWINDOWSTART= -2
STOPWINDOWEND= 2

filtercds4train <- .%>%subset(tag=='CCDS')

setpos<- .%>%{strand(.)<-'+';.}




string2onehot<-function(scol) vapply(c('A','C','T','G'),function(lev)as.numeric(scol==lev),rep(1,length(scol)))
dnaseq2onehot <- function(mat,pre){
	mat<-as.matrix(mat);
	lapply(1:ncol(mat),function(n) string2onehot(mat[,n])%>%set_colnames(paste0(pre,n,'.',colnames(.))))%>%purrr::reduce(cbind)
}


get_transcript_compartments<-function(topcdsexons,compartments){
	#add compartment to the mapped reads
	trcomps<-data_frame(chr=as.character(topcdsexons@seqnames),transcript_id=topcdsexons$transcript_id)%>%
		mutate(compartment=compartments[chr])%>%
		# distinct(transcript_id,compartment)%>%
		{setNames(CharacterList(as.list(.$compartment)),.$transcript_id)}%>%
		unlist
	trcomps
}

get_cdsread_trmap <- function(topcdsreads,topcdsexons,topcdsmap){
	stopifnot('compartment' %in% colnames(mcols(topcdsexons)))
	#reads as a gr
	topcdsreadsgr<-topcdsreads%>%as("GRanges")
	#get the shifts and store in hte gr
	topcdsreadsgr$length<-qwidth(topcdsreads)
	#now map these to transcript space
	cdsread_trmap<-topcdsreadsgr%>%subsetByOverlaps(topcdsexons,type='within')%>%splitmaptotranscripts(topcdsexons)

	# testwidde<-cdsread_trmap%>%subset(width>40)%>%.[1])
	cdsread_trmap$length<-topcdsreadsgr$length[cdsread_trmap$xHits]

	cdsread_trmap$compartment<-topcdsexons$compartment[cdsread_trmap$transcriptsHits]

	#select only those which mapped cleanly onto transcripts
	cdsread_trmap%<>%trim
	cdsread_trmap<-cdsread_trmap%>%subset(width==length)

	starcods_trmap <- topcdsmap%>%resize(3,ignore.strand=TRUE)%>%setNames(.,as.character(seqnames(.)))

	cdsread_trmap$phase <- ((start(cdsread_trmap) - start(starcods_trmap[seqnames(cdsread_trmap)])) %%3)


	assert_that({
		non3ttr<-coverage(cdsread_trmap)[topcdsmap]%>%lengths%>%map_dbl(`%%`,3)%>%.[.!=0]%>%names;
		length(non3ttr)==0#all the cds should be multiples of 3
	})

	assert_that(has_name(mcols(cdsread_trmap),c('length','compartment','phase')))

	cdsread_trmap

}







get_cdsread_trmap<-mymemoise(get_cdsread_trmap)

# get_cdsread_trmap<-function(topcdsreads,topcdsexons,startcods_trmap){
# 	#reads as a gr
# 	topcdsreadsgr<-topcdsreads%>%as("GRanges")
# 	#get the shifts and store in hte gr
# 	topcdsreadsgr$length<-qwidth(topcdsreads)
# 	#now map these to transcript space


# 	cdsread_trmap<-topcdsreadsgr%>%splitmaptotranscripts(topcdsexons)

# 	# testwidde<-cdsread_trmap%>%subset(width>40)%>%.[1])
# 	cdsread_trmap$length<-topcdsreadsgr$length[cdsread_trmap$xHits]

# 	cdsread_trmap$compartment<-topcdsexons$compartment[cdsread_trmap$transcriptsHits]

# 	#select only those which mapped cleanly onto transcripts
# 	cdsread_trmap%<>%trim
	
# 	cdsread_trmap<-cdsread_trmap%>%subset(width==length)

# 	cdsread_trmap$phase <- ((start(cdsread_trmap) - start(starcods_trmap[seqnames(cdsread_trmap)])) %%3 )

# 	non3ttr<-coverage(cdsread_trmap)[topcdsmap]%>%lengths%>%map_dbl(`%%`,3)%>%.[.!=0]%>%names
# 	stopifnot(length(non3ttr)==0)#all the cds should be multiples of 3

# 	cdsread_trmap
# }
# oldenv$cdsread_trmap
# newenv$cdsread_trmap

# intersect(oldenv$cdsread_trmap,newenv$cdsread_trmap)
# intersect(newenv$cdsread_trmap,oldenv$cdsread_trmap)

#for now let's not consider base pairs with no strong pwm


get_offsets<-function(cdsread_trmap,compartments){

	readsizes <- select_toplengths(cdsread_trmap$length)

	possibleoffsets <-seq(3,max(readsizes),by=1)%>%setNames(.,.)
	# possibleoffsets <-seq(6,max(readsizes),by=3)%>%setNames(.,.)

	stopifnot(!is.null(cdsread_trmap$compartment))
		
	phaselist<-list(0,1,2,0:2)%>%setNames(c(0,1,2,'all'))

	offset_cds_scores <-	lapply(possibleoffsets,function(offset){
			lapply(unique(compartments),function(compartment_i){
				lapply(readsizes,function(length_i){
					lapply(phaselist,function(phase_i){
						if((length_i - 6 - offset) < 0) return(data_frame(score=NA))
						
						cat('.')
						
						cdsread_trmap%>%
							subset(compartment==compartment_i)%>%
							subset(phase%in%phase_i)%>%
							subset(length==length_i)%>%
							resize(1,'start',ignore.strand=T)%>%
							shift(offset)%>%
							countOverlaps(topcdsmap)%>%
							`>`(0)%>%sum%>%data_frame(score=.)
					})%>%bind_rows(.id='phase')
				})%>%bind_rows(.id='length')
			})%>%bind_rows(.id='compartment')
		})%>%bind_rows(.id='offset')

	offset_cds_scores%<>% filter( ((as.numeric(offset) %%3) ==0) | (phase=='all') )
	offset_cds_scores%<>%mutate_at(vars(everything()),as.numeric)
	offset_cds_scores$compartment<-unique(compartments)[offset_cds_scores$compartment]


	bestscores<-offset_cds_scores%>%group_by(length,phase,compartment)%>%slice(which.max(score))

	list(bestscores,offset_cds_scores)
}


# get_offsets<-mymemoise(get_offsets)

#define compartments
get_compartments<-function(cds,circs=DEFAULT_CIRC_SEQS){
	compartments <- rep('nucl',length(seqlevels(cds)))%>%setNames(seqlevels(cds))
	circs_in_data <- intersect(circs,names(compartments))
	compartments[circs_in_data] <- circs_in_data
	compartments
}

#' - select the ones that have the best score for each length and phase
#We now have optimal offsets per readlength/phase


select_toplengths <- function(x) x %>%table%>%{./sum(.)}%>%keep(. > 0.05)%>%names%>%as.numeric%>%setNames(.,.)

#



apply_cds_offsets<- function(reads_tr,bestscores){
	compartmentcol <- if('compartment'%in%colnames(bestscores)) 'compartment' else NULL
	phasecol <- if('phase'%in%colnames(bestscores)) 'phase' else NULL
	joincols <- c('length',compartmentcol,phasecol)
	stopifnot(joincols %in% colnames(mcols(reads_tr)))

	reads_tr%>%
		mcols%>%
		as_tibble%>%
		safe_left_join(bestscores%>%
			ungroup%>%
			# select(-phase,-score),by=c('length','compartment'),
			select(offset,one_of(joincols)),by=joincols,
			allow_missing=TRUE
		)%>%
		.$offset
}

# reads_tr<-oldenv$cdsread_trmap
fp<-function(gr) ifelse(strand(gr)=='-',end(gr),start(gr))
tp<-function(gr) ifelse(strand(gr)=='-',start(gr),end(gr))
fpend<-function(x)resize(x,1,'start')
tpend<-function(x)resize(x,1,'end')

gr1<-GRanges(c('chr1:5-6:+'))
gr2<-GRanges(c('chr1:50-51:+','chr1:40-51:+'))

downstream_dist_till<-function(gr1,gr2){
	(fp(gr2)[precede(gr1,gr2)] - tp(gr1)) * ( ifelse(strand(gr1)=='+',1,-1))
}

upstream_dist_till<-function(gr1,gr2){
	(tp(gr2)[follow(gr1,gr2)] - fp(gr1)) * ( ifelse(strand(gr1)=='+',-1,1))
}

tpmost<-function(cds,groupvar='transcript_id'){
	
	ids<-data_frame(id=seq_along(cds),end=end(cds),strand=as.vector(strand(cds)),groupvar=mcols(cds)[[groupvar]])%>%group_by(groupvar)%>%slice(which.max(end*ifelse(strand=='-',-1,1)))%>%.$id
	cds[ids]
}
fpmost<-function(cds,groupvar='transcript_id'){
	ids<-data_frame(id=seq_along(cds),start=start(cds),strand=as.vector(strand(cds)),groupvar=mcols(cds)[[groupvar]])%>%group_by(groupvar)%>%slice(which.max(start*ifelse(strand=='-',1,-1)))%>%.$id
	cds[ids]
}


clip_start <- function(x,n) resize(x,width(x)-n,fix='end')
clip_end <- function(x,n) resize(x,width(x)-n,fix='start')
setstrand<-function(x) {strand(x)<-Rle('+') 
	x}

testthat::expect_equal(downstream_dist_till(GRanges(c('chr1:10:-')),GRanges(c('chr1:5:-','chr1:40-51:+'))),5)
testthat::expect_equal(downstream_dist_till(GRanges(c('chr1:10:+')),GRanges(c('chr1:14:+','chr1:40-51:+'))),4)

testthat::expect_equal(upstream_dist_till(GRanges(c('chr1:10:-')),GRanges(c('chr1:11:+','chr1:12:-','chr1:40-51:+'))),2)

testthat::expect_equal(upstream_dist_till(GRanges(c('chr1:10:+')),GRanges(c('chr1:11:+','chr1:6:+','chr1:40-51:+'))),4)
testthat::expect_equal(upstream_dist_till(GRanges(c('chr1:10:-')),GRanges(c('chr1:11:+','chr1:14:-','chr1:40-51:+'))),4)



shorten_ranges<-function(gr){gr <- shift(gr,- min(start(gr))) ;gr}

# cdsbak<-cds
# cds<-cdsbak%>%head(100e3)


# cds<-c(
# cds%>%subset(strand=='+')%>%split(.,.$transcript_id)%>%.[1],
# cds%>%subset(strand=='-')%>%split(.,.$transcript_id)%>%.[1]
# )
# cds%<>%unlist


get_stop_sites <- function(reads_tr,topcdsmap,STOPWINDOWSTART,STOPWINDOWEND){
	message('getting stop sites for trainng')
	assert_that(mcols(reads_tr)%>%has_name('cdsshift'))
	setpos<- .%>%{strand(.)<-'+';.}

	prestops <- topcdsmap%>%setpos%>%resize(3,'end')%>%resize(1,'start')%>%setNames(seqnames(.))
	starts <- topcdsmap%>%setpos%>%resize(1,'start')%>%setNames(seqnames(.))

	#all stop codons!
	# stopifnot(topExonseq[prestops%>%resize(3)%>%shift(3)]%>%translate%>%`==`('*'))

	#get the psites predicted to -+2 of the first nt of our pre-stop codons

	#get the distance to the 1 site of the prestop
	reads_tr$stopdist<-start(reads_tr)+reads_tr$cdsshift - start(prestops[seqnames(reads_tr)])
	reads_tr$startdist<-start(reads_tr)+reads_tr$cdsshift - start(starts[seqnames(reads_tr)])

	# reads_tr$stopdist<-start(reads_tr) - start(prestops[seqnames(reads_tr)])
	# reads_tr$startdist<-start(reads_tr) - start(starts[seqnames(reads_tr)])

	startpsites<-reads_tr%>%subset(between(startdist,STOPWINDOWSTART,0))%>%{.$dist<-factor(.$startdist);.}
	stoppsites<-reads_tr%>%subset(between(stopdist,0,STOPWINDOWEND))%>%{.$dist<-factor(.$stopdist);.}

	# c(startpsites,stoppsites)

	# startpsites downstream_dist_till
	# brower()

	# downstream_dist_till(fpend(startpsites),fpend(starts)) - startpsites$cdsshift

	# downstream_dist_till(fpend(stoppsites),fpend(prestops)) - stoppsites$cdsshift


	# stoppsites<-reads_tr%>%subset(between(stopdist,STOPWINDOWSTART,STOPWINDOWEND))%>%{.$dist<-.$stopdist;.}

	c(startpsites,stoppsites)
}

#TODO - length+comaprtment not just length
get_all_stop_sites <- function(bam,allcds,STOPWINDOWSTART,STOPWINDOWEND,bestscores){

	message('getting stop sites for training seq model - all stop sites')

	#work out start and stop codons from the cds


	startcodons <- fpmost(allcds,'transcript_id')%>%resize(3,'start')

	prestopcodons <- tpmost(allcds,'transcript_id')%>%resize(3,'end')


	starts_ends<-c(startcodons,prestopcodons)
	
	# browser()

	assert_that(startcodons%>%strandshift(0)%>%getSeq(FaFile(REF),.)%>%translate%>%`==`('M')%>%mean%T>%message%>%`>`(0.8))
	assert_that(prestopcodons%>%strandshift(3)%>%getSeq(FaFile(REF),.)%>%str_subset(neg=TRUE,'N')%>%DNAStringSet%>%translate%>%`==`('*')%>%mean%T>%message%>%`>`(0.8))


	reads_tr <- bam%>%BamFile(yieldSize=NA)%>%readGAlignments(param=ScanBamParam(which=starts_ends))
	reads_tr <- GRanges(reads_tr,length=qwidth(reads_tr))

# 	teststart<-GRanges(as.character(startcodons[1]))
# #	teststart<-GRanges(as.character(startcodons[10]))
# 	reads_tr %<>%subsetByOverlaps(teststart)
# 	reads_tr%<>%shift(-start(teststart)+100)
# 	teststart%<>%shift(-start(teststart)+100)



	reads_tr$startdist <- - downstream_dist_till(fpend(reads_tr),fpend(startcodons))	

	reads_tr$stopdist <- - downstream_dist_till(fpend(reads_tr),fpend(prestopcodons))	

	reads_tr$phase <- ifelse(abs(reads_tr$startdist)<abs(reads_tr$stopdist),
		(reads_tr$startdist) %%3,
		(reads_tr$stopdist)  %%3
	)


	maxcdsshift <- bestscores$offset%>%max

	reads_tr$stopdist[between(reads_tr$stopdist,-20,20)] %>% add(maxcdsshift)%>%table
	reads_tr$stopdist[between(reads_tr$stopdist,-20,20)] %>% add(maxcdsshift)%>%table


	reads_tr$compartment <- as(compartments,'List')[seqnames(reads_tr)]%>%unlist(use.names=F)

	reads_tr$cdsshift <- reads_tr%>%apply_cds_offsets(bestscores)

	table((reads_tr$startdist + reads_tr$cdsshift)%>%.[between(.,-20,20)])
	table((reads_tr$stopdist + reads_tr$cdsshift)%>%.[between(.,-20,20)])


	reads_tr$startdist <- reads_tr$startdist + reads_tr$cdsshift
	reads_tr$stopdist <- reads_tr$stopdist + reads_tr$cdsshift





	#So our cds shift is a negative factor
	#above means that things landing just upstrema of starts (our errors to correct, will not yet be)
	#negative enough - i.e. they will be positive
	#and stops that are shifted outside the CDS will be too negative

	# reads_tr$startdist%>%.[between(.,-20,20)]%>%na.omit%>%txtdensity
	# reads_tr$stopdist%>%.[between(.,-20,20)]%>%na.omit%>%txtdensity


	# reads_tr$startdist%>%.[order(abs(.))]%>%.[.!=0]%>%between(.,-STOPWINDOWSTART,0)%>%na.omit%>%any
	# #all stop codons!
	# # stopifnot(topExonseq[prestops%>%resize(3)%>%shift(3)]%>%translate%>%`==`('*'))

	# #get the psites predicted to -+2 of the first nt of our pre-stop codons

	# #get the distance to the 1 site of the prestop
	# reads_tr$stopdist<-start(reads_tr)+reads_tr$cdsshift - start(prestops[seqnames(reads_tr)])
	# reads_tr$startdist<-start(reads_tr)+reads_tr$cdsshift - start(starts[seqnames(reads_tr)])

	# # reads_tr$stopdist<-start(reads_tr) - start(prestops[seqnames(reads_tr)])
	# # reads_tr$startdist<-start(reads_tr) - start(starts[seqnames(reads_tr)])


	startpsites<-reads_tr%>%subset(between(startdist,0,STOPWINDOWEND))%>%{.$dist<-.$startdist;.}
	stoppsites<-reads_tr%>%subset(between(stopdist,STOPWINDOWSTART,0))%>%{.$dist<-.$stopdist;.}
 	

	startpsites$dist%>%table

	reads_tr%>%subset(between(stopdist,STOPWINDOWSTART,0))%>%.$stopdist%>%table

	out <- 	c(startpsites,stoppsites)
 	
 	stopifnot((STOPWINDOWSTART:STOPWINDOWEND)%in%out$dist)

 	out$dist %<>% as.factor

	c(out)

}

# get_stop_sites <- function(reads_tr,topcdsmap,STOPWINDOWSTART,STOPWINDOWEND){
# 	message('getting stop sites for trainng')
# 	assert_that(mcols(reads_tr)%>%has_name('cdsshift'))
# 	setpos<- .%>%{strand(.)<-'+';.}

# 	prestops <- topcdsmap%>%setpos%>%resize(3,'end')%>%resize(1,'start')%>%setNames(seqnames(.))
# 	starts <- topcdsmap%>%setpos%>%resize(1,'start')%>%setNames(seqnames(.))

# 	#all stop codons!
# 	# stopifnot(topExonseq[prestops%>%resize(3)%>%shift(3)]%>%translate%>%`==`('*'))

# 	#get the psites predicted to -+2 of the first nt of our pre-stop codons

# 	#get the distance to the 1 site of the prestop
# 	# reads_tr$stopdist<-start(reads_tr)+reads_tr$cdsshift - start(prestops[seqnames(reads_tr)])
# 	# reads_tr$startdist<-start(reads_tr)+reads_tr$cdsshift - start(starts[seqnames(reads_tr)])

# 	reads_tr$stopdist<-start(reads_tr) - start(prestops[seqnames(reads_tr)])
# 	reads_tr$startdist<-start(reads_tr) - start(starts[seqnames(reads_tr)])

# 	stoppsites<-reads_tr%>%subset(between(stopdist,0,STOPWINDOWEND))%>%{.$dist<-.$stopdist;.}
# 	startpsites<-reads_tr%>%subset(between(startdist,STOPWINDOWSTART,0))%>%{.$dist<-.$startdist;.}

# 	c(startpsites,stoppsites)

# 	# stoppsites<-reads_tr%>%subset(between(stopdist,STOPWINDOWSTART,STOPWINDOWEND))%>%{.$dist<-.$stopdist;.}

# 	# c(stoppsites)
# }

#This should be modified to take either a fastq file or the exonsequence
#And it should be okay with missing sequences

setMethod('seqinfo','DNAStringSet',function(x){
	x%>%{Seqinfo(names(.),nchar(.))}
})

get_seqforrest_traindata <- function(reads,seq,nbp=2){


	stopifnot('length' %in% colnames(mcols(reads)))
	stopifnot(all(seqnames(reads)%in%seqinfo(seq)@seqnames))


	suppressWarnings({fp<-reads[TRUE]%>%resize(nbp,'start')%>%resize(nbp*2,'end')})
	fp_inbounds <- !is_out_of_bounds(fp,seqinfo(seq))

	suppressWarnings({tp <- reads%>%resize(nbp,'end')%>%resize(nbp*2,'start')})
	tp_inbounds <- !is_out_of_bounds(tp,seqinfo(seq))

	inbounds <- fp_inbounds & tp_inbounds

	fp<-fp[inbounds]
	tp<-tp[inbounds]
	
	startseq <- getSeq(seq,fp)
	endseq <- getSeq(seq,tp)

	seqmat <- cbind( dnaseq2onehot(startseq,'fp.'),dnaseq2onehot(endseq,'tp.'))

	seqmat%<>%cbind(length=reads$length[inbounds])%>%as.data.frame

	if(has_name(mcols(reads),'dist')) seqmat %<>% cbind(dist=reads$dist[inbounds])


	seqmat

}

get_seqoffset_model <-  function(data4readshift){
	require(ranger)
	assert_that(has_name(data4readshift,c('length','fp.2.A','tp.2.A')))
		
	seqshiftmodel <- ranger::ranger(formula= dist ~ . ,
		data=data4readshift,
		importance='permutation'
	)


	
	seqshiftmodel

}

get_predictedshifts<-function(seqshiftmodel,data4readshift){
	outpred <- predict(seqshiftmodel,data=data4readshift)$prediction

	outpred <- outpred%>%as.character%>%as.numeric

	assert_that(all((outpred %% 1)==0))

	outpred
}



get_top_trs<-function(bam,exons,cds,startcods){
	message('getting top trs')
	#get counts over all cds
	cdsstren<-mymemoise(bamCount)(bam,cds)
	cds$count<-cdsstren
	trcounts<-cds$count%>%split(cds$transcript_id)%>%map_dbl(sum)


	#gene tr relationshiop df
	gtrdf<-exons%>%mcols%>%.[,c('gene_id','transcript_id')]%>%as.data.frame%>%distinct



	#' - For a given library, Now take the trancsript with the most overlaps, for the top 1k genes

	startcodcount <- startcods%>%.$transcript_id%>%table

	simpletrs<-startcodcount[startcodcount==1]%>%enframe%>%select(transcript_id=name)
	ccds <- cds%>%mcols%>%as.data.frame%>%select(transcript_id,tag)%>%filter(tag=='CCDS')%>%distinct
	stopifnot(nrow(ccds)>1e3)


	#get the tr with the highest signal per gene
	#get the tr with the highest signal per gene

	toptrs <- gtrdf%>%
		left_join(enframe(trcounts,'transcript_id','count'))%>%
		semi_join(simpletrs)%>%
		semi_join(ccds)%>%
		group_by(gene_id)%>%
		slice(which.max(count))%>%
		arrange(desc(count))%>%
		.$transcript_id%>%
		head(1e3)


	#' - Ensuring these have only 1 start codon, for simplicity

	topcds <- cds%>%
		subset(transcript_id %in% toptrs)%>%
		identity
		# head(1e3)

	topstartcods<-startcods[match(topcds%>%.$transcript_id%>%unique,startcods$transcript_id)]

	#get the exons for these
	topcdsexons <- exons%>%subset(transcript_id %in% toptrs)

	#' - Map our reads to these transcripts

	#mapped cds
	topcdsmap<-topcds%>%pmapToTranscripts(topcdsexons%>%split(.,.$transcript_id)%>%.[topcds$transcript_id])%>%reduce
	#also start codons
	starcods_trmap<-topstartcods%>%mapToTranscripts(topcdsexons%>%split(.$transcript_id))
	starcods_trmap%<>%setNames(as.character(seqnames(.)))
	return(list(topcdsexons,topcds,topcdsmap))

}

get_seqoffset_model_elim <-  function(traindata4readshift){
	require(ranger)
	stopifnot((-2:2) %in% traindata4readshift$dist)


	assert_that(is.factor(traindata4readshift$dist))

	assert_that(has_name(traindata4readshift,c('length','fp.2.A','tp.2.A')))

	# stopifnot(traindata4readshift$dist %in% c(STOPWINDOWEND:STOPWINDOWSTART))

	istestset = sample(((1:nrow(traindata4readshift)) %% 5) == 0)

	trainset <- traindata4readshift[!istestset,]
	testset <- traindata4readshift[istestset,]
	trainset<-trainset[sample(1:nrow(trainset))%>%head(10e3),]
	trainset%>%nrow
	
	tabs<-list()
	varstouse = colnames(trainset)%>%str_subset(neg=TRUE,'length|dist')

	# browser()
	# exportenv
	minerror <- NA
	while(length(varstouse)>0){
		seqshiftmodel <- ranger::holdoutRF(formula= dist ~ . ,
			# data=	traindata4readshift[,c(varstouse,'length','dist')]
			data=	trainset%>%select(one_of(varstouse),dist,length)
		# importance='permutation'
		)
		leastimportantvar <- seqshiftmodel$variable.importance%>%.[!names(.)%in%c('length','dist')]%>%abs%>%sort%>%head(1)%>%names
		varstouse <- setdiff(varstouse,	leastimportantvar)
		message('removed ');message(leastimportantvar)
		newerr <- seqshiftmodel$rf2$prediction.error
		if(isTRUE( newerr> (minerror+0.1))) break
		message(newerr)



		# tabs <- append(tabs,mean(testset$dist==predict(seqshiftmodel,testset)$prediction))
		tabs <- append(tabs,setNames(seqshiftmodel$rf2$prediction.error,leastimportantvar))
	}

	#So if we use 10k examples we can trim off the less important variables pretty easily.
	# keepvars <- (tabs%>%unlist < (tabs%>%unlist%>%min)+0.05)%>%.[!.]%>%names


	browser()

	#use ones we didn't eliminate
	keepvars <- colnames(trainset)%>%str_subset(neg=TRUE,'length|dist')%>%setdiff(names(tabs))

	seqshiftmodel <- ranger::ranger(formula= dist ~ . ,
			# data=	traindata4readshift[,c(varstouse,'length','dist')]
			data=	testset%>%select(one_of(keepvars),dist,length)
	)

	seqshiftmodel_allvars <- ranger::ranger(formula= dist ~ . ,
			# data=	traindata4readshift[,c(varstouse,'length','dist')]
			data=	testset%>%select(everything())
	)

	seqshiftmodel_len_only <- ranger::ranger(formula= dist ~ . ,
			# data=	traindata4readshift[,c(varstouse,'length','dist')]
			data=	testset%>%select(dist,length)
	)
	
	list(seqshiftmodel,seqshiftmodel_allvars,seqshiftmodel_len_only)

}
fp <-function(gr)ifelse(strand(gr)=='-',end(gr),start(gr))
tp <-function(gr)ifelse(strand(gr)=='-',start(gr),end(gr))
strandshift<-function(gr,shift) shift(gr , ifelse( strand(gr)=='-',- shift,shift))

readtopreads<-function(bam,top_cds=topcds){bam%>%BamFile(yieldSize=NA)%>%readGAlignments(param=ScanBamParam(which=top_cds))}

#testing funcs
cdscov <- function(reads2use,cds2use){
	out = coverage(reads2use)[cds2use]
	out <- suppressWarnings({unlist(out)})
	out<-out%>%as.vector
	stopifnot((length(out)%%3)==0)
	out
}

get_periodicity_scores<-function(cdstosamp,cdsread_trmap){
	message('.')
	assert_that(has_name(mcols(cdsread_trmap),c('seqshift','cdsshift')))
	psites2use<-cdsread_trmap%>%
		keepSeqlevels(cdstosamp,'coarse')%>%
		resize(1,'start')%>%
		shift(.$cdsshift)

	cds2use <- topcdsmap[cdstosamp]

	cds2use %<>% clip_start(30)%>%clip_end(30)

	periodicfrac<-cdscov(psites2use,cds2use)%>%get3bpfrac

	data_frame(
		no_seqshift = psites2use%>%cdscov(cds2use)%>%get3bpfrac,
		with_seqshift = psites2use %>% shift(-.$seqshift)%>%cdscov(cds2use)%>%get3bpfrac
	)
}
# get_periodicity_scores <- mymemoise(get_periodicity_scores)



get_frac_inds<-function(vect,l,r){
	lind <- floor(l*length(vect))
	rind <- floor(r*length(vect))
	vect[lind:rind]
}
# get_frac_inds(1:20,0.3,0.36)
get3bpfrac<-function(x) x %>%
	# add(c(1000000,0,0))%>%
	fft%>%
	abs%>%
	`^`(2)%>%
	{sum(get_frac_inds(.,0.3,0.36))/sum(.) }




# readtopreads<-mymemoise(readtopreads)


mimport<-mymemoise(function(...)rtracklayer::import(...))
mapToTranscripts<-mymemoise(GenomicFeatures::mapToTranscripts)
pmapToTranscripts<-mymemoise(GenomicFeatures::pmapToTranscripts)
readGAlignments<-mymemoise(GenomicAlignments::readGAlignments)
splitmaptotranscripts <- mymemoise(function(reads,transcripts) reads%>%split(.,ceiling(seq_along(.)/ 50e3))%>%lapply(.%>%mapToTranscripts(transcripts%>%split(.$transcript_id)))%>%Reduce(f=c))


#get exons
 if(!exists('gtf_gr')) gtf_gr<-mimport(con=gtf,format='gtf')
# exons <- gtf_gr%>%subset(type=='exon')

# if(!is('exons','GRanges')) exons <- gtf_gr%>%subset(type=='exon')
# if(!is('cds','GRanges')) cds <- gtf_gr%>%subset(type=='CDS')
# if(!exists('startcodsa')) startcods <- gtf_gr%>%subset(type=='start_codon')
exons <- gtf_gr%>%subset(type=='exon')
cds <- gtf_gr%>%subset(type=='CDS')
startcods <- gtf_gr%>%subset(type=='start_codon')

#First select some top exons for the offset finding process
# stopifnot(identical(bak,list(bam,exons,cds,startcods)))
# bak <- list(bam,exons,cds,startcods)


#get counts over all cds
cdsstren<-mymemoise(bamCount)(bam,cds)
cds$count<-cdsstren
trcounts<-cds$count%>%split(cds$transcript_id)%>%map_dbl(sum)

#gene tr relationshiop df
gtrdf<-exons%>%mcols%>%.[,c('gene_id','transcript_id')]%>%as.data.frame%>%distinct



#' - For a given library, Now take the trancsript with the most overlaps, for the top 1k genes

startcodcount <- startcods%>%.$transcript_id%>%table

simpletrs<-startcodcount[startcodcount==1]%>%enframe%>%select(transcript_id=name)
ccds <- cds%>%mcols%>%as.data.frame%>%select(transcript_id,tag)%>%filter(tag=='CCDS')%>%distinct
stopifnot(nrow(ccds)>1e3)
hasfputr <- gtf_gr%>%
	subset(
		type=='UTR')%>%subsetByOverlaps(gtf_gr%>%subset(type=='start_codon')%>%strandshift(-1)
	)%>%
	subset(width>40)%>%
	.$transcript_id%>%data_frame(transcript_id=.)
#get the tr with the highest signal per gene
#get the tr with the highest signal per gene
toptrs <- gtrdf%>%
	left_join(enframe(trcounts,'transcript_id','count'))%>%
	semi_join(simpletrs)%>%
	semi_join(ccds)%>%
	semi_join(hasfputr)%>%
	group_by(gene_id)%>%
	slice(which.max(count))%>%
	arrange(desc(count))%>%
	.$transcript_id%>%
	head(1e3)


#' - Ensuring these have only 1 start codon, for simplicity
topcds <- cds%>%
	subset(transcript_id %in% toptrs)%>%
	identity
	# head(1e3)

topstartcods<-startcods[match(topcds%>%.$transcript_id%>%unique,startcods$transcript_id)]

#get the exons for these
topcdsexons <- exons%>%subset(transcript_id %in% toptrs)

#' - Map our reads to these transcripts

#mapped cds
topcdsmap<-topcds%>%pmapToTranscripts(topcdsexons%>%split(.,.$transcript_id)%>%.[topcds$transcript_id])%>%reduce
topcdsmap%<>%setNames(.,as.character(seqnames(.)))

#also start codons
starcods_trmap<-topstartcods%>%mapToTranscripts(topcdsexons%>%split(.$transcript_id))
starcods_trmap%<>%setNames(as.character(seqnames(.)))


library(GenomicAlignments)
#get reads over them
if(!exists('topcdsreadsbak'))topcdsreadsbak <- bam %>% readtopreads
topcdsreads<-topcdsreadsbak


compartments <- get_compartments(cds,DEFAULT_CIRC_SEQS)

topcdsexons$compartment <- get_transcript_compartments(topcdsexons,compartments)

# cdsread_trmap <- get_cdsread_trmap(topcdsreads,topcdsexons,startcods_trmap)

cdsread_trmap <- get_cdsread_trmap(topcdsreads,topcdsexons,topcdsmap)
cdsread_trmap %<>% setpos


topExonseq <- getSeq(x=FaFile(REF),topcdsexons)%>%
	split(topcdsexons$transcript_id)%>%
	lapply(.%>%unlist)%>%DNAStringSet



# # # cdsread_trmap<-oldenv$cdsread_trmap%>%{mcols(.)$seqshift=NULL;mcols(.)$startdist=NULL;mcols(.)$stopdist=NULL;.}


#get the cds offsets, and export them
#okay so the get offsets procedure now seems to work fine when used with the old cdsread_trmap procedure
c(allbestscores,offset_cds_scores) %<-% mymemoise(get_offsets)(cdsread_trmap,compartments)

outfolder%>%dir.create(showWarnings=F,rec=T)
allbestscores%>%write_tsv(file.path(outfolder,'cdsmax_offsets.tsv'))




if(!USEPHASE){
	bestscores<-allbestscores%>%filter(is.na(phase))%>%ungroup%>%select(-phase)
}else{
	bestscores<-allbestscores%>%filter(!is.na(phase))
}
if(USERIBOSEQC){
	#Get riboqc_cutoffs
	sample = bam%>%dirname%>%basename 
	riboqccutoffs <- str_interp(here('pipeline/riboqc/data/${sample}/_P_sites_calcs'))%T>%{stopifnot(file.exists(.))}
	riboqcdf <- riboqccutoffs%>%fread%>%select(compartment=comp,length=read_length,offset=cutoff)
	bestscores<-riboqcdf%>%filter(length %in% allbestscores$length)
}



#get the cds shfit offsets for our data
cdsread_trmap$cdsshift <- apply_cds_offsets(cdsread_trmap,bestscores)


#get only those with a define shift, and for which this doesn't put them over the edge of tr
cdsread_trmap <- cdsread_trmap %>% subset(!is.na(.$cdsshift)) %>% intersect(.,trim(.))

topstoppsites <- get_stop_sites(cdsread_trmap,topcdsmap,STOPWINDOWSTART,STOPWINDOWEND)

cds4train <- cds%>%filtercds4train

allstoppsites <- mymemoise(get_all_stop_sites)(bam,cds4train,STOPWINDOWSTART,STOPWINDOWEND,bestscores)

stoppsites <- allstoppsites
allstoppsites$disttable
assert_that(allstoppsites%>%length%>%`==`(336076),msg='number of start/stop sites to train with change')

topstoppsites$dist%>%table%>%txtplot
allstoppsites$dist%>%table%>%txtplot
stoppsites$dist%>%table%>%txtplot

bestscores
stoppsites<-topstoppsites

(function(){


	if('chr1' %in% seqnames(stoppsites)){
	 	traindata4readshift <- get_seqforrest_traindata(stoppsites,FaFile(REF),nbp=2)
	}else{
		traindata4readshift <- get_seqforrest_traindata(stoppsites,topExonseq,nbp=2)

	}

	seqshiftmodel <- mymemoise(get_seqoffset_model)(traindata4readshift)

	c(seqshiftmodel,seqshiftmodel_allvars,seqshiftmodel_len_only) %<-% mymemoise(get_seqoffset_model_elim)(traindata4readshift)


	#we can't get sequence if their flanks extend out of the genome/tr

	cdsread_trmap_inbound <- cdsread_trmap[
		cdsread_trmap%>%resize(width(.)+4,'center')%>%is_out_of_bounds(.)%>%not
	]

	data4readshift <- get_seqforrest_traindata(cdsread_trmap_inbound,topExonseq)

	#it's the model thats actually sure

	cdsread_trmap_inbound$seqshift <-  get_predictedshifts(seqshiftmodel,data4readshift)




	#' ## Testing
	#' 
	#' We can see if our procedure worked by quickly looking at the fft (fast fourier transform) on the top1k genes
	#' binned randomly into 5 groups, before and after the application of our sequence specific cutoffs (to ALL reads),
	#' rather than just those around the stop codon

	seqshift_periodicities<- seqnames(cdsread_trmap_inbound)%>%
		unique%>%
		sample%>%
		# head(100)%>%
		as.character%>%
		split(seq_along(.)%%5)%>%
		mclapply(F=get_periodicity_scores,cdsread_trmap_inbound)%>%
		bind_rows

	seqshift_periodicities%>%unlist%>%txtplot

})()

stop()


#+ spectral_coefficient_strip_plot, fig.width =4,fig.height=4,out.width=400,out.height=450,dev='pdf',include=TRUE,eval=TRUE

# {
# 	seqshift_periodicities%>%
# 		gather(set,spectral_coefficient)%>%
# 		qplot(data=.,x=set,color=set,y=spectral_coefficient,geom=c('point'))+theme_bw()
# }

{

setpos<- .%>%{strand(.)<-'+';.}
prestops <- topcdsmap%>%setpos%>%resize(3,'end')%>%resize(1,'start')
prestarts <- topcdsmap%>%setpos%>%resize(3,'start')%>%resize(1,'start')
cdsread_trmap$startdist<-start(cdsread_trmap) - start(prestarts[seqnames(cdsread_trmap)])
cdsread_trmap$stopdist<-start(cdsread_trmap) - start(prestops[seqnames(cdsread_trmap)])

#Get riboqc_cutoffs
sample = bam%>%dirname%>%basename 
riboqccutoffs <- str_interp(here('pipeline/riboqc/data/${sample}/_P_sites_calcs'))%T>%{stopifnot(file.exists(.))}
riboqcdf <- riboqccutoffs%>%fread%>%select(length=read_length,riboqc_shift=cutoff)

starts <- topcdsmap%>%setpos%>%resize(1,'start')%>%setNames(seqnames(.))
cdsread_trmap_inbound$startdist <- start(cdsread_trmap_inbound) - start(starts[cdsread_trmap_inbound@seqnames])

stops <- topcdsmap%>%setpos%>%resize(3,'end')%>%fpend%>%setNames(seqnames(.))
cdsread_trmap_inbound$stopdist <- start(cdsread_trmap_inbound) - start(stops[cdsread_trmap_inbound@seqnames])


fcods=20
lflank=(fcods*3)
rflank=(fcods*3)+2
mpoint=((fcods*2)+1)*3
epoint=((fcods*4)+2)*3

# cdsread_trmap%>%as.data.frame%>%mutate(d=startdist+stopdist)%>%.$d%>%hist
# cdsread_trmap%>%as.data.frame%>%mutate(d=startdist+stopdist)%>%.$d%>%between(-lflank,rflank)%>%table
metaplotlabs<-c(paste0('-',lflank/2),'AUG',lflank/2,paste0('mid -',lflank/2),'mid',paste0('mid +',lflank/2),paste0('end -',lflank/2),'stop',paste0('end +',lflank/2))
metaplotbreaks<-c(-lflank/2,0,lflank/2,mpoint-lflank/2,mpoint,mpoint+lflank/2,epoint-(lflank/2),epoint,epoint+(lflank/2))
#

seqshiftfuncs <- list(
	Riboqc = .%>% safe_left_join(riboqcdf,by=c('length'))%>% mutate(startdist=startdist+riboqc_shift,stopdist=stopdist+riboqc_shift),
	CDSmax_shift = .%>% mutate(startdist=startdist+cdsshift,stopdist=stopdist+cdsshift),
	seqshift = .%>% mutate(startdist=startdist+cdsshift,stopdist=stopdist+cdsshift)%>%mutate(startdist=startdist-seqshift,stopdist=stopdist-seqshift)
)

#

#' We can also simply look at a metaplot of the frames, (comparing riboqc and cdsmax shift)
plist=lapply(seqshiftfuncs,mymemoise(function(seqshiftfunc){
	#
	p=cdsread_trmap_inbound%>%
		# sample(100e3)%>%
		# shift(.$cdsshift)%>%
		mcols()%>%as_tibble	%>%
		select(startdist,stopdist,length,cdsshift,seqshift)%>%
		# filter(between(abs(startdist-stopdist),-lflank,rflank))%>%
		seqshiftfunc%>%
		mutate(
			mdist = startdist - ((floor(((startdist-stopdist+1)/3)/2))*3),
			instart=between(startdist,-lflank,rflank),
			inend=between(stopdist,-lflank,rflank),
			inmiddle=between(mdist,-lflank,rflank)
		)%>%
		filter((instart|inend|inmiddle) & (!(instart&inend&inmiddle)))%>%
		# filter(instart)%>%
		mutate(x = NA,
			x=ifelse(instart,startdist,x),
			x=ifelse(inmiddle,mdist+mpoint,x),
			x=ifelse(inend,stopdist+epoint,x) ) %>%
		mutate(length=as.character(length))%>%
		bind_rows(.,mutate(.,length='all'))%>%
		group_by(x,length)%>%tally%>%
		# filter(length=='all')%>%
		mutate(phase = factor(x%%3))%>%
		# arrange(desc(n))///
		ggplot(aes(x=x,ymax=n,ymin=0,color=phase))+geom_linerange()+
		facet_grid(scale='free',length~.)+
		geom_vline(linetype=2,xintercept=c(0,epoint+2),alpha=I(0.3))+
		scale_x_continuous(labels=metaplotlabs,breaks=metaplotbreaks)+
		# scale_y_continuous(limits=c(0,300))+
		# coord_cartesian(xlim=c(-lflank,rflank))+
		geom_rect(aes(xmin=0,xmax=epoint,ymin=0,ymax=Inf,alpha=I(0.01)),fill='grey',color='white')+
		scale_fill_discrete(guide=F)+
		theme_minimal()	
	p
}))

#+ riboshift_comparison plot, fig.width =12,fig.height=12,out.width=1200,out.height=1250,dev='pdf',include=TRUE,eval=TRUE

if(isknitr) {
	ggpubr::ggarrange(ncol=3,plotlist=plist,labels=names(seqshiftfuncs))
}else{
	shiftplotfile <- 'shift_methods_comp.pdf'
	shiftplotfile <- file.path(outfolder,shiftplotfile)
	pdf(shiftplotfile,w=14,h=10)
	print(ggpubr::ggarrange(ncol=3,plotlist=plist,labels=names(seqshiftfuncs)))
	dev.off()
	message(normalizePath(shiftplotfile))
}


}

#' .... All of which seems to indicate our procedure is doing a better job of recovering actual P-site positions. Note that it's p
#' site positions are consistent between reads - no dramatic difference between 28 and the others - and manage to recover periodicity
#' in specific read lengths where it was missing. 
#' 
#' Okay so the strand flip error in the stoppsites function is fixed now. We have consstent improvement, however the behavior of the cds shift function is somewhat erratic...
#' 
#' Okay so using the riboseq cutoffs does badly, the 

if((!isknitr) & (interactive())) {
	'/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/R/offsets3.R'%>%
	{rmarkdown::render(knitr::spin(knit=F,.),output_file=here('reports',basename(str_replace(.,'.R$','.html'))))}
}


inenv<-new.env()
testout<-load('pipeline/riboqc/data/RPI12_PolyP0_2/_for_SaTAnn')
