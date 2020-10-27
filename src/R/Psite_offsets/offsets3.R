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
suppressMessages({library(assertthat)})
suppressMessages({library(multitaper)})

reduce <- GenomicRanges::reduce

MAPQTHRESH <- 50
USEPHASE <- TRUE
USERIBOSEQC <- FALSE
base::source(here('src/R/Rprofile.R'))



for(fname in lsf.str('package:GenomicRanges')) assign(fname,get(fname,'package:GenomicRanges'))
for(fname in lsf.str('package:dplyr')) assign(fname,get(fname,'package:dplyr'))


argv <- c(
	bam = here('pipeline/star/data/E175_ribo_1/E175_ribo_1.bam')%T>%{stopifnot(file.exists(.))},
	gtf = here('pipeline/my_gencode.vM12.annotation.gtf'),
	REF = here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'),
	outfolder = here('pipeline/seqshift_reads/data/E175_ribo_1/')
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

string2onehot<-function(scol) vapply(c('A','C','T','G'),function(lev)as.numeric(scol==lev),rep(1,length(scol)))%>%matrix(ncol=4,dimnames=list(NULL,c('A','C','T','G')))

dnaseq2onehot <- function(mat,pre){
	mat<-as.matrix(mat);
	lapply(1:ncol(mat),function(n) string2onehot(mat[,n,drop=FALSE])%>%set_colnames(paste0(pre,n,'.',colnames(.))))%>%purrr::reduce(cbind)
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

	assert_that(all(has_name(mcols(cdsread_trmap),c('length','compartment','phase'))))

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
	assert_that(all(has_name(mcols(reads_tr),joincols)))

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
	assert_that(all(mcols(reads_tr)%>%has_name('cdsshift')))
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


get_seqforrest_traindata <- function(trainreads,seq,topcdsmap=NULL,nbp=2,trim=TRUE){

	stopifnot('length' %in% colnames(mcols(trainreads)))
	stopifnot(all(seqnames(trainreads)%in%seqinfo(seq)@seqnames))


	fp<-trainreads[TRUE]%>%resize(nbp,'start')%>%resize(nbp*2,'end')
	fp_inbounds <- !is_out_of_bounds(fp,seqinfo(seq))

	tp <- trainreads%>%resize(nbp,'end')%>%resize(nbp*2,'start')
	tp_inbounds <- !is_out_of_bounds(tp,seqinfo(seq))
	
	inbounds <- fp_inbounds & tp_inbounds
	inbounds%>%table
	
	if(trim){
		fp<-fp[inbounds]
		tp<-tp[inbounds]
	}else{
		assert_that(all(fp_inbounds))
		assert_that(all(tp_inbounds))
	}
	startseq <- getSeq(seq,fp)
	endseq <- getSeq(seq,tp)


	seqmat <- cbind( dnaseq2onehot(startseq,'fp.'),dnaseq2onehot(endseq,'tp.'))

	seqmat%<>%cbind(length=trainreads$length[inbounds])%>%as.data.frame

	if(! is.null(topcdsmap)){
		dist <- (start(trainreads) - start(topcdsmap[seqnames(trainreads)]))
		seqmat$phase = as.factor(abs(-(dist%%3)))
	}

	if(has_name(mcols(trainreads),'dist')){
		seqmat %<>% cbind(dist=trainreads$dist[inbounds])
		#note that we need to negate the phase to get it from the 'dist', so e.g. -1 means the 2nd base of the stop,
		#while 1 means the 3rd base of the one before AUG
		seqmat$phase = as.factor(abs(-(as.numeric(as.character((seqmat)$dist))%%3)))
	}
	
	seqmat

}

get_seqoffset_model <-  function(data4readshift){
	require(ranger)
	assert_that(all(has_name(data4readshift,c('length','fp.2.A','tp.2.A'))))
		
	seqshiftmodel <- ranger::ranger(formula= dist ~ . ,
		data=data4readshift,
		importance='permutation'
	)


	
	seqshiftmodel

}

get_predictedshifts<-function(seqshiftmodel,data4readshift,probcutoff=0.5){
	require(ranger)
	outpred <- predict(seqshiftmodel,data=data4readshift)$prediction

	if(is.matrix(outpred)){
		maxcol <- max.col(outpred)
		is_highprob <- outpred[matrix(c(1:nrow(outpred),maxcol),ncol=2)] > probcutoff
		outpred=ifelse(is_highprob,colnames(outpred)[maxcol],0)

	}

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

#TODO - length+comaprtment not just length
get_offset_training_sites <- function(bam,allcds){
	# browser()
	message('getting stop sites for training seq model - all stop sites')

	#work out start and stop codons from the cds

	startcodons = allcds%>%sort_grl_st%>%resize_grl(3,'start')%>%unlist
	prestopcodons = allcds%>%sort_grl_st%>%resize_grl(3,'end')%>%unlist

	starts_ends<-c(startcodons,prestopcodons)
	
	assert_that(startcodons%>%sample(100)%>%strandshift(0)%>%getSeq(FaFile(REF),.)%>%translate%>%`==`('M')%>%mean%T>%message%>%`>`(0.8))
	assert_that(prestopcodons%>%sample(100)%>%strandshift(3)%>%getSeq(FaFile(REF),.)%>%str_subset(neg=TRUE,'N')%>%DNAStringSet%>%translate%>%`==`('*')%>%mean%T>%message%>%`>`(0.8))

	reads_tr <- bam%>%BamFile(yieldSize=NA)%>%readGAlignments(param=ScanBamParam(which=starts_ends))
	reads_tr <- GRanges(reads_tr,length=qwidth(reads_tr))


	reads_tr$startdist <- downstream_dist_till(fpend(reads_tr),fpend(startcodons))	
	reads_tr$stopdist <- downstream_dist_till(fpend(reads_tr),fpend(prestopcodons))	

	reads_tr$phase <- ifelse(abs(reads_tr$startdist)<abs(reads_tr$stopdist),
		(3-reads_tr$startdist) %%3,
		(3-reads_tr$stopdist)  %%3
	)
	
	reads_tr$dist <- ifelse(abs(reads_tr$startdist)<abs(reads_tr$stopdist),
		(reads_tr$startdist),
		(reads_tr$stopdist)
	)

	reads_tr%<>%subset(dist >= 4)%>%subset(dist <= (length-4))

	reads_tr
}


get_seqoffset_model_elim <-  function(traindata4readshift,usephase=F,doelim=T){


	require(ranger)

	# stopifnot((-2:2) %in% traindata4readshift$dist)


	assert_that(is.factor(traindata4readshift$dist))

	assert_that(all(has_name(traindata4readshift,c('length','fp.2.A','tp.2.A'))))


	if(usephase){
		#note that we need to negate the phase to get it from the 'dist', so e.g. -1 means the 2nd base of the stop,
		#while 1 means the 3rd base of the one before AUG
		traindata4readshift$phase <- as.factor(abs(-(as.numeric(as.character((traindata4readshift)$dist))%%3)))
	}else{
		traindata4readshift$phase <- NULL
	}
	# stopifnot(traindata4readshift$dist %in% c(STOPWINDOWEND:STOPWINDOWSTART))

	istestset = sample(((1:nrow(traindata4readshift)) %% 5) == 0)

	# trainset <- traindata4readshift[!istestset,]
	# testset <- traindata4readshift[istestset,]
	allset <- traindata4readshift
	
	tabs<-list()
	varstouse = colnames(allset)%>%str_subset(neg=TRUE,'length|dist|phase')

	# exportenv
	minerror <- NA



	if(doelim){
		while(length(varstouse)>0){

			holdoutRF <- mymemoise(ranger::holdoutRF)

			seqshiftmodel <- holdoutRF(formula= dist ~ . ,
				data=	allset%>%select(one_of(varstouse),dist,length)
			)
			#
			leastimportantvar <- seqshiftmodel$variable.importance%>%.[!names(.)%in%c('length','dist','phase')]%>%abs%>%sort%>%head(1)%>%names

			newerr <- seqshiftmodel$rf2$prediction.error
			if(is.na(minerror)) minerror <- newerr#set it the first time

			message(newerr)
			if(isTRUE( newerr > (minerror*1.1))) break


			varstouse <- setdiff(varstouse,	leastimportantvar)

			message('removed ');message(leastimportantvar)
			

			oldseqshiftmodel<-seqshiftmodel
			# tabs <- append(tabs,mean(testset$dist==predict(seqshiftmodel,testset)$prediction))
			tabs <- append(tabs,setNames(seqshiftmodel$rf2$prediction.error,leastimportantvar))
		}
			#use ones we didn't eliminate
		keepvars <- colnames(traindata4readshift)%>%str_subset(neg=TRUE,'length|dist|phase')%>%setdiff(names(tabs))

		if(usephase) keepvars %<>% append('phase')


		

	}else{
		seqshiftmodel=NULL
	}

	#So if we use 10k examples we can trim off the less important variables pretty easily.
	# keepvars <- (tabs%>%unlist < (tabs%>%unlist%>%min)+0.05)%>%.[!.]%>%names

	
	
	allseqvars <-  colnames(allset)%>%str_subset(neg=TRUE,'length|dist|phase')

	seqshiftmodel_allvars <- ranger::ranger(formula= dist ~ . ,
			# data=	traindata4readshift[,c(varstouse,'length','dist')]
			data=	allset%>%select(everything(),dist,length)%>%mutate(length=as_factor(length)),
			probability=FALSE,num.threads=20
	)
	# seqshiftmodel_allvars

	# seqshiftmodel_GAonly <- ranger::ranger(formula= dist ~ . ,
	# 		data=	traindata4readshift[,c(varstouse,'length','dist')]
	# 		data=	allset%>%select(matches('[ft]p\\.\\d\\.[GA]'),dist,length),probability=TRUE,num.threads=20
	# )
	seqshiftmodel_GAonly<-NULL
	
	list(seqshiftmodel,seqshiftmodel_allvars,	seqshiftmodel_GAonly)

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

get_periodicity_scores<-function(cdstosamp,cdsread_trmap,seqshiftisneg=FALSE){
	message('.')

	assert_that(all(has_name(mcols(cdsread_trmap),c('seqshift','cdsshift'))))
	psites2use<-cdsread_trmap%>%
		keepSeqlevels(cdstosamp,'coarse')%>%
		resize(1,'start')%>%
		shift(.$cdsshift)

	cds2use <- topcdsmap[cdstosamp]

	cds2use %<>% clip_start(30)%>%clip_end(30)

	periodicfrac<-cdscov(psites2use,cds2use)%>%get3bpfrac

	data_frame(
		no_seqshift = psites2use%>%cdscov(cds2use)%>%get3bpfrac,
		with_seqshift = psites2use %>% shift(ifelse(seqshiftisneg,-1,1) * .$seqshift)%>%cdscov(cds2use)%>%get3bpfrac
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

qwidth(topcdsreadsbak)%>%table%>%sort%>%cumsum%>%{./last(.)}%>%`>`(0.05)
0

compartments <- get_compartments(cds,DEFAULT_CIRC_SEQS)

topcdsexons$compartment <- get_transcript_compartments(topcdsexons,compartments)

# cdsread_trmap <- get_cdsread_trmap(topcdsreads,topcdsexons,startcods_trmap)

# cdsread_trmap <- mymemoise(get_cdsread_trmap)(topcdsreads,topcdsexons,topcdsmap)
# cdsread_trmap %<>% setpos


topExonseq <- mymemoise(getSeq)(x=FaFile(REF),topcdsexons)%>%
	split(topcdsexons$transcript_id)%>%
	lapply(.%>%unlist)%>%DNAStringSet



# # # cdsread_trmap<-oldenv$cdsread_trmap%>%{mcols(.)$seqshift=NULL;mcols(.)$startdist=NULL;mcols(.)$stopdist=NULL;.}


#get the cds offsets, and export them
#okay so the get offsets procedure now seems to work fine when used with the old cdsread_trmap procedure
# c(allbestscores,offset_cds_scores) %<-% mymemoise(get_offsets)(cdsread_trmap,compartments)


allbestscores <- read_tsv('ext_data/offsets_manual.tsv')
if(USEPHASE){
	bestscores = allbestscores%>%rowwise%>%mutate(phase = list(0:2))%>%unnest(phase)
}

# outfolder%>%dir.create(showWarnings=F,rec=T)
# allbestscores%>%write_tsv(file.path(outfolder,'cdsmax_offsets.tsv'))

# if(!USEPHASE){
# 	bestscores<-allbestscores%>%filter(is.na(phase))%>%ungroup%>%select(-phase)
# }else{
# 	bestscores<-allbestscores%>%filter(!is.na(phase))
# }

if(USERIBOSEQC){
	#Get riboqc_cutoffs
	sample = bam%>%dirname%>%basename 
	riboqccutoffs <- str_interp(here('pipeline/riboqc/data/${sample}/_P_sites_calcs'))%T>%{stopifnot(file.exists(.))}
	riboqcdf <- riboqccutoffs%>%fread%>%select(compartment=comp,length=read_length,offset=cutoff)
	bestscores<-riboqcdf%>%filter(length %in% allbestscores$length)
}

metainfo = read_tsv('data/metainfo.tsv')
cdsspl <- cds%>%split(.,.$protein_id)
bestcds = metainfo%>%filter(isbest)%>%.$protein_id%>%cdsspl[.]
# #get the cds shfit offsets for our data
# cdsread_trmap$cdsshift <- apply_cds_offsets(cdsread_trmap,bestscores)
# #get only those with a define shift, and for which this doesn't put them over the edge of tr
# cdsread_trmap <- cdsread_trmap %>% subset(!is.na(.$cdsshift)) %>% .[!is_out_of_bounds(.)]

# topstoppsites <- get_stop_sites(cdsread_trmap,topcdsmap,STOPWINDOWSTART,STOPWINDOWEND)

alltrainreads <- mymemoise(get_offset_training_sites)(bam,bestcds)
# assert_that(allstoppsites%>%length%>%`-`(336076),msg='number of start/stop sites to train with change')
alltrainreads%<>%subset(length %in% bestscores$length)
alltrainreads$dist%<>%as.factor

# (function(){

if('chr1' %in% seqnames(alltrainreads)){
 	traindata4readshift <- get_seqforrest_traindata(alltrainreads,FaFile(REF),nbp=2)
}else{
	#in case I go back to using the old sites...
	traindata4readshift <- get_seqforrest_traindata(alltrainsites,topExonseq,nbp=2)
}

get_seqoffset_model_elim <- mymemoise(get_seqoffset_model_elim)

c(seqshiftmodel,seqshiftmodel_allvars,	seqshiftmodel_GAonly) %<-% get_seqoffset_model_elim(
	traindata4readshift,
	usephase=USEPHASE,
	doelim=FALSE)

stop()

predreadshift <- get_predictedshifts(seqshiftmodel_allvars,traindata4readshift)

nf2n <- .%>%as.character%>%as.numeric

txtplot(table(nf2n(traindata4readshift$dist) - predreadshift))

table(nf2n(traindata4readshift$dist) - predreadshift)%>%.[which.max(.)]




# #we can't get sequence if their flanks extend out of the genome/tr
# cdsread_trmap_inbound <- cdsread_trmap[
# 	cdsread_trmap%>%resize(width(.)+4,'center')%>%is_out_of_bounds(.)%>%not
# ]

# data4readshift <- get_seqforrest_traindata(cdsread_trmap_inbound,topExonseq,topcdsmap)

# #it's the model thats actually sure
# # cdsread_trmap_inbound$seqshift <-  get_predictedshifts(seqshiftmodel,data4readshift)
# # cdsread_trmap_inbound$seqshift <-  get_predictedshifts(seqshiftmodel_allvars,data4readshift,prob=0.25)

# cdsread_trmap_inbound$seqshift <-  get_predictedshifts(seqshiftmodel_allvars,data4readshift,prob=0.25)

# Psite_model$new%>%args

psite_model <- Psite_model$new(bestscores,seqshiftmodel_allvars,compartments,referencefasta=FaFile(REF))
psite_model%>%saveRDS(file.path(outfolder,'seqshiftmodel.rds'))
# save.image()
# save.image()

{

is3nt = bestcds%>%width%>%sum%>%mod(3)%>%`==`(0)

cds4test = bestcds[is3nt]%>%head(100)

bcdspsites = psites%>%subsetByOverlaps(cds4test)

bestcdsreadsbak  <- bam%>%BamFile(yieldSize=NA)%>%readGAlignments(param=ScanBamParam(which=cds4test%>%unlist))
bestcdsreads=bestcdsreadsbak
bestcdsreads = alltrainreads

ov = findOverlaps(resize(as(bestcdsreads,'GRanges'),1),bestcds%>%sort_grl_st%>%resize_grl(.,sum(width(.))+42,'end'),select='arbitrary')
ov = findOverlaps(resize(as(bestcdsreads,'GRanges'),1),bestcds,select='arbitrary')
# ov = ov%>%as.data.table%>%group_by(queryHits)%>%sample_n(1)
ov=ov%>%as.data.table
isneg = (strand(bestcdsreads)=='-')%>%as.vector
phaseoffsets = unlist(bestcds)$phase[ov[[1]]]
mcols(bestcdsreads)$phase=NA
mcols(bestcdsreads)$phase[!isneg] = ((start(bestcdsreads)[!isneg] - start(unlist(bestcds))[ov[[1]][!isneg]])+phaseoffsets[!isneg])%%3
mcols(bestcdsreads)$phase[isneg] = (end(unlist(bestcds))[ov[[1]][isneg]] - end(bestcdsreads)[isneg] +phaseoffsets[isneg])%%3
mcols(bestcdsreads)$phase%>%table
mcols(bestcdsreads)$phase%>%is.na%>%table%>%{./sum(.)}
}

if(!is(bestcdsreads,"GRanges")) mcols(bestcdsreads)$length <- qwidth(bestcdsreads)
offsets = psite_model$get_offsets(bestcdsreads%>%subset(!is.na(phase))%>%subset(length%in%traindata4readshift$length)%>%head(100000))

seqadjustedpsites = psite_model$get_psites(bestcdsreads%>%subsetByOverlaps(cds4test)%>%subset(length%in%traindata4readshift$length))

# seqadjustedpsites = psite_model$get_psites(alltrainreads)
seqadjustedpsites

cds4test = cds4test%>%sort_grl_st%>%resize_grl(sum(width(.))+3,'end')

seqadjustedpsites%>%countOverlaps(cdswinds)
seqadjustedpsites%>%mapToTranscripts(cdswinds)



#not getting any overlap here
psitecovseqadj <- seqadjustedpsites%>%as("GRanges")%>%mapToTranscripts(cdswinds)%>%
	        coverage

psitecovseqadj%>%.[sum(runLength(.))>100]%>%lapply(head,30)%>%sapply(as.vector)%>%log1p(.)%>%rowMeans%>%txtplot
psitecovseqadj%>%.[sum(runLength(.))>100]%>%lapply(tail,30)%>%sapply(as.vector)%>%log1p(.)%>%rowMeans%>%txtplot
psitecov%>%.[sum(runLength(.))>100]%>%lapply(tail,30)%>%sapply(as.vector)%>%log1p(.)%>%rowMeans%>%txtplot

psitecovseqadj%>%head(4000)%>%.[sum(runLength(.))>100]%>%lapply(head,30)%>%sapply(as.vector)%>%log1p(.)%>%rowMeans%>%txtplot
psitecov%>%head(4000)%>%.[sum(runLength(.))>100]%>%lapply(head,30)%>%sapply(as.vector)%>%log1p(.)%>%rowMeans%>%txtplot

get_sitedf<-function(testpsitecovs,codposdfspl,stop_codons=stopcodons){
		sitedf <- testpsitecovs%>%lapply(as.vector)%>%stack%>%set_colnames(c('count','gene'))
		
		testpsitecodons = codposdfspl%>%.[names(testpsitecovs)]%>%bind_rows

		sitedf%<>%group_by(gene)%>%mutate(phase=as_factor(((1:length(count))-1)%%3) )
		sitedf$codon=NA
		sitedf$codon[seq(1,nrow(sitedf),by=3)] = testpsitecodons$x
		sitedf$codon[seq(2,nrow(sitedf),by=3)] = testpsitecodons$x
		sitedf$codon[seq(3,nrow(sitedf),by=3)] = testpsitecodons$x
		acod = testpsitecodons%>%group_by(protein_id)%>%mutate(acod=lead(x))%>%.$acod
		sitedf$a_codon = NA
		sitedf$a_codon[seq(1,nrow(sitedf),by=3)] = acod
		sitedf$a_codon[seq(2,nrow(sitedf),by=3)] = acod
		sitedf$a_codon[seq(3,nrow(sitedf),by=3)] = acod
		sitedf <- sitedf%>%filter(!codon %in% stopcodons,!a_codon %in% stop_codons)
		sitedf
}
codposdf<-readRDS('data/codposdf.rds')	

codposdfspl=codposdf%>%split(.,.$protein_id)
stopcodons <-  (GENETIC_CODE=='*')%>%.[.]%>%names
sitedf = get_sitedf(psitecovseqadj,codposdfspl)
sitedfnoadj = get_sitedf(psitecov[names(psitecovseqadj)],codposdfspl)

txtplot(
	sitedfnoadj%>%group_by(codon)%>%summarise(mcount=sum(count))%>%.$mcount,
	sitedf%>%group_by(codon)%>%summarise(mcount=sum(count))%>%.$mcount
)

sitedf%>%group_by(a_codon)%>%summarise(mcount=sum(count))%>%filter(mcount>7.5e3)

library(MatrixModels)
library(MASS)
mtheta=0.32

glmfit = glm4(count ~ 0 + gene + codon +a_codon + phase, data=sitedf,family=negative.binomial(mtheta),MXITER=400,doFit=T, sparse=T, verbose=T)
glmfitnoadj = glm4(count ~ 0 + gene + codon + a_codon + phase, data=sitedfnoadj,family=negative.binomial(mtheta),MXITER=400,doFit=T, sparse=T, verbose=T)


glmfit_p_codon <-  list(seqadj = glmfit,noseqadj = glmfitnoadj)%>%
	map_df(.id='sample', ~ coef(.)%>%enframe)%>%
	filter(name%>%str_detect('^codon'))%>%
	mutate(p_codon = name%>%str_replace('codon',''))%>%
	filter(translate(DNAStringSet(p_codon))!='*')%>%
	select(sample,p_codon,codon_dt_glm=value)

glmfit_p_codon%>%spread(sample,codon_dt_glm)%>%{txtplot(.[[2]],.[[3]])}



# ################################################################################
# ########Test the quality of the model
# ################################################################################
	


# 	#' ## Testing
# 	#' 
# 	#' We can see if our procedure worked by quickly looking at the fft (fast fourier transform) on the top1k genes
# 	#' binned randomly into 5 groups, before and after the application of our sequence specific cutoffs (to ALL reads),
# 	#' rather than just those around the stop codon

# 	seqshift_periodicities<- seqnames(cdsread_trmap_inbound)%>%
# 		unique%>%
# 		sample%>%
# 		# head(300)%>%
# 		as.character%>%
# 		split(seq_along(.)%%5)%>%
# 		mclapply(F=get_periodicity_scores,cdsread_trmap_inbound,seqshiftisneg=FALSE)%>%
# 		bind_rows

# 	seqshift_periodicities%>%unlist%>%txtplot

# 	stop()







# # })()

# # stop()


# #+ spectral_coefficient_strip_plot, fig.width =4,fig.height=4,out.width=400,out.height=450,dev='pdf',include=TRUE,eval=TRUE

# # {
# # 	seqshift_periodicities%>%
# # 		gather(set,spectral_coefficient)%>%
# # 		qplot(data=.,x=set,color=set,y=spectral_coefficient,geom=c('point'))+theme_bw()
# # }

# {

# setpos<- .%>%{strand(.)<-'+';.}
# prestops <- topcdsmap%>%setpos%>%resize(3,'end')%>%resize(1,'start')
# prestarts <- topcdsmap%>%setpos%>%resize(3,'start')%>%resize(1,'start')
# cdsread_trmap$startdist<-start(cdsread_trmap) - start(prestarts[seqnames(cdsread_trmap)])
# cdsread_trmap$stopdist<-start(cdsread_trmap) - start(prestops[seqnames(cdsread_trmap)])

# #Get riboqc_cutoffs
# sample = bam%>%dirname%>%basename 
# riboqccutoffs <- str_interp(here('pipeline/riboqc/data/${sample}/_P_sites_calcs'))%T>%{stopifnot(file.exists(.))}
# riboqcdf <- riboqccutoffs%>%fread%>%select(length=read_length,riboqc_shift=cutoff)

# starts <- topcdsmap%>%setpos%>%resize(1,'start')%>%setNames(seqnames(.))
# cdsread_trmap_inbound$startdist <- start(cdsread_trmap_inbound) - start(starts[cdsread_trmap_inbound@seqnames])

# stops <- topcdsmap%>%setpos%>%resize(3,'end')%>%fpend%>%setNames(seqnames(.))
# cdsread_trmap_inbound$stopdist <- start(cdsread_trmap_inbound) - start(stops[cdsread_trmap_inbound@seqnames])


# fcods=20
# lflank=(fcods*3)
# rflank=(fcods*3)+2
# mpoint=((fcods*2)+1)*3
# epoint=((fcods*4)+2)*3

# # cdsread_trmap%>%as.data.frame%>%mutate(d=startdist+stopdist)%>%.$d%>%hist
# # cdsread_trmap%>%as.data.frame%>%mutate(d=startdist+stopdist)%>%.$d%>%between(-lflank,rflank)%>%table
# metaplotlabs<-c(paste0('-',lflank/2),'AUG',lflank/2,paste0('mid -',lflank/2),'mid',paste0('mid +',lflank/2),paste0('end -',lflank/2),'stop',paste0('end +',lflank/2))
# metaplotbreaks<-c(-lflank/2,0,lflank/2,mpoint-lflank/2,mpoint,mpoint+lflank/2,epoint-(lflank/2),epoint,epoint+(lflank/2))
# #

# seqshiftfuncs <- list(
# 	Riboqc = .%>% safe_left_join(riboqcdf,by=c('length'))%>% mutate(startdist=startdist+riboqc_shift,stopdist=stopdist+riboqc_shift),
# 	CDSmax_shift = .%>% mutate(startdist=startdist+cdsshift,stopdist=stopdist+cdsshift),
# 	seqshift = .%>% mutate(startdist=startdist+cdsshift,stopdist=stopdist+cdsshift)%>%mutate(startdist=startdist-seqshift,stopdist=stopdist-seqshift)
# )

# #

# #' We can also simply look at a metaplot of the frames, (comparing riboqc and cdsmax shift)
# plist=lapply(seqshiftfuncs,mymemoise(function(seqshiftfunc){
# 	#
# 	p=cdsread_trmap_inbound%>%
# 		# sample(100e3)%>%
# 		# shift(.$cdsshift)%>%
# 		mcols()%>%as_tibble	%>%
# 		select(startdist,stopdist,length,cdsshift,seqshift)%>%
# 		# filter(between(abs(startdist-stopdist),-lflank,rflank))%>%
# 		seqshiftfunc%>%
# 		mutate(
# 			mdist = startdist - ((floor(((startdist-stopdist+1)/3)/2))*3),
# 			instart=between(startdist,-lflank,rflank),
# 			inend=between(stopdist,-lflank,rflank),
# 			inmiddle=between(mdist,-lflank,rflank)
# 		)%>%
# 		filter((instart|inend|inmiddle) & (!(instart&inend&inmiddle)))%>%
# 		# filter(instart)%>%
# 		mutate(x = NA,
# 			x=ifelse(instart,startdist,x),
# 			x=ifelse(inmiddle,mdist+mpoint,x),
# 			x=ifelse(inend,stopdist+epoint,x) ) %>%
# 		mutate(length=as.character(length))%>%
# 		bind_rows(.,mutate(.,length='all'))%>%
# 		group_by(x,length)%>%tally%>%
# 		# filter(length=='all')%>%
# 		mutate(phase = factor(x%%3))%>%
# 		# arrange(desc(n))///
# 		ggplot(aes(x=x,ymax=n,ymin=0,color=phase))+geom_linerange()+
# 		facet_grid(scale='free',length~.)+
# 		geom_vline(linetype=2,xintercept=c(0,epoint+2),alpha=I(0.3))+
# 		scale_x_continuous(labels=metaplotlabs,breaks=metaplotbreaks)+
# 		# scale_y_continuous(limits=c(0,300))+
# 		# coord_cartesian(xlim=c(-lflank,rflank))+
# 		geom_rect(aes(xmin=0,xmax=epoint,ymin=0,ymax=Inf,alpha=I(0.01)),fill='grey',color='white')+
# 		scale_fill_discrete(guide=F)+
# 		theme_minimal()	
# 	p
# }))

# maxn<-plist%>%map_dbl(~max(.$data$n))%>%max
# for(i in seq_along(plist)) plist[[i]] <- plist[[i]]+scale_y_continuous(limits=c(0,maxn/2))

# #+ riboshift_comparison plot, fig.width =12,fig.height=12,out.width=1200,out.height=1250,dev='pdf',include=TRUE,eval=TRUE

# if(isknitr) {
# 	ggpubr::ggarrange(ncol=3,plotlist=plist,labels=names(seqshiftfuncs))
# }else{
# 	suppressWarnings(dir.create(outfolder,rec=TRUE))
# 	shiftplotfile <- 'shift_methods_comp.pdf'
# 	shiftplotfile <- file.path(outfolder,shiftplotfile)
# 	pdf(shiftplotfile,w=14,h=10)
# 	print(ggpubr::ggarrange(ncol=3,plotlist=plist,labels=names(seqshiftfuncs)))
# 	dev.off()
# 	message(normalizePath(shiftplotfile))
# }


# }







# #







# data4readshift%>%as.data.frame%>%select(matches('^[ft]p'))%>%
# 	apply(2,table)%>%.['1',TRUE]%>%
# 	enframe%>%separate(name,into=c('end','pos','base'),sep='\\.')%>%
# 	group_by(end,pos)%>%mutate(value = value/sum(value))%>%
# 	spread(base,value)





# #Now let's check up/down stream of codons





# shiftplotfile <- 'shift_methods_comp.pdf'
# shiftplotfile <- file.path(outfolder,shiftplotfile)
# pdf(shiftplotfile,w=14,h=10)
# print(ggpubr::ggarrange(ncol=3,plotlist=plist,labels=names(seqshiftfuncs)))
# dev.off()
# message(normalizePath(shiftplotfile))


# library(Biostrings)

# cdsseq <- topExonseq[topcdsmap]

# for(codon in names(GENETIC_CODE)){
# 	codonlocs <- vmatchPattern(codon,cdsseq)
# 	codonmatchgr<-unlist(codonlocs)%>%{GRanges(names(.),.)}%>%setNames(NULL)

# 	tmp1<-codonmatchgr%>%resize(30,'end')%>%mapFromTranscripts(topcdsmap)%>%head
# 	seqnames(tmp1)%in%names(topExonseq)
# 	seqinfo(tmp1)
# 	topExonseq[tmp1[,NULL]]

# 	# vmatchPDict(names(GENETIC_CODE)[1:2],cdsseq[1:10])
# }

# cdsread_trmap%>%.$length%>%table

# makeTxDbFromGRanges(gtf_gr)

	
# forget(readGAlignments)
# bam%>%file.exists
# bamfile<-BamFile(bam,yield=)



# ##Let's make fake, positive control reads

# #for each of our 

# allstoppsites%>%head%>%as.data.frame



# #Simulate Riboseq reads
# #length distribution
# lengthdist <- c(`26` = 51599L, `27` = 82461L, `28` = 57037L, `29` = 43928L)%>%{./sum(.)}
# #initial_frame_dist 
# initial_frame_dist <- c(0.6,0.3,0.1)
# testexons <- topcdsexons[1:100]
# testcds <- topcdsmap[topcdsexons$transcript_id%>%unique]

# sample<-base::sample

# starts_ends <- testcds%>%head(2)%>%width%>%divide_by(3)%>%map(seq_len)%>%map(sample,1e3,replace=TRUE)%>%
# 	map(~ .+ sample(c(0,1,2),p=initial_frame_dist,rep=TRUE,size=length(.)))%>%
# 	map(~data_frame(start=.,end= . + as.numeric(names(lengthdist))[sample(seq_along(lengthdist),rep=TRUE,size=length(.))]))%>%


# ####Look more modelishly at reads

# ################################################################################
# ######### Distribution around codons
# ################################################################################
	

# prestops <- topcdsmap%>%resize(3,'end')%>%resize(1,'start')
# stopreads <- subsetByOverlaps(cdsread_trmap,prestops)
# leftseg <- fpend(stopreads)%>%downstream_dist_till(prestops)
# rightseg <- tpend(stopreads)%>%upstream_dist_till(prestops)

# table(stopreads$length)

# #so yeah very much left seg influence by right....
# table(leftseg,rightseg)%>%apply(1,function(x) x / sum(x))
# table(leftseg,rightseg)%>%apply(1,function(x) x / sum(x))%>%.[c(7,8,9),]


# codons <- names(GENETIC_CODE)

# codondisttabs<-lapply(codons,function(codon){

# 	codonlocs <- topExonseq[topcdsmap%>%setpos]%>%
# 		# vmatchPattern(codon,.)%>%
# 		vmatchPattern(codon,.)%>%
# 		unlist%>%
# 		subset((end%%3)==0)%>%
# 		GRanges(names(.),.)

# 	codonlocs%<>%mapFromTranscripts(topcdsmap)

	

# 	lendisttab<-codonlocs%>%mergeByOverlaps(cdsread_trmap)%>%{tibble(dist=start(.$.)-start(.$cdsread_trmap),length=.$length)}

# 	lendisttab
# })

# codondisttabs%<>%setNames(codons)
# codondisttabs%<>%bind_rows(.id='codon')

# codistplot <- here('plots/codondistplots/codondistbycodon.pdf')
# codistplot%>%dirname%>%dir.create(showWarnings=FALSE)
# pdf(codistplot)
# codondisttabs%>%group_by(codon,dist,length)%>%tally%>%
# 	filter(dist>5)%>%filter(dist<length-5)%>%
# 	ungroup%>%
# 	mutate(dist=floor(dist/3))%>%
# 	group_by(codon,length,dist)%>%summarise(n=sum(n))%>%
# 	group_by(codon,length)%>%
# 	mutate(n=n/min(n))%>%
# 	ggplot(aes(y=n,x=dist,color=codon))+facet_grid(scale='free',length~.)+geom_line()+
# 	scale_x_continuous(breaks=1:20)+
# 	scale_color_discrete(guide=FALSE)+
# 	theme_bw()
# dev.off()
# normalizePath(codistplot)%>%message

# codistplot <- here('plots/codondistplots/codondist.pdf')
# pdf(codistplot)
# codondisttabs%>%group_by(codon,dist,length)%>%tally%>%
# 	filter(dist>5)%>%filter(dist<length-5)%>%
# 	group_by(codon,length)%>%
# 	mutate(n=n/min(n))%>%
# 	ggplot(aes(y=n,x=dist,color=codon))+facet_grid(scale='free',length~.)+geom_line()+
# 	scale_x_continuous(breaks=seq(1,30,3),minor_breaks=1:30,name="Codon bp1 - 5' Read End")+
# 		scale_color_discrete(guide=FALSE)+
# 	theme_minimal()+
# 	theme(panel.grid.minor = element_line( size=0.5))
# dev.off()
# normalizePath(codistplot)%>%message



# codonsmagvect<-codondisttabs%>%group_by(codon)%>%summarise(signal=sum(between(dist,6,16)))%>%arrange(desc(signal))%>%.$codon


# stopcodons<-codons%>%DNAStringSet%>%translate%>%as.character%>%is_in(c('*'))%>%codons[.]
# stopcodonsoratg<-c('ATG',stopcodons)

# codistplot <- here('plots/codondistplots/codondist.pdf')
# pdf(codistplot)
# codondisttabs%>%group_by(codon,dist,length)%>%tally%>%
# 	filter(codon %in% c(codonsmagvect[1:4],rev(codonsmagvect)[1:4]))%>%
# 	filter(!codon %in% stopcodonsoratg)%>%
# 	filter(dist>3)%>%filter(dist<length-3)%>%
# 	group_by(codon,length)%>%
# 	mutate(n=n/min(n))%>%
# 	ggplot(aes(y=n,x=dist,color=codon))+facet_grid(scale='free',length~.)+geom_line()+
# 	scale_x_continuous(breaks=seq(1,30,3),minor_breaks=1:30,name="Codon bp1 - 5' Read End")+
# 		scale_color_discrete()+
# 	theme_minimal()+
# 	theme(panel.grid.minor = element_line( size=0.5))
# dev.off()
# normalizePath(codistplot)%>%message


# codondisttabs[[9]]%>%.$dist%>%table%>%
# 	.[order(as.numeric(names(.)))]%>%txtplot(y=.,x=as.numeric(names(.)))


# codondisttabs%>%setNames(codons)%>%bind_rows(.id='codon')
# 	filter(dist==5,length==29)%>%group_by(codon,dist)%>%tally

# codonlocs%>%mergeByOverlaps(cdsread_trmap)%>%{tibble(dist=start(.$.)-start(.$cdsread_trmap),length=.$length)}%>%table%>%
# 	.[order(as.numeric(names(.)))]%>%txtplot(y=.,x=as.numeric(names(.)))




#TODO - test this
#' Why am I not getting great periodicity?
#' What am I not getting 

#length
#total psite count
#psite count at the edges






#TODO - count just outside the CDS here - both to get our edge scores and to allow access to seqs

#and in the middle

#Now derive various states for our ORFs.



#' .... All of which seems to indicate our procedure is doing a better job of recovering actual P-site positions. Note that it's p
#' site positions are consistent between reads - no dramatic difference between 28 and the others - and manage to recover periodicity
#' in specific read lengths where it was missing. 
#' 
#' Okay so the strand flip error in the stoppsites function is fixed now. We have consstent improvement, however the behavior of the cds shift function is somewhat erratic...
#' 
#' Okay, found problem, frozen randomeness in the model training due to the caching... take NOTE!
#' 
#' Now I"m fucked again, can't figure out why, probability trees apparently hurt??

 # if((!isknitr) & (interactive())) {'/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/R/offsets3.R'%>% {rmarkdown::render(knitr::spin(knit=F,.),output_file=here('reports',basename(str_replace(.,'.R$','.html'))))} } 

 # inenv<-new.env() testout<-load('pipeline/riboqc/data/RPI12_PolyP0_2/_for_SaTAnn')