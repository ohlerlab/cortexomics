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
USERIBOSEQC <- FALSE
source(here('src/R/Rprofile.R'))



for(fname in lsf.str('package:GenomicRanges')) assign(fname,get(fname,'package:GenomicRanges'))
for(fname in lsf.str('package:dplyr')) assign(fname,get(fname,'package:dplyr'))


argv <- c(
	bam = here('pipeline/star/data/E13_ribo_2/E13_ribo_2.bam')%T>%{stopifnot(file.exists(.))},
	gtf = here('pipeline/my_gencode.vM12.annotation.gtf'),
	REF = here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'),
	outfolder = here('pipeline/seqshift_reads/data/E13_ribo_2/')
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




get_cdsread_trmap<-memoise(get_cdsread_trmap)




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




#offsetreads[c(1,length(offsetreads))]%>%apply_psite_offset(12)

get_offsets<-function(offsetreads,compartments,cdstotest = topcds ){

	readsizes <- select_toplengths(qwidth(offsetreads))

	possibleoffsets <-seq(3,max(readsizes),by=1)%>%setNames(.,.)
	# possibleoffsets <-seq(6,max(readsizes),by=3)%>%setNames(.,.)

		
	if(is.null(mcols(offsetreads)$compartment)){
		mcols(offsetreads)$compartment <- compartments%>%as('SimpleList')%>%.[seqnames(offsetreads)]%>%unlist
	}
	readcompartments <- unique(mcols(offsetreads)$compartment)

	phaselist<-list(0,1,2,0:2)%>%setNames(c(0,1,2,'all'))

	nophase = (!'phase' %in% colnames(mcols(offsetreads))) | (n_distinct(mcols(offsetreads)$phase)==1)

	if(nophase){
		mcols(offsetreads)$phase = NA
		phaselist = list('all'=NA)
		message('Phase info not in reads, doing phaseless offsets')
	}

	if(!'length' %in% colnames(mcols(offsetreads))){
		mcols(offsetreads)$length = qwidth(offsetreads)
	}

	offset_cds_scores <-	mclapply(possibleoffsets,function(offset){
			lapply(unique(compartments)%>%intersect(readcompartments),function(compartment_i){
				lapply(readsizes,function(length_i){
					lapply(phaselist,function(phase_i){
						if((length_i - 6 - offset) < 0) return(NULL)
						
						cat('.')

						score <- offsetreads%>%
							subset(compartment==compartment_i)%>%
							subset(phase%in%phase_i)%>%
							subset(length==length_i)%>%
							apply_psite_offset(offset)%>%
							countOverlaps(cdstotest)%>%
							`>`(0)%>%sum

						if(score!=0) 
						stopifnot(score!=0)
					
						data.table(score)

					})%>%bind_rows(.id='phase')
				})%>%bind_rows(.id='length')
			})%>%bind_rows(.id='compartment')
		})%>%bind_rows(.id='offset')

	offset_cds_scores%<>% filter( ((as.numeric(offset) %%3) ==0) | (phase=='all') )
	offset_cds_scores%<>%mutate_at(vars(everything()),as.numeric)
	offset_cds_scores$compartment<-unique(compartments)[offset_cds_scores$compartment]

	bestscores<-offset_cds_scores%>%group_by(length,phase,compartment)%>%slice(which.max(score))

	if(any(bestscores$score==0)) browser()

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


# reads_tr<-topcdsreads




shorten_ranges<-function(gr){gr <- shift(gr,- min(start(gr))) ;gr}

# cdsbak<-cds
# cds<-cdsbak%>%head(100e3)


# cds<-c(
# cds%>%subset(strand=='+')%>%split(.,.$transcript_id)%>%.[1],
# cds%>%subset(strand=='-')%>%split(.,.$transcript_id)%>%.[1]
# )
# cds%<>%unlist

# reads_tr<-topcdsreads
# topcdsmap<-topcds

#This should be modified to take either a fastq file or the exonsequence
#And it should be okay with missing sequences

setMethod('seqinfo','DNAStringSet',function(x){
	x%>%{Seqinfo(names(.),nchar(.))}
})

#trainreads<-topcdsreads
#seq=FaFile(REF)
#nbp=2
setMethods('resize',c('GAlignments','numeric','character'),function(x,width,where='start',...){

	if(where=='start'){
		return(qnarrow(x,start=1,width=width,...))
	}else if(where=='end'){
		return(qnarrow(x,end=qwidth(reads),width=width,...))
	}else{
		stop('dont know how to centre on that')
	}
})



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
get_seqshift_trainsites <- function(bam,allcds,STOPWINDOWSTART,STOPWINDOWEND,bestscores,startonly=TRUE){
	# browser()
	message('getting stop sites for training seq model - all stop sites')

	#work out start and stop codons from the cds


	startcodons <- allcds[isfpmost(allcds,'transcript_id')]%>%resize(3,'start')

	starts_ends<-c(startcodons)
	
	assert_that(startcodons%>%strandshift(0)%>%getSeq(FaFile(REF),.)%>%translate%>%`==`('M')%>%mean%T>%message%>%`>`(0.8))

	reads <- bam%>%BamFile(yieldSize=NA)%>%readGAlignments(param=ScanBamParam(which=starts_ends))
	reads <- GRanges(reads,length=qwidth(reads))
	reads <- reads[overlapsAny(resize(reads,width(reads)-3,'center'),starts_ends)]

	mcols(reads)$cdsshift <- reads%>%get_cds_offsets(bestscores,compartments)
	reads%<>%subset(!is.na(cdsshift))
	psites <- reads%>%resize(1)%>% strandshift(mcols(.)$cdsshift)
	
	names(startcodons)<-paste0('startcodon_',seq_along(startcodons))

	#map tot he codon (plus the 2bp other side of it)
	startpsites <- psites%>%mapToTranscripts(startcodons%>%resize(1,'start')%>%resize(1+STOPWINDOWEND-STOPWINDOWSTART,'center'))

	startstoppsites <- 	c(startpsites)		
	
	startstoppsites%<>%sample%>%.[!duplicated(.$xHits)]

 	reads <- reads[startstoppsites$xHits]
 	reads$dist <- as.factor(1-STOPWINDOWSTART-start(startstoppsites))



	c(reads)

}

get_seqoffset_model_elim <-  function(traindata4readshift,usephase=F,doelim=T){


	require(ranger)

	stopifnot((-2:2) %in% traindata4readshift$dist)

	traindata4readshift$dist%<>%as.factor
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

	if(doelim){while(length(varstouse)>0){

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
	}}

	#So if we use 10k examples we can trim off the less important variables pretty easily.
	# keepvars <- (tabs%>%unlist < (tabs%>%unlist%>%min)+0.05)%>%.[!.]%>%names

	#use ones we didn't eliminate
	keepvars <- colnames(traindata4readshift)%>%str_subset(neg=TRUE,'length|dist|phase')%>%setdiff(names(tabs))

	if(usephase) keepvars %<>% append('phase')

	seqshiftmodel  <- ranger::ranger(formula= dist ~ . ,
			# data=	traindata4readshift[,c(varstouse,'length','dist')]
			data=	allset%>%select(keepvars,dist,length),num.threads=20
	)
	
	allseqvars <-  colnames(allset)%>%str_subset(neg=TRUE,'length|dist|phase')

	seqshiftmodel_allvars <- ranger::ranger(formula= dist ~ . ,
			# data=	traindata4readshift[,c(varstouse,'length','dist')]
			data=	allset%>%select(everything(),dist,length),
			probability=TRUE,num.threads=20
	)

	seqshiftmodel_allvars

	seqshiftmodel_GAonly <- ranger::ranger(formula= dist ~ . ,
			# data=	traindata4readshift[,c(varstouse,'length','dist')]
			data=	allset%>%select(matches('[ft]p\\.\\d\\.[GA]'),dist,length),num.threads=20
	)
	
	list(seqshiftmodel,seqshiftmodel_allvars,	seqshiftmodel_GAonly)

}




get_seqoffset_model_elim_rl <- function(traindata4readshift,...){
	lens <- traindata4readshift$length%>%unique%>%setNames(.,.)
	stopifnot(length(lens)>0)
	lapply(lens,function(sublen){
		get_seqoffset_model_elim(traindata4readshift[traindata4readshift$length==sublen,],...)
	})
}
# get_seqoffset_model_elim_rl(traindata4readshift)


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

topcdsmap<-topcds%>%pmapToTranscripts(topcdsexons%>%split(.,.$transcript_id)%>%.[topcds$transcript_id])%>%reduce
names(topcdsmap) <- seqnames(topcdsmap)


#' - Map our reads to these transcripts
  
#get reads over them
if(!exists('topcdsreadsbak'))topcdsreadsbak <- bam %>% readtopreads
#TODO
topcdsreads<-topcdsreadsbak

qwidth(topcdsreadsbak)%>%table%>%sort%>%cumsum%>%{./last(.)}%>%`>`(0.05)
0

compartments <- get_compartments(cds,DEFAULT_CIRC_SEQS)

topcdsexons$compartment <- get_transcript_compartments(topcdsexons,compartments)

# cdsread_trmap <- get_cdsread_trmap(topcdsreads,topcdsexons,startcods_trmap)

add_offset_metacols <- function(offsetreads,compartment){	
	if(is.null(mcols(offsetreads)$compartment)){
		mcols(offsetreads)$compartment <- compartments%>%as('SimpleList')%>%.[seqnames(offsetreads)]%>%unlist
	}

	if(!'phase' %in% colnames(mcols(offsetreads))){
		mcols(offsetreads)$phase = NA
		phaselist = list('all'=NA)
		message('Phase info not in reads, doing phaseless offsets')
	}

	if(!'length' %in% colnames(mcols(offsetreads))){
		mcols(offsetreads)$length = qwidth(offsetreads)
	}

	offsetreads
}
topcdsreads <- topcdsreads%>%add_offset_metacols




################################################################################
########Now run the psite determination procedure
################################################################################


#get the cds offsets, and export them
#okay so the get offsets procedure now seems to work fine when used with the old cdsread_trmap procedure
STOPWINDOWSTART = -5
STOPWINDOWEND = 5
LENGTHSUBSET = c(28)

get_offsets <- mymemoise(get_offsets)
c(allbestscores,offset_cds_scores) %<-% get_offsets(topcdsreads,compartments)

outfolder%>%dir.create(showWarnings=F,rec=T)
allbestscores%>%write_tsv(file.path(outfolder,'cdsmax_offsets.tsv'))

bestscores<-allbestscores%>%filter(is.na(phase))%>%ungroup%>%select(-phase)

if(USERIBOSEQC){
	#Get riboqc_cutoffs
	sample = bam%>%dirname%>%basename
	riboqccutoffs <- str_interp(here('pipeline/riboqc/data/${sample}/_P_sites_calcs'))%T>%{stopifnot(file.exists(.))}
	riboqcdf <- riboqccutoffs%>%fread%>%select(compartment=comp,length=read_length,offset=cutoff)
	bestscores<-riboqcdf%>%filter(length %in% allbestscores$length)
}



#get the cds shfit offsets for our data
mcols(topcdsreads)$cdsshift <- get_cds_offsets(topcdsreads,bestscores)
topcdsreads%<>%subset(!is.na(cdsshift))

topcdsreads%<>%subset(mcols(.)$length %in% LENGTHSUBSET)

# #get only those with a define shift, and for which this doesn't put them over the edge of tr
# cdsread_trmap <- cdsread_trmap %>% subset(!is.na(.$cdsshift)) %>% intersect(.,trim(.))

# topstoppsites <- get_stop_sites(topcdsreads,topcds,STOPWINDOWSTART,STOPWINDOWEND)

cds4train <- cds%>%filtercds4train

seqoffseqtrainreads <- mymemoise(get_seqshift_trainsites)(bam,cds4train,STOPWINDOWSTART,STOPWINDOWEND,bestscores)

seqoffseqtrainreads%<>% subset(length %in% LENGTHSUBSET)

library(testthat)



# test_that("the trainin  site function works the way I think I think it does",{

# 	testread1 <- mymemoise(get_seqshift_trainsites)(bam,cds4train%>%subset(protein_id%in%sample(protein_id,100)),STOPWINDOWSTART,STOPWINDOWEND,bestscores,startonly=TRUE)

# 	(testread1$cdsshift+as.numeric(as.character(testread1$dist))) == abs(fp(testread1)-fp(cds4train[resize(testread1,1)%>%precede(cds4train)]))

# })


trainsites <- seqoffseqtrainreads

# assert_that(seqoffseqtrainreads%>%length%>%`-`(336076),msg='number of start/stop sites to train with change')

# toptrainsites$dist%>%table%>%txtplot
#seqoffseqtrainreads$dist%>%table%>%txtplot
trainsites$dist%>%table%>%txtplot


traindata4readshift <- get_seqforrest_data(trainsites,FaFile(REF),nbp=2)

# get_seqoffset_model_elim <- mymemoise(get_seqoffset_model_elim)
get_seqoffset_model_elim <- mymemoise(get_seqoffset_model_elim_rl)

# c(seqshiftmodel,seqshiftmodel_allvars,	seqshiftmodel_GAonly) %<-% get_seqoffset_model_elim(
c(seqshiftmodel,seqshiftmodel_allvars,	seqshiftmodel_GAonly) %<-% get_seqoffset_model_elim_rl(
	traindata4readshift,
	doelim=FALSE)

psite_model <- Psite_model$new(bestscores,seqshiftmodel_allvars,referencefasta=REF,compartments)

psite_model%>%saveRDS(file.path(outfolder,'seqshiftmodel.rds'))
save.image('offsetsgenomicspace.RData')