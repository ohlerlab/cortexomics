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
suppressMessages({library(bamsignals)})
suppressMessages({library(zeallot)})
suppressMessages({library(stringr)})
suppressMessages({library(here)})
suppressMessages({library(assertthat)})
suppressMessages({library(multitaper)})
suppressMessages({library(GenomicAlignments)})
suppressMessages({library(GenomicFeatures)})

library(parallel)

reduce <- GenomicRanges::reduce

MAPQTHRESH <- 200
USEPHASE <- FALSE
USERIBOSEQC <- FALSE
source(here('src/R/Rprofile.R'))



for(fname in lsf.str('package:GenomicRanges')) assign(fname,get(fname,'package:GenomicRanges'))
for(fname in lsf.str('package:dplyr')) assign(fname,get(fname,'package:dplyr'))


argv <- c(
	bam = here('pipeline/star/data/P0_total_1/P0_total_1.bam')%T>%{stopifnot(file.exists(.))},
	gtf = here('pipeline/my_gencode.vM12.annotation.gtf'),
	REF = here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'),
	shiftmodel = 'pipeline/seqshift_reads/data/P0_total_1/seqshiftmodel.rds',
	outfolder = 'riboseq_quant/data/P0_total_1/'
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

segment_ribocount<-function(exonsubset,cdssubset,bam,psite_model,REF,
	STARTCLIP = 15,
	ENDCLIP = 5,
	nbp=2,
	startflank=9,
	endflank=9,
	do_seqshift=TRUE,
	revstrand=FALSE){

	# cat(names(exonsubset)[1])
	
	stopifnot(!is.null(names(exonsubset)))
	assert_that(all(names(exonsubset)==names(cdssubset)))

	exons2use

	#to ensure that our reads don't overlap two ranges
	readwindow = cdssubset%>%unlist%>%reduce%>%c(.,gaps(.)%>%subset(width<400))%>%reduce

	#get only reads which actually overlap cds
	riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=MAPQTHRESH,which=readwindow)

	#read in reas
	reads <- readGAlignments(bam,param=riboparam)
	#invert stran if needed, use only the ones on right strand
	table(strand(reads))
	if(revstrand) reads%<>%invertStrand
	reads%<>%subsetByOverlaps(cdssubset)
	mcols(reads)$length <- qwidth(reads)

	stopifnot(max(psite_model$offsets$length)<100)
	if(!is.null(psite_model)) mcols(reads)$cdsshift <- get_cds_offsets(reads,psite_model$offsets,psite_model$compartments)

	if(do_seqshift){
		seqdata <- get_seqforrest_data(reads,FaFile(REF),nbp=2)

		mcols(topcdsreads)$seqshift <-  get_predictedshifts(psite_model$seqshiftmodel,seqdata,prob=0.5)

		psitesmapped <- reads%>%
					subset(!is.na(cdsshift))%>%
					apply_psite_offset(c('cdsshift','seqshift'))%>%as("GRanges")%>%
					mapToTranscripts(exonsubset)
	}else if(!is.null(psite_model)){
		psitesmapped <- reads%>%
					subset(!is.na(cdsshift))%>%
					apply_psite_offset(c('cdsshift'))%>%as("GRanges")%>%
					mapToTranscripts(exonsubset)
	}else{
		psitesmapped <- reads%>%qnarrow(start=floor(qwidth(.)/2),width=1)%>%as("GRanges")%>%mapToTranscripts(exonsubset)
	}

	cdsmap <- cdssubset%>%unlist%>%pmapToTranscripts(exonsubset[names(.)])%>%reduce

	assert_that(length(cdsmap)==length(exonsubset))

	assert_that(all(table(seqnames(cdsmap))==1),msg='CDS shoudl be one unbroken stretch of the transcript')
	assert_that(all(names(exonsubset) %in% seqlevels(cdsmap)),msg='CDS shoudl be one unbroken stretch of the transcript')

	startclipbp<-	(3*STARTCLIP) + startflank
	stopclipbp<-	(3*ENDCLIP) + endflank

	covvects <- psitesmapped%>%coverage

	expcds <- cdsmap %>% resize(width(.)+startflank+endflank,'center')

	#make sure our expansions worked
	stopifnot(! any (is_out_of_bounds(expcds)))

	seqlengths(expcds) <- sum(width(exonsubset))
		
	# expcds <- trim(expcds)


	covvects <- covvects[expcds]
	covvects <- lapply(covvects,as.vector)
	#Get, periodicity scores
	vectftests<-covvects%>%lapply(ftestvect)%>%simplify2array

	cat('..')
	# message('memory in use ',gigsused())

	data.frame(
		spec_coef = vectftests[1,],
		spec_pval = vectftests[2,],
		length = covvects%>%map_dbl(length),
		total = covvects%>%map_dbl(sum),
		startcount = covvects%>%map_dbl(possibly(~sum(.[1:startclipbp]),NA)),
		enddcount = covvects%>%map_dbl(possibly(~sum(.[(length(.)-stopclipbp+1):length(.)]),NA)),
		centercount = covvects%>%map_dbl(possibly(~sum(.[(startclipbp+1):(length(.)-stopclipbp)]),NA))
	)%>%rownames_to_column('protein_id')

}


################################################################################
	











# cds2use < cds%>%split(.,.$protein_id)%>%sample(50)
cds2use <- cds%>%split(.,.$protein_id)
p2trdf<-with(cds,data_frame(protein_id,transcript_id))%>%distinct
cdstrs<-data_frame(protein_id=names(cds2use))%>%safe_left_join(p2trdf)%>%.$transcript_id
exons2use <- exons%>%split(.,.$transcript_id)%>%.[cdstrs]
names(exons2use)<-names(cds2use)


exons2use <- sort(exons2use)
exons2use <- exons2use%>%resize_grl(9,'end')
exons2use <- exons2use%>%resize_grl(9,'start')



#also expand the exons so they always extend at least flank base pairs up and downstream



# #
# test<-'ENSMUSP00000110057'
# exons2use%<>%.[test]
# cds2use%<>%unlist(use.names=F)%>%subset(protein_id==test)%>%split(.,.$protein_id)



# warning('test mode')
# testp <- 'ENSMUSP00000026425'
# splitids <- test
# splitids = 'ENSMUSP00000000001'
message('producing read counts')

splitids <- cds$protein_id%>%unique%>%sort%>%split(floor((1:length(.))/300))

bams <- argv['bam']

if(interactive())bams <- here('pipeline/star/data/*/*.bam')%>%Sys.glob%>%str_subset(neg=TRUE,'transcript')%>%str_subset('_ribo_|total')
# bams<-bams[3]
for(bam in bams){
	message(bam)
	library(glue)
	thissamp<-bam%>%dirname%>%basename
	shiftmodel = here(glue('pipeline/seqshift_reads/data/{thissamp}/seqshiftmodel.rds'))
	outfolder = here(glue('pipeline/riboseq_quant/data/{thissamp}/'))




	rdsfile<-here(glue('{outfolder}/allgenes.rds'))
	segcoutdffile<-file.path(outfolder,'segment_counts_df.tsv')
	# if(file.exists(segcoutdffile)) next
	# allgenesrds<-str_interp('pipeline/ribosegcount/${thissamp}/allgenes.rds')
	psite_model <- shiftmodel%>%possibly(readRDS,NULL)(.)

	isrev <- fread(here('src/sample_parameter.csv'))%>%filter(sample_id==thissamp)%>%.$protocol%>%`==`('reverse')

	detectCores()
	message(getwd())
	message(outfolder)

	message(paste0('counting ', length(unlist(splitids)),' genes' ))

	# ns_all_genes<-splitids%>%mclapply(mc.cores=detectCores(),function(idchunk){
	# splitids[1]->idchunk
	ns_all_genes<-splitids%>%lapply(function(idchunk){
		# safely(segment_ribocount)(
		segment_ribocount(
			exons2use[idchunk],
			cds2use[idchunk],
			bam,
			psite_model,
			REF,
			STARTCLIP = 15,
			ENDCLIP = 5,
			startflank = 9,
			endflank = 9 ,
			do_seqshift=FALSE,
			revstrand = isrev
		)
	})

	quantfailed <- ns_all_genes%>%map('error')%>%map_lgl(is.null)
	message(paste0(sum(quantfailed), ' of ',length(quantfailed), ' genes batches of genes were quantified successfully'))

	require(glue)
	rdsfile<-here(glue('{outfolder}/allgenes.rds'))

	rdsfile%>%dirname%>%dir.create(rec=TRUE,showWarn=FALSE,.)
	# ns_all_genes<-readRDS(rdsfile)
	message('saving object to')
	message(rdsfile)

	ns_all_genes%>%saveRDS(rdsfile)

	#readRDS(rdsfile)->ns_all_genes


	if('result' %in% names(ns_all_genes[[1]])){
		segcountdf<-ns_all_genes%>%map('result')%>%bind_rows%>%select(protein_id=matches('id'),everything())
	}else{
		segcountdf<-ns_all_genes%>%bind_rows%>%select(protein_id=matches('id'),everything())
	}

	segcoutdffile<-file.path(outfolder,'segment_counts_df.tsv')

	message('writing table to')
	message(segcoutdffile)

	segcountclasses<-list(protein_id = "character", spec_coef = "numeric", spec_pval = "numeric",
	    length = "numeric", total = "numeric", startcount = "numeric",
    	enddcount = "numeric", centercount = "numeric")

	stopifnot(identical(dput(map(segcountdf,class)),segcountclasses))
	stopifnot(nrow(segcountdf)>0)
	segcountdf%>%write_tsv(segcoutdffile)

}







if(interactive()){

ns_all_genes = here('pipeline/riboseq_quant/data/*//segment_counts_df.tsv')%>%Sys.glob%>%str_subset('ribo|total')%>%setNames(.,.)%>%map_df(fread,.id='file')
# #why so diff
# library(testthat)
# ns_all_genes%>%bind_rows
# #Note that this test is pretty impefect since it fails to filter by read quality, or size.
test_that("Signal from this matches signal from bamcounts fairly closely",{

	test = ns_all_genes$protein_id[1]
	testcds<-cds%>%
		#subset(gene_name=='Satb2')%>%
	subset(protein_id==test)
	testcount<-testcds%>%bamsignals::bamCount(bam,.)%>%sum
	segcounts <- ns_all_genes%>%bind_rows%>%subset(protein_id==test)
	segcount<-ns_all_genes%>%bind_rows%>%subset(protein_id==test)%>%filter(protein_id==test)%>%{.$startcount+.$enddcount+.$centercount}
	expect_true( (abs(testcount-segcount)/testcount) < 0.2)
})


test_that("psite models should not exist for the total RNA libraries, and they should look fairly reasonable for the ones that are ribo",{
	shiftmodelfiles = Sys.glob(glue('pipeline/seqshift_reads/data/*/seqshiftmodel.rds'))
	totalsamples <- 	sampleparams%>%filter(assay=='total')%>%.$sample_id
	nototalpsitemodels <- shiftmodelfiles%>%str_detect(paste(totalsamples,collapse='|'))%>%any
	expect_false(nototalpsitemodels)

	 
	
	offsetcompdf <- Sys.glob(glue('pipeline/seqshift_reads/data/*ribo*/cdsmax*'))%>%setNames(.,.)%>%map_df(.id='file',~ fread(.))%>%mutate(file=basename(dirname(file)))%>%filter(compartment!='chrM')%>%
		group_by(file,length,offset)%>%filter(is.na(phase))%>%slice(which.max(na.omit(score)))%>%
		distinct(file,length,offset)%>%spread(file,offset)%>%as.data.frame
	libshave3offsets <- offsetcompdf[,-1]%>%apply(1,is.na)%>%`!`%>%rowSums%>%`>`(3)
	expect_true(all(libshave3offsets))

})

test_that("The situation with extending the flanks for each gene has been dealt with appropriately",{
	#DEAL WITH above by checks
	#make sure no link between strand and periodicity
	strandspeccoefs<-ns_all_genes%>%
		mutate(neg = protein_id %in% (cds%>%subset(strand=='-')%>%.$protein_id%>%unique))%>%
		group_by(neg)%>%summarise(spec_coef=median(na.omit(spec_coef)))
	expect_true(strandspeccoefs$spec_coef%>%{abs(.[1]-.[2]) / mean(.)} < 0.5)

})
test_that('subseting the grls works like I think',{
	expect_true(identical(unlist(tmpgrl,use.names=F)[lag(cuminds,1,0)+1] , tmpgrl%>%lapply(head,1)%>%GRangesList%>%unlist(use.names=F)))
})


}