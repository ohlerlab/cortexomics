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
suppressMessages({library(GenomicAlignments)})
suppressMessages({library(GenomicFeatures)})
suppressMessages({library(GenomicFiles)})
suppressMessages({library(bamsignals)})
suppressMessages({library(zeallot)})
suppressMessages({library(stringr)})
suppressMessages({library(here)})
suppressMessages({library(assertthat)})
suppressMessages({library(multitaper)})

library(parallel)

reduce <- GenomicRanges::reduce

MAPQTHRESH <- 200
USEPHASE <- FALSE
USERIBOSEQC <- FALSE
source(here('src/R/Rprofile.R'))



for(fname in lsf.str('package:GenomicRanges')) assign(fname,get(fname,'package:GenomicRanges'))
for(fname in lsf.str('package:dplyr')) assign(fname,get(fname,'package:dplyr'))


argv <- c(
	bam = here('pipeline/star/data/E13_ribo_2/E13_ribo_2.bam')%T>%{stopifnot(file.exists(.))},
	gtf = here('pipeline/my_gencode.vM12.annotation.gtf'),
	REF = here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'),
	shiftmodel = 'pipeline/seqshift_reads/data/E13_ribo_2/seqshiftmodel.rds',
	outfolder = 'riboseq_quant/data/E13_ribo_2/'
)


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
segment_ribocount<-function(exonsubset,cdssubset,bam,psite_model,REF,
	STARTCLIP = 15,
	ENDCLIP = 5,
	startflank = 9,
	endflank = 9 ,
	nbp=2,
	do_seqshift=TRUE,
	revstrand=FALSE){

	cat(names(exonsubset)[1])
	
	stopifnot(!is.null(names(exonsubset)))
	assert_that(all(names(exonsubset)==names(cdssubset)))

	require(multitaper)
	require(GenomicFeatures)
	require(GenomicAlignments)


	#get only reads which actually overlap cds
	riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=MAPQTHRESH,which=unlist(cdssubset))
	#read in reas
	reads <- readGAlignments(bam,param=riboparam)
	#invert stran if needed, use only the ones on right strand
	table(strand(reads))
	if(revstrand) reads%<>%invertStrand
	reads%<>%subsetByOverlaps(cdssubset)
	mcols(reads)$length <- qwidth(reads)

	if(!is.null(psite_model)) mcols(reads)$cdsshift <- get_cds_offsets(reads,psite_model$offsets,psite_model$compartments)


	if(do_seqshift){
	

		seqdata <- get_seqforrest_data(reads,FaFile(REF),nbp=2)

		mcols(topcdsreads)$seqshift <-  get_predictedshifts(psite_model$seqshiftmodel,seqdata,prob=0.5)

		psites <- reads%>%
					subset(!is.na(cdsshift))%>%
					apply_psite_offset(c('cdsshift','seqshift'))%>%as("GRanges")%>%
					mapToTranscripts(exonsubset)
	}else if(!is.null(psite_model)){
		psites <- reads%>%
					subset(!is.na(cdsshift))%>%
					apply_psite_offset(c('cdsshift'))%>%as("GRanges")%>%
					mapToTranscripts(exonsubset)
	}else{
		psites <- reads%>%qnarrow(start=floor(qwidth(.)/2),width=1)
	}

	cdsmap <- cdssubset%>%unlist%>%pmapToTranscripts(exonsubset[names(.)])%>%reduce

	assert_that(length(cdsmap)==length(exonsubset))

	assert_that(all(table(seqnames(cdsmap))==1),msg='CDS shoudl be one unbroken stretch of the transcript')
	assert_that(all(names(exonsubset) %in% seqlevels(cdsmap)),msg='CDS shoudl be one unbroken stretch of the transcript')

	startclipbp<-(3*STARTCLIP) + startflank
	stopclipbp<-(3*ENDCLIP) + endflank

	covvects <- psites%>%coverage

	expcds<-cdsmap%>%resize(width(.)+startflank+endflank,'center')
	seqlengths(expcds) <- sum(width(exonsubset))
	expcds <- trim(expcds)

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

# test<-'ENSMUST00000000188'
# exons2use%<>%.[test]
# cds2use%<>%unlist(use.names=F)%>%subset(transcript_id==test)%>%split(.,.$transcript_id)


splitids <- cds$protein_id%>%unique%>%sort%>%split(floor((1:length(.))/300))

# warning('test mode')
# testp <- 'ENSMUSP00000026425'
# splitids <- testp
# splitids = 'ENSMUSP00000000001'
message('producing read counts')

bams <- here('pipeline/star/data/*/*.bam')%>%Sys.glob%>%str_subset(neg=TRUE,'transcript')%>%str_subset('_ribo_')



Sys.glob(here('pipeline/riboseq_quant/data/*'))
segcoutdffile<-file.path(,'segment_counts_df.tsv')

for(bam in bams){
	library(glue)
	thissamp<-bam%>%dirname%>%basename
	shiftmodel = here(glue('pipeline/seqshift_reads/data/{thissamp}/seqshiftmodel.rds'))
	outfolder = here(glue('pipeline/riboseq_quant/data/{thissamp}/'))

	rdsfile<-here(glue('{outfolder}/allgenes.rds'))
	segcoutdffile<-file.path(outfolder,'segment_counts_df.tsv')
	if(file.exists(segcoutdffile)) next
	# allgenesrds<-str_interp('pipeline/ribosegcount/${thissamp}/allgenes.rds')


	isrev <- fread(here('src/sample_parameter.csv'))%>%filter(sample_id==thissamp)%>%.$protocol%>%`==`('reverse')

	detectCores()
	message(getwd())
	message(outfolder)

	message(paste0('counting ', length(unlist(splitids)),' genes' ))
	ns_all_genes<-splitids%>%lapply(mc.cores=detectCores(),function(idchunk){
	# ns_all_genes<-splitids%>%lapply(function(idchunk){
		safely(segment_ribocount)(
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

	require(glue)
	rdsfile<-here(glue('{outfolder}/allgenes.rds'))

	rdsfile%>%dirname%>%dir.create(rec=TRUE,showWarn=FALSE,.)
	# ns_all_genes<-readRDS(rdsfile)
	message('saving object to')
	message(rdsfile)

	ns_all_genes%>%saveRDS(rdsfile)

	#readRDS(rdsfile)->ns_all_genes


	segcountdf<-ns_all_genes%>%map('result')%>%bind_rows%>%select(protein_id=matches('id'),everything())

	segcountdf%>%filter(protein_id=='ENSMUSP00000000001')

	segcoutdffile<-file.path(outfolder,'segment_counts_df.tsv')

	message('writing table to')
	message(segcoutdffile)
	segcountdf%>%write_tsv(segcoutdffile)

}




# #why so diff
# library(testthat)
# ns_all_genes%>%bind_rows
# #Note that this test is pretty impefect since it fails to filter by read quality, or size.
# test_that("Signal from this matches signal from bamcounts fairly closely",{

# 	testcds<-cds%>%subset(gene_name=='Pa2g4')%>%subset(protein_id==testp)
# 	testcount<-testcds%>%bamsignals::bamCount(bam,.)%>%sum
# 	segcounts <- ns_all_genes%>%map('result')%>%bind_rows%>%subset(protein_id==testp)%>%filter(protein_id==testp)
# 	segcount<-ns_all_genes%>%map('result')%>%bind_rows%>%subset(protein_id==testp)%>%filter(protein_id==testp)%>%{.$startcount+.$enddcount+.$centercount}
# 	expect_true( (abs(testcount-segcount)/testcount) < 0.4)

# })

