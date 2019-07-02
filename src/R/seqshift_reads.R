



testfunc<-function(exonsubset,cdssubset,bam,psite_model,REF,
	STARTCLIP = 15,
	ENDCLIP = 5,
	startflank = 9,
	endflank = 9 ,
	nbp=2,
	do_seqshift=TRUE){

	# browser()

	assert_that(all(names(exonsubset)==names(cdssubset)))
	require(GenomicFeatures)
	require(GenomicAlignments)

	#get only reads which actually overlap cds
	riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=200,which=unlist(cdssubset))
	reads <- readGAlignments(bam,param=riboparam)

	#seperate out the junction reads
	jreads<-reads[unlist(njunc(reads)>0)]

	reads<-reads[njunc(reads)==0]



	#expand the exons, so that we can always get the sequence for reads on the endge
	exonsubset%<>%unlist%>%setNames(NULL)
	maxrl<-psite_model$offsets$length%>%max
	exonsubset[isfpmost(exonsubset,'transcript_id')]%<>%resize(width(.)+startflank+maxrl+nbp,'end')
	exonsubset[istpmost(exonsubset,'transcript_id')]%<>%resize(width(.)+endflank+maxrl+nbp,'start')
	exonsubset%<>%split(.,.$transcript_id)

	#first carry out mappign of the spliced reads
	#This is a bit convoluted - the more obvious approach with pmaptotrancripts
	#was somehow creating a vector too big for memory, so we chunk the junction reads, unlist
	#them intot heir start sites, map to transcripts, then resize them once mapped
	jmapped <- jreads%>%split(.,ceiling(seq_along(.)/ 5e3))%>%as.list%>%lapply(function(jreads){

		ov <- findOverlaps(jreads,exonsubset)
		ovenc<-encodeOverlaps(as(jreads[queryHits(ov)],"GRangesList"),exonsubset[subjectHits(ov)])
		ov<-ov[isCompatibleWithSplicing(ovenc)]
		if(length(ov)==0) return(NULL)
		jreads <- jreads[queryHits(ov)]
		qwidths<-qwidth(jreads)
			

		jmapped <- jreads %>%narrow(start=1,width=1)  %>% as("GRanges")%>%pmapToTranscripts(exonsubset[subjectHits(ov)])
		end(jmapped) <- end(jmapped)-1+unlist(qwidths)
		jmapped
	})

	#anc 
	suppressWarnings(jmapped%<>%Reduce(f=c))


	exonsubset%>%width%>%cumsum
	cdssubset%>%width

	#now carry out mapping of unsliced reads
	reads <- as(reads,'GRanges')
	reads <- reads[overlapsAny(reads%>%resize(width(.)+nbp),exonsubset,type='within')]
	mapped <- suppressWarnings({mapToTranscripts(reads,exonsubset)})

	#and combine
	mapped <- append(mapped,jmapped)

	mapped$length <- width(mapped)

	exonsubsetchrs<-seqnames(exonsubset)%>%runValue

	chrs<-unlist(exonsubsetchrs[match(seqnames(mapped),names(exonsubset))])%>%as.vector
	mapped$compartment<-psite_model$compartments[chrs]

	mapped$cdsoffset <- mapped%>%apply_cds_offsets(psite_model$offsets)

	mapped%<>%subset(!is.na(cdsoffset))

	exonseq<-extractTranscriptSeqs(FaFile(REF),exonsubset)

	suppressWarnings(
		mapped <- mapped[!is_out_of_bounds(mapped,seqinfo(exonseq))]
		)

	# mapped%>%seqlevels
	# seqinfo(exonseq)%>%seqlevels

	# mapped%<>%head(5)



	 # reads%>%resize(2,'start')%>%resize(4,'end')%>%head(5)%>%getSeq(x=FaFile(REF),.)

	#this is necessary since our exonseqs are already reverse complemented for the negative strands
	#and mapped reads are always positive.
	strand(mapped)<-'+'

	# reads%>%resize(2,'start')%>%resize(4,'end')%>%head(5)%>%getSeq(x=FaFile(REF),.)
	
	# get_seqforrest_traindata(mapped,exonseq,trim=FALSE)%>%head(5)



	#turn the mapped reads into psites

	mapped <- mapped[,c('length','cdsoffset')]%>%resize(1,'start')%>%shift(.$cdsoffset)
	
	if(do_seqshift){
		mapped$seqoffset <- get_seqforrest_traindata(mapped,exonseq,trim=FALSE)%>%
			predict(psite_model$seqshiftmodel,.)%>%
			.$prediction%>%as.character%>%as.numeric

		mapped %<>%shift(.$seqoffset)
	}

	cdsmap <- cdssubset%>%unlist(use.names=F)%>%pmapToTranscripts(exonsubset[.$transcript_id])%>%reduce
	


	assert_that(length(cdsmap)==length(exonsubset))
	assert_that(all(table(seqnames(cdsmap))==1),msg='CDS shoudl be one unbroken stretch of the transcript')
	assert_that(all(names(exonsubset) %in% seqlevels(cdsmap)),msg='CDS shoudl be one unbroken stretch of the transcript')

	startclipbp<-(3*STARTCLIP) + startflank
	stopclipbp<-(3*ENDCLIP) + endflank

	covvects <- mapped%>%coverage

	expcds<-cdsmap%>%resize(width(.)+startflank+endflank,'center')

	cdsmap%>%resize(width(.)+startflank+endflank,'center')

	
	covvects <- covvects[expcds]
	covvects <- lapply(covvects,as.vector)
	#Get, periodicity scores
	vectftests<-covvects%>%lapply(ftestvect)%>%simplify2array

	data.frame(
		spec_coef = vectftests[1,],
		spec_pval = vectftests[2,],
		length = covvects%>%map_dbl(length),
		total = covvects%>%map_dbl(sum),
		startcount = covvects%>%map_dbl(possibly(~sum(.[1:startclipbp]),NA)),
		enddcount = covvects%>%map_dbl(possibly(~sum(.[(length(.)-stopclipbp+1):length(.)]),NA)),
		centercount = covvects%>%map_dbl(possibly(~sum(.[(startclipbp+1):(length(.)-stopclipbp)]),NA))
	)%>%rownames_to_column('transcript_id')

}


cds2use <- cds%>%split(.,.$protein_id)%>%head(100)

cds2use<-topcds%>%split(.,.$protein_id)%>%head(100)

exons2use <- exons%>%subset(transcript_id %in% unlist(cds2use)$transcript_id)
exons2use%<>%split(.,.$transcript_id)

cds2use%<>%unlist
# cds2use<-cds2use[isfpmost(cds2use,'transcript_id')]%>%resize(1,'start')%>%shift(-3)
cds2use%<>%split(.,.$transcript_id)


# test<-'ENSMUST00000000188'
# exons2use%<>%.[test]
# cds2use%<>%unlist(use.names=F)%>%subset(transcript_id==test)%>%split(.,.$transcript_id)
# 0

testfunc(
	exons2use,
	cds2use,
	bam,
	psite_model,
	REF,
	STARTCLIP = 15,
	ENDCLIP = 5,
	startflank = 9,
	endflank = 9 ,
	do_seqshift=FALSE
)




# testfunc(
# 	exons2use,
# 	cds2use,
# 	bam,
# 	psite_model,
# 	REF,
# 	STARTCLIP = 15,
# 	ENDCLIP = 5,
# 	startflank = 9,
# 	endflank = 9 ,
# 	do_seqshift=TRUE)

# mapped<-mapped%>%head(5)
# reads<-reads%>%head(5)

# 	strand(mapped)<-'+'

# mapped$seqoffset <- get_seqforrest_traindata(mapped,exonseq,trim=FALSE)%>%
# 	predict(psite_model$seqshiftmodel,.)%>%
# 	.$prediction%>%
# 	as.character%>%
# 	as.numeric



# mapFromTranscripts(mapped,exons2use)
# reads




####
#alternative approaches - try looking at independance of end distances
#Try looking at artifical positive data
#look at read length distribution around pos and neg







