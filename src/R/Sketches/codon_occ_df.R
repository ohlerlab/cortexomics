bam = ''
get_genomic_psites <- function(bam,windows,offsets,mapqthresh=200,comps=c('chrM'='chrM')) {
  require(GenomicAlignments)
  require(hashmap)
  riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=windows)
  reads <- readGAlignments(bam,param=riboparam)
  mcols(reads)$length <- qwidth(reads)
  #
  #if no 'offsets' data then just take read midpoints (e.g. for RNAseq)
  if(is.null(offsets)){
    mcols(reads)$offset <- floor(qwidth(reads)/2)
  }else{
    #else check we have offsets
  stopifnot('length' %in% colnames(offsets))
  stopifnot('offset' %in% colnames(offsets))
  stopifnot('comp' %in% colnames(offsets))
    #
    #define our compartment for all of the seqnames in the object
    comps = comps[names(comps)%in%seqnames(reads)]
    useqnms <- as.character(unique(seqnames(reads)))%>%setdiff(names(comps))
    compmap = safe_hashmap(c(useqnms,names(comps)),c(rep('nucl',length(useqnms)),comps))
    stopifnot(all(compmap$values()%in%offsets$comp))
    #now fetch an offset for all of our 
    mcols(reads)$offset <-
      data.frame(length=mcols(reads)$length,
        comp=compmap[[as.character(seqnames(reads))]])%>%
      safe_left_join(offsets%>%select(offset,length,comp),allow_missing=TRUE,by=c('length','comp'))%>%.$offset
  }
  #get rid of the reads that have no offset
  reads <- reads%>%subset(!is.na(mcols(reads)$offset))
  psites <- apply_psite_offset(reads,c('offset'))%>%as("GRanges")
  #annotate psites with length
  mcols(psites)$length <- mcols(reads)$length
  psites
}


bestcdsseq = cdsgrl%>%extractTranscriptSeqs(x=FaFile(REF),.)

codposdf = lapply(bestcdsseq,function(cdsseq){
	codonmat = codons(cdsseq)%>%{cbind(pos = .@ranges@start,as.data.frame(.))}%>%
		identity
})

sitedf <- testpsitecovs%>%lapply(as.vector)%>%stack%>%set_colnames(c('count','gene'))
