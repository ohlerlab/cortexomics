
spl_mapFromTranscripts <- function(trspacegr,exons_grl){

  exons_tr<-exons_grl%>%unlist%>%mapToTranscripts(exons_grl)%>%.[names(.)==seqnames(.)]
  ov <- findOverlaps(trspacegr,exons_tr)

  trspacegr_spl <- suppressWarnings({trspacegr[queryHits(ov)]%>%pintersect(exons_tr[subjectHits(ov)])})
  genomic_trspacegr <- mapFromTranscripts(
  trspacegr_spl,
  # exons_tr[subjectHits(ov)]%>%split(.,seqnames(.))
  exons_grl
  )
  genomic_trspacegr$xHits <- queryHits(ov)[genomic_trspacegr$xHits]
  genomic_trspacegr
}

resize_grl_startfix<-function(grl,width){
  #what follows is some slightly black magic using S4 vectors
  #Integerlist which showings how much we'd need to trim that exon to get to to the desired transcript length
  trim =  cumsum(width(grl)) - width 
  #Where trim is greater than the exon width, we drop it
  drop = trim >=  width(grl)
  grl = grl[!drop]
  #vector showing location of the new 3' end of each transcript
  newends = cumsum(elementNROWS(grl))
  #vector with the amount we need to trim each new 3' end by
  endtrims=trim[IntegerList(as.list(elementNROWS(grl)))]@unlistData
  #finally, use these to trim
  grl@unlistData[newends] <- resize(grl@unlistData[newends], width(grl@unlistData[newends]) - endtrims  )
  grl
  
}

str_order_grl<-function(grl){order( start(grl)*(((strand(grl)!='-')+1)*2 -3) )}
sort_grl_st <- function(grl)grl[str_order_grl(grl),]
resize_grl_endfix <- function(grl,width){
  grl = invertStrand(grl)%>%sort_grl_st
  
  grl = resize_grl_startfix(grl,width)
  invertStrand(grl)%>%sort_grl_st
}
resize_grl <- function(grl,width,fix='start',check=TRUE){
  stopifnot(all(width>0))
  assert_that(all(all(diff(str_order_grl(grl))==1) ),msg = "grl needs to be 5'->3' sorted")
  if(fix=='start'){
    grl = resize_grl_startfix(grl,width)
  }else if(fix=='end'){
    grl = resize_grl_endfix(grl,width)
  }else if(fix=='center'){
    grlwidths = sum(width(grl)) 
    diffs = (width - grlwidths)
    
    grl = resize_grl_startfix(grl,grlwidths + ceiling(diffs/2))
    grl = resize_grl_endfix(grl,grlwidths + diffs)
    
  }
  if(check){
    startstoolow <- any(start(grl)<=0)
    if(any(startstoolow)){
      stop(str_interp("${sum(startstoolow)} ranges extended below 1 .. e.g. ${head(which(startstoolow,1))}"))
    }
    grlseqs <- as.vector(unlist(use.names=F,seqnames(grl)[IntegerList(as.list(rep(1,length(grl))))]))
    endstoohigh <- any((end(grl)>seqlengths(grl)[grlseqs])%in%TRUE)
    if(any(endstoohigh)){
      stop(str_interp("${sum(endstoohigh)} ranges extended below above seqlength .. e.g. ${head(which(endstoohigh,1))}"))
    }
  }
  grl
}
trim_grl <- function(grl,bp,end='tp'){
  if(end=='tp'){
    resize_grl(grl,sum(width(grl)) - bp,fix='start')
  }else if(end=='fp'){
    resize_grl(grl,sum(width(grl)) - bp,fix='end')
  }else {
    stop("end should be fp or tp")
  }
}

fmcols <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)[start(grl@partitioning)]
}
fmcols_List <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)%>%split(grl@partitioning)
}
# base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
try(silent=T,{library(colorout)})
library(tidyverse)
suppressMessages(library(magrittr))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(Rsamtools))
suppressMessages(library(here))
select <- dplyr::select

library(GenomicFeatures)
windback <- 60
windforw <- 60
outprefix = 'ribotrans'

# gtf = here('pipeline/my_gencode.vM12.annotation.gtf')
# fafile = Rsamtools::FaFile('pipeline/my_GRCm38.p5.genome.chr_scaff.fa')

args = commandArgs(trailingOnly=TRUE)
gtf = args[1]
fafile = args[2]
outprefix = args[3]

# rm(gtf_gr,exons,cds)
if(!exists('gtf_gr')) gtf_gr<-rtracklayer::import(con=gtf,format='gtf')
if(!exists('exons',where=baseenv())) exons <- gtf_gr%>%subset(type=='exon')
if(!exists('cds',where=baseenv())) cds <- gtf_gr%>%subset(type=='CDS')


cdsgrl = cds%>%split(.,.$transcript_id)
exonsgrl = exons%>%split(.,.$transcript_id)


#get expanded exons - so that we never hit the edge of the exons when expanding
#our cds by the necessary amount

cdsseq = extractTranscriptSeqs(x=FaFile(fafile),cdsgrl)
endcods = translate(subseq(cdsseq,-3,-1))%>%table
cdshasstop = sort(endcods)%>%tail(1)%>%names%>%`==`('*')
if(!cdshasstop){
	library(assertthat)
	cdsgrlext = cdsgrl%>%resize_grl(sum(width(.))+3,'start')
	cdsseq = extractTranscriptSeqs(x=FaFile(fafile),filtcdsgrl)
}else{
	cdsgrlext = cdsgrl%>%resize_grl('start')
	cdsgrl = cdsgrl %>%resize_grl(sum(width(.))-3,'start')
}
is3bpmult = cdsgrl%>%width%>%sum%>%`%%`(3)%>%`==`(0)
filtcdsgrl = cdsgrl[is3bpmult]
filtcdsgrl = filtcdsgrl%>%sort_grl_st
message(str_interp('filtered out ${sum(!is3bpmult)} ORFs for not being multiples of 3bp long'))
cdsseq = extractTranscriptSeqs(x=FaFile(fafile),cdsgrlext)
validORF = cdsseq%>%translate(if.fuzzy.codon='X')%>%str_detect('^M[^*X]+\\*$')%>%setNames(names(cdsgrlext))
filtcdsgrl = filtcdsgrl[validORF[names(filtcdsgrl)]]
message(str_interp('filtered out ${sum(!validORF)} ORFs for not starting with M and ending with *'))
#
length(filtcdsgrl)
message(str_interp('${length(filtcdsgrl)} ORFs left'))
#
#get exons for our cds
cdsexonsgrl <- exonsgrl[fmcols(filtcdsgrl,transcript_id)]
cdsexonsgrl%<>%sort_grl_st
#get an object representing the CDS In transript space
cdstrspace = pmapToTranscripts(filtcdsgrl,cdsexonsgrl[fmcols(filtcdsgrl,transcript_id)])
#ensure all cds map cleanly to the exons
stopifnot(cdstrspace%>%elementNROWS%>%`==`(1))

cdsstartpos = setNames(start(cdstrspace@unlistData),names(cdstrspace))
endpos = sum(width(cdsexonsgrl))-end(cdstrspace@unlistData)

#expand our first exon when needed
startposexpansion = pmax(0,windback - cdsstartpos + 1)
#expand/trim the 5' end of the exons
cdsexonsgrl@unlistData[start(cdsexonsgrl@partitioning)]%<>%resize(width(.)+startposexpansion,'end')
#expand or trim the last exon when needed
endposexpansion = pmax(0,windforw - endpos)
cdsexonsgrl@unlistData[cdsexonsgrl@partitioning@end]%<>%resize(width(.)+endposexpansion,'start')

{
#now map our cds to that
cds_exptrspc = pmapToTranscripts(filtcdsgrl,cdsexonsgrl)
stopifnot(cds_exptrspc%>%elementNROWS%>%`==`(1))

expcds_exptrspc=cds_exptrspc
stopifnot(!any(expcds_exptrspc%>%elementNROWS%>%`>`(1)))
expcds_exptrspc%<>%unlist
#expand our cds exons
expcds_exptrspc%<>%resize(width(.)+windback,'end',ignore.strand=TRUE)
#and expand the 3' ends
expcds_exptrspc%<>%resize(width(.)+windforw,'start',ignore.strand=TRUE)
#now back to genome space
expcdsgenspace = spl_mapFromTranscripts(expcds_exptrspc,cdsexonsgrl)
expcdsgenspace = split(expcdsgenspace,names(expcdsgenspace))
#get the sequences
expcdsgenspaceseq <- 
	expcdsgenspace%>%
	sort_grl_st%>%
	extractTranscriptSeqs(x=FaFile(fafile))

}

{
fastanames <- paste(sep='|',
	fmcols(cdsexonsgrl,transcript_id),
	fmcols(cdsexonsgrl,gene_id),
	fmcols(cdsexonsgrl,havana_gene),
	fmcols(cdsexonsgrl,havana_transcript),
	fmcols(cdsexonsgrl,transcript_name),
	fmcols(cdsexonsgrl,gene_name),
	sum(width(cdsexonsgrl)),
	paste0('UTR5:1-',start(cds_exptrspc)-1),
	paste0('CDS:',start(cds_exptrspc),'-',end(cds_exptrspc)),
	paste0('UTR3:',1+end(cds_exptrspc),'-',sum(width(cdsexonsgrl))),
	'|')

names(expcdsgenspaceseq) <- fastanames

library(rtracklayer)

as.data.frame(unlist(cds_exptrspc))%>%select(seqnames,start,end)%>%write_tsv(paste0(outprefix,'_trcds.tsv'))
}

{
#also write our cds coordinates to disk in the new trspace
new_trspc_anno <- c(
	GRanges(as.character(seqnames(cds_exptrspc)),IRanges(1,sum(width(cdsexonsgrl))))%>%{.$type='exon';.},
	cds_exptrspc%>%unlist%>%{.$type='CDS';.}
)
new_trspc_anno%>%export(paste0(outprefix,'_trspaceanno.gtf'))
#write the expanded exons sequences to disk
writeXStringSet(expcdsgenspaceseq,paste0(outprefix,'.fa'))
}

if(FALSE){
# cds4salmonseq%>%writeXStringSet('pipeline/')
expcdsgenspaceseq%>%writeXStringSet(here('pipeline/gencode.vM12.trimtrs.fa'))

expcdata = getSeq(x=Rsamtools::FaFile('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.pc_transcripts_filter.fa'))

exheader = expcdata%>%names%>%.[1000]
exttr = exheader%>%str_extract('ENSMUST\\d+')
cdsgrl[exttr]%>%pmapToTranscripts(exonsgrl[exttr])
exonsgrl[exttr]%>%width%>%sum

testr = 'ENST00000340870.6'
which(cdstrspace%>%elementNROWS%>%`!=`(1))


}
# current strategy - quantify with salmon and dp
# but this means that we transcripts that hardly have any UTR
# This makes the counts reflect differences in UTR lenght rather than just TE
# I should quantify TE using the 




 