#load the gtf
library(rtracklayer)
library(here)
library(tidyverse)
gtf = here('annotation/gencode.vM12.annotation.gtf')
STARTTRIMCODS <- 15
ENDTRIMCODS <- 5
if(!exists('gtf_gr')) gtf_gr<-import(con=gtf,format='gtf')
# gtf_grbak <- gtf_gr
# gtf_gr<-gtf_grbak

# gtf_gr <- gtf_grbak[gtf_grbak$transcript_id%in%unique(gtf_grbak$transcript_id)[666:676],]
# gtf_gr<-sort(gtf_gr)

exons <- gtf_gr%>%subset(type=='exon')
cds <- gtf_gr%>%subset(type=='CDS')



trimlen = (STARTTRIMCODS*3) + (STARTTRIMCODS*3)
#define longe enough ids
longenough_ids = cds%>%split(.,.$protein_id)%>%width%>%sum%>%.[.>trimlen]%>%names
#now subset
cds %<>% subset(protein_id %in% longenough_ids)
exons %<>% subset(transcript_id %in% unique(cds$transcript_id))
genes <- gtf_gr %>% subset(type=='gene') %>% subset(gene_id %in% exons$gene_id)

#select some cds regions and corresponding transciropts
testcds <- cds%>%split(.,.$protein_id)
testtrs <- exons%>%subset(transcript_id%in%unlist(testcds)$transcript_id)%>%split(.,.$transcript_id)
#get start codons
startcods <- testcds%>%revElements(.,any(strand(.)=='-'))%>%lapply('[',1)%>%GRangesList%>%unlist%>%resize(1,'start')
endcods <- testcds%>%revElements(.,any(strand(.)=='+'))%>%lapply('[',1)%>%GRangesList%>%unlist%>%resize(1,'end')
#map our start codons to transcript space
library(GenomicFeatures)
matching_trs <- testtrs[startcods$transcript_id]
startcods_trspace <- pmapToTranscripts(startcods,matching_trs)
endcods_trspace <- pmapToTranscripts(endcods,matching_trs)
#this function maps granges in transcript space back to the genome, poentially creating more
#than n ranges for n in put ranges, if they cross exon boundaries
spl_mapFromTranscripts<-function(trspacegr,exons_grl){
  exons_tr<-exons_grl%>%unlist%>%setNames(paste0('exon_',seq_along(.)))%>%mapToTranscripts(exons_grl)
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

cds_trspace <- startcods_trspace
end(cds_trspace) <- end(endcods_trspace)

start(cds_trspace) <- start(cds_trspace) + (STARTTRIMCODS*3)
end(cds_trspace) <- end(cds_trspace) - (ENDTRIMCODS*3)

trimmed_cds <- cds_trspace%>%
  spl_mapFromTranscripts(testtrs)%>%#map back to the genome
  split(.,.$xHits)%>%
  unlist(use.names=TRUE)

mcols(trimmed_cds) <- testcds[trimmed_cds$xHits]%>%.[IntegerList(as.list(rep(1,length(.))))]%>%unlist%>%
  {mcols(.)[,c('gene_id','transcript_id','protein_id','gene_name','transcript_name','type','phase','tag','gene_type')]}
names(trimmed_cds)<-NULL

trimmed_cds%>%export(gtf%>%str_replace('.gtf$','.trimmed.gtf'))

grl = GRangesList(list(c(
  GRanges('a:3-6:+',foo=1),
  GRanges('a:8-10:+'),
  GRanges('a:13-15:+')
),  c(GRanges('a:3-6:-'),
  GRanges('a:8-10:-'),
  GRanges('a:13-15:-',bar=2)
)))


data.frame(p=names(unlist(grl@partitioning)),start=grl@unlistData@ranges@start,w=grl@unlistData@ranges@width,strand=grl@unlistData@strand)%>%
  mutate(oorder=1:n())%>%
  group_by(p)%>%
  arrange(p,ifelse(strand=='-',-1,1)*start)%>%
  mutate(cw = lag(cumsum(w),1,0))%>%
  # mutate(cw = lead(cumsum(w),1,0))%>%
  filter(cw < )

#first eliminate everything for which cw (the lenght of prior exons) is equal to or greater than the nbp
#Then, for the _last_ one of each group, set width equal to nbp - cw


  # lapply(c(1,3,5,10,15,20),function(nbps){
  nbps= if(length(nbps)==1) rep(nbps,length(grl)) else nbps

  resizeddt <- data.frame(
    p=names(unlist(grl@partitioning)),
    start=grl@unlistData@ranges@start,
    w=grl@unlistData@ranges@width,
    strand=grl@unlistData@strand
  )%>%
  group_by(p)%>%
  mutate(oorder=1:n())%>%
  arrange(p,ifelse(strand=='-',-1,1)*start)%>%
  # mutate(cw = lead(cumsum(w),1,0))%>%
  mutate(nbp = nbps[p])%>%
  filter( lag(cumsum(w),1,0) < nbp)%>%
  # mutate(nw = (1:n() == n()))
  mutate(w = ifelse((1:n()) == n(), w + (nbp - cumsum(w)),w))

  resizeddt%>%arrange(p,oorder)%>%mutate(pind = )

# })
spl_mapFromTranscripts<-function(trspacegr,exons_grl){
  exons_tr<-exons_grl%>%unlist%>%setNames(paste0('exon_',seq_along(.)))%>%mapToTranscripts(exons_grl)
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
nbps = 7
fix = 'start'
resize_grl <- function(grl,nbps,fix='start'){
  require(magrittr)
  require(GenomicFeatures)
  stopifnot(is(grl,'GRangesList'))
  stopifnot(identical(grl,sort(grl)))

  leftinds <- IntegerList(as.list(rep(1,length(grl))))
  rightinds <- IntegerList(as.list(elementNROWS(grl)))

  leftexons <- grl[leftinds]
  rightexons <- grl[rightinds]


  #make unlisted grl object
  onames <- names(grl)
  names(grl) <- seq_along(grl)
  ulgrl <- unlist(grl,use.names=TRUE)

  #numbers showing distance from left or right in original grl
  ulgrl$left_enum <- seq_along(ulgrl) - lag(cumsum(elementNROWS(grl)),1,0)[names(ulgrl)]
  ulgrl$right_enum <- cumsum(elementNROWS(grl))[names(ulgrl)] - seq_along(ulgrl) + 1

  erows = elementNROWS(grl)

  #get the fp coordinate of each transcript
  grlstrand <- (strand(ulgrl[ulgrl$left_enum==1]))
  fpends <- ulgrl[ifelse(grlstrand=='-',
    which(ulgrl$right_enum==1),
    which(ulgrl$left_enum==1)
  )]%>%resize(1)

  start(ulgrl)[ulgrl$left_enum==1] <- 1
  end(ulgrl)[ulgrl$right_enum==1] <- 1e9

  tspaceinds <- IntegerList(split(1:sum(erows),names(ulgrl)))%>%revElements(grlstrand=='-')

  ulgrl%<>%split(.,names(.))

  expspacefps <- pmapToTranscripts(fpends,ulgrl)

  expspaceranges <- expspacefps[rep(seq_along(grl),elementNROWS(grl))]

  ends <- width(grl)%>%revElements(grlstrand=='-')%>%cumsum%>%as('IntegerList')
  starts <- ends
  starts[leftinds] <- 0
  starts[-1*leftinds] = ends[(-1)*rightinds] 
  ends%<>%subtract(1)
  starts%<>%unlist
  ends%<>%unlist

  expspaceexons <- GRanges(
    rep(names(ulgrl),elementNROWS(grl)),
    IRanges(
        rep(start(expspacefps),elementNROWS(grl)) + unlist(starts),
        rep(start(expspacefps),elementNROWS(grl)) + unlist(ends)
    ),strand=rep(grlstrand,erows),
    originalind = unlist(tspaceinds)
  )
  

  expspaceexonsmerge <- GenomicRanges::reduce(expspaceexons) %>%resize(nbps,fix,ignore.strand=TRUE)%>%rep(erows)

  start(expspaceexons)[starts==0] %<>% pmin(.,start(expspaceexonsmerge)[starts==0])
  end(expspaceexons)[lead(starts==0,1,TRUE)] %<>% pmax(.,end(expspaceexonsmerge)[lead(starts==0,1,TRUE)])


  expspaceexons <- GenomicRanges::pintersect(expspaceexons,expspaceexonsmerge)%>%
    subset(hit)%>%
    subset(start<=end)



  gspacered<-spl_mapFromTranscripts(expspaceexons,ulgrl)
  grlorder <- expspaceexons$originalind[gspacered$xHits]
  mcols(gspacered) <- mcols(unlist(grl))[grlorder,]
  gspacered[order(grlorder),]%>%split(rep(onames,erows)[grlorder])


}
cds_grl <- cds%>%split(.$protein_id)%>%sort

resize_grl(cds_grl,100)

resize_grl <- function(grl,nbps,fix='start'){
  stopifnot(identical(grl,sort(grl)))

  leftinds <- IntegerList(as.list(rep(1,length(grl))))
  rightinds <- IntegerList(as.list(elementNROWS(grl)))

  leftexons <- grl[leftinds]
  rightexons <- grl[rightinds]


  #make unlisted grl object
  onames <- names(grl)
  names(grl) <- seq_along(grl)
  ulgrl <- unlist(grl,use.names=TRUE)

  grlisneg <- strand(grl[leftinds])%>%`==`('-')%>%setNames(NULL)%>%unlist
  width(grl)

  #numbers showing distance from left or right in original grl
  ulgrl$left_enum <- seq_along(ulgrl) - lag(cumsum(elementNROWS(grl)),1,0)[names(ulgrl)]
  ulgrl$right_enum <- cumsum(elementNROWS(grl))[names(ulgrl)] - seq_along(ulgrl) + 1

  erows = elementNROWS(grl)

  #When changing the tp end, we first eliminate those ranges whose start is more than the desired distance from the
  #fp end of the transcript, then resize the final region so that it's width is the desired distance
  if(fix=='start'){
    data.frame(s=start(ulgrl))
  }
  #reverse procedure for the 

  #get the fp coordinate of each transcript
  grlstrand <- (strand(ulgrl[ulgrl$left_enum==1]))
  fpends <- ulgrl[ifelse(grlstrand=='-',
    which(ulgrl$right_enum==1),
    which(ulgrl$left_enum==1)
  )]%>%resize(1)

  start(ulgrl)[ulgrl$left_enum==1] <- 1
  end(ulgrl)[ulgrl$right_enum==1] <- 1e9

  tspaceinds <- IntegerList(split(1:sum(erows),names(ulgrl)))%>%revElements(grlstrand=='-')

  ulgrl%<>%split(.,names(.))

  expspacefps <- pmapToTranscripts(fpends,ulgrl)

  expspaceranges <- expspacefps[rep(seq_along(grl),elementNROWS(grl))]

  ends <- width(grl)%>%revElements(grlstrand=='-')%>%cumsum%>%as('IntegerList')
  starts <- ends
  starts[leftinds] <- 0
  starts[-1*leftinds] = ends[(-1)*rightinds] 
  ends%<>%subtract(1)
  starts%<>%unlist
  ends%<>%unlist

  expspaceexons <- GRanges(
    rep(names(ulgrl),elementNROWS(grl)),
    IRanges(
        rep(start(expspacefps),elementNROWS(grl)) + unlist(starts),
        rep(start(expspacefps),elementNROWS(grl)) + unlist(ends)
    ),strand=rep(grlstrand,erows),
    originalind = unlist(tspaceinds)
  )
  
  data.frame(p=names(unlist(grl@partitioning)),grl@unlistData@ranges@start,w=grl@unlistData@ranges@width,strand=grl@unlistData@strand)%>%group_by(p)%>%mutate(wcumsum=ifelse(strand=='-',sum(w)-cumsum(w),cumsum(w)))%>%mutate(oorder=1:n())

  expspaceexonsmerge <- GenomicRanges::reduce(expspaceexons) %>%resize(nbps,fix,ignore.strand=TRUE)%>%rep(erows)

  start(expspaceexons)[starts==0] %<>% pmin(.,start(expspaceexonsmerge)[starts==0])
  end(expspaceexons)[lead(starts==0,1,TRUE)] %<>% pmax(.,end(expspaceexonsmerge)[lead(starts==0,1,TRUE)])


  expspaceexons <- GenomicRanges::pintersect(expspaceexons,expspaceexonsmerge)%>%
    subset(hit)%>%
    subset(start<=end)



  gspacered<-spl_mapFromTranscripts(expspaceexons,ulgrl)
  grlorder <- expspaceexons$originalind[gspacered$xHits]
  mcols(gspacered) <- mcols(unlist(grl))[grlorder,]
  gspacered[order(grlorder),]%>%split(rep(onames,erows)[grlorder])


}

#lapply(1:20,function(i){resize_grl(grl,i,fix='start')%>%width})
# cds_grl<-gtf_gr%>%{split(.,.$protein_id)}
resize_grl(sort(cdsgrl), 3)