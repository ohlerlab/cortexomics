library(Gviz)
library(Rsamtools) 
library(tidyverse)
library(rtracklayer)
library(stringr)
library(biomaRt)

strandedBamImport <- function (file, selection){
  if (!file.exists(paste(file, "bai", sep = "."))){ 
    stop("Unable to find index for BAM file '", file, "'. You can   build an index using the following command:\n\t", 
         "library(Rsamtools)\n\tindexBam(\"", file, "\")")
  }
sinfo <- scanBamHeader(file)[[1]] 
res <- if (!as.character(seqnames(selection)[1]) %in% names(sinfo$targets)) { 
  mcols(selection) <- DataFrame(score = 0) 
  selection 
}else { 

  param <- ScanBamParam(what = c("pos", "qwidth", "strand"), 
                        which = selection, flag = 
                          scanBamFlag(isUnmappedQuery = FALSE)) 
  x <- scanBam(file, param = param)[[1]] 

    rwidth = x[['qwidth']]
  
    gr <- GRanges(strand=x[["strand"]], ranges=IRanges(x[["pos"]], 
                                                     width = rwidth), seqnames=seqnames(selection)[1]) 
  grs <- split(gr, strand(gr)) 
  cov <- lapply(grs[c("+", "-")], function(y) coverage(ranges(y), 
                                                       width=end(selection))) 

  if(length(pos)==0){ 
    mcols(selection) <- DataFrame(plus=0, minus=0) 
    selection 
  }else{ 
    cov = cov%>%map(.%>%list%>%as('SimpleRleList')%>%setNames(seqnames(selection)[1])%>%as("GRanges")%>%subset(score!=0))
    colnames(mcols(cov[['+']])) = 'plus'
    colnames(mcols(cov[['-']])) = 'minus'
    cov[['+']]$minus = 0
    cov[['-']]$plus = 0
    cov[['-']]$minus = -1 * cov[['-']]$minus
    mcols(cov[['-']]) <- mcols(cov[['-']])[c('plus','minus')]   
    c(cov[['+']],cov[['-']])
    # GRanges(seqnames = seqnames(selection)[1], 
    #         ranges=IRanges(start=head(pos, -1), end=tail(pos, -1)), 
    #         plus=as.numeric(cov[["+"]][head(pos, -1)]), 
    #         minus=-as.numeric(cov[["-"]][head(pos, -1)])) 
  } 
} 
return(res) 
} 

datafile <- 'data/star/data/P0_ribo_2/P0_ribo_2.bam'
transcriptfile = '~/projects/cortexomics/data/transcripts.gff3'

# #biomart import gtrack code
# mart = biomaRt::useMart('ENSEMBL_MART_MOUSE')
# listDatasets(mart)
# mart = biomaRt::useMart('ENSEMBL_MART_MOUSE',dataset = 'mc57bl6nj_gene_ensembl')
# biomTrack <- BiomartGeneRegionTrack(mart,
#                                     chromosome = as.character(gchr), start = gstart, end = gend,
#                                     name = "ENSEMBL")
# 

########
transcriptfile <- '~/projects/cortexomics/data/static_local/gencode.vM12.annotation.gtf'
transcripts <- transcriptfile%>% {rtracklayer::import(.)}
transcripts = transcriptfile %>% import
vgene = transcripts%>%subset(gene_name%>%str_detect('Pa2g4') & (type=='gene'))
gchr = vgene%>%seqnames
gstart = (vgene%>%start)
gend = (vgene%>%end)
#create the gene track for gviz
gtrack = transcripts%>%{
  transcripts = .
  transsubset=transcripts%>%subsetByOverlaps(vgene)
  transsubset$feature %<>% as.character
  transsubset%>%subset(feature=='UTR')
  transsubset<-transsubset[ transsubset$feature %in% c("CDS","UTR")]
  transsubset$exon = transsubset$exon_id
  transsubset$gene = transsubset$gene_name
  transsubset$gene = transsubset$transcript_name
  transsubset$symbol = transsubset$transcript_name
  transsubset$feature %<>% as.character
  transsubset$transcript= transsubset$transcript_id
  # transsubset$group = transsubset$transcript_id 
  gtrack = GeneRegionTrack(transsubset,thinBoxFeature=c("UTR"),showId=TRUE,geneSymbol=TRUE)
}

#create the datatrack for gviz
datafile = "~/projects/cortexomics/data/star/data/E13_ribo_1/E13_ribo_1.bam"
dTrack <- DataTrack(range=datafile, genome="mm10", name="Coverage", 
                    window=-1, chromosome=gchr, importFunction=strandedBamImport, 
                    stream=TRUE) 

plotTracks(list(dTrack,gtrack),chr = as.character(gchr), from = gend - 500, to = gend, col=c("red", "blue"), 
           groups=c("+", "-"), type="hist", col.histogram=NA)





##Code dump

strandedBamImport <- function (file, selection){
  if (!file.exists(paste(file, "bai", sep = "."))){ 
    stop("Unable to find index for BAM file '", file, "'. You can   build an index using the following command:\n\t", 
         "library(Rsamtools)\n\tindexBam(\"", file, "\")")
  }
  sinfo <- scanBamHeader(file)[[1]] 
  res <- if (!as.character(seqnames(selection)[1]) %in% names(sinfo$targets)) { 
    mcols(selection) <- DataFrame(score = 0) 
    selection 
  }else { 
    
    param <- ScanBamParam(what = c("pos", "qwidth", "strand"), 
                          which = selection, flag = 
                            scanBamFlag(isUnmappedQuery = FALSE)) 
    x <- scanBam(file, param = param)[[1]] 
    
    rwidth = x[['qwidth']]
    
    gr <- GRanges(strand=x[["strand"]], ranges=IRanges(x[["pos"]], 
                                                       width = rwidth), seqnames=seqnames(selection)[1]) 
    grs <- split(gr, strand(gr)) 
    cov <- lapply(grs[c("+", "-")], function(y) coverage(ranges(y), 
                                                         width=end(selection))) 
    
    if(length(pos)==0){ 
      mcols(selection) <- DataFrame(plus=0, minus=0) 
      selection 
    }else{ 
      browser()
      cov = cov%>%map(.%>%list%>%as('SimpleRleList')%>%setNames(seqnames(selection)[1])%>%as("GRanges")%>%subset(score!=0))
      colnames(mcols(cov[['+']])) = 'plus'
      colnames(mcols(cov[['-']])) = 'minus'
      cov[['+']]$minus = 0
      cov[['-']]$plus = 0
      cov[['-']]$minus = -1 * cov[['-']]$minus
      mcols(cov[['-']]) <- mcols(cov[['-']])[c('plus','minus')]   
      c(cov[['+']],cov[['-']])
      # GRanges(seqnames = seqnames(selection)[1], 
      #         ranges=IRanges(start=head(pos, -1), end=tail(pos, -1)), 
      #         plus=as.numeric(cov[["+"]][head(pos, -1)]), 
      #         minus=-as.numeric(cov[["-"]][head(pos, -1)])) 
    } 
  } 
  return(res) 
} 

datafile <- 'data/star/data/P0_ribo_2/P0_ribo_2.bam'
transcriptfile = '~/projects/cortexomics/data/transcripts.gff3'

# #biomart import gtrack code
# mart = biomaRt::useMart('ENSEMBL_MART_MOUSE')
# listDatasets(mart)
# mart = biomaRt::useMart('ENSEMBL_MART_MOUSE',dataset = 'mc57bl6nj_gene_ensembl')
# biomTrack <- BiomartGeneRegionTrack(mart,
#                                     chromosome = as.character(gchr), start = gstart, end = gend,
#                                     name = "ENSEMBL")
# 

########
transcriptfile <- '~/projects/cortexomics/data/static_local/gencode.vM12.annotation.gtf'
transcripts <- transcriptfile%>% {rtracklayer::import(.)}
transcripts = transcriptfile %>% import
vgene = transcripts%>%subset(gene_name%>%str_detect('Pa2g4') & (type=='gene'))
gchr = vgene%>%seqnames
gstart = (vgene%>%start)
gend = (vgene%>%end)
#create the gene track for gviz
gtrack = transcripts%>%{
  transcripts = .
  transsubset=transcripts%>%subsetByOverlaps(vgene)
  transsubset$feature %<>% as.character
  transsubset%>%subset(feature=='UTR')
  transsubset<-transsubset[ transsubset$feature %in% c("CDS","UTR")]
  transsubset$exon = transsubset$exon_id
  transsubset$gene = transsubset$gene_name
  transsubset$gene = transsubset$transcript_name
  transsubset$symbol = transsubset$transcript_name
  transsubset$feature %<>% as.character
  transsubset$transcript= transsubset$transcript_id
  # transsubset$group = transsubset$transcript_id 
  gtrack = GeneRegionTrack(transsubset,thinBoxFeature=c("UTR"),showId=TRUE,geneSymbol=TRUE)
}

#create the datatrack for gviz
datafile = "~/projects/cortexomics/data/star/data/E13_ribo_1/E13_ribo_1.bam"
dTrack <- DataTrack(range=datafile, genome="mm10", name="Coverage", 
                    window=-1, chromosome=gchr, importFunction=strandedBamImport, 
                    stream=TRUE) 

plotTracks(list(dTrack,gtrack),chr = as.character(gchr), from = gend - 500, to = gend, col=c("red", "blue"), 
           groups=c("+", "-"), type="hist", col.histogram=NA)

