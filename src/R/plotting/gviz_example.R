
library(purrr)
library(here)
library(ggExtra)
library(ggpubr)
library(tidyverse)
library(Rcpp)
library(doMC)
library(rstan)
#BiocManager::install('rstan')
library(tidyverse)
library(magrittr)
library(data.table)
library(stringr)
library(magrittr)
library(splines)
library(parallel)
#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(assertthat))
suppressMessages(library(limma))
message('...done')

#define an annotation gtf
transcriptfile <- here('pipeline/my_gencode.vM12.annotation.gtf')
#get the bame files
# datafiles <- Sys.glob(here::here('data/mergedbigwigs/*/*/*'))%>%
datafiles <- Sys.glob(here::here('data/bigwigs/*/*/*'))%>%
  #filter out the ones with transcript in the name
  grep(v=TRUE,inv=TRUE,patt='transcript')
#subset the bam files, for speed
datafiles <- datafiles[c(1:4)]


message('loading transcripts')
transcripts <- transcriptfile%>% {rtracklayer::import(.)}
transcripts = transcriptfile %>% import

message('getting target gene')

vgene = transcripts%>%subset(gene_name%>%str_detect('Pa2g4') & (type=='gene'))
gchr = vgene%>%seqnames
gstart = (vgene%>%start)
gend = (vgene%>%end)

#create the gene track for gviz
gtrack = transcripts%>%{
  transcripts = .
  transsubset=transcripts%>%subsetByOverlaps(vgene)
  transsubset$feature = transsubset$type %>% as.character
  transsubset<-transsubset[ transsubset$feature %in% c("CDS","UTR")]
  transsubset$exon = transsubset$exon_id
  transsubset$gene = transsubset$gene_name
  transsubset$gene = transsubset$transcript_name
  transsubset$symbol = transsubset$transcript_name
  transsubset$feature = transsubset$feature %>% as.character
  transsubset$transcript= transsubset$transcript_id
  # transsubset$group = transsubset$transcript_id 
  gtrack = Gviz::GeneRegionTrack(transsubset,thinBoxFeature=c("UTR"),showId=TRUE,chromosome=gchr,geneSymbol=TRUE)
}

#function
make_track<-function(file,trackname,isneg=FALSE,ylims=NULL,...){
  require(Gviz)
  if(str_detect(trackname,'total')) isneg=TRUE
  if(isneg){ 
    importfunction = import.bw.neg
  }else{
    importfunction= function(file,selection) import.bw(file,sel=selection)
  }
  DataTrack(file,name = trackname,chromosome=gchr,stream=TRUE,importFunction = importfunction,ylim = ylims)
}


#now make track objects from each bigwig
bwTracks <- 
  datafiles%>%
  data_frame(file=.)%>%
  mutate(trackname = file%>%basename%>%str_replace('(.transcript)?.bw','') )%>%
  # filter(trackname%>%str_detect('ribo.*pos|(total.*neg)'))%>%
  {map2(.$file,.$trackname,.f = make_track)}

# atracks <-map(bam_files[c(1:4,17:20)],.f = .%>%AlignmentsTrack(.,isPaired=FALSE))
# 
tracks <- c(bwTracks,gtrack,codtrack,Gviz::GenomeAxisTrack())
plotTracks(tracks,chr = as.character(gchr), from = gend - 210, to = gend, type="h")


plotTracks(list(AlignmentsTrack(bam_files%>%str_subset('E13_ribo_1')),gtrack),chromosome = gchr, from = gend-300, to = gend,main='E13_ribo_1')
plotTracks(AlignmentsTrack(bam_files%>%str_subset('E16_ribo_1'),gtrack),chromosome = gchr, from = gend-300, to = gend,main='E16_ribo_1')

plotTracks(list(AlignmentsTrack(bam_files%>%str_subset('P0_ribo_1')),gtrack),chromosome = gchr, from = gend-300, to = gend,main='P0_ribo_1')
plotTracks(list(AlignmentsTrack(bam_files%>%str_subset('P0_ribo_2')),gtrack),chromosome = gchr, from = gend-300, to = gend,main='P0_ribo_2')
plotTracks(AlignmentsTrack(bam_files[17]),chromosome = gchr, from = gend-300, to = gend)


