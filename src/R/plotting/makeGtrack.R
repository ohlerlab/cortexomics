
transcriptfile <- '~/projects/cortexomics/data/static_local/gencode.vM12.annotation.gtf'
transcripts <- transcriptfile%>% {rtracklayer::import(.)}
transcripts = transcriptfile %>% import
message('getting gene track')
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
  gtrack = Gviz::GeneRegionTrack(transsubset,thinBoxFeature=c("UTR"),showId=TRUE,geneSymbol=TRUE)
}
