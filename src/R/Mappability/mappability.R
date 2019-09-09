library(rtracklayer)
library(Rsamtools)
library(GenomicFeatures)
library(here)
library(tidyverse)
library(magrittr)
#fasta files
fafile <- FaFile('my_GRCm38.p5.genome.chr_scaff.fa')

#load the cds
cds <- import(here('pipeline/my_gencode.vM12.annotation.cds.gtf'))
exons <- import(here('pipeline/my_gencode.vM12.annotation.gtf'))%>%subset(type=='exon')
#we need the cdsd to be 
cdsmap<-cds%>%setNames(.,as.character(.))%>%mapToTranscripts(exons%>%split(.$transcript_id))
mcols(cdsmap) <- mcols(cds)[cdsmap$xHits,]
#load the tracks
mappatrack <- import(here('pipeline/mappability/mappability_27.bedgraph'))
#add CDS
seqinfo(mappatrack) <- seqinfo(cdsmap)[
	c(seqinfo(mappatrack)@seqnames,setdiff(seqinfo(cdsmap)@seqnames,seqinfo(mappatrack)@seqnames)),]

mappatrack%<>%coverage(weight='score')

cdssamp <- cdsmap%>%sample(1000)

stop()

cdsmap$nomapbases <- mappatrack[cdsmap]%>%mean
cdsmap$cds_id = paste0(seqnames(cdsmap),'_',start(cdsmap),'_',end(cdsmap))

data.frame(seqnames(cdsmap),names(cdsmap),cdsmap$nomapbases)


mcols(cds)[,c('nomapbases','cds_id')] %>% as.data.frame%>%write_tsv('pipeline/mappability/cds')

stop('mappability - now save it')

if(interactive()){


cdsmap_df <- cdsmap%>%as.data.frame%>%
	select(transcript_id,nomapbases,width)%>%
	group_by(gene_id,transcript_id)%>%
	summarise(width = sum(width),nomapbases = sum(nomapbases))


cdsmap_df%<>%mutate(nomap_frac = nomapbases / width)

pdfexpr<-function(file,expr,...){
	dir.create(dirname(file))
	pdf(file,...)
	expr
	dev.off()
	message(normalizePath(file))
}
mappahistplotfile <- here('plots/mapppbility/mappability.pdf')
pdfexpr(mappahistplotfile,hist((floor(cdsmap_df$nomap_frac/0.1)*0.1),50))

}
