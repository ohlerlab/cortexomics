library(rtracklayer)
library(Rsamtools)
library(GenomicFeatures)
fafile <- FaFile('my_GRCm38.p5.genome.chr_scaff.fa')

mappatrack <- import('mappability/mappability_27.bedgraph')

cds <- import('my_gencode.vM12.annotation.cds.gtf')
exons <- import('my_gencode.vM12.annotation.gtf')%>%subset(type=='exon')

cds<-cds%>%mapToTranscripts(exons%>%split(.$transcript_id))

seqinfo(mappatrack) <- seqinfo(cds)[
	c(seqinfo(mappatrack)@seqnames,setdiff(seqinfo(cds)@seqnames,seqinfo(mappatrack)@seqnames)),]

mappatrack%<>%coverage(weight='score')

#cds$nomapbases <- 


cdssamp <- cds%>%sample(1000)

cut_number( mappatrack[cdssamp]%>%sum / width(cdssamp) ,5)%>%table


cdsmap_df <- cds%>%as.data.frame%>%
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

pdfexpr(mappahistplotfile,hist((floor(cdsmap_df$nomap_frac/0.1)*0.1),50))


cdsmap_df$nomap_frac%>%`>`(0.05)%>%mean


cds%>%subsetByOverlaps(GRanges('chr1:195132153-195146569'))
cds%>%subsetByOverlaps(GRanges('chr1:195132153-195146569'))%>%width%>%sum
25/503

mappability[]