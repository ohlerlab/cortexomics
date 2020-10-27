library(rtracklayer)
library(Rsamtools)
library(GenomicFeatures)
library(here)
library(tidyverse)
library(magrittr)
#fasta files
fafile <- FaFile(here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'))

#load the cds
cds <- import(here('pipeline/my_gencode.vM12.annotation.cds.gtf'))
exons <- import(here('pipeline/my_gencode.vM12.annotation.gtf'))%>%subset(type=='exon')
#we need the cdsd to be 
#load the tracks
unmappatrack <- import(here('pipeline/mappability/mappability_27.bedgraph'))
#add CDS
seqinfo(unmappatrack) <- seqinfo(cdsmap)[
	c(seqinfo(unmappatrack)@seqnames,setdiff(seqinfo(cdsmap)@seqnames,seqinfo(unmappatrack)@seqnames)),]
unmappatrack%<>%coverage(weight='score')

bestcds<-cds%>%split(.,.$protein_id)%>%.[best_protein_ids]%>%unlist
bcdstr<-bestcds%>%{pmapToTranscripts(.,split(exons,exons$transcript_id)[.$transcript_id])}
bcdstr <- GenomicRanges::reduce(bcdstr,with.revmap=TRUE)
mcols(bcdstr) <- mcols(bestcds)[bcdstr$revmap%>%{.@unlistData[.@partitioning@end]},]
names(bcdstr) <- bcdstr$protein_id
unmappacovtracks <- unmappatrack[bcdstr]%>%setNames(bcdstr$protein_id)
bestcdsmappab<-unmappacovtracks%>%mean

unmapgr <- unmappacovtracks[is3nt]%>%as("GRanges")
end(unmapgr) = ceiling(end(unmapgr)/3)*3
start(unmapgr) = ifelse(start(unmapgr)==1,1,(lag(end(unmapgr),1,0)+1))
unmappacovtracks =unmapgr%>%coverage(weight='score')


stop()

cdssamp <- cdsmap%>%sample(1000)


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


#ideally what we want to do is include whole codons, if and only if the first basepair of the
#codon is mappable.

n3pid = 'ENSMUSP00000114761'
bestcds%>%.[names(.)==n3pid]%>%
bestcds%>%.[names(.)==n3pid]%>%width
bestcds%>%.[names(.)==n3pid]%>%getSeq(x=fafile)%>%translate%>%as.character
bestcds%>%subset(protein_id==n3pid)%>%{pmapToTranscripts(.,split(exons,exons$transcript_id)[.$transcript_id])}%>%width
bcdstr%>%subset(protein_id==n3pid)%>%width%>%`%%`(3)
is3nt = psitecov%>%runLength%>%sum%>%`%%`(3)%>%`==`(0)







################################################################################
########Plot mappability vs Entropy
################################################################################
	
fread('pipeline/mappability/mappability_27.bedgraph')
library(LSD)
base::source('https://raw.githubusercontent.com/stineb/LSD/master/R/LSD.heatscatter.R')
#now plot
plotfile<- here(paste0('plots/','pme_vs_mappability','.pdf'))
pdf(plotfile)
	heatscatter(pmes,bestcdsmappab,ggplot=TRUE)+
	ggtitle(paste0(scattertitle))+
	scale_x_continuous(paste0('% Maximum Entropy'))+
	scale_y_continuous(paste0('% Unmappable Bases'))+
	ggtitle(paste0('PME vs Mappability'))+
	theme_bw()
dev.off()
normalizePath(plotfile)

bestcdsmappab

