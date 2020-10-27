
#function that takes all the granges in a region, and compresses them down so introns
#are all a fixed distance
library(Gviz)
library(rtracklayer)
library(tidyverse)
library(magrittr)
source(here::here('src/R/Rprofile.R'))

mymem<-projmemoise(function(){rnorm(1)})



mymem()
mycache$keys()
projmemoise

rtracklayer::import('/fast/groups/ag_ohler/work/dharnet_m/cortexomics/pipeline/my_gencode.vM12.annotation.gtf')

stop()

options(ucscChromosomeNames=FALSE) 

#toy ranges to shift
exrangestoshift <- c(GRanges('a',IRanges(4,8)),GRanges('a',IRanges(12,13)),GRanges('a',IRanges(20,30)),GRanges('a',IRanges(40,42),strand='-'))
#

#this is the function that shifts a bunch of ranges to minimize the space between them, down to a ceiling maxwidth
maxwidth=1
get_comp_shift<-function(rangestoshift,maxwidth=1){
	rangegaps <- rangestoshift%>%{strand(.)<-'*';.}%>%gaps%>%.[-1]
	rangegaps%>%width
	#maxwidth <- 4
	toshift<-rep(0,length(rangestoshift))
	i=3
	for(i in seq_along(rangegaps)){
		adj <- width(rangegaps[i])-maxwidth
		if(adj>0){
			rightofgap <- start(rangestoshift)>end(rangegaps[i])
			toshift[rightofgap] <- toshift[rightofgap] + adj
		}
	}
	toshift %<>% {. + min(start(rangestoshift))-1}
	#rangestoshift<-GenomicRanges::shift(rangestoshift,-rangestoshift$toshift)
	toshift
}
compress_gaps<-function(gr,maxwidth){GenomicRanges::shift(rangestoshift,-get_comp_shift(rangestoshift,maxwidth))}
rangestoshift

compress_gaps(rangestoshift,1)

rangestoshift%>%gaps%>%width
compress_gaps(rangestoshift,2)%>%{strand(.)<-'*';.}%>%gaps%>%width

shiftedranges = compress_gaps(rangestoshift,1)

#generate fake data
cov = coverage(exrangestoshift)
covch = cov[['a']]
covch[covch!=0] %<>% {.*rpois(length(.), as.vector(.) * 32)}
cov[['a']]=covch
fakecov = as(cov,'GRanges')

#Now we need to shift this over to the compressed space
ov = findOverlaps(fakecov,exrangestoshift)
fakecov = shift(fakecov[ov@from],-shiftedranges$toshift[ov@to])





########
library(rtracklayer)

fread('head -n 100 pipeline/my_gencode.vM12.annotation.gtf')%>%
	as.matrix%>%
	apply(1,paste,collapse='\t')%>%
	lapply(paste,collapse='\n')%>%
	{import(text=.,format='gtf')}

assert_that(file.exists(annofile))
# annotation <- annofile%>% {rtracklayer::import(.)}


#get counts for the exons
bam=dexseqbams[[1]]
windows = annotation

#get offset values
offsets<-read_tsv(here('ext_data/offsets_manual.tsv'))

#get the bams for each library
bams <- Sys.glob(here('pipeline/star/data/*/*.bam'))%>%str_subset(negate=TRUE,'transcript')
#
dexseqbams<-bams%>%str_subset('_(ribo)')
dexseqbamsrna<-bams%>%str_subset('_(total)')

getwd()
annotation <- fread('grep -ie "satb2" /fast/groups/ag_ohler/work/dharnet_m/cortexomics/pipeline/my_gencode.vM12.annotation.gtf')%>%
	as.matrix%>%apply(1,paste0,collapse='\t')%>%paste0(collapse='\n')%>%import(text=.,format='gtf')%>%
	subset(gene_name=='Satb2')

annotation$toshift <- get_comp_shift(annotation)


exontrack =
  annotation%>%
  shift(.,-.$toshift)%>%
  # {.$transcript_id %<>% str_replace('\\.[0-9]+$','');.}%>%
  subset(type%in%c('UTR','CDS'))%>%
  .[,'type']%>%
  {.$feature<-as.character(.$type);.}%>%
  Gviz::GeneRegionTrack(.,thinBoxFeature=c("UTR"))
  
library(Gviz)

bams=dexseqbams

get_psite_datatrack_comp <- function(bam,annotation,offsets){
	#get the reads
	browser()
	psites = get_genomic_psites(bam,annotation,offsets)
	#collapse the reads into a psite coverage track
	psites = table(start(psites))%>%{GRanges(seqnames(psites)[1],IRanges(as.numeric(names(.)),w=1),score=as.vector(.))}
	#	
	ov = findOverlaps(psites,annotation)
	psites = shift(psites[ov@from],-psites$toshift[ov@to])
	#
  	st = as.character(strand(annotation[1]))
  	#
  	psites%>%subset(strand==st)
}

#now plot
options(ucscChromosomeNames=FALSE)
plotfile<- here(paste0('plots/','gviztest','.pdf'))
pdf(plotfile,w=21,h=7)
plotTracks(list(
	map(bams, .f = get_psite_datatrack_comp,annotation,offsets),
 	exontrack
))
dev.off()
normalizePath(plotfile)


get_genomic_psites(dexseqbams[[1]],annotation,offsets)%>%subset(strand=='-')%>%coverage%>%as('GRanges')%>%subset(score!=0)%>%coverage%>%sum%>%sum
