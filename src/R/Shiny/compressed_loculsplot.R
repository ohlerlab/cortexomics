
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
rangestoshift <- c(GRanges('a',IRanges(4,8)),GRanges('a',IRanges(12,13)),GRanges('a',IRanges(20,30)),GRanges('a',IRanges(40,42),strand='-'))
#

#this is the function that shifts a bunch of ranges to minimize the space between them, down to a ceiling maxwidth
maxwidth=1
narrow_gap<-function(rangestoshift,maxwidth=1){
	rangegaps <- rangestoshift%>%{strand(.)<-'*';.}%>%gaps%>%.[-1]
	rangegaps%>%width
	#maxwidth <- 4
	rangestoshift$toshift<-0
	i=3
	for(i in seq_along(rangegaps)){
		adj <- width(rangegaps[i])-maxwidth
		if(adj>0){
			rightofgap <- start(rangestoshift)>end(rangegaps[i])
			rangestoshift$toshift[rightofgap] <- rangestoshift$toshift[rightofgap] + adj
		}
	}
	rangestoshift$toshift %<>% {. + min(start(rangestoshift))-1}
	rangestoshift<-GenomicRanges::shift(rangestoshift,-rangestoshift$toshift)
  rangestoshift
}

rangestoshift

narrow_gap(rangestoshift,1)

rangestoshift%>%gaps%>%width
narrow_gap(rangestoshift,2)%>%{strand(.)<-'*';.}%>%gaps%>%width

shiftedranges = narrow_gap(rangestoshift,1)

#generate fake data
cov = coverage(exrangestoshift)
covch = cov[['a']]
covch[covch!=0] %<>% {.*rpois(length(.), as.vector(.) * 32)}
cov[['a']]=covch
fakecov = as(cov,'GRanges')

#Now we need to shift this over to the compressed space
ov = findOverlaps(fakecov,exrangestoshift)
fakecov = shift(fakecov[ov@from],-shiftedranges$toshift[ov@to])

memoise
plotTracks(list(
  DataTrack(fakecov,type='hist'),
  Gviz::AnnotationTrack(shiftedranges)
))

