
#function that takes all the granges in a region, and compresses them down so introns
#are all a fixed distance




rangestoshift <- c(GRanges('a',IRanges(4,8)),GRanges('a',IRanges(12,13)),GRanges('a',IRanges(20,30)),GRanges('a',IRanges(40,42),strand='-'))
#

narrow_gaps<-function(rangestoshift,maxwidth=1){
	rangegaps <- rangestoshift%>%{strand(.)<-'*';.}%>%gaps
	maxwidth <- 4
	rangestoshift$toshift<-0
	for(i in seq_along(rangegaps)){
		adj <- width(rangegaps[i])-maxwidth
		if(adj>0){
			rightofgap <- start(rangestoshift)>end(rangegaps[i])
			rangestoshift$toshift[rightofgap] <- rangestoshift$toshift[rightofgap] - adj
		}
	}
	rangestoshift<-GenomicRanges::shift(rangestoshift,rangestoshift$toshift)
	# rangestoshift%>%{strand(.)<-'*';.}%>%gaps%>%width
	rangestoshift
}

               