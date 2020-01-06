
library(Gviz)
options(ucscChromosomeNames=FALSE)
select<-dplyr::select


displayPars(orftrack) <- list(height=6)
displayPars(peptide_track) <- list(height=6)


mypdfdev <- function(filename,...){
	if(tools::file_ext(filename)!='pdf'){
		filename <- paste0(filename,'.pdf')
	}
	dir.create(dirname(filename),showWarnings=FALSE)
	
	pdf(filename,...)
	message(normalizePath(filename))
}
mysvglitedev <- function(filename,...){
	if(tools::file_ext(filename)!='svg'){
		filename <- paste0(filename,'.svg')
	}
	dir.create(dirname(filename),showWarnings=FALSE)
	message(normalizePath(filename))
	svglite(filename,...)
}



get_riboproftrack<- function(exons_tr,bigwigpair,rangename){

	stopifnot(c('+','-') %in% names(bigwigpair))

	# for(i in ls()) assign(i,get(i),envir=.GlobalEnv)
	profilegrange <- 
		suppressWarnings(lapply(bigwigpair,import,which = unlist(exons_tr)))%>%
		# suppressWarnings(lapply(bigwigpair,import,which = unlist(testwindow)))
		{for(strandi in names(.)) strand(.[[strandi]]) <- strandi;.}%>%
		{suppressWarnings(Reduce(f=c,.))}%>%
		# subsetByOverlaps(testwindow)%>%
		subset(score>0)
	# profilegrange<<-profilegrange
	#now map our prifle data to the exons space
	gr <- suppressWarnings(mapToTranscripts(profilegrange,exons_tr))

	# (mapToTranscripts(testwindow,exons_tr))
	if(length(gr)==1) gr <- c(gr,gr)
	gr$score<-profilegrange$score[gr$xHits];
	# gr%>%subset(score>30)
	# profilegrange%>%subset(score>30)

	datname = bigwigpair%>%extract_id%>%.[[1]]
	
	scoremat <- rep(0,length(gr)*3)%>%matrix(ncol=3)
	scoremat[matrix( c(1:length(gr),(start(gr)%%3)+1 ) ,ncol=2 ) ] <- gr$score
	mcols(gr) = scoremat
	rgbvect <- c('red','green','blue')%>%sort

	if(length(gr)==0) {groupvect <- NULL}else{
		# groupvect <- groupvect[(start(gr)%%3)+1]
		groupvect <- paste0('frame ',1:3)
	} 
	DataTrack(gr,name=datname,chr=rangename,groups=groupvect,col = rgbvect, fill=rgbvect, cex.title=0.3,legend=FALSE,genome='transcript')
}


transcript_orfplot <- function(trname,trcellname,
	exons,
	orftrack, # made with the orf dt
	peptide_track,
	bigwigs2plot #list of pairs (strands) of wigs in genome space
	){
	range = exons[[trname]]


	trname%<>%str_replace_all('[\\[\\]]','')

	plottitle <- str_interp("Riboseq Read Profile for:\n${trname} = ${range}\n")		

	
	#construct track showing our exons
	exontrack<-
		GRanges(trname,IRanges(c(1,cumsum(width(range))[-length(range)]+1),cumsum(width(range))))%>%
		{AnnotationTrack(.,name='exons',col='black',fill='yellow',strand='*',genome=genome,id=paste0('exon_',seq_along(.)),
			showFeatureId=F,
			chr=trname)}


	#name the plot file
	transprofplotfile <- str_interp('./plots/loci_riboprofiles/${trname}_${trcellname}_prof')
#		transprofplotfile%<>%str_replace_all('[^a-zA-z0-9_\\./]+','')

	orfs <- orftrack@range%>%subset(seqnames%in%trname)%>%GR2DT
	
	#calculate plot borders
	plotstart <- orfs$start%>%min
	plotend <- orfs$end%>%max
	orfswidth <- plotend-plotstart
	plotstart = max(0,plotstart - (orfswidth))
	rangelength = sum(width(range))
	plotend = min(rangelength,plotend + (orfswidth))
	plotgenomewindow <- GRanges(trname,IRanges(plotstart,plotend))%>%split(.,seqnames(.))

	#get the profile track in transcript space
	ribotracks <- bigwigs2plot%>%
		map(~get_riboproftrack(setNames(GRangesList(range),trname),.,trname))

	max_y <- ribotracks%>%map(~.@data[,overlapsAny(.@range,orftrack@range,ign=TRUE)]%>%max)%>%unlist%>%max%>%divide_by(10)%>%ceiling%>%multiply_by(10)
	max_y=log2(max_y)

	for(i in seq_along(ribotracks)){
		displayPars(ribotracks[[i]])$ylim <- c(0,max_y)

	}
	displayPars(peptide_track)$col.id <- 'black'
	# add_legendtrack <- function(ribotracks){
	# 	l<-length(ribotracks)

	# 	for(i in seq_along(ribotracks))displayPars(ribotracks[[i]]) <- list(height=1)
	# 	displayPars(ribotracks[[l]]) <- list(height=1.2)
	# 	displayPars(ribotracks[[l]]) <- list(legend=TRUE,lineheight.legend=0.2)

	# 	ribotracks
	# }
	# ribotracks %<>% add_legendtrack
	chromosome(orftrack)<-trname
	chromosome(peptide_track)<-trname
	
	#open pdf file	
	mypdfdev(transprofplotfile)
	
	plotTracks(main=plottitle,cex.main=0.3,
		from=plotstart,to=plotend,#zoomed in on the orf in question
		c(
			ribotracks, # plot the riboseq signal
			exontrack,
			orftrack,
			peptide_track,
			GenomeAxisTrack(range=GRanges(trname,IRanges(0,rangelength)),id=trname)
		),
		type='hist',
		col.labels='black',
		chr=seqnames(range),
		transformation = function(x) {log2(x+1)},
		NULL
	)
	dev.off()

}

#seqinfo object for orftrack
transcript_seqinfo <- exons[unique(peptide_hits_df$transcript_id)]%>%
	width%>%sum%>%{Seqinfo(seqnames=names(.),seqlength=.)}

#track with the orf's in transcript space
orftrack<- peptide_hits_df%>%
	select(seqnames=transcript_id,start=orfstart,end=orfend,ORF_id_tr)%>%
	mutate(strand='+')%>%
	makeGRangesFromDataFrame(keep=T)%>%
#	{x=makeGRangesFromDataFrame(.);mcols(x)<-cbind(mcols(x),select(.,-seqnames,-start,-end));x}%>%
	{seqinfo(.)<-transcript_seqinfo[seqlevels(.)];.}%>%
	unique%>%
	# {mcols(.)$contains_peptide <- ifelse(.$contains_peptide,'Peptide+','-');.} %>%
	{Gviz::AnnotationTrack(
		name='ORFs',.,
		feature=sort(c('red','green','blue'))[(start(.)%%3)+1],red='red',green='green',blue='blue',
		id=.$ORF_id_tr,showFeatureId=TRUE,cex=0.3,fontcolor.item="black")}

#track with the orf's in transcript space
peptide_track<- peptide_hits_df%>%
	select(seqnames=transcript_id,start=start,end=end,peptide)%>%
	mutate(strand='+')%>%
	makeGRangesFromDataFrame(keep=T)%>%
#	{x=makeGRangesFromDataFrame(.);mcols(x)<-cbind(mcols(x),select(.,-seqnames,-start,-end));x}%>%
	{seqinfo(.)<-transcript_seqinfo[seqlevels(.)];.}%>%
	# {mcols(.)$contains_peptide <- ifelse(.$contains_peptide,'Peptide+','-');.} %>%
	unique%>%
	{Gviz::AnnotationTrack(
		name='Peptide',.,
#		feature=sort(c('red','green','blue'))[(start(.)%%3)+1],red='red',green='green',blue='blue',
		id=.$peptide,showFeatureId=TRUE,cex=0.3,fontcolor.item="black")}



peptide_trs <- peptide_hits_df$transcript_id%>%unique

# for (trname in addtr){
for (trname in peptide_trs){
	trcellnames <- peptide_hits_df%>%
		filter(transcript_id==trname)%>%.$sample%>%
		extract_oneof(cellnames)%>%unique

	trcellnames%<>%str_subset('OD5P')
	
	for(trcellname in trcellnames){
		#select which tracks to show

		bigwigs2plot <- bigwigpairlist%>%.[names(.)%>%str_detect(trcellname)]


		transcript_orfplot(trname,trcellname,exons,orftrack,peptide_track,bigwigs2plot)
	}
}



