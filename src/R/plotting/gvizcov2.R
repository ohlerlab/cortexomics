library(Gviz)
options(ucscChromosomeNames=FALSE)
mypdfdev <- function(filename,...){
	if(tools::file_ext(filename)!='pdf'){
		filename <- paste0(filename,'.pdf')
	}
	dir.create(dirname(filename),showWarnings=FALSE)
	message(normalizePath(filename))
	pdf(filename,...)
}
mysvglitedev <- function(filename,...){
	if(tools::file_ext(filename)!='svg'){
		filename <- paste0(filename,'.svg')
	}
	dir.create(dirname(filename),showWarnings=FALSE)
	message(normalizePath(filename))
	svglite(filename,...)
}

# datafiles <- allbigwigs

# make_track<-function(file, trackname, sizefactor, vtr,isneg=FALSE,...){

#   selection =   GRanges(vtr$transcript_id,IRanges(1,sum(width(vtrexons))))

#   if(str_detect(trackname,'total')) isneg=TRUE
#   if(isneg){ 
#     importfunction = import.bw.neg
#     # ylims = c(-30,0)
#   }else{
#     importfunction= function(file,selection) import.bw(file,sel=selection)
#     # ylims = c(0,10)
#   }
#   DataTrack(file,name = trackname,chromosome=vtr$transcript_id,stream=TRUE,importFunction = importfunction)

# }

# ribotracks <- datafiles%>%
#   data_frame(file=.)%>%
#   mutate(trackname =extract_id(file,sampleids),
#   	strand=ifelse(str_detect(file,'neg'),'-','+')
#   )%>%
#   split(.,.$strand)%>%
#   map(.%>%{map2(.$file,.$trackname,.f = make_track)}%>%setNames(.,.$trackname))


options(ucscChromosomeNames=F)



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
	DataTrack(gr,name=	datname,chr=rangename,groups=groupvect,col = rgbvect, fill=rgbvect, cex.title=0.3,legend=FALSE,genome='transcript')
}

cellnames <- allbigwigs%>%dirname%>%basename%>%str_extract('[^_]+')%>%unique
#this is ugly but the samples are badly named, no time to fix
cellnames%<>%grep(v=TRUE,inv=TRUE,patt='OMM475')



genome='mine'
genomefile='pipeline/my_hg38.fa'


# all_linc_orfs <- satannorfs%>%map(.%>%.$ORFs_tx%>%subset(clipid(seqnames) %in% clipid(linctranscripts$transcript_id)))
all_linc_orfs <- satannorfs%>%map(.%>%.$ORFs_tx%>%subset((seqnames) %in% (linctranscripts$transcript_id)))



vmatchPattern('SYLRRHLDF',allprotseqs)




#seqinfo object for orftrack
transcript_seqinfo <- exons[unique(linc_orfs_dt_with_linctable$transcript_id)]%>%
	width%>%sum%>%{Seqinfo(seqnames=names(.),seqlength=.)}

#track with the orf's in transcript space
orftrack<- linc_orfs_dt_with_linctable%>%
	DT2GR(seqinf=transcript_seqinfo)%>%
	{mcols(.)$contains_peptide <- ifelse(.$contains_peptide,'Peptide+','-');.} %>%
	unique%>%
	{Gviz::AnnotationTrack(
		name='ORFs',.,
		feature=sort(c('red','green','blue'))[(start(.)%%3)+1],red='red',green='green',blue='blue',
		id=.$contains_peptide,showFeatureId=TRUE)}

displayPars(orftrack) <- list(height=2)




transcript_orfplot <- function(ranges_with_goodorf,
	orftrack, # made with the orf dt
	bigwigpairlist #list of pairs (strands) of wigs in genome space
	){}


for (rangename in names(ranges_with_goodorf)){
	#get the rnage we need
	range = ranges_with_goodorf[[rangename]]


	rangename%<>%str_replace_all('[\\[\\]]','')

	plottitle <- str_interp("Riboseq Read Profile for:\n${rangename} = ${range}\n")		

	
	#construct track showing our exons
	exontrack<-
		GRanges(rangename,IRanges(c(1,cumsum(width(range))[-length(range)]+1),cumsum(width(range))))%>%
		{AnnotationTrack(name='exons',col='black',fill='yellow',strand='*',genome=genome,id=paste0('exon_',seq_along(.)),showFeatureId=TRUE,
			chr=rangename)}


	#name the plot file
	transprofplotfile <- str_interp('plots/loci_riboprofiles/${rangename}_prof')
#		transprofplotfile%<>%str_replace_all('[^a-zA-z0-9_\\./]+','')

	orfs <- orftrack@range%>%subset(seqnames%in%rangename)
	
	#calculate plot borders
	plotstart <- orfs$start%>%min
	plotend <- orfs$end%>%max
	orfswidth <- plotend-plotstart
	plotstart = max(0,plotstart - (orfswidth))
	rangelength = sum(width(range))
	plotend = min(rangelength,plotend + (orfswidth))
	plotgenomewindow <- GRanges(rangename,IRanges(plotstart,plotend))%>%split(.,seqnames(.))

	#select which tracks to show
	whichwigs <- bigwigpairlist%>%names%>%
		extract_oneof(cellnames)%>%
		is_in(extract_oneof(orfs$sample,cellnames))

	#get the profile track in transcript space
	ribotracks <- bigwigpairlist[whichwigs]%>%
		map(~get_riboproftrack(plotgenomewindow,.,rangename))

	
	# add_legendtrack <- function(ribotracks){
	# 	l<-length(ribotracks)

	# 	for(i in seq_along(ribotracks))displayPars(ribotracks[[i]]) <- list(height=1)
	# 	displayPars(ribotracks[[l]]) <- list(height=1.2)
	# 	displayPars(ribotracks[[l]]) <- list(legend=TRUE,lineheight.legend=0.2)

	# 	ribotracks
	# }
	# ribotracks %<>% add_legendtrack

	#open pdf file	
	mypdfdev(transprofplotfile)
	
	plotTracks(main=plottitle,cex.main=0.5,
		from=plotstart,to=plotend,#zoomed in on the orf in question
		c(
			ribotracks, # plot the riboseq signal
			exontrack,
			orftrack,
			GenomeAxisTrack(range=GRanges(rangename,IRanges(0,rangelength)),id=rangename)
		),
		type='hist',
		chr=seqnames(range)
	)
	dev.off()


}