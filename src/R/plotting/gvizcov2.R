
library(Gviz)
library(here)
library(rtracklayer)
options(ucscChromosomeNames=FALSE)
select<-dplyr::select



make_gtrack<-function(annotationgr){
	nctrs <- annotationgr$transcript_type%>%split(annotationgr$transcript_id)%>%map(unique)%>%keep(Negate(identical),'protein_coding')%>%names
    
    annotationgr$feature = annotationgr$type %>% as.character
    annotationgr<-annotationgr[ (annotationgr$feature %in% c("CDS","UTR")) | annotationgr$transcript_id %in% 	nctrs]

    annotationgr$exon = annotationgr$exon_id
    annotationgr$gene = annotationgr$gene_name
    annotationgr$gene = annotationgr$transcript_name
    annotationgr$symbol = annotationgr$transcript_name
    annotationgr$feature = annotationgr$feature %>% as.character
    annotationgr$transcript= annotationgr$transcript_id
    annotationgr$group = annotationgr$transcript_id 
    nccolor<-ifelse(annotationgr$transcript_id %in% 	nctrs,'blue','orange')
     Gviz::GeneRegionTrack(annotationgr,thinBoxFeature=c("UTR"),showId=TRUE,geneSymbol=TRUE,fill=nccolor)
}
gtrack =transcripts%>%subsetByOverlaps(vgene)%>%make_gtrack





library(GenomicFeatures)


if(!exists('allanno')) allanno<-here('pipeline/my_gencode.vM12.annotation.gtf')%>%import
bams <- Sys.glob(here('pipeline/star/data/*/*.bam'))%>%str_subset(neg=T,'transcript|mappability|accepted|unmapped')


bamname<-'foo'
bamfile=bams%>%str_subset('E13_ribo_1')
genegr=vexonsred
genegr%>%width%>%sum
genegr%>%width

#so pa2g4 is 
#rapidly gets a covergae profile using bamsignals. 
get_Profile_gr<-function(genegr,bamfile,signame=NULL,rangename=NULL,nolift = FALSE, score_negative=FALSE,invstrand=FALSE,logtrans=TRUE){
	require(bamsignals)
	require(tools)

	# if(is.null(names(genegr))) rangename = ''
	if(is.null(signame)) signame = bamfile%>%file_path_sans_ext%>%basename

	stopifnot(genegr ==sort(genegr))
	stopifnot(genegr ==reduce(genegr))
	stopifnot(n_distinct(unique(seqnames(genegr)))==1)
	stopifnot(n_distinct(unique(strand(genegr)))==1)
	
	# strand(genegr) <-'-'
	if(invstrand) genegr <- invertStrand(genegr)

	isneg <- '-' %in% strand(genegr)
	


	# signal <-  bamsignals::bamProfile(bamfile,genegr%>%tail(2)%>%head(1)%>%resize(4,'start'),ss=TRUE)%>%as.list%>%{if(isneg) rev(.) else . }%>%Reduce(f=cbind)%>%.['sense',]%>%Rle
	signal <-  bamsignals::bamProfile(bamfile,genegr,ss=TRUE)%>%
		as.list%>%
		{if(isneg) rev(.) else . }%>%
		Reduce(f=cbind)%>%
		.['sense',]%>%
		Rle

	if(invstrand) signal = rev(signal)

	nzsig <- which(signal!=0)
	genomecoords<-GRanges(rangename,IRanges(nzsig,w=1),score=signal[nzsig])

	if(logtrans) genomecoords$score = sign(genomecoords$score) * log2(abs(genomecoords$score))
	if(score_negative) genomecoords$score = genomecoords$score * -1

	# if(logtrans) genomecoords$score = log2(abs(genomecoords$score)+1)
		# browser()
	if(nolift) return(DataTrack(genomecoords,name=signame))

	genomecoords <- mapFromTranscripts(genomecoords,GRangesList(genegr)%>%setNames(rangename))
	
	stopifnot(genomecoords$xHits == seq_along(nzsig))
	mcols(genomecoords) <- DataFrame(score=signal[nzsig])

	DataTrack(genomecoords,name=signame)

}

# gname <- vgene$gene_name
gname <- 'Orc3'

vgene <- allanno%>%subset(gene_name==gname)%>%subset(type=='gene')




#Get the annoation track lifted over

subanno<-subsetByOverlaps(allanno,vgene,ignore.strand=T)
vexonsred <-subanno%>%subset(gene_name==gname)%>%subset(type=='exon')%>%reduce
liftanno <- GenomicFeatures::mapToTranscripts(subanno,GRangesList(list(vexonsred))%>%setNames(gname))
strand(liftanno)<-'+'
stopifnot(liftanno$xHits == 1:length(liftanno))
mcols(liftanno) <- mcols(subanno)[liftanno$xHits,]
gtrack <- liftanno%>%
	subset(!type%in%c('gene','transcript'))%>%
	make_gtrack





# subsetByOverlaps(liftanno,GRanges(gname,500))

#  liftanno%>%
# 	subset(!type%in%c('gene','transcript'))%>%
# 	subsetByOverlaps(GRanges(gname,500))



testanno<-gtrack@range%>%
	subsetByOverlaps(GRanges(gname,3750))
testanno
liftanno%>%subset(transcript_id%in%testanno$transcript)
subanno%>%subset(transcript_id%in%testanno$transcript)

# subanno%>%subset(type=='exon')%>%gaps%>%width


bamdf<-fread(here('pipeline','sample_parameter.csv'))%>%
	filter(!fraction %in% 'input')%>%
	filter(!time%in%c('E145','E175'))%>%
	filter(!str_detect(group,'test'))%>%
	arrange(time,-is.na(fraction),str_extract(sample_id,'_\\d'),assay=='total')

bamlist<-bamdf$sample_id%>%setNames(.,.)%>%
	map_chr(~paste0('pipeline/star/data/',.,'/',.,'.bam' ))%T>%
	{stopifnot(file.exists(.))}

bamtracks <- bamlist%>%imap(function(x,nm) {

	if(bamdf$assay[bamdf$sample_id==nm]=='total'){
		get_Profile_gr(x,genegr=vexonsred,nolift=TRUE,rangename=gname,logtrans=T,score_negative=TRUE,invstrand=TRUE)
	}else{
		get_Profile_gr(x,genegr=vexonsred,nolift=TRUE,rangename=gname,logtrans=T)

	}
})

tracks<-bamtracks

normalize_tracks<-function(tracks,sizefacts){
	for(i in seq_along(tracks)){
		tracks[[i]]@data%<>%divide_by(sizefacts[[i]])
	}
	tracks
}

transform_abs_vals<-function(tracks){
	for(i in seq_along(tracks)){
		tracks[[i]]@data%<>%{sign(.) * log2(abs(.))}
	}
	tracks
}



equal_y_axes<-function(tracks){
	datarange <- tracks%>%map(~.@data%>%range)

	max_y <- datarange%>%unlist%>%abs%>%max%>%divide_by(10)%>%ceiling%>%multiply_by(10)
	stopifnot(max_y>0)
	
	datarange_has_zero <- datarange%>%map_lgl(~ 0 %in% .)
	stopifnot(datarange_has_zero)
	isnegtrack <- datarange%>%map_lgl(~ min(.) < 0)
	
	for(i in seq_along(tracks)){
		if(isnegtrack[i]){
			displayPars(tracks[[i]])$ylim <- c(-max_y,0)
		}else{
			displayPars(tracks[[i]])$ylim <- c(0,max_y)
		}
	}
	tracks
}

sizefactors <- read_tsv(here('pipeline','sizefactors.tsv'))
sizefactors <- sizefactors$sizefactor%>%setNames(sizefactors$sample_id)

# bamtracks <- normalize_tracks(bamtracks,sizefacts=sizefactors[names(bamtracks)])
# bamtracks <- transform_abs_vals(bamtracks)
bamtracks <- equal_y_axes(bamtracks)

tracks <- c(bamtracks,gtrack,	GenomeAxisTrack(GRanges(gname,IRanges(1,width(vgene)))))



# pfile <- here('plots','coverage','tmp.pdf')
pfile <- here('plots','coverage',paste0(gname,'.pdf'))
pdf(pfile,w=14*2,h=14*2)
plotTracks(main=gname,cex.main=1,
	from=1,to=sum(width(vexonsred)),#zoomed in on the orf in question
	tracks,
	type='hist',
	col.labels='black',
	chr=gname
	# transformation = function(x) {log2(x+1)}
)
dev.off()
message(normalizePath(pfile))
'/fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/coverage/tmp.pdf'
'/fast/work/groups/ag_ohler/dharnet_m/cortexomics/plots/coverage/Orc3.pdf'

