library(Gviz)
library(GenomicAlignments)
library(GenomicFiles)

options(ucscChromosomeNames=FALSE)

annofile <- 'my_gencode.vM12.annotation.gtf'
txdb <- makeTxDbFromGFF(annofile)
generegiontrack <- GeneRegionTrack(txdb,geneSymbol=TRUE)

anno <- rtracklayer::import(annofile,which)

genes<-rtracklayer::readGFF(annofile, tags=c('gene_name','gene_id'), filter=list(type='gene'))%>%
	makeGRangesFromDataFrame( keep.extra.columns=TRUE)
names(genes)<-genes$gene_id	

exons<-rtracklayer::readGFF(annofile, tags=c('gene_name','gene_id','transcript_id','transcript_name'), filter=list(type='exon'))%>%
	makeGRangesFromDataFrame( keep.extra.columns=TRUE)

transcripts<-rtracklayer::readGFF(annofile, tags=c('gene_name','gene_id','transcript_id','transcript_name'), filter=list(type='transcript'))%>%
	makeGRangesFromDataFrame( keep.extra.columns=TRUE)
transcripts%<>%setNames(.,.$transcript_id)


monobams <- 'star/data/RPI6_80SE16_*/RPI6_80SE16_*.bam'%>%Sys.glob%>%{stopifnot(length(.)>0);.}
polybams <- 'star/data/RPI8_PolyE16_*/RPI8_PolyE16_*.bam'%>%Sys.glob%>%{stopifnot(length(.)>0);.}


govmonobam <- GenomicAlignments::summarizeOverlaps(ignore.strand=F,exons%>%split(.,.$transcript_id),BamFile(monobams))
govpolybam <- GenomicAlignments::summarizeOverlaps(ignore.strand=F,exons%>%split(.,.$transcript_id),BamFile(polybams))
govmonobam%<>%assay
govpolybam%<>%assay


countdf <- cbind(govmonobam,govpolybam)%>%as.data.frame%>%rownames_to_column('transcript_id')
countdf %<>% filter(RPI6_80SE16_2.bam > 3)

countdf%>%colnames%>%dput
c("RPI6_80SE16_2.bam", "RPI8_PolyE16_2.bam")

countdf %<>% mutate(monosomal = RPI6_80SE16_2.bam/RPI8_PolyE16_2.bam > 1 &(RPI6_80SE16_2.bam>20))

pdf('tmp.pdf')
countdf%>%qplot(data=.,x=log10(RPI8_PolyE16_2.bam+1),y=log10(RPI6_80SE16_2.bam+1),color=monosomal,geom='point')+
# scale_x_log10()+scale_y_log10()+
# geom_smooth(method='lm')
scale_color_discrete(name='')+
ggtitle('monosomal vs poly')
dev.off()


alltracks <- c(
	map(monobams,~ AlignmentsTrack(.,name=basename(.))),
	map(polybams,~ AlignmentsTrack(.,name=basename(.))),
	list(generegiontrack),
	IdeogramTrack(genome= 'mm10'),
	GenomeAxisTrack()
)



{
plotwindow<-transcripts['ENSMUST00000142608']
plotwindow <- countdf%>%filter(monosomal)%>%.$transcript_id%>%{transcripts[.]}%>%subset(width>10e3)%>%sample(1)
#generegiontrack <- GeneRegionTrack(txdb,geneSymbol=TRUE,chromosome=seqnames(plotwindow)[1],from=start(plotwindow)[1],to=end(plotwindow)[1])

#pdf('tmp.pdf'%>%normalizePath%>%{message(.);.},h=10,w=7)
plotTracks(chromosome=seqnames(plotwindow)[1],from=start(plotwindow)[1] - 1e3,to=end(plotwindow)[1]+1e3,strand=strand(plotwindow)[1],
	c(alltracks,AnnotationTrack(plotwindow)),sizes=c(2,2,2,2,1,1,1/2,0.5),
	type='coverage'
)
# dev.off()
message(plotwindow%>%names)
countdf%>%filter(transcript_id==plotwindow$transcript_id)
}

#ENSMUSG00000031448 - weird spike on one exon
#ENSMUST00000181301 - too damn big
#ENSMUST00000031627 - looks pretty even....
#ENSMUST00000159166 - this one looks funny, but actually probably just spike due to low counts
#ENSMUST00000168704 - promising maybe
#ENSMUST00000142327 - spike but in polysomal stuff
#ENSMUST00000054697 - another odd spike in last exon...
#ENSMUST00000085206 - Tracks pretty well... kiiind of a monosome spike at the start
#ENSMUST00000137553 - spike at the termination point
#ENSMUST00000019069 - meh
#start going for only 20+ mono her
#ENSMUST00000142608 - another spike at the termination point
#ENSMUST00000092834 - termination point
#ENSMUST00000106729 - meh - just at start
#ENSMUST00000067468 - just at start
#ENSMUST00000019069 - just spikes at start not real?
############

genes <- genes(txdb,'TXNAME')

genes <- genes(txdb,'AME')

fread%>%args

gtfdat <- fread(annofile,nrows=20)

metadata<-gtfdat[[9]]%>%str_extract_all('.*?;')%>%enframe%>%unnest

metadata%>%.$value%>%str_split_fixed('\\"',n=2)%>%as_data_frame%>%str_


#now let's write a function that compresses down the genome view of a track so we can see it on the 


trtracks <- alltracks[c(3,5)] 
trtracks[[1]]

plotTracks(chromosome=seqnames(plotwindow)[1],from=start(plotwindow)[1] - 1e3,to=end(plotwindow)[1]+1e3,strand=strand(plotwindow)[1],
	c(),
	type='coverage'
)

function (file, selection)
{
    indNames <- c(sub("\\.bam$", ".bai", file), paste(file, "bai",
        sep = "."))
    index <- NULL
    for (i in indNames) {
        if (file.exists(i)) {
            index <- i
            break
        }
    }
    if (is.null(index))
        stop("Unable to find index for BAM file '", file, "'. You can build an index using the following command:\n\t",
            "library(Rsamtools)\n\tindexBam(\"", file, "\")")
    pairedEnd <- parent.env(environment())[["._isPaired"]]
    if (is.null(pairedEnd))
        pairedEnd <- TRUE
    bf <- BamFile(file, index = index, asMates = pairedEnd)
    param <- ScanBamParam(which = selection, what = scanBamWhat(),
        tag = "MD", flag = scanBamFlag(isUnmappedQuery = FALSE))
    reads <- if (as.character(seqnames(selection)[1]) %in% names(scanBamHeader(bf)$targets))
        scanBam(bf, param = param)[[1]]
    else list()
    md <- if (is.null(reads$tag$MD))
        rep(as.character(NA), length(reads$pos))
    else reads$tag$MD
    if (length(reads$pos)) {
        layed_seq <- sequenceLayer(reads$seq, reads$cigar)
        region <- unlist(bamWhich(param), use.names = FALSE)
        ans <- stackStrings(layed_seq, start(region), end(region),
            shift = reads$pos - 1L, Lpadding.letter = "+", Rpadding.letter = "+")
        names(ans) <- seq_along(reads$qname)
    }
    else {
        ans <- DNAStringSet()
    }
    return(GRanges(seqnames = if (is.null(reads$rname)) character() else reads$rname,
        strand = if (is.null(reads$strand)) character() else reads$strand,
        ranges = IRanges(start = reads$pos, width = reads$qwidth),
        id = reads$qname, cigar = reads$cigar, mapq = reads$mapq,
        flag = reads$flag, md = md, seq = ans, isize = reads$isize,
        groupid = if (pairedEnd) reads$groupid else seq_along(reads$pos),
        status = if (pairedEnd) reads$mate_status else rep(factor("unmated",
            levels = c("mated", "ambiguous", "unmated")), length(reads$pos))))
}


myimportfun <- function (file, selection) {
    indNames <- c(sub("\\.bam$", ".bai", file), paste(file, "bai",
        sep = "."))
    index <- NULL
    for (i in indNames) {
        if (file.exists(i)) {
            index <- i
            break
        }
    }
    if (is.null(index))
        stop("Unable to find index for BAM file '", file, "'. You can build an index using the following command:\n\t",
            "library(Rsamtools)\n\tindexBam(\"", file, "\")")
    pairedEnd <- parent.env(environment())[["._isPaired"]]
    if (is.null(pairedEnd))
        pairedEnd <- TRUE
    bf <- BamFile(file, index = index, asMates = pairedEnd)
    param <- ScanBamParam(which = selection, what = scanBamWhat(),
        tag = "MD", flag = scanBamFlag(isUnmappedQuery = FALSE))
    reads <- if (as.character(seqnames(selection)[1]) %in% names(scanBamHeader(bf)$targets))
        scanBam(bf, param = param)[[1]]
    else list()
    md <- if (is.null(reads$tag$MD))
        rep(as.character(NA), length(reads$pos))
    else reads$tag$MD
    if (length(reads$pos)) {
        layed_seq <- sequenceLayer(reads$seq, reads$cigar)
        region <- unlist(bamWhich(param), use.names = FALSE)
        ans <- stackStrings(layed_seq, start(region), end(region),
            shift = reads$pos - 1L, Lpadding.letter = "+", Rpadding.letter = "+")
        names(ans) <- seq_along(reads$qname)
    }
    else {
        ans <- DNAStringSet()
    }
    return(GRanges(seqnames = if (is.null(reads$rname)) character() else reads$rname,
        strand = if (is.null(reads$strand)) character() else reads$strand,
        ranges = IRanges(start = reads$pos, width = reads$qwidth),
        id = reads$qname, cigar = reads$cigar, mapq = reads$mapq,
        flag = reads$flag, md = md, seq = ans, isize = reads$isize,
        groupid = if (pairedEnd) reads$groupid else seq_along(reads$pos),
        status = if (pairedEnd) reads$mate_status else rep(factor("unmated",
            levels = c("mated", "ambiguous", "unmated")), length(reads$pos))))
}


readgr <- myimportfun(polybams[1],plotwindow)


splitcigars <- readgr$cigar%>%str_extract_all('\\d+(M|N)')

