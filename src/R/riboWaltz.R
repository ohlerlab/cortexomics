suppressMessages({library(svglite)})
suppressMessages({library(readr)})
suppressMessages({library(Biostrings)})
suppressMessages({library(Rsamtools)})
suppressMessages({library(rtracklayer)})
suppressMessages({library(GenomicFeatures)})
suppressMessages({library(stringr)})
suppressMessages({library(data.table)})
suppressMessages({library(assertthat)})
suppressMessages({library(parallel)})
suppressMessages({library(Biostrings)})
suppressMessages({library(dplyr)})
suppressMessages({library(riboWaltz)})
suppressMessages({library(purrr)})

for(fname in lsf.str('package:dplyr')) assign(fname,get(fname,'package:dplyr'))

source('../src/R/Rprofile.R')

argv <- c(
	transcriptbam = 'riboWaltz/RPI8_PolyE16_2/RPI8_PolyE16_2.star_transcript.bam',
	gtf = 'my_gencode.vM12.annotation.gtf',
	outfolder = 'riboWaltz/RPI8_PolyE16_2_RUST/'
)


argv[] <- commandArgs(trailing=TRUE)

for (nm in names(argv)) assign(nm,argv[[nm]])

RUST <- str_detect(outfolder,'RUST')


message(getwd())
message(argv)
# save.image()
# stop()
# on.exit(save.image())

i=1
message(i);i=i+1
sampnames <- basename(dirname(transcriptbam))

riboWaltzanno <- create_annotation(gtf)

newtranscriptbam <- file.path(outfolder,basename(transcriptbam))

file.link(transcriptbam,newtranscriptbam)
message(i);i=i+1

reads_list <- bamtolist(bamfolder = dirname(transcriptbam), list_name = sampnames,annotation = riboWaltzanno, transcript_align=TRUE ,rm_version=T)
0

if(RUST) {
	nureads_list<- reads_list-
	reads_list<-reads_list %>% map(~ distinct(.,transcript,end5,end3,.keep_all=TRUE)%>%setindex('stop_pos'))
	setindex(ureads_list[[1]],'stop_pos')
}
message(i);i=i+1

0
example_length_dist_zoom <- rlength_distr(reads_list, sample = sampnames, cl = 99)

# filtered_list <- length_filter(data = reads_list, length_filter_mode = "custom",
#  				length_filter_vector = 14:30)

# pfiltered_list <- length_filter(data = reads_list, length_filter_mode = "periodicity",
#  				periodicity_threshold = 50)






message(i);i=i+1
psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")


message(i);i=i+1

reads_psite_list <- psite_info(reads_list, psite_offset)

0
message(i);i=i+1
example_psite_region <- region_psite(reads_psite_list, riboWaltzanno, sample = sampnames)


message(i);i=i+1
psite_cds <- psite_per_cds(reads_psite_list, riboWaltzanno)

message(i);i=i+1
example_frames_stratified <- frame_psite_length(reads_psite_list, sample = sampnames[[1]],
                                                region = "all", cl = 90)

message(i);i=i+1
example_frames <- frame_psite(reads_psite_list, sample = sampnames[[1]], region = "all")

message(i);i=i+1



example_metaprofile <- metaprofile_psite(reads_psite_list, mm81cdna, sample = sampnames[[1]],
                                         utr5l = 20, cdsl = 40, utr3l = 20,
                                         plot_title = "auto")


message(i);i=i+1
readlens <- reads_psite_list[[1]]$length%>%table%>%keep(~.>1e3)%>%names%>%sort

example_metaprofile_i<- rep(NA,length(readlens))%>%setNames(readlens)
example_metaheatmap_compi <- rep(NA,length(readlens))%>%setNames(readlens)








example_metaheatmap_compi%>%names

message(i);i=i+1
mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.5)),
    colhead = list(fg_params=list(cex = 0.5)),
    rowhead = list(fg_params=list(cex = 0.5)))


ribowaltzpdf <- paste0(outfolder,'/',sampnames[[1]],'_ribowaltplots.pdf')

0
getwd()
ribowaltzpdf%>%dirname%>%dir.create

message(i);i=i+1
pdf(ribowaltzpdf,w=14,h=7)
example_length_dist_zoom[["plot"]]
grid::grid.newpage()
# example_ends_heatmap[["plot"]]
gridExtra::grid.table(psite_offset, theme = mytheme,rows=NULL)
example_psite_region[["plot"]]
example_frames_stratified[["plot"]]
example_metaprofile[["plot"]]
example_frames[["plot"]]

metaprofile_psite%>%debug
metaprofile_psite%>%undebug

for(readlen in readlens%>%keep(~.%in%c(16,21,27,40))){
	try({
	message(paste0('comparison plots for readlength:',readlen))

	reads_psite_readlen <- reads_psite_list[[sampnames[[1]]]][length == as.numeric(readlen)]

	example_metaprofile_i <- metaprofile_psite(setNames(list(reads_psite_readlen),sampnames[[1]]), riboWaltzanno, sample = sampnames[[1]],
	                                            length_range = readlen, utr5l = 20, cdsl = 60, 
	                                            transcripts = reads_psite_readlen$transcript%>%unique,
	                                            utr3l = 20, plot_title = "auto")
	print(example_metaprofile_i[['plot']])

	comparison_dt <- list()
	
	
	comparison_dt[[paste0("subsample_",readlen,"nt")]] <- reads_psite_readlen
	comparison_dt[["whole_sample"]] <- reads_psite_list[[sampnames[[1]]]]

	names_list <- list( paste0("subsample_",readlen,"nt"),"whole_sample" )%>%setNames(c(paste0("Only_",readlen),'All'))

	scale_facts <- comparison_dt%>%map_dbl(~ 1e6 / nrow(.))%>%setNames(names(comparison_dt))
	example_metaheatmap_compi <- metaheatmap_psite(comparison_dt, riboWaltzanno, sample = names_list,
	                                         utr5l = 20, cdsl = 40, utr3l = 20, log = F, scale_factors=scale_facts)
	print(example_metaheatmap_compi[['plot']])

	example_frames <- frame_psite(reads_psite_list%>%map(~.[length==as.numeric(readlen)]), sample = sampnames[[1]], region = "all")
	print(example_frames[["plot"]])
	})
}

dev.off()

normalizePath(ribowaltzpdf)

suppressMessages({library(svglite)})
suppressMessages({library(readr)})
suppressMessages({library(Biostrings)})
suppressMessages({library(Rsamtools)})
suppressMessages({library(rtracklayer)})
suppressMessages({library(GenomicFeatures)})
suppressMessages({library(GenomicAlignments)})
suppressMessages({library(stringr)})
suppressMessages({library(data.table)})
suppressMessages({library(assertthat)})
suppressMessages({library(parallel)})
suppressMessages({library(Biostrings)})
suppressMessages({library(dplyr)})
suppressMessages({library(riboWaltz)})
suppressMessages({library(purrr)})

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


testbam<-'star/data/RPI8_PolyE16_2/RPI8_PolyE16_2.bam'
testbam<-'star/data/E13_total_1/E13_total_1.bam'
testbam <- '../../Ribo_Lausanne/pipeline/star/data/OD5P_ctrl_1/OD5P_ctrl_1.bam'

subexons <- exons%>%subset(strand=='+')%>%subset(seqnames%in%c('chr1','chr2'))
testreadsnu <- readGAlignments(testbam,param=ScanBamParam(which=subexons))
testreads <- testreadsnu%>%as("GRanges")%>%unique
genome<-FaFile('my_GRCm38.p5.genome.chr_scaff.fa')

readstartseqs <- 	testreads%>%subset(strand=='+')%>%as("GenomicRanges")%>%resize(2,'start')%>%resize(4,'end')%>%getSeq(x=genome)
readendseqs <- 		testreads%>%subset(strand=='+')%>%as("GenomicRanges")%>%resize(2,'end')%>%resize(4,'start')%>%getSeq(x=genome)

startseqdf<-list(start=readstartseqs,end=readendseqs)%>%lapply(.%>%as.matrix%>%apply(2,table)%>%sweep(.,2,colSums(.),FUN='/')%>%round(2)%>%set_colnames(1:4)%>%as.data.frame%>%rownames_to_column('base')%>%
	gather(position,frequency,-base))%>%bind_rows(.id='start_or_end')%>%
	mutate(start_or_end=factor(start_or_end,unique(start_or_end)))

startseqdf$position%<>%as.numeric
cutlabels = c('cut -2','cut -1','cut +1','cut +2')

cutplottitle<-str_interp('Read start/end base frequencies\nBamfile ${testbam}')

startseqdf%>%
	ggplot(aes(color=base,x=as.numeric(position),y=frequency))+
	scale_x_continuous(labels=cutlabels)+
	geom_line()+facet_grid(start_or_end ~.)+
	ggtitle(cutplottitle)


acut<-readstartseqs%>%str_detect('.A..')%>%mean
atwopos<-readstartseqs%>%str_detect('A...')%>%mean

readstartseqs%>%str_detect('AA..')%>%mean
acut*atwopos


readstarts <- testreads %>% as("GenomicRanges")%>%resize(1,'start')
readstarts$readlength <- width(testreads)




readstartstr <- mapToTranscripts(readstarts,exons%>%split(.,.$transcript_id))
readstartstr$readlength <- readstarts$readlength[readstartstr$xHits]

#now put scores in 
readlens = 26:31
maincovvect <- readstartstr%>%subset(readlength==27)%>%coverage
subcovvect <- readstartstr%>%subset(readlength==29)%>%coverage

trstouse<- readstartstr@seqnames%>%table%>%sort%>%tail(1000)%>%names


vlen<-length(allexprvect)

offsetcors <- lapply(-20:20,function(i){
	cor(
		log(allexprvect+1)%>%.[100:(vlen-100)],
		log(subexprvect+1)%>%.[(100+i):((vlen-100)+i)]
)
})
offsetcors%>%unlist%>%txtplot
offsetcors%>%setNames(-20:20)%>%unlist%>%sort%>%round(2)
sum(maincovvect)%>%{names(.[.>10])}
#now create coverage vectors for the biggest read lengths.

#now, create a big 
splitinds <- sample(1:1000)%>%split(.,ceiling(./(1000/40)))

offsetcors <- (-10:2) %>% setNames(.,.)%>% lapply(function(i){
	lapply(splitinds,function(splitind){
		allexprvect <- maincovvect[trstouse[splitind]]%>%unlist
		subexprvect <- subcovvect[trstouse[splitind]]%>%unlist
		vlen = length(allexprvect)
		cor(
			log(allexprvect+1)%>%.[100:(vlen-100)],
			log(subexprvect+1)%>%.[i(100+i):((vlen-100)+i)]
		)
	})
})

offsetcordf <-offsetcors%>%map(unlist)%>%setNames(-10:10)%>%enframe%>%unnest

offsetcordf%>%group_by(name)%>%summarise(mean=mean(value),lc=quantile(value,0.25),hc=quantile(value,0.75))%>%arrange(desc(mean))

%>%{txtplot(as.numeric(.$name),.$value)}


