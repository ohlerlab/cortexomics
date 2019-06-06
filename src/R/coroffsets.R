
################################################################################
########My own correlation based procedure
################################################################################


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
# testbam <- '../../Ribo_Lausanne/pipeline/star/data/OD5P_ctrl_1/OD5P_ctrl_1.bam'

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





# vlen<-length(allexprvect)

# offsetcors <- lapply(-20:20,function(i){
# 	cor(
# 		log(allexprvect+1)%>%.[100:(vlen-100)],
# 		log(subexprvect+1)%>%.[(100+i):((vlen-100)+i)]
# )
# })

# offsetcors%>%unlist%>%txtplot
# offsetcors%>%setNames(-20:20)%>%unlist%>%sort%>%round(2)
# sum(maincovvect)%>%{names(.[.>10])}
#now create coverage vectors for the biggest read lengths.

#now, create a big 
splitinds <- sample(1:1000)%>%split(.,ceiling(./(1000/40)))

offsets = -6:6


i= -6
splitind=splitinds[[1]]

offsetcors <- (offsets) %>% setNames(.,.)%>% lapply(function(i){
	lapply(splitinds,function(splitind){
		allexprvect <- maincovvect[trstouse[splitind]]%>%unlist
		subexprvect <- subcovvect[trstouse[splitind]]%>%unlist
		vlen = length(allexprvect)
		cor(
			log(allexprvect+1)%>%.[100:(vlen-100)],
			log(subexprvect+1)%>%.[(100+i):((vlen-100)+i)]
		)
	})
})

offsetcordf <-offsetcors%>%map(unlist)%>%setNames(offsets)%>%enframe%>%unnest

offsetcordf%>%group_by(name)%>%summarise(mean=mean(value),lc=quantile(value,0.25),hc=quantile(value,0.75))%>%arrange(desc(mean))

# %>%{txtplot(as.numeric(.$name),.$value)}


#Figure out the most expressed transcripts with Riboseq
#Get say the 1000 most expressed transcxripts
#Now from these, produce vectors of read length, sequence surrounding cut sites, 
#Now see if we can maximize the efficiency of the procedure to align these tracks
mostexpressed <- get_most_expressed()

reads_on_expressed <- get_reads(mostexpressed)

reads_on_expressed <- annotate_reads_flankseq(reads_on_expressed)

# read_vect <- 






