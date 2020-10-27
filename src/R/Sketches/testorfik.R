library(ORFik)
library(GenomicFiles)
library(GenomicFeatures)
library(GenomicAlignments)
filter<-dplyr::filter
select<-dplyr::select
slice<-dplyr::slice
filter<-dplyr::filter
txdb <- makeTxDbFromGFF('pipeline/my_gencode.vM12.annotation.gtf')
footprints <- BamFile

footprints <-'pipeline/star/data/E13_ribo_1/E13_ribo_1.bam'%>%
	BamFile(file=.,yieldSize=NA)%>%
	readGAlignments

footprints%>%qwidth%>%table

footprints %<>% .[qwidth(.)<32]

Sys.time()
shifts <- detectRibosomeShifts(footprints, txdb, stop = TRUE)
Sys.time()

filtshifts <- shifts%>%filter(fragment_length<30)


p_sites <- ORFik::shiftFootprints(footprints,shifts$fragment_length,shifts$offsets_start)

txNames <- txNamesWithLeaders(txdb)
windows <- getStartStopWindows(txdb,txNames)

hitMapStart <- metaWindow(footprints, windows$start)
hitMapStop <- metaWindow(footprints, windows$stop)


hitMapStart <- metaWindow(p_sites, windows$start)
hitMapStop <- metaWindow(p_sites, windows$stop)

p<-ggplot(hitMapStart, aes(x = factor(position), y = avg_counts, fill = factor(frame))) +
	geom_bar(stat = "identity") +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
	labs(title = paste0("Length 29 over START of canonical CDS")) +
	xlab("\nshift from first START nucleotide [bp]") +
	ylab("Averaged counts") +
	guides(fill = FALSE)
)

 p2<- ggplot(hitMapStop, aes(x = factor(position), y = avg_counts, fill = factor(frame))) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = paste0("Length 29 over STOP of canonical CDS")) +
    xlab("\nshift from last STOP nucleotide [bp]") +
    ylab("Averaged counts") +
    guides(fill = FALSE)
pdf('tmp.pdf');print(p);print(p2);dev.off()


