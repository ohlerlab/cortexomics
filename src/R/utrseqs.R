library(Rsamtools)
library(GenomicFeatures)
library(rtracklayer)
gtf <- 'pipeline/my_gencode.vM12.annotation.gtf'
fafile <- 'pipeline/my_GRCm38.p5.genome.chr_scaff.fa'
tputrfile <- 'tputrs.fa'
fputrfile <- 'fputrs.fa'



stopifnot(file.exists(fafile))
stopifnot(file.exists(gtf))
fafileob = Rsamtools::FaFile(fafile)
Rsamtools::IndexFa(fafile)
#get utrs
#load the gtf as a granges
gtf_gr<-rtracklayer::import(con=gtf,format='gtf')
#makea. transcript db object
gtftxdb <- makeTxDbFromGRanges(gtf_gr)
#get utrs from this as grangeslists
tputrs <- threeUTRsByTranscript(gtftxdb, use.names = TRUE)
fputrs <- fiveUTRsByTranscript(gtftxdb, use.names = TRUE)

#extract sequence and print to a file
 GenomicFeatures::extractTranscriptSeqs(tputrs,x=fafileob)%>%
 	writeXStringSet(tputrfile)
 GenomicFeatures::extractTranscriptSeqs(fputrs,x=fafileob)%>%
 	writeXStringSet(fputrfile)
