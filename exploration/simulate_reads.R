

library(rtracklayer)
library(GenomicFiles)
library(Rsamtools)
library(Biostrings)

args <- c('my_gencode.vM12.annotation.gff3','my_GRCm38.p5.genome.chr_scaff.fa','tRNAs.gtf')
args <- commandArgs(trailingOnly = TRUE) 

annot <- args[1] %>% import

genomefile <- args[2]

trnas <- args[3] %>% import %>% unique

rrnas <- annot%>%subset(gene_type=='rRNA')%>%subset(type=='exon')

ebp1<-annot%>%subset(type=='transcript')%>%subset(transcript_name=='Pa2g4-001')
ebp1short<-ebp1%>%resize(28+9)
randomcoding<-annot%>%subset(type=='CDS')%>%subset(width>28)%>%sample(1000)

#takes a length, a region, a genome, optional scoring later?

regions <- c(GRanges(c('chr1:10-40:+')),trnas[1:10,NULL])
regions <- rrnas

strandshift<-function(gr,shift) {shift = ifelse(strand(gr)=='-',-shift,shift);shift(gr,shift)}
simulate_reads <- function(regions, readlength, genome, fqfile){

	dir.create(dirname(fqfile),showWarn=F,rec=T)

	n_pos = width(regions) - readlength +1
	stopifnot(all(n_pos>0))
	reads <- rep(regions,n_pos)%>%resize(readlength,'start')
	shifts <- lapply(n_pos-1,seq,from=0)%>%unlist
	reads <- strandshift(reads[,NULL],shifts)
	names(reads) <- paste0(names(reads),'_',as.character(reads))
	#
	readfile <- tempfile(fileext='.bed')
	export(reads,readfile)
	message(str_interp('generated ${length(reads)} reads'))
	#run perl script to go from fasta to reads with perfect quality 
	cmd <- str_interp('bedtools getfasta -s -fi ${genomefile} -bed ${readfile}  -name | perl -lanpe \'s/>/@/ ; s/(^[^@]+)/\\1\n+/; if($1){$a="I" x length($1); s/\\+/+\n$a/}\' | gzip > ${fqfile}')
	system(cmd)
}

simulate_reads(ebp1,readlength=28,genome,'input/test_ebp1/test_ebp1.fastq.gz')

namegr<-function(x,name) x%>%setNames(rep(name,length(.)))
c(rrnas[,NULL]%>%namegr('rRNA'),ebp1short[,NULL]%>%namegr('ebp1'))%>%
	simulate_reads(readlength=28,genome,'input/test_rRNA/test_rRNA.fastq.gz')

c(trnas[,NULL]%>%namegr('tRNA'),ebp1short[,NULL]%>%namegr('ebp1'))%>%
	simulate_reads(readlength=28,genome,'input/test_tRNA/test_tRNA.fastq.gz')

simulate_reads(randomcoding,readlength=16,genome,'input/test_16bpRNA/test_16bpRNA.fastq.gz')


#call it on our regions
###Now test this
testumapread <- GRanges('chr1:171101906-171101933:-')

testumapread%>%mergeByOverlaps(trnas)