

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

randomcoding<-annot%>%subset(type=='transcript')%>%subset(transcript_name=='Pa2g4-001')

#takes a length, a region, a genome, optional scoring later?

regions <- c(GRanges(c('chr1:10-40:+')),trnas[1:10,NULL])
regions <- rrnas

simulate_reads <- function(regions, readlength, genome, fqfile){

	n_pos = width(regions) - readlength +1
	reads <- rep(regions,n_pos)%>%resize(readlength,'start')
	shifts <- lapply(n_pos-1,seq,from=0)%>%unlist
	reads <- shift(reads[,NULL],shifts)
	names(reads) <- as.character(reads)
	#
	readfile <- tempfile(fileext='.bed')
	export(reads,readfile)
	message(str_interp('generated ${length(reads)} reads'))
	#run perl script to go from fasta to reads with perfect quality 
	cmd <- str_interp('bedtools getfasta -fi ${genomefile} -bed ${readfile}  -name | perl -lanpe \'s/>/@/ ; s/(^[^@]+)/\\1\n+/; if($1){$a="I" x length($1); s/\\+/+\n$a/}\' | gzip > ${fqfile}')
	system(cmd)
}
'input/test_ebp1/test_ebp1.fastq.gz'%>%dirname%>%dir.create

simulate_reads(randomcoding,readlength=28,genome,'input/test_ebp1/test_ebp1.fastq.gz')

simulate_reads(c(rrnas,randomcoding),readlength=28,genome,'input/test_rRNA/test_rRNA.fastq.gz')

simulate_reads(c(trnas[,NULL],randomcoding[,NULL]),readlength=28,genome,'input/test_tRNA/test_tRNA.fastq.gz')



#call it on our regions