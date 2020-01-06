#this script is going to look at abias in cut sites. It needs a bunch of kmer counts, assumes they are equally spaced around the cut sites

library(data.table)
library(Rsamtools)
library(rtracklayer)
library(tidyverse)
library(stringr)
library(magrittr)
Sys.glob('cutsequences/*/*.txt.gz')

# kmerfiles<-Sys.glob('pipeline/cutsequences/*/cutseqs.txt.gz*')
# for(kmerfile in kmerfiles){
# 	# try({
# 	sample <- kmerfile %>%dirname%>%basename
# 	cuttype = kmerfile %>%str_extract('right')%>%replace(.,is.na(.),'left')

# 	kmercounts<-kmerfile%>%read_table(col_names=F)

# 	nbases <- kmercounts[[2]][[1]]%>%str_length
# 	basemat<-kmercounts[[2]]%>%str_split_fixed('',nbases)
# 	bases<-unique(basemat[,1])

# 	basefreqs<-sapply(1:ncol(basemat),function(colnum){
# 		sapply(bases,function(base){
# 			sum(kmercounts[[1]][basemat[,colnum] == base])
# 		})
# 	}) 

# 	basefreqdf<-basefreqs%>%as.data.frame%>%rownames_to_column('base')%>%gather(position,count,-base)%>%mutate(position=str_replace(position,'V','')%>%as.numeric%>%{.-(max(.)/2)})
# 	basefreqdf%<>%group_by(position)%>%mutate(freq = count/sum(count))
# 	nreads = sum(kmercounts[[1]])
# 	freqplot<-basefreqdf%>%ggplot(aes(x=position,y=freq,color=base,group=base))+geom_line()+theme_bw()+scale_y_continuous(limits=c(0,1))+ggtitle(paste0('Base Frequencies around\n',nreads,' ',cuttype,' cut sites \n',sample))


# 	dir.create('./plots/')
# 	ggsave(freqplot,file=paste0('./plots/',sample,'_',cuttype,'_','basefreqs.pdf')%>%normalizePath%T>%message)
# 	# })
# }

#  file = 'pipeline/star/data/test/test.bam'
# library(Rsamtools)
# getfaseqs <- function (file, selection){
#   if (!file.exists(paste(file, "bai", sep = "."))){ 
#     stop("Unable to find index for BAM file '", file, "'. You can   build an index using the following command:\n\t", 
#          "library(Rsamtools)\n\tindexBam(\"", file, "\")")
#   }
# sinfo <- scanBamHeader(file)[[1]] 
# res <- if (!as.character(seqnames(selection)[1]) %in% names(sinfo$targets)) { 
#   mcols(selection) <- DataFrame(score = 0) 
#   selection 
# }else { 

#   param <- ScanBamParam(what = c("pos", "qwidth", "strand"), 
#                         which = selection, flag = 
#                           scanBamFlag(isUnmappedQuery = FALSE)) 
#   x <- scanBam(file, param = param)[[1]] 

#     rwidth = x[['qwidth']]
  
#     gr <- GRanges(strand=x[["strand"]], ranges=IRanges(x[["pos"]], 
#                                                      width = rwidth), seqnames=seqnames(selection)[1]) 
#   grs <- split(gr, strand(gr)) 
#   cov <- lapply(grs[c("+", "-")], function(y) coverage(ranges(y), 
#                                                        width=end(selection))) 

#   if(length(pos)==0){ 
#     mcols(selection) <- DataFrame(plus=0, minus=0) 
#     selection 
#   }else{ 
#     cov = cov%>%map(.%>%list%>%as('SimpleRleList')%>%setNames(seqnames(selection)[1])%>%as("GRanges")%>%subset(score!=0))
#     colnames(mcols(cov[['+']])) = 'plus'
#     colnames(mcols(cov[['-']])) = 'minus'
#     cov[['+']]$minus = 0
#     cov[['-']]$plus = 0
#     cov[['-']]$minus = -1 * cov[['-']]$minus
#     mcols(cov[['-']]) <- mcols(cov[['-']])[c('plus','minus')]   
#     c(cov[['+']],cov[['-']])
#     # GRanges(seqnames = seqnames(selection)[1], 
#     #         ranges=IRanges(start=head(pos, -1), end=tail(pos, -1)), 
#     #         plus=as.numeric(cov[["+"]][head(pos, -1)]), 
#     #         minus=-as.numeric(cov[["-"]][head(pos, -1)])) 
#   } 
# } 
# return(res) 
# }




strandshift <-function(gr,k){GenomicRanges::shift(gr,ifelse(as.logical(strand(gr)=='-'),-k,k))}

bamfile = 'pipeline/star/data//E13_ribo_1/E13_ribo_1.bam'
fastafile = 'pipeline/my_GRCm38.p5.genome.chr_scaff.fa'
cds4filt<-import('pipeline/my_gencode.vM12.annotation.cdsfilt.gtf')

param <- ScanBamParam(what = c("rname","pos", "qwidth", "strand"),which=cds4filt%>%head(100)%>%{strand(.)<-'*';.} ,flag = 
                      scanBamFlag(isUnmappedQuery = FALSE))
bam <- GenomicAlignments::readGappedReads(bamfile, param = param)
bam %<>% as('GRanges')
bam <- bam[between(bam$qwidth,25,29)]
bam <- bam[strand(bam)=='+']
mcols(bam)<-NULL

resizebam <- bam%>% resize(1,'start')
resizebam %<>% strandshift(0)
seqs<-scanFa(file = fastafile, resizebam)
as.character(seqs)%>%table%>%{./sum(.)}

	str_split_fixed('',6)%>%
	apply(2,table) %>%
	map(~.[c('A','C','G','T')])%>%
	simplify2array%>%
	apply(2,function(x) x/sum(x,na.rm=TRUE) )

