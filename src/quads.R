
message('loading libraries')
library(pqsfinder)
library(rtracklayer)
library(magrittr)
library(data.table)
library(rtracklayer)
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(assertthat))
suppressMessages(library(tidyverse))
library(GenomicRanges)

message('...done')

#load arguments
getwd()
args <- c(
	seqfile = 'motseqs/five_prime_UTR/xtailtranslregged/xtailtranslregged.fa',
	seqfileback = 'motseqs/five_prime_UTR/nontransregged.fa',
  outputfolder = 'gquadsearch'
)
args <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
for(i in names(args)) assign(i,args[i])
dir.create(outputfolder,showWarn=T)

seqsetname <- seqfile%>%dirname%>%basename
seqbackgroundname <- paste0(seqfileback%>%dirname%>%basename,'_background')
seq <- rtracklayer::import(seqfile,format='fasta')
backseq <- rtracklayer::import(seqfileback,format='fasta')

gqadseqob <- mclapply(mc.cores=12,F=pqsfinder,seq,strand='+')
backgqadseqob <- mclapply(mc.cores=12,F=pqsfinder,backseq,strand='+')


gquadscores <- bind_rows(
	gqadseqob%>%
	# sample(100)%>%
		map_df(~data_frame(score=c(0,.@elementMetadata$score),nm=c(0,.@elementMetadata$nm)),.id='gene')%>%
		mutate(set=seqsetname),
	backgqadseqob%>%
	# sample(100)%>%
		map_df(~data_frame(score=c(0,.@elementMetadata$score),nm=c(0,.@elementMetadata$nm)),.id='gene')%>%
		mutate(set=seqbackgroundname)
)

gquadscores%>%group_by(gene)%>%dplyr::slice(which.max(score))%>%{split(.$score,.$set)}%>%{wilcox.test(.[[1]],.[[2]])}
gquadscores%>%filter(nm<2)%>%group_by(gene)%>%dplyr::slice(which.max(score))%>%{split(.$score,.$set)}%>%{wilcox.test(.[[1]],.[[2]])}

gquadscores%>%{split(.$score,.$set)}%>%{wilcox.test(.[[1]],.[[2]])}

gquadscores%>%group_by(set,gene)%>%summarize(score=any(score>0))%>%{(table(.$score,.$set))}
gquadscores%>%group_by(set,gene)%>%summarize(score=any(score>0))%>%{fisher.test(table(.$score,.$set))}
gquadscores%>%group_by(set,gene)%>%summarize(score=any(score>35))%>%{fisher.test(table(.$score,.$set))}
gquadscores%>%group_by(set,gene)%>%summarize(score=any(score>50))%>%{fisher.test(table(.$score,.$set))}

gquadscores%>%filter(nm<2)%>%{split(.$score,.$set)}%>%{wilcox.test(.[[1]],.[[2]])}

gquadscores%>%filter(nm<2)%>%group_by(set,gene)%>%summarize(score=any(score>50))%>%{fisher.test(table(.$score,.$set))}


gquadscores%>%group_by(set,gene)%>%summarize(score=any(score>0))%>%{fisher.test(table(.$score,.$set))}


scoredistplotfile = normalizePath(file.path(outputfolder,'gquadscoredists.pdf')
pdf(scoredistplotfile); 
print(ggplot(gquadscores%>%group_by(gene)%>%dplyr::slice(which.max(score)),aes(x= score))+
	geom_histogram(binwidth=1)+
	coord_cartesian(xlim=c(0,150))+
	scale_x_continuous('Pqsfinder Gquad scores')+
	facet_grid(set ~ .,scale='free')
)
print(ggplot(gquadscores%>%group_by(gene)%>%dplyr::slice(which.max(score)),aes(x= score,color=set))+
	geom_density()+
	coord_cartesian(xlim=c(0,150))+
	scale_x_continuous('Pqsfinder Gquad scores')

)
print(gquadscores%>%group_by(set,gene)%>%summarize(score=any(score>50))%>%ggplot(aes(y=as.numeric(score),x=set))+stat_summary(fun.data=mean_cl_boot,geom='errorbar')+stat_summary(fun.data=mean_cl_boot,geom='bar',alpha=I(0.5))+
	scale_y_continuous('Fraction with strong (>50 pqsfinder) Quadruplexes')
	)
dev.off()
message(scoredistplotfile)
