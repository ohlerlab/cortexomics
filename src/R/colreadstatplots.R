library(parallel)

readlengths <- 
	Sys.glob('collapse_reads/*/*collreadstats.txt')%>%
	grep(v=T,inv=T,patt='test')%>%
	# str_subset('Poly|80S')%>%
	map(.%>%{
		message(.)
		lines <- readLines(.)
		tabstart = lines%>%str_detect('>read_lengths_unique')%>%which
		tab <- fread(.,skip=tabstart+2)
		tab$sample <- dirname(.)%>%basename
		tab%<>%set_colnames(c('readlength','count','sample'))
		tab
	})%>%bind_rows

readlengths$readlength = readlengths$readlength - 8


getfastqlengths<-function(fastqfile,catbin='zcat'){
	if(str_detect(fastqfile,'collapse')) catbin <- 'cat'
	out<-fread(
		str_interp("${catbin} ${fastqfile} | awk 'NR%4==2{a[length($0)] += 1} END { for (key in a) { print key \"\\t\" a[key] } }' ")
		)
	message('.')
	out
}



getreadlengthcounts <- function(gtf,bamfile){
	fread(str_interp("tail -n +4 ${gtf} | awk '{print $1,$4,$5,$6,$8,\"-\"}'  | samtools view -L - ${bamfile} | awk '{a[length($10)] += 1} END { for (key in a) { print key \"\\t\" a[key] } }' "))
}

fastqfiles<-
	c(
		# Sys.glob("collapse_reads/*/*.fastq.gz"),
		Sys.glob("input/*/*.fastq.gz"),
		Sys.glob("cutadapt_reads/*/*.fastq.gz"),
		Sys.glob("collapse_reads/*/*.fastq.gz"),
		Sys.glob("trim_reads/*/*.fastq.gz"),
		Sys.glob("processed_reads/*/*.fastq.gz")
	)


readlengths <- fastqfiles%>%
	setNames(.,.)%>%
	mclapply(safely(getfastqlengths),mc.cores=20)
readlengths <- readlengths%>%map('result')%>%	bind_rows(.id='fastq')%>%
	mutate(processing_stage = basename(dirname(dirname(fastq))),fastq = basename(dirname(fastq)))%>%
	set_colnames(c('sample','readlength','count','processing_stage'))

readlengths%>%write_tsv('readlengthstats.tsv')

frac80sbams<-Sys.glob("star/data/*/*.bam")%>%grep(v=T,inv=T,patt='transcript')

cdsgtf<-"my_gencode.vM12.annotation.cds.gtf"

rlfrac80sbams<- frac80sbams %>% 
	setNames(.,.)%>%
	mclapply(function(.) safely(getreadlengthcounts)(cdsgtf,.),mc.cores=20)
rlfrac80sbams<- rlfrac80sbams%>%map('result')	%>%bind_rows(.id='fastq')%>%
	mutate(processing_stage = 'aligned_to_cds',fastq = basename(dirname(fastq)))%>%
	set_colnames(c('sample','readlength','count','processing_stage'))

readlengths <- rbind(readlengths,rlfrac80sbams)

readlengths$processing_stage%<>%factor(.,unique(.))
readlengths%>%filter(processing_stage=='aligned_to_cds')
#
pdf('../plots/80S_readlendists_prcstages.pdf'%>%normalizePath%T>%message,w=12,h=12)
print(
readlengths%>%
	# filter(sample%>%str_detect('ribo'))%>%
	filter(sample%>%str_detect('80S'))%>%
	group_by(sample,processing_stage)%>%
	mutate(count_frac=count/sum(count))%>%
	ggplot(aes(x=readlength,y=count_frac,color=sample))+
		facet_grid(processing_stage~.,scale='free')+
		geom_line()+
		scale_x_continuous(breaks=seq_len(max(readlengths$readlength)))+
		theme_minimal()+
		theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank())
		)
dev.off()

pdf('../plots/Poly_readlendists_prcstages.pdf'%>%normalizePath%T>%message,w=12,h=12)
print(
readlengths%>%
	# filter(sample%>%str_detect('ribo'))%>%
	filter(sample%>%str_detect('Poly'))%>%
	group_by(sample,processing_stage)%>%
	mutate(count_frac=count/sum(count))%>%
	ggplot(aes(x=readlength,y=count_frac,color=sample))+
		facet_grid(processing_stage~.,scale='free')+
		geom_line()+
		scale_x_continuous(breaks=seq_len(max(readlengths$readlength)))+
		theme_minimal()+
		theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank())
		)
dev.off()

pdf('../plots/Ribo_readlendists_prcstages.pdf'%>%normalizePath%T>%message,w=12,h=12)
print(
readlengths%>%
	# filter(sample%>%str_detect('ribo'))%>%
	filter(sample%>%str_detect('ribo'))%>%
	group_by(sample,processing_stage)%>%
	mutate(count_frac=count/sum(count))%>%
	ggplot(aes(x=readlength,y=count_frac,color=sample))+
		facet_grid(processing_stage~.,scale='free')+
		geom_line()+
		scale_x_continuous(breaks=seq_len(max(readlengths$readlength)))+
		theme_minimal()+
		theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank())
		)
dev.off()

#cludge - get rid of data form some small test files
readlengths%<>%group_by(sample,readlength,processing_stage)%>%slice(which.max(count))

stop()

getbamlengths<-function(fastqfile){
	fread(
		str_interp("samtools view -F 4 star/data/RPI5_80SE16_1/RPI5_80SE16_1.bam | awk '{a[length($10)] += 1} END { for (key in a) { print key \"\\t\" a[key] } }' ")
		)
}

getbamlengths

bamTagFilter(object)




#####FOrdifferent sample groups

pdf(w=14,h=3,'../plots/readlendists_Poly_80S.pdf')
readlengths%>%
	# filter(sample%>%str_detect('ribo'))%>%
	filter(sample%>%str_detect('80S|Poly'))%>%
	group_by(sample)%>%
	mutate(count_frac=count/sum(count))%>%
	ggplot(aes(x=readlength,y=count_frac,color=sample))+
		geom_line(aes(alpha=sample=='RP2_80SE13_2'))+
		scale_x_continuous(breaks=seq_len(max(readlengths$readlength)))+
		theme_minimal()+
		theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank())
dev.off()
normalizePath('../plots/readlendists_Poly_80S.pdf')


pdf(w=14,h=3,'../plots/readlendists_totalribo.pdf')
readlengths%>%
	filter(sample%>%str_detect('ribo'))%>%
	# filter(!sample%>%str_detect('80S|Poly'))%>%
	group_by(sample)%>%
	mutate(count_frac=count/sum(count))%>%
	ggplot(aes(x=readlength,y=count_frac,color=sample))+
		geom_line()+
		scale_x_continuous(breaks=seq_len(max(readlengths$readlength)))+
		theme_minimal()+
		theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank())
dev.off()
normalizePath('../plots/readlendists_totalribo.pdf')


pdf(w=14,h=3,'../plots/readlendists_all.pdf')
readlengths%>%
	filter(sample%>%str_detect('ribo|80S|Poly'))%>%
	# filter(!sample%>%str_detect('80S|Poly'))%>%
	group_by(sample)%>%
	mutate(count_frac=count/sum(count))%>%
	ggplot(aes(x=readlength,y=count_frac,color=sample))+
		geom_line()+
		scale_x_continuous(breaks=seq_len(max(readlengths$readlength)))+
		theme_minimal()+
		theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank())
dev.off()
normalizePath('../plots/readlendists_all.pdf')


readlengths%>%group_by(sample)%>%summarize(max_rl = readlength[which.max(count)])%>%filter(!str_detect(sample,'total'))%>%arrange(desc(max_rl))


pdf(w=14,h=3,'../plots/readlendists_80SE13.pdf')
readlengths%>%
	filter(sample%>%str_detect('80SE13'))%>%
	# filter(!sample%>%str_detect('80S|Poly'))%>%
	group_by(sample)%>%
	mutate(count_frac=count/sum(count))%>%
	ggplot(aes(x=readlength,y=count_frac,color=sample))+
		geom_line()+
		scale_x_continuous(breaks=seq_len(max(readlengths$readlength)))+
		theme_minimal()+
		theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank())
dev.off()
normalizePath('../plots/readlendists_80SE13.pdf')


# readlengths%>%group_by(sample)%>%summarize(max_rl = readlength[which.max(count)])%>%filter(!str_detect(sample,'total'))%>%arrange(desc(max_rl))i