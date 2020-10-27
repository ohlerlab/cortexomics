library(svglite)
library(readr)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(stringr)
library(data.table)
library(assertthat)
library(parallel)

filter<-dplyr::filter
slice<-dplyr::slice
select<-dplyr::select
summarise<-dplyr::summarise

source('../src/R/Rprofile.R')

#initiaal form when testing
#let's look at the scores for our
satannfiles <- Sys.glob('SaTAnn/*/*Final_ORFs*')%T>%{stopifnot(length(.)>0)}
satannorfs <- 
	# Sys.glob('SaTAnn/*/*Final_ORFs*')%>%
	satannfiles%>%
	setNames(.,basename(dirname(.)))%>%
	mclapply(load_objs)

#now filter out the weird GRanges columns, for now, and aggregate into one data table 
all_orfs <- satannorfs%>%map(.%>%.$ORFs_tx)
all_orfs_gen <- satannorfs%>%map(.%>%.$ORFs_gen)%>%GRangesList%>%unique

for (i in seq_along(all_orfs))all_orfs[[i]]$Protein%<>%as.character

orfs_dt <- all_orfs%>%map(.%>%{
	issimplelist <- mcols(.)%>%vapply(is,TRUE,'atomic')
	.[,issimplelist]%>%GR2DT
})%>%bind_rows(.id='sample')



##Outputing merged files for uORFs etc.

#get the canonical annotation
annocds <- read_compressed_gfile('my_gencode.vM12.annotation.gtf','CDS')
exons <- read_compressed_gfile('my_gencode.vM12.annotation.gtf','exon')
exons %<>% split(.,.$transcript_id)

#genomic seqinfo object
gseqinfo<-fread('cat *.chromsizes')%>%set_colnames(c('seqnames','seqlengths'))%>%do.call(Seqinfo,.)

#transcript seqinfo object
trseqinfo<-exons%>%width%>%sum%>%enframe%>%set_colnames(c('seqnames','seqlengths'))%>%do.call(Seqinfo,.)

#the genomic locations of all ORFs found by the pipeline
orfgendt<-all_orfs_gen%>%as.list%>%map(.%>%{.$ORF_id_tr<-names(.);GR2DT(.)})%>%bind_rows%>%unique

#filter the table of ORFs to uORFs, then join them 
uorfdt <- orfs_dt %>% 
	filter(ORF_category_Tx=='uORF')  %>% 
	select(-seqnames,-start,-end,-width,-strand)%>%
	left_join(orfgendt)

#now export to a file 
uorfdt_cdsfilt <- uorfdt%>%
	DT2GR(gseqinfo)%>%
	subsetByOverlaps(annocds,invert=TRUE)

#export cdsfilt
uorfdt_cdsfilt %>% export('SaTAnn/uORFs.gtf')

#also get counts
bamfiles<-Sys.glob("star/data/*/*.bam")%>%grep(v=T,inv=T,patt='transcript')%>%str_subset('(ribo_|total_)..bam')
uroffcounts<- bamfiles %>% mclapply(mc.cores=20,function(bamfile)
	featureCounts(annot.ext=uorfgtf,isGTFAnnotationFile=TRUE,files=bamfile,GTF.featureType='sequence_feature',nthreads=1,strandSpecific=1,juncCounts=1)
)

#get countmat
uorfcountmat <- uroffcounts%>%
	map('counts')%>%
	#extract colnames, filenames with slashes as dots
	map( ~ set_colnames(.,colnames(.)%>%str_split('\\.')%>%.[[1]]%>%tail(2)%>%head(1)))%>%
	do.call(cbind,.)

#now export this count matrix as a file for xtail to use
uorfcountmat%>%
	as_tibble(rownames='feature_id')%>%
	mutate(feature_id = paste('uORF_',feature_id))%>%
	write_tsv('SaTAnn/uORFs.feature_counts')
######################Number of hits etc.


# orfs_dt%<>%group_by(sample)%>%mutate(fdr =  p.adjust(pval, method='fdr'))
# orfs_dt %<>% dplyr::filter(fdr<0.05)%>%ungroup


#TODO - eveyrthing in linc

n_genes_translated <- orfs_dt %>% group_by(sample)%>%dplyr::summarise(n_genes=n_distinct(gene_id))

readthroughs <- satannorfs%>%map_df(.%>%.$ORFs_readthroughs%>%length)%>%stack%>%set_colnames(c('n_readthroughs','sample'))
readthroughs$sample%<>%as.character
readthroughs$n_readthroughs%<>%as.numeric

genetypehits <- satannorfs%>%lapply(
	.%>%.$ORFs_tx%>%.[T,c('gene_biotype','gene_id')]%>%mcols%>%as.data.frame%>%distinct%>%group_by(gene_biotype)%>%tally%>%
	mutate(gene_biotype = ifelse(gene_biotype%>%str_detect('pseudogene'),'pseudogene',gene_biotype))%>%
	group_by(gene_biotype)%>%summarise(n=sum(n))%>%
	dplyr::filter(gene_biotype %in% c('antisense','lincRNA','protein_coding','pseudogene')))%>%
	bind_rows(.id='sample')%>%
	spread(gene_biotype,n)
#but does our uorf 


satann_summary_table <- genetypehits%>%left_join(readthroughs)

satann_summary_table%>%write_tsv('satann_summary.tsv'%T>%{message(normalizePath(.))})

orfs_dt%>%distinct(gene_id,ORF_category_Tx,.keep_all=T)%>%.$ORF_category_Tx%>%table%>%stack

satannorfs[[1]]$ORFs_tx[,c('gene_biotype','gene_id')]%>%mcols%>%as.data.frame

uorfgtf<-"SaTAnn/uORFs.gtf"
cdsgtf<-"my_gencode.vM12.annotation.cds.gtf"
bamfiles<-"star/data/RPI2_80SE13_2/RPI2_80SE13_2.bam"

bamfiles<-Sys.glob("star/data/*/*.bam")%>%grep(v=T,inv=T,patt='transcript')%>%str_subset('ribo_..bam')



getreadlengthcounts <- function(gtf,bamfile){
	fread(str_interp("tail -n +4 ${gtf} | awk '{print $1,$4,$5,$6,$8,\"-\"}'  | samtools view -L - ${bamfile} | awk '{a[length($10)] += 1} END { for (key in a) { print key \"\\t\" a[key] } }' "))
}

uorfreadlengths <- mclapply(bamfiles,function(.) getreadlengthcounts(uorfgtf,.))
uorfreadlengths%<>%setNames(bamfiles)%>%bind_rows(.id='bamfile')%>%mutate(region='uORF',sample=basename(dirname(bamfile)))%>%select(readlength=V1,count=V2,region,sample)
cdsreadlengths <- mclapply(bamfiles,function(.) getreadlengthcounts(cdsgtf,.))
cdsreadlengths%<>%setNames(bamfiles)%>%bind_rows(.id='bamfile')%>%mutate(region='CDS',sample=basename(dirname(bamfile)))%>%select(readlength=V1,count=V2,region,sample)

uorfcdsrls<-rbind(
	uorfreadlengths,
	cdsreadlengths
)



#
for (i in 1:2){
	if(i==1) pdf('../plots/uorf_vs_cds_readlengths.pdf'%>%normalizePath%T>%message,w=12,h=12)
	if(i==2) svglite('../plots/uorf_vs_cds_readlengths.svg'%>%normalizePath%T>%message,w=12,h=12)
	print(
	uorfcdsrls%>%
		# filter(sample%>%str_detect('ribo'))%>%
		# filter(sample%>%str_detect('80S|Poly'))%>%
		group_by(sample,region)%>%
		mutate(count_frac=count/sum(count))%>%
		ggplot(aes(x=readlength,y=count_frac,color=sample))+
			facet_grid(region~.,scale='free')+
			geom_line()+
			scale_x_continuous(breaks=seq_len(max(uorfcdsrls$readlength)))+
			theme_minimal()+
			theme(panel.grid.minor.x=element_blank(),panel.grid.minor.y=element_blank())
			)
	dev.off()
}

