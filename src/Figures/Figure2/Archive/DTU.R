
#https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html#drimseq
################################################################################
################################################################################
base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("src/Figures/Figure0/0_load_annotation.R")
}
library(GenomicFeatures)
if(!'DRIMSeq'%in%installed.packages()) BiocManager::install("DRIMSeq")
iso_tx_countdata <- readRDS(here('data/iso_tx_countdata.rds'))
library(tximport)
library(tidyverse)
#get the transript level files
rnasalmonfiles = Sys.glob(here('pipeline/salmon/data/*/quant.sf'))
dpoutfiles = Sys.glob(here('pipeline/deepshapeprime/fakesalmonfiles/*ribo*/*'))
#naem them
allquantfiles = c(rnasalmonfiles,dpoutfiles)
names(allquantfiles) <- allquantfiles%>%dirname%>%basename
#load only the shared transcripts
dptrs = dpoutfiles[[1]]%>%fread%>%.$Name
salmontrs = allquantfiles[[1]]%>%fread%>%.$Name%>%str_extract('ENSMUST\\w+')
trs = intersect(dptrs,salmontrs)
#Now cget our transcript-level counts
tx2genemap = mcols(cds)%>%as.data.frame%>%distinct(transcript_id,gene_id)


################################################################################
########3' UTR analysis
################################################################################
#get 3' UTRs

#get table of these per transcript

#using our abundances, get the estimated 3' UTR length per gene at each timepoint

#get the change in average 3' UTR length per gene

#compare this to xtail changes in gene-TE


gtftxdb = makeTxDbFromGRanges(gtf_gr)

tputrs = threeUTRsByTranscript(gtftxdb,use.names=TRUE)

genetputrlens = iso_tx_countdata$abundance%>%
	as.data.frame%>%rownames_to_column('transcript_id')%>%
	pivot_longer(-transcript_id,names_to='dataset',values_to='TPM')%>%
	separate(dataset,c('time','assay','rep'))%>%
	left_join(sum(width(tputrs))%>%enframe('transcript_id','tputrlen'))%>%
	mutate(tputrlen = replace_na(tputrlen,100))%>%
	left_join(tx2genemap)%>%
	filter(assay=='total')%>%
	group_by(gene_id,time,rep)%>%
	summarise(tputrlen=tputrlen*(TPM/sum(TPM)))%>%
	group_by(gene_id,time)%>%
	summarise(tputrlen=mean(tputrlen))

genetputrlenchange = genetputrlens%>%group_by(gene_id)%>%mutate(tp_utr_change = tputrlen - tputrlen[time=="E13"])


genetputrlenchange%>%filter(time=='P0')%>%filter(between(tp_utr_change,-1e3,1e3))%>%.$tp_utr_change%>%txtdensity

# xtailfoldchange%>%
# 	select(gene_id,gene_name,time,te_log2fc=log2fc)%>%
# 	inner_join(genetputrlenchange)%>%
# 	filter(time=='E175')%>%
# 	filter(tp_utr_change!=0)%>%
# 	filter(between(tp_utr_change,-1e3,1e3))%>%
# 	{quicktest(.$tp_utr_change,.$te_log2fc)}

# fread('tables/xtailTEchange.tsv')%>%
# 	inner_join(genetputrlenchange)%>%
# 	filter(time=='E175')%>%
# 	{fisher.test(table(.$tp_utr_change>10,.$up==1))}

inclusiontable(tputrs$gene_id , cds$gene_id)

# exons%>%setdiff(cds,rev=TRUE)

#create our Drimseq objects
cts <- iso_tx_countdata$counts
cts <- cts[rowSums(cts) > 0,]

tx2genemap%<>%set_colnames(c('tr_id','g_id'))

txdf = tx2genemap%>%select(GENEID=g_id,TXNAME=tr_id)%>%group_by(GENEID)%>%mutate(ntx=n())

txdf <- txdf[match(rownames(cts),txdf$TXNAME),]
all(rownames(cts) == txdf$TXNAME)

counts <- data.frame(gene_id=txdf$GENEID,
                     feature_id=txdf$TXNAME,
                     cts)
library(DRIMSeq)
allcountdesign = colnames(iso_tx_countdata$counts)%>%data.frame(sample=.)%>%separate(sample,into=c('time','assay','rep'),remove=F)
allcountdesign = allcountdesign%>%arrange(assay=='ribo')%>%mutate(assay=as_factor(assay))%>%as.data.frame%>%set_rownames(.$sample)

allcountdesign%<>%mutate(sample_id=sample)
d <- dmDSdata(counts=counts, samples=allcountdesign%>%as.data.frame)
n.small = allcountdesign%>%group_by(time,assay)%>%tally%>%.$n%>%min
n = nrow(allcountdesign)
d <- dmFilter(d,
              min_samps_feature_expr=n.small, min_feature_expr=10,
              min_samps_feature_prop=n.small, min_feature_prop=0.1,
              min_samps_gene_expr=n, min_gene_expr=10)
design_full <- model.matrix(~assay*time, data=DRIMSeq::samples(d))




#Run drimseq
d = d
d <- dmPrecision(d, design=design_full)
d <- dmFit(d, design=design_full)
d <- dmTest(d)
d %>% saveRDS(here('data/d.rds'))


