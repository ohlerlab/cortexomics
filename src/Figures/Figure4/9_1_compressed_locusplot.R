
#function that takes all the granges in a region, and compresses them down so introns
#are all a fixed distance
library(Gviz)
library(rtracklayer)
library(tidyverse)
library(magrittr)
library(Rsamtools)
library(GenomicAlignments)

mymem<-projmemoise(function(){rnorm(1)})

if(!exists("cdsgrl")) {
 # base::source(here::here('src/R/Rprofile.R'))
  base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
}
tpcols <- c(E13='#214098',E145='#2AA9DF',E16='#F17E22',E175='#D14E28',P0='#ED3124')


mymem()
mycache$keys()
projmemoise



options(ucscChromosomeNames=FALSE) 

#toy ranges to shift
exrangestoshift <- c(GRanges('a',IRanges(4,8)),GRanges('a',IRanges(12,13)),GRanges('a',IRanges(20,30)),GRanges('a',IRanges(40,42),strand='-'))
#

#this is the function that shifts a bunch of ranges to minimize the space between them, down to a ceiling maxwidth
maxwidth=1
get_comp_shift<-function(rangestoshift,maxwidth=1){
	rangegaps <- rangestoshift%>%{strand(.)<-'*';.}%>%gaps%>%.[-1]
	rangegaps%>%width
	#maxwidth <- 4
	toshift<-rep(0,length(rangestoshift))
	i=3
	for(i in seq_along(rangegaps)){
		adj <- width(rangegaps[i])-maxwidth
		if(adj>0){
			rightofgap <- start(rangestoshift)>end(rangegaps[i])
			toshift[rightofgap] <- toshift[rightofgap] + adj
		}
	}
	toshift %<>% {. + min(start(rangestoshift))-1}
	#rangestoshift<-GenomicRanges::shift(rangestoshift,-rangestoshift$toshift)
	toshift
}
# get_comp_shift(igeneanno,maxwidth=Inf)
# compress_gaps<-function(gr,maxwidth){GenomicRanges::shift(gr,-get_comp_shift(gr,maxwidth))}

# gr=c(GRanges('a:10-12'),GRanges('a:22-24'),GRanges('a:39-42'))
# compress_gaps(gr,2)




# # shiftedranges = compress_gaps(rangestoshift,1)

# #generate fake data
# cov = coverage(exrangestoshift)
# covch = cov[['a']]
# covch[covch!=0] %<>% {.*rpois(length(.), as.vector(.) * 32)}
# cov[['a']]=covch
# fakecov = as(cov,'GRanges')

# #Now we need to shift this over to the compressed space
# ov = findOverlaps(fakecov,exrangestoshift)
# fakecov = shift(fakecov[ov@from],-shiftedranges$toshift[ov@to])





########
library(rtracklayer)

# fread('head -n 100 pipeline/my_gencode.vM12.annotation.gtf')%>%
# 	as.matrix%>%
# 	apply(1,paste,collapse='\t')%>%
# 	lapply(paste,collapse='\n')%>%
# 	{import(text=.,format='gtf')}

# assert_that(file.exists(annofile))
# annotation <- annofile%>% {rtracklayer::import(.)}


#get counts for the exons
# bam=dexseqbams[[1]]
# windows = annotation

#get offset values
offsets<-read_tsv(here('ext_data/offsets_manual.tsv'))

#get the bams for each library
# bams <- Sys.glob(here('pipeline/star/data/*/*.bam'))%>%str_subset(negate=TRUE,'transcript')
#
# dexseqbams<-bams%>%str_subset('_(ribo)')
# dexseqbamsrna<-bams%>%str_subset('_(total)')

# bams=dexseqbams

# annotation%<>%GenomicRanges::shift(-annotation$toshift)
totbams = Sys.glob('pipeline/star/data/*total*/*[12].bam')
ribobams = Sys.glob('pipeline/star/data/*ribo*/*[12].bam')
ribobams%<>%setNames(.,basename(dirname(.)))
totbams%<>%setNames(.,basename(dirname(.)))

if(!exists('peptidemsdata')){
	peptidemsdata = fread('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/MS_Data_New/Ages_Brain_PEP_summ.txt')%>%
	rename('gene_name'=Gene.names)

	peptidemsdata%>%colnames
	normpepmat = peptidemsdata%>%select(matches('Intensity.*input'))%>%as.matrix%>%
		proDA::median_normalization(.)

	normpepmat = colnames(normpepmat)%>%str_extract('(?<=\\.).*?(?=_input)')%>%
		str_replace('p','')%>%
		{split(seq_along(.),.)}%>%
		lapply(.,function(colgrp){
			normpepmat[,colgrp]%>%rowMedians(na.rm=T)
		})%>%
		simplify2array
	rownames(normpepmat) <- peptidemsdata$Sequence



	# 'AAPAETDQR',proteinsequences['ENSMUST00000042857']
	# vmatchPattern('AAPAETDQR',translate(proteinsequences['ENSMUST00000042857']))
	# vmatchPattern('AAPAETDQR',translate(proteinsequences['ENSMUST00000114415']))

	extrseq = function(gr)GenomicFeatures::extractTranscriptSeqs(gr,x=FaFile('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'))

	# vmatchPattern('MERRSESPCLRDSPD',translate(
	# 	extrseq(cdsgrl[['ENSMUST00000042857']][1]%>%split(.,.$transcript_id))
	# ))%>%GRanges%>%mapFromTranscripts(cdsgrl)

	# vmatchPattern('MERRSESPCLRDSPD',translate(
	# 	extrseq(cdsgrl[['ENSMUST00000114415']][1]%>%split(.,.$transcript_id))
	# ))%>%GRanges%>%mapFromTranscripts(cdsgrl)

	# cdsgrl%<>%sort

	libstg = c(ribobams,totbams)%>%dirname%>%basename%>%str_split_fixed('_',3)%>%.[,1]
	libstg%<>%setNames(c(ribobams,totbams)%>%dirname%>%basename)
	allvoom<-readRDS('data/allvoom.rds')
	libsizes = allvoom@.Data[[1]]%>%{setNames(.$lib.size,rownames(.))}
	rbamcov = allvoom@.Data[[1]]$lib.size
	offsets%<>%mutate(comp='nucl')

	cdsgrl = cdsgrl[str_order_grl(cdsgrl)]
	proteinsequences = GenomicFeatures::extractTranscriptSeqs(cdsgrl,x=FaFile('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'))
	aaproteinsequences = translate(proteinsequences)
	trlist = ids_nrgname%>%distinct(gene_name,transcript_id)%>%{split(.[[2]],.[[1]])}
	iso_tx_countdata <- readRDS(here('data/iso_tx_countdata.rds'))

}

igene='Flna'

COMPRESS=T
for(COMPRESS in c(TRUE,FALSE)){
for(igene in c('Satb2','Flna','Nes','Bcl11b','Tle4'))
# for(igene in c('Satb2','Bcl11b'))
{
{

types = c('UTR','CDS')
if(!COMPRESS) types = c(types,'gene')


igeneanno <- fread(str_interp('grep -ie ${igene} /fast/groups/ag_ohler/work/dharnet_m/cortexomics/ext_data/gencode.vM12.annotation.gtf'))%>%
	as.matrix%>%apply(1,paste0,collapse='\t')%>%paste0(collapse='\n')%>%import(text=.,format='gtf')%>%
	subset(gene_name==igene)%>%
	subset(type%in%types)
stopifnot(length(igeneanno)>0)
if(COMPRESS){
	igeneanno$toshift <- get_comp_shift(igeneanno,maxwidth=10)
}else{
	igeneanno$toshift <- 0
}
igeneanno$transcript_id%<>%trimids


exontrack =
  igeneanno%>%
  shift(.,-.$toshift)%>%
  # {.$transcript_id %<>% str_replace('\\.[0-9]+$','');.}%>%
  subset(type%in%setdiff(types,'gene'))%>%
  .[,c('type','transcript_id')]%>%
  {.$feature<-as.character(.$type);.}%>%
  {.$transcript<-as.character(.$transcript_id);.}%>%
  Gviz::GeneRegionTrack(.,thinBoxFeature=c("UTR"))
  
library(Gviz)
}


get_cov_track <- function(igeneanno,bams=ribobams,offsets_=offsets){
	#get the reads
	# igeneanno = igeneanno%>%subset(gene_name=igene)
	st = as.character(strand(igeneanno[1]))
	ichr = as.character(seqnames(igeneanno))[[1]]
	# psites = lapply(ribobams,function(bam)get_genomic_psites(bam,igeneanno,offsets_,comps=NULL))
	riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=200,which=igeneanno)
  	psites <- lapply(bams,function(bam)readGAlignments(bam,param=riboparam))

	psites%<>%lapply(coverage)
	normcov = map2(psites,libsizes[names(psites)],~ (.x / .y)*1e6)
	tpnormcov = normcov%>%{split(.,libstg[names(.)])} %>%lapply(function(trackpair) purrr::reduce(.f=`+`,trackpair))
	tpnormcov = tpnormcov%>%map(GRanges)%>%GRangesList%>%unlist%>%subset(score!=0)
	bp1tpnormcov=tpnormcov@ranges%>%as('IntegerList')%>%setNames(.,seq_along(.))%>%stack%>%as.data.frame
	bp1tpnormcov$stage = names(tpnormcov)[as.integer(bp1tpnormcov$name)]
	bp1tpnormcov$score = tpnormcov$score[as.integer(bp1tpnormcov$name)]
	bp1tpnormcov = bp1tpnormcov%>%select(start=value,stage,score)%>%
		spread(stage,score)%>%
		mutate(seqnames=ichr,end=start)%>%
		GRanges(.)
	

	ov = findOverlaps(bp1tpnormcov,igeneanno)
	ov = ov[!duplicated(queryHits(ov)),]
	bp1tpnormcov = GenomicRanges::shift(bp1tpnormcov[ov@from],-igeneanno$toshift[ov@to])
	strand(bp1tpnormcov) = st
	# bp1tpnormcov
	stopifnot(all(tps%in%colnames(mcols(bp1tpnormcov))))
	# for(tp in tps) mcols(bp1tpnormcov)[[tp]]%<>%{log2(.+(.5*min(.,na.rm=T)))}
	DataTrack(bp1tpnormcov,type='hist',col = tpcols,groups=tps)
}
# get_cov_track(igeneanno)

 
peptidemsdata_=peptidemsdata;cdsgrl_=cdsgrl;aaproteinsequences_ = aaproteinsequences
get_gene_pepgr = function(igene,peptidemsdata_=peptidemsdata,cdsgrl_=cdsgrl,aaproteinsequences_ = aaproteinsequences){
	gpep_df = peptidemsdata_%>%filter(gene_name==igene)
	peptides=gpep_df$Sequence
	# peptides = 'AAPAETDQR'
	ig_seqs = aaproteinsequences_[trlist[[igene]]]
	#
	peptidegr = peptides%>%setNames(.,.)%>%map_df(.id='peptide',function(peptide){
		vmatchPattern(peptide,ig_seqs)%>%
		as.data.frame%>%
		transmute(start,end,width,seqnames=names(ig_seqs)[group])
	})%>%GRanges

	gpeptidegr = peptidegr%>%{
		x=.
		x = resize(x,1,'start',ignore.strand=TRUE)
		end(x)  = ((start(x)-1)*3)+1
		start(x)  = end(x)
		# end(x)  = start(x)
		out = mapFromTranscripts(x,cdsgrl_)
		out$transcript_id = names(cdsgrl_)[out$transcriptsHits]
		out$peptide = x$peptide[out$xHits]
		out
		}
	# gpep	
	gpeptidegr$xHits=NULL
	gpeptidegr$transcriptsHits=NULL
	mcols(gpeptidegr)%<>%cbind(normpepmat[match(gpeptidegr$peptide,rownames(normpepmat)),])
	gpeptidegr
}
# get_gene_pepgr('Satb2')

get_peptide_track<-function(igeneanno){
	igeneanno%<>%subset(type!='gene')
	igene = igeneanno$gene_name[1]
	pepgr = get_gene_pepgr(igene)
	pepgr = pepgr[!duplicated(pepgr$peptide),]
	# igeneanno = annotation%>%subset(gene_name=igene)
	ov = findOverlaps(pepgr,igeneanno)
	ov = ov[pepgr$transcript_id[queryHits(ov)]== igeneanno$transcript_id[subjectHits(ov)],]
	ov = ov[!duplicated(queryHits(ov)),]
	pepgr = GenomicRanges::shift(pepgr[ov@from],-igeneanno$toshift[ov@to])
	# start(pepgr)%>%txtplot
	pepgr$transcript_id=NULL
	pepgr = unique(pepgr)
	pepgr = sort(pepgr)
	# pepgr[60:70]
	# pepgr = pepgr[70:80]
	DataTrack(pepgr[,tps],cols = tpcols,groups=tps)	
	# pepgr
}

motifpattern='TGTANATA'

getmotiftrack = function(motifpattern,igeneanno){
		igeneanno%<>%subset(type!='gene')
spligeneanno = igeneanno%>%split(.,.$transcript_id)
spligeneanno%<>%sort_grl_st
trs = names(spligeneanno)
igeneannoseq = spligeneanno%>%extrseq

motifgr = vmatchPattern(motifpattern,fixed=F,igeneannoseq)%>%
	as.data.frame%>%
	transmute(start,end,width,seqnames=names(igeneannoseq)[group])%>%
	GRanges%>%
	mapFromTranscripts(spligeneanno)

ov = findOverlaps(motifgr,igeneanno)
ov = ov[trs[motifgr$transcriptsHits][queryHits(ov)]== igeneanno$transcript_id[subjectHits(ov)],]
ov = ov[!duplicated(queryHits(ov)),]
motifgr = GenomicRanges::shift(motifgr[ov@from],-igeneanno$toshift[ov@to])
motifgr = unique(motifgr)
AnnotationTrack(motifgr)
}


peaktrack = function(gr,igeneanno){		
	ov = findOverlaps(gr,igeneanno,select='first')
	gr = GenomicRanges::shift(gr,-igeneanno$toshift[ov])
	gr = unique(gr)
	AnnotationTrack(gr)
}
getstart_track = function(igeneanno){
	spligeneanno = igeneanno%>%subset(type=='CDS')%>%split(.,.$transcript_id)%>%.[str_order_grl(.)]
	spligeneanno = spligeneanno%>%resize_grl(1,'start')%>%unlist
	spligeneanno %<>% GenomicRanges::shift(-spligeneanno$toshift)
	AnnotationTrack(spligeneanno,background.panel = "grey")
}

# get_peptide_track('Satb2',annotation)

# ranges('Satb2',annotation)
tps = names(tpcols)

zhangetal_pum2clip = GRanges('chr1:56794013-56794101:-')

{
nametrack = function(x,tnm) x%>%{.@name=tnm;.}
#now plot
options(ucscChromosomeNames=FALSE)
compressstr = if(COMPRESS) '' else '_NONcompressed_'
plotfile<- here(paste0('plots/locusplot_',compressstr,igene,'.pdf'))
pdf(plotfile,w=24,h=12)
plotTracks(list(
	get_cov_track(igeneanno,totbams)%>%nametrack('RNA-Seq (RPM)'),
	get_cov_track(igeneanno,ribobams)%>%nametrack('Ribo-Seq (RPM)'),
	get_peptide_track(igeneanno)%>%nametrack('log2(Normalized Intensity)'),
 	exontrack%>%nametrack('Transcripts'),
 	getstart_track(igeneanno)%>%nametrack('AUG'),
 	getmotiftrack('TGTANATA',igeneanno)%>%nametrack('Pum2 Motifs'),
 	peaktrack(zhangetal_pum2clip,igeneanno)%>%nametrack('Zhang et al Clip')
 	
),col = tpcols)
dev.off()
message(normalizePath(plotfile))
}

}
}

#checking if satb2 pum2s overlap the clip for pum2
if(F){
	satb2trs = igeneanno$transcript_id%>%unique
	satb2exons = exonsgrl[satb2trs]
	satb2isseq = satb2exons%>%sort_grl_st%>%extrseq

	motifgr = vmatchPattern(motifpattern,fixed=F,satb2isseq)%>%
		as.data.frame%>%
		transmute(start,end,width,seqnames=names(satb2isseq)[group])%>%
		GRanges%>%
		mapFromTranscripts(exonsgrl[satb2trs])
	motifgr%<>%unique

	motseqs = motifgr%>%split(.,seq_along(.))%>%extrseq
	stopifnot(all(as.character(motseqs)%>%str_detect('TGTA.ATA')))
	#http://genesdev.cshlp.org/content/31/13/1354.full
	#from table s1
	zhangetal_pum2clip = GRanges('chr1:56794013-56794101:-')
	motifgr%>%overlapsAny(zhangetal_pum2clip)
	motifgr%>%distance(zhangetal_pum2clip)

}
