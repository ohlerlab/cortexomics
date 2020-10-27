#!/usr/bin/env Rscript
	# facet_grid(scale='free',ifelse(model%in%names(assay2model),'Seq Data','MS') ~ . )+
message('loading libraries')
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(DESeq2))
suppressMessages(library(assertthat))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(here))
suppressMessages(library(biomaRt))
suppressMessages(library(testthat))
library(conflicted)
library(zeallot)
library(splines)
library(limma)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)



source('src/R/Rprofile.R')
#get realistic distribution of site strengths from things with only one transcriopt
bamfile = '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/star/data/E13_ribo_1/E13_ribo_1.bam'

#are there enough highly expressed genes with only one transcript?
#transcripts <- 
gtf_gr%>%as.data.frame%>%head(10e3)%>%group_by(gene_id)%>%summarise(n=n_distinct(transcript_id))%>%.$n%>%table
#better yet, just coding regions that only overlap a single transcript.
cdsorig <- cds%>%subset(seqnames%in%c('chr1','chr2'))%>%subset(strand=='+')
riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=200,which=cdsorig)
reads <- readGAlignments(bamfile,param=riboparam)

cdscounts <- cdsorig%>%split(.,.$protein_id)%>%countOverlaps(reads)%>%enframe('pid','count')
cds_orig$cdscount <- cdscounts[match(cds_orig$protein_id,names(cdscounts))]

fpextcds <- cdsorig %>%sort%>%split(.,.$protein_id)%>%ext_grl(12,'end')%>%unlist
singletrcds <- fpextcds%>%.[countOverlaps(.,gtf_gr%>%subset(type=='transcript'))==1]
singletrcds%<>%resize(width(.)-.$phase,'end')

#Now determine how many of these guys are so well covered that I can use them as a fairly 
cov = reads%>%coverage
widths=25:31

singletrcds$readdensity <- singletrcds%>%countOverlaps(reads)%>%divide_by(width(singletrcds))

#thisi s pretty bimoal - lots at 0 and near 1
singletrcds$bases_covered = mean(cov[singletrcds]>0)
singletrcds$bases_fpcov = mean(fpcov[singletrcds]>0)
(bases_fpcov>0.5)%>%table
(bases_covered>0.5)%>%table
(bases_covered==1)%>%table


highcovcds%>%length

widths <- 25:31
width_i=27
readlenfpnum <- coverage(reads%>%subset(width==width_i)%>%GRanges%>%resize(1))[highcovcds]
fpnum <- fpcov[highcovcds]
readlenfpnum[284]/sum(fpnum[284])
readlenfpprob <- readlenfpnum/sum(fpnum)
library(multitaper)

readlenfpprob%>%purrr::reduce(c)%>%as.vector%>%head(200:300)%>%ftestvect

c(readlenfpprob%>%purrr::reduce(c)%>%.[seq(1,length(.),by=3)]%>%sum,
readlenfpprob%>%purrr::reduce(c)%>%.[seq(2,length(.),by=3)]%>%sum,
readlenfpprob%>%purrr::reduce(c)%>%.[seq(3,length(.),by=3)]%>%sum)%>%{./sum(.)}


c(fpcov[highcovcds]%>%purrr::reduce(c)%>%.[seq(1,length(.),by=3)]%>%sum,
fpcov[highcovcds]%>%purrr::reduce(c)%>%.[seq(2,length(.),by=3)]%>%sum,
fpcov[highcovcds]%>%purrr::reduce(c)%>%.[seq(3,length(.),by=3)]%>%sum)%>%{./sum(.)}


probvect <- readlenfpprob%>%purrr::reduce(c)

get_readfp_probvect <- function(reads,widths,highcovcds){
	fpcov = reads%>%as("GRanges")%>%subset(width%in%widths)%>%resize(1)%>%coverage

	fpnum <- fpcov[highcovcds]

	lapply(widths,function(width_i){
		readlencov<-reads%>%subset(width==width_i)%>%GRanges%>%resize(1)
		readlenfpnum <- coverage(readlencov)[highcovcds]
		readlenfpprob <- readlenfpnum/sum(fpnum)
		probvect <- readlenfpprob%>%purrr::reduce(c)
		probvect
	})%>%setNames(paste0('rl',widths))%>%RleList
}

highcovcds <- singletrcds%>%subset(bases_fpcov>1/10)
medcovcdsd <- singletrcds%>%subset(bases_fpcov<1/10)%>%subset(bases_fpcov> 1/20)
readlfpcovvects  <- get_readfp_probvect(reads,widths,highcovcds)

# trlens = c(100,400,1000)%>%setNames(letters[seq_along(.)])
# trcounts = c(20,50,80)

generate_fakecov<-function(readlfpcovvects ,trlens,trcounts){
	#I want to preserve autocorrelative structure of hte Riboseq, so I'll sample from this
	#in continous tracts.
	#That also needs to be done for a set of read lengths
	pvectlen <- readlfpcovvects [[1]]%>%length
	#
	#
	trlens<-unlist(trlens)
	stopifnot(!is.null(names(trlens)))
	stopifnot(length(trlens)==length(trcounts))
	stopifnot(all(trlens<pvectlen))
	randomstarts <- floor(runif(length(trlens),1,pvectlen-trlens))
	pvectsegs <- map2(randomstarts,randomstarts+trlens-1,.f=IRanges)%>%IRangesList
	readlenprops <- readlfpcovvects%>%sum%>%{./sum(.)}
	#
	#
	# readlenprop<-readlenprops[[1]]
	#
	rlnames<-names(readlfpcovvects)
	widths<-rlnames%>%str_replace('rl','')%>%as.numeric
	out <- mapply(widths,readlenprops,FUN=function(width_i,readlenprop){
		rlname <- paste0('rl',width_i)
		segprobs <- RleList('a'=readlfpcovvects [[rlname]])[pvectsegs%>%setNames(rep('a',length(.)))]
		segprobs%<>%setNames(names(trlens))
		segprobs[segprobs!=0] <- (segprobs[segprobs!=0] / sum(segprobs)) * readlenprop * trcounts

	# })
	
	# out%>%lapply(sum)%>%simplify2array%>%rowSums

		segints = segprobs - segprobs 
		segints[segprobs!=0] <-  segprobs[segprobs!=0]%>%lapply(as.vector)%>%map(~rpois(length(.),.x))%>%map(Rle)%>%RleList
		#now use this integer vector as a fake coverage vector to make the granges

		as(segints,"GRanges")%>%subset(score!=0)%>%{suppressMessages({resize(.,width_i)})}%>%rep(.,.$score)%>%{.$score<-NULL;.}
	})%>%{suppressMessages({purrr::reduce(.,c)})}
	out$length <- width(out)
	resize(out,1)
}

highcovcdswhole <- cds%>%subset(protein_id %in% highcovcds$protein_id)
highcovcdswhole%<>%split(.,.$protein_id)
highcovcdswhole

fakefpcov<-generate_fakecov(
		readlfpcovvects,
		width(medcovcdsd)%>%setNames(medcovcdsd$exon_id),
		medcovcdsd$readdensity*width(medcovcdsd)
	)


conflict_prefer('rowMeans','Matrix')

exons <- gtf_gr%>%split(.,.$transcript_id)
medcovcdsdtrs <- medcovcdsd$transcript_id%>%unique


fakefpcov <- fakefpcov%>%resize(1)%>%mapFromTranscripts(medcovcdsd%>%{setNames(.,.$exon_id)})%>%{seqinfo(.)<-seqinfo(FaFile(REF));.}%>%coverage

# fpcov[medcovcdsd]%>%lapply(as.vector)%>%mclapply(ftestvect)%>%simplify2array%>%.[1,]%>%na.omit%>%log10%>%txtdensity
# fakefpcov[medcovcdsd]%>%lapply(as.vector)%>%mclapply(ftestvect)%>%simplify2array%>%.[1,]%>%na.omit%>%log10%>%txtdensity

#but for this to be valid count magnitudes must be similiar
fpcov[medcovcdsd]%>%sum%>%add(1)%>%log10%>%hist(plot=F)%>%{txtplot(.$mids,.$density,ylim=c(0,2))}
fakefpcov[medcovcdsd]%>%sum%>%add(1)%>%log10%>%hist(plot=F)%>%{txtplot(.$mids,.$density,ylim=c(0,2))}

#teests on the high coverage cdsd
realftests<-fpcov[medcovcdsd]%>%lapply(as.vector)%>%mclapply(ftestvect)
fakeftests<-fakefpcov[medcovcdsd]%>%lapply(as.vector)%>%mclapply(ftestvect)

#significance the same?
realftests%>%simplify2array%>%.[2,]%>%na.omit%>%hist(plot=F)%>%{txtplot(.$mids,.$density,ylim=c(0,2))}
fakeftests%>%simplify2array%>%.[2,]%>%na.omit%>%hist(plot=F)%>%{txtplot(.$mids,.$density,ylim=c(0,2))}

realftests%>%simplify2array%>%.[2,]%>%{.<0.05}%>%table
fakeftests%>%simplify2array%>%.[2,]%>%{.<0.05}%>%table

#what about the fold difference of density vs. width?

fakeseqinfo<-lengths(fakefpcov)%>%enframe%>%as.data.frame%>%{Seqinfo(.[[1]],.[[2]])}
is_out_of_bounds(medcovcdsd,fakeseqinfo)%>%which
is_out_of_bounds(medcovcdsd,seqinfo(exons))%>%which

medcovcdsd[162]

max(end(medcovcdsd%>%subset(seqnames=='chr2')))
max(end(exons%>%unlist%>%subset(seqnames=='chr2')))

seqinfo(exons)
fakefpcov%>%lengths


csd


################################################################################
########Now with kallisto 
################################################################################

kallistofile <- 'pipeline/kallisto_old//data/E13_total_2/out.d/abundance.tsv'
kallistores<-kallistofile%>%fread

highcovcds
fakecovreaadsd <- kallistores%>%
	filter(target_id %in% highcovcds$transcript_id)%>%
	{
	generate_fakecov(
		readlfpcovvects,
		setNames(.$length,.$target_id),
		setNames(.$est_counts,.$target_id)
	)
}



#We also need to 



#Now armed without our function that gives us a set of probabilities for each readlength, of n bp






#for each of our singletr cds
#Get a track of read fp sites emmision props, per read length, beginning in phase 0


#read length?

#checking with a bam file export

GR(c('a:4-5','a:7-8'))%>%split(.,seqnames(.))%>%countOverlaps(GR('a:5-8'))


#realistic distribution of transcript abundances



#Now given what's known about the transcript levels form the RNA 
#Use it to generate several test sets - 1 scenario where all isoforms are equally expressed
#Scenario where the same amount of diversity exists at transcript and translatome level.
#Run SaTAnn
#Look at results.


#Note - I think this is probably a bit conservative, because i think that the variance in Riboseq sites in the lowly

#Note - how do I realistically distribute Riboseq signal?
#I can take really expressed things and then assume this is an accurate picture of 
#the variance within the same CDS
#Then for each of these (I start a bit upstream of each start codon)
#I need amodel I can input an abundance into that outputs the distribution of riboseq along a sequence.
#I guess this will depend on a whole bunhc of things - like the particular libar
	#This is actually an interesting question though.
	
#I need to be dealing with 5' peaks? Actually no, not ith E13 data
#But I should calculate the densities 
#expressed genes is just going to be higher than in the highly translated, housekeeping genes.
#Still, let's try