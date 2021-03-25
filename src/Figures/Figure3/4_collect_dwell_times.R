

library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(GenomicAlignments)
require(txtplot)
require(Rsamtools)
require(rlang)
#coverageplots

library(data.table)

base::source(here::here('src/R/Rprofile.R'))
if(!exists('cdsgrl')) base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
if(!exists('iso_tx_countdata')) load('data/1_integrate_countdata.R')
# base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))

displaystagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
stageconv = names(displaystagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
stagecols <- displaystagecols%>%setNames(names(stageconv))

CODOONSOFINTEREST<-DNAStringSet(c('TCA','TCG','TCC','TCT'))
codons2scan <- c(CODOONSOFINTEREST,reverseComplement(CODOONSOFINTEREST))
# codons2scan <- c(CODOONSOFINTEREST,(CODOONSOFINTEREST))
ctrlcodons<-c('TTT','GGG','AAA','CCC')
codons2scan <- c(CODOONSOFINTEREST,ctrlcodons)
#codons2scan <- CODOONSOFINTEREST

FLANKCODS<-15

ribocovtrs<-readRDS(here('data/ribocovtrs.rds'))

cds2use <- cdsgrl[ribocovtrs]

ref<-Rsamtools::FaFile(fafile)

cdsseq <- cds2use%>%extractTranscriptSeqs(x=fafile)



allcodons=getGeneticCode()
#get code for mitcondrial coding genes

bamfiles<-Sys.glob(here::here('pipeline/star/data/*/*_*bam'))
library(here)
bams <- here('pipeline/star/data/*/*.bam')%>%Sys.glob%>%
	str_subset(neg=TRUE,'transcript')%>%
	str_subset('_ribo_|Poly')

i=1
cdscodons <- 	lapply(seq_along(allcodons),function(i){
	#	
	codon=names(allcodons)[[i]]
	message(codon)
	codmatches<-vmatchPattern(pattern=codon,cdsseq%>%subseq(4,Inf))#exclude the start ccodon
	#
	nmatches = 	codmatches%>%elementNROWS 
	#
	strands<-strand(cds2use[names(codmatches)])[as(rep(1,length(cds2use)),'IntegerList')]%>%unlist
	strands <- rep(as.vector(strands),nmatches)
	matchgr<-codmatches%>%unlist%>%GRanges(names(.),.,strand=strands)
	matchgr%<>%IRanges::shift(3)
	matchgr%<>%subset(start %%3 == 1)
	matchgr = matchgr[start(matchgr)>FLANKCODS]
	trends = cds2use%>%width%>%.[seqnames(matchgr)]%>%sum
	matchgr = matchgr[end(matchgr)<trends - FLANKCODS]
	#
	codmatchonchr<-matchgr%>%mapFromTranscripts(cds2use)
	codmatchonchr%<>%subset(width==3)
	#
	codmatchwindows<-codmatchonchr%>%resize(width(.)+(2*(3*FLANKCODS)),'center')
	#
	# codmatchwindows%>%sample(1e3)%>%getSeq(x=ref)
	codmatchwindows
})

cdscodons%<>%setNames(names(allcodons))

codons2scan<-allcodons%>%names%>%
	# str_subset('TTC|GTC|CAC|AAC|ATG')%>%
	identity
codons2scan%<>%setNames(.,.)
#check this worked

stopifnot(all( c('A','T','G') == (cdscodons[['ATG']]%>%sample(1e3)%>%getSeq(x=ref)%>%as.matrix%>%apply(2,table)%>%simplify2array%>%head((1+FLANKCODS)*3)%>%tail(3)%>%unlist)%>%names))

bams2scan<-bams%>%
	identity

widths <- (25:31)%>%setNames(.,paste0('rl',.))
bam=bams[1]
getreadfps <- function(bam,cds2use,mapqthresh=200){
	  riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=cds2use%>%unlist)
	  reads <- readGAlignments(bam,param=riboparam)
	  reads%<>%as('GRanges')
	  reads$length <- width(reads)
	  reads%<>%subset(length %in% widths)
	  reads %<>% resize(1,'start')
	  reads
}
fpsitelist <- mclapply(bams2scan,mymemoise(getreadfps),cds2use)
fpsitelist%<>%setNames(bams2scan%>%basename%>%str_replace('\\.\\w+$',''))

################################################################################
########Now collect density over codons
################################################################################

#Now run through the bam file and get a profile for each codon
{
fprofilemats_unproc<-mclapply(mc.cores=1,fpsitelist,mymemoise(function(psites){
	#
	cdscounts <- countOverlaps(cds2use,psites)
	zerocds <- names(cdscounts)[cdscounts < 32]
	#
	#
	cdscounts_densitys <- cdscounts%>%divide_by(sum(width(cds2use)))%>%setNames(names(cds2use))
	#
	stopifnot(all(names(cdscounts_densitys)==names(cds2use)))

	mclapply(mc.cores=20,codons2scan,function(codon){
		#
		message(codon)
		#
		codmatchwindows <- cdscodons[[codon]]
		#
		lapply(widths,function(length_i){
			#
			pospsites<-psites%>%subset(length==length_i)%>%subset(strand=='+')
			negpsites<-psites%>%subset(length==length_i)%>%subset(strand=='-')
			#
			poswinds<-codmatchwindows%>%subset(strand=='+')
			negwinds<-codmatchwindows%>%subset(strand=='-')
			profiles <- 
			c(
				pospsites%>%coverage%>%.[poswinds],
				negpsites%>%coverage%>%.[negwinds],
				NULL
			)
			nposwinds<-poswinds%>%length
			profiles[-(1:nposwinds)]%<>%revElements
			windnames<-c(names(poswinds),names(negwinds))	
			#now rather than covert to a matrix, we can do this to 
			profiles <- profiles[!windnames%in%zerocds]
			normdensities <- cdscounts_densitys[windnames[!windnames%in%zerocds]]
			
			np<-length(profiles)
			n_col<-length(profiles[[1]])

			logsumpvect<-sapply(1:n_col,function(i){
				profiles[as(rep(i,np),'IntegerList')]%>%sum%>%
					mean(na.rm=TRUE)
			})
			if(any(!is.finite(logsumpvect)))
			invisible(logsumpvect)
			enframe(logsumpvect,'position','signal')
		})%>%bind_rows(.id='readlen')
	})%>%bind_rows(.id='codon')
}))%>%bind_rows(.id='sample')
}


################################################################################
########PRocess these
################################################################################

if(all(fprofilemats_unproc$sample %in% 1:100)) fprofilemats_unproc$sample %<>% {names(fpsitelist)[as.numeric(.)]}

codonprofiles<-fprofilemats_unproc
codonprofiles%<>%filter(!codon %in% c('TAG','TAA','TGA'))
codonprofiles%<>%mutate(position = position - 1 - (FLANKCODS*3))
codonprofiles%<>%group_by(readlen)%>%filter(any(signal!=0))
codonprofiles%<>%group_by(readlen,codon,sample)%>%
	mutate(occ_nonorm=signal)%>%
	mutate(occupancy = signal / median(signal))
codonprofiles%<>%select(-signal)
stopifnot(codonprofiles$occupancy%>%is.finite%>%all)
codonprofiles %>%saveRDS('data/codonprofiles.rds')