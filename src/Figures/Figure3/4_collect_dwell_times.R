

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
codonprofiles <- readRDS('data/codonprofiles.rds')



trna_ab_df_samp = allcodsigmean_isomerge[c("fraction", "time", "sample", "anticodon", "abundance", "codon",
"weightedusage", "availability","rep")]

codon_data <- trna_ab_df_samp%>%
	filter(fraction=='Total')%>%select(-fraction,-sample)%>%
    select(time,rep,codon,abundance,availability)%>%
    group_by(time,codon)%>%
    summarise_at(vars(one_of(c('abundance','availability'))),list(mean))


lexp = 3+3 #include positions for positions corresponding to bigger offsets
rexp = 3#include positions for positions corresponding to smaller offsets
# codonprofiles%>%
codondata <- codonprofiledat%>%
# codondata <- rustprofiledat%>%
	group_by(sample)%>%
	mutate(readlen=paste0('rl',readlen))%>%
	safe_left_join(offsets%>%select(readlen,offset))%>%
    # filter(position== -offset-3)%>%
    # mutate(-offset-3,-offset)%>%head%>%as.data.frame
	# group_by(sample,readlen,codon)%>%group_slice(1)%>%
    filter(position <= -(offset-rexp))%>%
    filter(position >= -(offset+lexp))%>%
    separate(sample,c('time','assay','rep'))%>%
    group_by(time,codon,rep)%>%
    summarise(dwell_time = sum(occ_nonorm))%>%
    left_join(codon_data)
	
#now plot
plotfile<- here(paste0('plots/','stage_dtdist','.pdf'))
pdf(plotfile)
codondata%>%
	ggplot(.,aes(fill=time,x=dwell_time))+
	geom_density(alpha=I(0.5))+
	scale_fill_manual(name='stage',values=stagecols)+
	scale_x_continuous(paste0('Broad Dwell Time'))+
	# scale_y_continuous(paste0(''))+
	ggtitle(paste0('Dwell Time Dist'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))

	
#now plot
repsumcodondata<-codondata%>%
	group_by(time,codon)%>%summarise_at(vars(abundance,availability,dwell_time),mean)
plotfile<- here(paste0('plots/','bDT_ab_vs_dt','.pdf'))
pdf(h=15,w=5,plotfile)
col1=sym('dwell_time')
col2=sym('abundance')
corlabel = repsumcodondata%>%filter(is.finite(!!col1),is.finite(!!col2))%>%
		group_by(time)%>%
		summarise(tidy(cor.test(!!col1, !!col2)))
corlabel = corlabel%>%
	mutate(
		pformat=format(p.value,format='e',digits=4),
		pvalstring = ifelse(p.value > 0.001,round(p.value,4),pformat),
		labl=paste0('rho = ',round(estimate,3),'\n','pval = ',pvalstring))
nlabel=repsumcodondata%>%group_by(time)%>%summarise(labl=paste0('N=',n()))
repsumcodondata%>%
	group_by(time,codon)%>%summarise_at(vars(abundance,availability,dwell_time),mean)%>%
	ggplot(.,aes(y=abundance,x=dwell_time))+
	geom_point(alpha=I(0.5))+
	geom_smooth(method='lm')+
	facet_grid(time~.)+
	scale_fill_manual(name='stage',values=stagecols)+
	scale_x_continuous(paste0('Broad Dwell Time'))+
		geom_text(show.legend=F,data=corlabel,
			hjust=1,vjust=1,x= Inf,y=Inf,aes(color=NULL,label=labl))+
	geom_text(show.legend=F,data=nlabel,
			hjust=0,vjust=1,x= -Inf,y=Inf,aes(color=NULL,label=labl))
	# scale_y_continuous(paste0(''))+
	ggtitle(paste0('Dwell Time Dist'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))	

#now plot
repsumcodondata<-codondata%>%
	group_by(time,codon)%>%summarise_at(vars(abundance,availability,dwell_time),mean)
plotfile<- here(paste0('plots/','bDT_av_vs_dt','.pdf'))
pdf(h=15,w=5,plotfile)
col1=sym('dwell_time')
col2=sym('availability')
corlabel = repsumcodondata%>%filter(is.finite(!!col1),is.finite(!!col2))%>%
		group_by(time)%>%
		summarise(tidy(cor.test(!!col1, !!col2)))
corlabel = corlabel%>%
	mutate(
		pformat=format(p.value,format='e',digits=4),
		pvalstring = ifelse(p.value > 0.001,round(p.value,4),pformat),
		labl=paste0('rho = ',round(estimate,3),'\n','pval = ',pvalstring))
nlabel=repsumcodondata%>%group_by(time)%>%summarise(labl=paste0('N=',n()))
repsumcodondata%>%
	ggplot(.,aes(y=availability,x=dwell_time))+
	geom_point(alpha=I(0.5))+
	geom_smooth(method='lm')+
	facet_grid(time~.)+
	scale_fill_manual(name='stage',values=stagecols)+
	scale_x_continuous(paste0('Broad Dwell Time'))+
		geom_text(show.legend=F,data=corlabel,
			hjust=1,vjust=1,x= Inf,y=Inf,aes(color=NULL,label=labl))+
	geom_text(show.legend=F,data=nlabel,
			hjust=0,vjust=1,x= -Inf,y=Inf,aes(color=NULL,label=labl))
	# scale_y_continuous(paste0(''))+
	ggtitle(paste0('Dwell Time Dist'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))


#now plot

repsumcodondata<-codondata%>%
	group_by(time,codon)%>%summarise_at(vars(abundance,availability,dwell_time,aacor_dwell_time),mean)
plotfile<- here(paste0('plots/','bDT_av_vs_dt','.pdf'))
pdf(h=15,w=5,plotfile)
col1=sym('aacor_dwell_time')
col2=sym('availability')
corlabel = repsumcodondata%>%filter(is.finite(!!col1),is.finite(!!col2))%>%
		group_by(time)%>%
		summarise(tidy(cor.test(!!col1, !!col2)))
corlabel = corlabel%>%
	mutate(
		pformat=format(p.value,format='e',digits=4),
		pvalstring = ifelse(p.value > 0.001,round(p.value,4),pformat),
		labl=paste0('rho = ',round(estimate,3),'\n','pval = ',pvalstring))
nlabel=repsumcodondata%>%group_by(time)%>%summarise(labl=paste0('N=',n()))
repsumcodondata%>%
	ggplot(.,aes(y=availability,x=aacor_dwell_time))+
	geom_point(alpha=I(0.5))+
	geom_smooth(method='lm')+
	facet_grid(time~.)+
	scale_fill_manual(name='stage',values=stagecols)+
	scale_x_continuous(paste0('Broad Dwell Time'))+
		geom_text(show.legend=F,data=corlabel,
			hjust=1,vjust=1,x= Inf,y=Inf,aes(color=NULL,label=labl))+
	geom_text(show.legend=F,data=nlabel,
			hjust=0,vjust=1,x= -Inf,y=Inf,aes(color=NULL,label=labl))
	# scale_y_continuous(paste0(''))+
	ggtitle(paste0('Dwell Time Dist'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))


