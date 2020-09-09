


library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(GenomicAlignments)
require(txtplot)
require(Rsamtools)
require(rlang)


FLANKCODS<-15
cdsids2use <- names(expcdsgrl)
cds2use <- expcdsgrl
REF <- 'yeast_test/Yeast.saccer3.fa'
ref<-Rsamtools::FaFile(REF)
cdsseq <- cds2use%>%extractTranscriptSeqs(x=ref)
allcodons=getGeneticCode()
#get code for mitcondrial coding genes


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
	#
	codmatchonchr<-matchgr%>%mapFromTranscripts(cds2use)
	codmatchonchr%<>%subset(width==3)
	#
	codmatchwindows<-codmatchonchr%>%resize(width(.)+(2*(3*FLANKCODS)),'center')
	#
	# codmatchwindows%>%sample(1e3)%>%getSeq(x=ref)
	seqlengths(codmatchwindows) <- seqlengths(ref)[unique(seqnames(codmatchwindows))]
	codmatchwindows<-codmatchwindows[!is_out_of_bounds(codmatchwindows)]
	codmatchwindows
})



cdscodons%<>%setNames(names(allcodons))

codons2scan<-allcodons%>%names%>%
	# str_subset('TTC|GTC|CAC|AAC|ATG')%>%
	identity
codons2scan%<>%setNames(.,.)
#check this worked

stopifnot(all( c('A','T','G') == (cdscodons[['ATG']]%>%sample(1e3)%>%getSeq(x=ref)%>%as.matrix%>%apply(2,table)%>%simplify2array%>%head((1+FLANKCODS)*3)%>%tail(3)%>%unlist)%>%names))
bams2scan <- bamfile
# psitelist <- mclapply(bams2scan,get_genomic_psites,cds2use%>%unlist)
# psitelist%<>%setNames(bams2scan%>%basename%>%str_replace('\\.\\w+$',''))

widths <- (25:31)%>%setNames(.,paste0('rl',.))
getreadfps <- function(bam,mapqthresh=200){
	  riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=cds2use%>%unlist)
	  reads <- readGAlignments(bam,param=riboparam)
	  reads%<>%as('GRanges')
	  reads$length <- width(reads)
	  reads%<>%subset(length %in% widths)
	  reads %<>% resize(1,'start')
	  reads
}
fpsitelist <- mclapply(bams2scan,mymemoise(getreadfps))
fpsitelist%<>%setNames(bams2scan%>%basename%>%str_replace('\\.\\w+$',''))

#codons2scan='GTC'%>%setNames(.,.)

tmp<-sort(cds2use[c(1,666)])%>%lapply(.%>%head(2)%>%.[,NULL])%>%GRangesList
ext_grl(tmp,1,fixend='start')

cds2useext <- cds2use

#####
################################################################################
########
################################################################################

#Now run through the bam file and get a profile for each codon

if(!file.exists(here('data/yst_fpprofilemats_unproc.rds'))){
	
codon=codons2scan[[1]]
psites=fpsitelist[[1]]
length_i='27'
yst_fpprofilemats_unproc<-mclapply(mc.cores=1,fpsitelist,mymemoise(function(psites){

# profilemats_unproc<-lapply(psitelist,function(psites){
	#
	cdscounts <- countOverlaps(cds2useext,psites)
	zerocds <- names(cdscounts)[cdscounts < 32]
	#
	#
	cdscounts_densitys <- cdscounts%>%divide_by(sum(width(cds2useext)))%>%setNames(names(cds2useext))
	#
	stopifnot(all(names(cdscounts_densitys)==names(cds2useext)))

	mclapply(mc.cores=1,codons2scan,function(codon){
		#
		message(codon)
		#
		codmatchwindows <- cdscodons[[codon]]%>%head
		#
		lapply(widths,function(length_i){
			#
			pospsites<-psites%>%subset(length==length_i)%>%subset(strand=='+')
			negpsites<-psites%>%subset(length==length_i)%>%subset(strand=='-')
			#
			poswinds<-codmatchwindows%>%subset(strand=='+')#%>%resize(.,width(.-((FLANKCODS-3)*3)),'center')
			negwinds<-codmatchwindows%>%subset(strand=='-')#%>%resize(.,width(.-((FLANKCODS-3)*3)),'center')
	
			profiles <- 
			c(
				pospsites%>%coverage%>%.[poswinds],
				negpsites%>%coverage%>%.[negwinds],
				NULL
			)

			nposwinds<-poswinds%>%length
			profiles[-(1:nposwinds)]%<>%revElements

			#get the profiles as a matrix, of nbp cols, and nrow codons
			# profilemat<-profiles%>%setNames(paste0('w',seq_along(.)))%>%as('List')%>%as("DataFrame")%>%as("Matrix")

			# profilemat<-profiles%>%simplify2array%>%t
			windnames<-c(names(poswinds),names(negwinds))
			#now rather than covert to a matrix, we can do this to 
			profiles <- profiles[!windnames%in%zerocds]
			normdensities <- cdscounts_densitys[windnames[!windnames%in%zerocds]]
			
			np<-length(profiles)
			n_col<-length(profiles[[1]])

			logsumpvect<-sapply(1:n_col,function(i){
				# i=1
				profiles[as(rep(i,np),'IntegerList')]%>%sum%>%divide_by(normdensities)%>%mean(na.rm=TRUE)
				# profiles[as(rep(i,np),'IntegerList')]%>%sum%>%sum
			})
			if(any(!is.finite(logsumpvect)))

			invisible(logsumpvect)
			enframe(logsumpvect,'position','signal')
		})%>%bind_rows(.id='readlen')
	})%>%bind_rows(.id='codon')
}))%>%bind_rows(.id='sample')
}

	saveRDS(yst_fpprofilemats_unproc,here('data/yst_fpprofilemats_unproc.rds'))
}else{
	yst_fpprofilemats_unproc<-readRDS(here('data/yst_fpprofilemats_unproc.rds'))
}

if(all(yst_fpprofilemats_unproc$sample %in% 1:100)) yst_fpprofilemats_unproc$sample %<>% {names(fpsitelist)[as.numeric(.)]}

# codonprofiles <- yst_fpprofilemats_unproc%>%map_depth(4,enframe,'position','signal')%>%map_depth(3,1)%>%map_depth(2,bind_rows,.id='readlen')%>%map_depth(1,bind_rows,.id='codon')%>%bind_rows(.id='sample')
codonprofiles<-yst_fpprofilemats_unproc
codonprofiles%<>%filter(!codon %in% c('TAG','TAA','TGA'))
codonprofiles%<>%mutate(position = position - 1 - (FLANKCODS*3))
codonprofiles%<>%group_by(readlen)%>%filter(any(signal!=0))
# if(!str_detect(codonprofiles$readlen,'rl')) codonprofiles$readlen%<>%as.numeric(.)%>%names(widths)[.]
codonprofiles%<>%group_by(readlen,codon,sample)%>%
	mutate(occ_nonorm=signal)%>%
	mutate(occupancy = signal / median(signal,na.rm=T))
codonprofiles%<>%select(-signal)
stopifnot(codonprofiles$occupancy%>%is.finite%>%all)


codonoccs<-codonprofiles%>%
	inner_join(offsets%>%select(readlen=length,offset)%>%mutate(readlen=paste0('rl',readlen)),by='readlen')%>%
	filter(sample%>%str_detect('ribo'))%>%
	filter(position <= -offset-3, position >= -offset-5)%>%
	# filter(position <= -offset, position >= -offset-5)%>%
	# filter(position <= -offset-3, position >= -offset-5)%>%
	# filter(readlen%in%c('rl29'))%>%
	ungroup%>%
	mutate(time=str_extract(sample,'[^_]+'))%>%
	group_by(time,codon)%>%
	# group_slice(2)
	summarise(occupancy=mean(occ_nonorm ,na.rm=T))


wbsalmondf%>%
	left_join(codonoccs)

profvarpca <- codonprofiles%>%
	split(.,.$sample)%>%
	map_df(.id='sample',.%>%
		split(.,list(.$readlen))%>%
		map_df( .id='readlen',.%>%
			mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
			filter(position> -numreadlen+6,position < -6)%>%
			# filter(sample=='E13_ribo_1')%>%
			ungroup%>%
			select(-numreadlen,-occ_nonorm,-readlen,-sample)%>%
			group_by(codon)%>%
			spread(position,occupancy)%>%
			{set_rownames(.[,-1],.$codon)}%>%
			princomp%>%{.$loadings[,1]}%>%{./.[which.max(abs(.))]}%>%enframe('position','pca1')
			# identity
		)
	)
profvarpca%<>%select(sample,readlen,position,pca1)

offsetspca<-offsets%>%mutate(readlen=paste0('rl',length))%>%filter(comp=='nucl')


plotfile<-'plots/yeast_fppos_vs_codon_pcascore.pdf'

grDevices::pdf(plotfile,w=12,h=12)
profvarpca%>%slice_by(sample,c(1,2,3,4,5,6))%>%
ggplot(data=.,aes(y=pca1,x=as.numeric(position)))+geom_point()+
	facet_grid(readlen~sample)+
		geom_vline(data=offsetspca,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=offsetspca,aes(xintercept= -offset-5),color=I('green'),linetype=2)
dev.off()
normalizePath(plotfile)

codonocccompdf <- codonoccs%>%left_join(
		fread('ext_data/weinberg_etal_2016_S2.tsv')%>%rename('codon':=Codon)
	)%>%left_join(
		fread('yeast_test/weinberg_yeast_ribo_skout/codons.csv')%>%rename('codon':='codon')
	)
pdf<-grDevices::pdf

#now plot
scattertitle = codonocccompdf%>%{cor.test(use='complete',(.[['occupancy']]),.[['RiboDensity at A-site']])}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n ${.$conf.low} - ${.$conf.high}')}
plotfile<- here(paste0('plots/','metaplotocc_vs_wbRibodens','.pdf'))
pdf(plotfile)
	ggplot(data=codonocccompdf,aes(x=occupancy,y=`RiboDensity at A-site`))+
	geom_point()+
	scale_x_continuous(paste0('codon density - metaplotocc'))+
	scale_y_continuous(paste0('codon density - weinberg et al'))+
	ggtitle(paste0('GLM gene term vs Riboseq Density'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)

