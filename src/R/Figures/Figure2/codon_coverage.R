
conflict_prefer("which", "Matrix")
conflict_prefer("rowSums", "BiocGenerics")
conflict_prefer("colMeans", "BiocGenerics")
conflict_prefer("setdiff", "dplyr")

ms_id2protein_id
best_uprotein_ids

library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(GenomicAlignments)
require(txtplot)
require(Rsamtools)
#coverageplots
offsets <- read_tsv(here('ext_data/offsets_manual.tsv'))

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


cdsids2use <- ms_id2protein_id%>%filter(uprotein_id %in% best_uprotein_ids)%>%.$protein_id

cds2use <- cds%>%subset(protein_id %in% cdsids2use)%>%split(.,.$protein_id)

REF <- 'pipeline/my_GRCm38.p5.genome.chr_scaff.fa'
ref<-Rsamtools::FaFile(REF)

cdsseq <- cds2use%>%extractTranscriptSeqs(x=ref)


allcodons=getGeneticCode()
#get code for mitcondrial coding genes

bamfiles<-Sys.glob(here::here('pipeline/star/data/*/*_*bam'))
library(here)
bams <- here('pipeline/star/data/*/*.bam')%>%Sys.glob%>%
	str_subset(neg=TRUE,'transcript')%>%
	str_subset('_ribo_|Poly')

cdscodons <- 	lapply(seq_along(allcodons),function(i){
	
	codon=names(allcodons)[[i]]
	message(codon)
	codmatches<-vmatchPattern(pattern=codon,cdsseq%>%subseq(4,Inf))#exclude the start ccodon

	nmatches = 	codmatches%>%elementNROWS 

	strands<-strand(cds2use[names(codmatches)])[as(rep(1,length(cds2use)),'IntegerList')]%>%unlist
	strands <- rep(as.vector(strands),nmatches)
	matchgr<-codmatches%>%unlist%>%GRanges(names(.),.,strand=strands)

	matchgr%<>%subset(start %%3 == 1)

	codmatchonchr<-matchgr%>%mapFromTranscripts(cds2use)
	codmatchonchr%<>%subset(width==3)

	codmatchwindows<-codmatchonchr%>%resize(width(.)+(2*(3*FLANKCODS)),'center')

	codmatchwindows

})

cdscodons%<>%setNames(names(allcodons))

codons2scan<-allcodons%>%names%>%
	# str_subset('TTC|GTC|CAC|AAC|ATG')%>%
	identity
codons2scan%<>%setNames(.,.)
#check this worked

stopifnot(all( c('A','T','G') == (cdscodons[['ATG']]%>%getSeq(x=ref)%>%as.matrix%>%apply(2,table)%>%simplify2array%>%head((1+FLANKCODS)*3)%>%tail(3)%>%unlist)%>%names))

bams2scan<-bams%>%
	# .[c(1:2,length(.)-1,length(.))] %>%
	identity
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
conflict_prefer('lag','dplyr')

tmp<-sort(cds2use[c(1,666)])%>%lapply(.%>%head(2)%>%.[,NULL])%>%GRangesList
ext_grl(tmp,1,fixend='start')

cds2useext <- ext_grl(sort(cds2use),45,fixend='start')%>%ext_grl(45,fixend='end')

#####
################################################################################
########
################################################################################

#Now run through the bam file and get a profile for each codon
{

fprofilemats_unproc<-mclapply(mc.cores=1,fpsitelist,mymemoise(function(psites){

# profilemats_unproc<-lapply(psitelist,function(psites){
	#
	cdscounts <- countOverlaps(cds2useext,psites)
	zerocds <- names(cdscounts)[cdscounts < 32]
	#
	#
	cdscounts_densitys <- cdscounts%>%divide_by(sum(width(cds2useext)))%>%setNames(names(cds2useext))
	#
	stopifnot(all(names(cdscounts_densitys)==names(cds2useext)))

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
			poswinds<-codmatchwindows%>%subset(strand=='+')#%>%resize(.,width(.-((FLANKCODS-3)*3)),'center')
			negwinds<-codmatchwindows%>%subset(strand=='-')#%>%resize(.,width(.-((FLANKCODS-3)*3)),'center')
			#
			# pospsites%>%coverage%>%.[poswinds]%>%length%>%divide_by(1e3)
			# negpsites%>%coverage%>%.[negwinds]%>%length%>%divide_by(1e3)
			#
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
			# profilemat%<>%sweep(1,STATS=cdscounts_densitys[windnames],FUN='/')
			# profilemat%<>%as.matrix
			# profilemat%<>%t
			# profilemat <- profilemat[	!profilemat[,1]%>%map_lgl(is.nan),]
			# profilemat <- profilemat[	profilemat%>%apply(1,function(x) all((is.finite(x))) ) ,]
			
		
			# sapply(i:n_col,function(i) profiles[as(rep(i,np),'IntegerList')]%>%sum%>%add(1)%>%log)
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

			# logsumpvect%T>%txtplot(xlab=as.character((-3*FLANKCODS):(2+3*FLANKCODS)),width=100)

			# rustsumvect<-sapply(1:n_col,function(i){
			# 	profiles[as(rep(i,np),'IntegerList')]%>%sum%>%{.>0}%>%sum
			# })%T>%txtplot(xlab=as.character((-3*FLANKCODS):(2+3*FLANKCODS)),width=100)

			# sumvect<-sapply(1:n_col,function(i){
			# 	profiles[as(rep(i,np),'IntegerList')]%>%sum%>%sum
			# })%T>%txtplot(xlab=as.character((-3*FLANKCODS):(2+3*FLANKCODS)),width=100)
					# invisible(list())
			invisible(logsumpvect)
			enframe(logsumpvect,'position','signal')
		})%>%bind_rows(.id='readlen')
	})%>%bind_rows(.id='codon')
}))%>%bind_rows(.id='sample')
}

fprofilemats_unproc%>%slice_by(sample,1)%>%group_by(codon)%>%summarise(m=mean(signal))%>%.$m%>%txtdensity

# fprofilemats_unproc2$signal%>%txtdensity
################################################################################
########PRocess these
################################################################################



if(all(fprofilemats_unproc$sample %in% 1:100)) fprofilemats_unproc$sample %<>% {names(fpsitelist)[as.numeric(.)]}

fprofilemats_unproc$codon%>%n_distinct
codons2scan%>%n_distinct
fprofilemats_unproc$position%>%unique

# codonprofiles <- fprofilemats_unproc%>%map_depth(4,enframe,'position','signal')%>%map_depth(3,1)%>%map_depth(2,bind_rows,.id='readlen')%>%map_depth(1,bind_rows,.id='codon')%>%bind_rows(.id='sample')
codonprofiles<-fprofilemats_unproc
codonprofiles%<>%filter(!codon %in% c('TAG','TAA','TGA'))
codonprofiles%<>%mutate(position = position - 1 - (FLANKCODS*3))
codonprofiles%<>%group_by(readlen)%>%filter(any(signal!=0))
if(!str_detect(codonprofiles$readlen,'rl')) codonprofiles$readlen%<>%as.numeric(.)%>%names(widths)[.]
codonprofiles%<>%group_by(readlen,codon,sample)%>%mutate(occ_nonorm=signal)%>%mutate(signal = signal / median(signal))
stopifnot(codonprofiles$signal%>%is.finite%>%all)

codonproftppos<-codonprofiles%>%distinct(readlen,codon)%>%mutate(tppos = 1-as.numeric(str_replace(readlen,'rl','')))
codonproftppos<-codonprofiles%>%distinct(readlen,codon)%>%mutate(tppos = 1-as.numeric(str_replace(readlen,'rl','')))
# offsets<- psite_model$offsets%>%filter(compartment=='nucl')%>%mutate(readlen=paste0('rl',length))%>%select(readlen,offset)

# offsets%>%write_tsv('ext_data/offsets_manual.tsv')