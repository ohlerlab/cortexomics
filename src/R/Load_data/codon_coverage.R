
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

shiftmodel <- 'pipeline/seqshift_reads/data/P0_ribo_1/seqshiftmodel.rds'
if(!file.exists(shiftmodel)){
	message('no psite model found')
	psite_model <- NULL
}else{
	psite_model <- readRDS(shiftmodel)
}

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

bams <- here('pipeline/star/data/*/*.bam')%>%Sys.glob%>%
	str_subset(neg=TRUE,'transcript')%>%
	str_subset('_ribo_')

bamfile=bams[1]


cdscodons <- 	lapply(seq_along(allcodons),function(i){
	
	codon=names(allcodons)[[i]]
	message(codon)
	codmatches<-vmatchPattern(pattern=codon,cdsseq)

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
# profilemats_unproc<-lapply(psitelist,function(psites){
fprofilemats_unproc<-mclapply(mc.cores=1,fpsitelist,mymemoise(function(psites){
	#
	cdscounts <- countOverlaps(cds2useext,psites)
	zerocds <- names(cdscounts)[cdscounts < 32]
	#
	stopifnot(all(names(cdscounts_densitys)==names(cds2use)))
	#
	cdscounts_densitys <- cdscounts%>%divide_by(sum(width(cds2use)))%>%setNames(names(cds2use))
	#

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

################################################################################
########PRocess these
################################################################################
	
stopifnot(all(is.finite(fprofilemats_unproc$signal)))

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
codonprofiles%<>%group_by(readlen,codon,sample)%>%mutate(signal = signal / median(signal))
stopifnot(codonprofiles$signal%>%is.finite%>%all)

################################################################################
########Plot occupancies at codon level
################################################################################
	
#plotting variance amongst codons at each point.
codonprofiles%>%
	ungroup%>%
	filter(sample==unique(sample)[1])%>%
	group_by(sample,readlen,position)%>%
	# filter(signal==0)%>%.$codon%>%unique
	# filter(signal!=0)%>%
	group_by(sample,position,readlen)%>%
	# filter(position==1)%>%
	filter(!is.nan(signal))%>%
	summarise(sdsig=sd(signal,na.rm=T)/mean(signal,na.rm=T))%>%
	group_by(sample,readlen)%>%
	filter(readlen=='rl27')%>%
	identity%>%
	# .$position%>%unique
	filter(between(position,-24,0))%>%
	arrange(position)%>%
	{txtplot(x=.$position,y=.$sdsig)}

	{txtplot(x=.$position,y=.$sdsig)}



codonproftppos<-codonprofiles%>%distinct(readlen,codon)%>%mutate(tppos = 1-as.numeric(str_replace(readlen,'rl','')))
codonproftppos<-codonprofiles%>%distinct(readlen,codon)%>%mutate(tppos = 1-as.numeric(str_replace(readlen,'rl','')))
offsets<- psite_model$offsets%>%filter(compartment=='nucl')%>%mutate(readlen=paste0('rl',length))%>%select(readlen,offset)

codonprofiles$codon%>%unique
codonprofiles$codon%>%str_subset('GTC')%>%unique
codons2scan%>%str_subset('GTC')
fprofilemats_unproc$codon%>%str_subset('GTC')%>%unique
fprofilemats_unproc$codon%>%n_distinct
codons2scan%>%n_distinct


samplestage <- unique(codonprofiles$sample)%>%{setNames(stageconv[str_replace(.,'_.*?_.*?$','')],.)}
pdf('tmp.pdf',w=24,h=14)
codonprofiles%>%
	# filter(codon%>%str_detect(c('GTC|AAC|ATG')))%>%
	filter(codon%>%str_detect(c('Glu-TTC|GTC|Val-CAC|AAC|ATG')))%>%
	ungroup%>%
	mutate(codon = as_factor(codon))%>%
	# filter(sample%>%str_detect(c('ribo')))%>%
	{print(ggplot(.,aes(position,signal,group=sample,
		color=samplestage[sample]))+
		scale_color_manual(values=displaystagecols)+
		scale_x_continuous(minor_breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=3),breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=9) )+
		facet_grid(codon~readlen)+
		geom_line(aes())+
		geom_vline(xintercept=0,linetype=2)+
		geom_vline(data=filter(codonproftppos,codon%in%.$codon),aes(xintercept=tppos),linetype=2)+
		coord_cartesian(xlim=c(-39,12))+
		geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=offsets,aes(xintercept= -offset-2-6),color=I('green'),linetype=2)+
		# geom_rect(data=offsets,color=I('black'),alpha = I(0.1),aes(x=NULL,y=NULL,xmin= -8-offset, xmax = -offset, ymin = 0 ,ymax = Inf))+
		theme_bw())}
dev.off()
normalizePath('tmp.pdf')


pdf('tmp.pdf',w=24,h=14)
codonprofiles%>%
	inner_join(offsets)%>%
	mutate(position = position+offset)%>%
	group_by(sample,codon,position)%>%summarise(signal=sum(signal))%>%
	# filter(codon%>%str_detect(c('GTC|AAC|ATG')))%>%
	# filter(codon%>%str_detect(c('Glu-TTC|GTC|Val-CAC|AAC|ATG')))%>%
	# filter(sample%>%str_detect(c('ribo')))%>%
	{print(ggplot(.,aes(position,signal,group=sample,
		color=samplestage[sample]))+
		geom_rect(color=I('black'),alpha = I(0.3),aes(xmin= -6, xmax = 3, ymin = 0 ,ymax = Inf))+
		scale_color_manual(values=displaystagecols)+
		scale_x_continuous(minor_breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=3),breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=9) )+
		facet_grid(codon~.)+
		geom_line(aes())+
		geom_vline(xintercept=0,linetype=2)+
		# geom_vline(data=codonproftppos,aes(xintercept=tppos),linetype=2)+
		coord_cartesian(xlim=c(-39,12))+
		# geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		theme_bw())}
dev.off()
normalizePath('tmp.pdf')

profilemats_unproc%>%saveRDS('data/fprofilemats_unproc.rds')

profilemats <- profilemats_unproc%>%setNames(.,bams[1:length(.)])

cdswidths<-cds%>%split(.,.$protein_id)%>%width%>%sum
longestcdsids<-enframe(cdswidths,'protein_id','length')%>%left_join(cds%>%mcols%>%as_tibble%>%distinct(protein_id,gene_id))%>%group_by(gene_id)%>%slice(which.max(length))%>%.$protein_id
cds_noov<-cds%>%subset(protein_id%in%longestcdsids)
cds_noov%<>%split(.,.$protein_id)
# cds_noov<-cds_noov%>%subset(strand=='+')%>%split(.,.$protein_id)%>%head(10)%>%lapply(head,1)%>%GRangesList%>%unlist%>%resize(3)

oligonucleotideFrequency(getSeq(FaFile(REF),cds_noov),3,step=3)


################################################################################
########Codon optimality scores based on codon frequencies, gene expression levels
################################################################################
	
codonfreqs<-oligonucleotideFrequency(extractTranscriptSeqs(FaFile(REF),cds_noov),3,step=3)
rownames(codonfreqs)<-names(cds_noov)
overallcodonfreqs<-codonfreqs%>%colSums
bestcds<- cds%>%subset(protein_id%in%best_protein_ids)%>%split(.,.$protein_id)
usedcodonfreqs<-oligonucleotideFrequency(extractTranscriptSeqs(FaFile(REF),bestcds),3,step=3)%>%set_rownames(names(bestcds))
codon_is_optimal <- overallcodonfreqs%>%enframe('codon','freq')%>%
	mutate(AA=as.character(translate(DNAStringSet(codon))))%>%
	group_by(AA)%>%mutate(optimal = freq==max(freq))%>%
	{setNames(.$optimal,.$codon)}

optimalcodons <- names(codon_is_optimal)[codon_is_optimal]

pc_optcodons <- ((usedcodonfreqs[,optimalcodons]%>%rowSums)/(usedcodonfreqs%>%rowSums))%>%setNames(rownames(usedcodonfreqs))%>%enframe('protein_id','pc_opt')

  exprs(countexprdata)

usedcodonfreqs[1:3,]%>%as.data.frame%>%rownames_to_column("protein_id")%>%gather(codon,protein_id)%>%
	left_join(tRNA_occ_df%>%distinct(codon,time,tRNA_expr=signal,occupancy))%>%
	group_by(time,protein_id)%>%
	summarise(occupancy=mean(occupancy),tRNA_expr=mean(tRNA_expr))

tRNA_occ_df
codons4occ <- tRNA_occ_df$codon%>%unique%>%setdiff(c('TAG','TAA','TGA'))

cdsexprvals <- (2^bestmscountvoom$E[,T])%>%{rownames(.)%<>%str_replace('_\\d+$','');.}

#length norm them
cdsexprvals <- cdsexprvals / fData(countexprdata)$length[match(rownames(cdsexprvals),fData(countexprdata)$protein_id)]


cdsexprdf<- cdsexprvals[,1:20] %>% as.data.frame%>% rownames_to_column('protein_id')%>%gather(sample,signal,-protein_id)%>%
	separate(sample,c('time','assay','rep'))%>%group_by(protein_id,time,assay)%>%summarise(signal=mean(signal))%>%
	filter(assay=='ribo')%>%
	spread(time,signal)%>%
	arrange(match(protein_id,rownames(usedcodonfreqs)))

weighted_codon_usage <- lapply(times,function(itime)(usedcodonfreqs[,]*cdsexprdf[[itime]])%>%colSums%>%{./sum(.)}%>%enframe('codon','weightedusage'))%>%bind_rows(.id='time')

tRNA_abchange <- tRNA_occ_df%>%
	group_by(codon)%>%
	filter(!is.na(signal))%>%
	nest%>%
	mutate(abundancechange = map_dbl(data,~{lm(data=.,signal ~ seq_along(time))$coef[2]}))

usage_v_abundance_df <- weighted_codon_usage%>%filter(time=='E13')%>%ungroup%>%select(-time)%>%left_join(tRNA_abchange%>%select(codon,abundancechange))%>%
	filter(!is.na(abundancechange),!is.na(weightedusage))

wus_tab_cors<-usage_v_abundance_df%>%
	filter(!codon%in%c('TAG','TAA','TGA'))%>%
	nest%>%mutate(cor = map_dbl(data,~ cor(.$abundancechange,.$weightedusage)))

wus_tab_cors$pval <- usage_v_abundance_df%>%	filter(!codon%in%c('TAG','TAA','TGA'))%>%{cor.test(.$abundancechange,.$weightedusage)}%>%.$p.value

pdf('tmp.pdf')
usage_v_abundance_df%>%
	filter(!codon%in%c('TAG','TAA','TGA'))%>%
	{
		print(
			qplot(data=.,x=.$abundancechange,y=.$weightedusage,label=codon,geom='blank')+
				scale_x_continuous('tRNA Abundance Signal Slope E13 - P0')+
				scale_y_continuous('RiboExpression Weighted E13 ')+
				geom_smooth(method='lm')+
				# facet_grid(tRNA_time~.)+
				geom_text()+
				geom_text(data=wus_tab_cors,aes(x=0,y=0,label=paste0('r = ',round(cor,2))))+
				theme_bw()
			)
	}
dev.off()
normalizePath('tmp.pdf')






weighted_codon_usage%>%filter(time=='E13')%>%ungroup%>%select(-time)%>%left_join(tRNA_occ_df%>%filter(time=='P0')%>%ungroup%>%select(-time))%>%{txtplot(.$signal,.$weighted_codon_usage)}
tRNA_occ_df%>%filter(time=='P0')%>%.$signal

timeoccscores <- map_df(.id='time',times,function(itime){
	timeoccvect <- tRNA_occ_df%>%filter(time==itime)%>%mutate(occupancy=occupancy-median(na.rm=T,occupancy))%>%{setNames(.$occupancy,.$codon)[codons4occ]}
	(t(usedcodonfreqs[,codons4occ])*(timeoccvect))%>%colMeans%>%enframe('protein_id','occ_score')
})

##Score for cds at particular time points
times<-tRNA_occ_df$time%>%unique%>%setNames(.,.)
itime <- times[1]



timeoccscores <- map_df(.id='time',times,function(itime){
	timeoccvect <- tRNA_occ_df%>%filter(time==itime)%>%mutate(occupancy=occupancy-median(na.rm=T,occupancy))%>%{setNames(.$occupancy,.$codon)[codons4occ]}
	(t(usedcodonfreqs[,codons4occ])*(timeoccvect))%>%colMeans%>%enframe('protein_id','occ_score')
})

timeSigscores <- map_df(.id='time',times,function(itime){
	timeoccvect <- tRNA_occ_df%>%filter(time==itime)%>%mutate(signal=signal-median(na.rm=T,signal))%>%{setNames(.$signal,.$codon)[codons4occ]}
	(t(usedcodonfreqs[,codons4occ])*(timeoccvect))%>%colMeans(na.rm=T)%>%enframe('protein_id','tRNA_expr_score')
})

score_techange_df<-timeoccscores%>%left_join(timeSigscores)%>%left_join(ms_id2protein_id%>%distinct(protein_id,gene_name))%>%
	left_join(allTEchangedf)

occchange_vs_te_df <- score_techange_df%>%
	group_by(up,down,protein_id)%>%
	# filter(time!='P0')%>%
	# group_slice(1:2)%>%
	nest%>%
	mutate(occhange = map_dbl(data,~{lm(data=.,occ_score ~ seq_along(time))$coef[2]}))%>%
	mutate(tps = map_dbl(data,nrow))%>%
	select(-data)

tRNAchange_vs_te_df <- score_techange_df%>%
	group_by(up,down,protein_id)%>%
	# filter(time!='P0')%>%
	# group_slice(1:2)%>%
	nest%>%
	mutate(tRNA_score_change = map_dbl(data,~{lm(data=.,tRNA_expr_score ~ seq_along(time))$coef[2]}))%>%
	mutate(tps = map_dbl(data,nrow))%>%
	select(-data)

pdf('plots/figures/figure2/te_change_vs_occscorechange.pdf')
occchange_vs_te_df%>%
	filter(!(up&down))%>%
	mutate(TEchange_class = case_when(
		up==1 ~ 'up',
		down==1 ~ 'down',
		TRUE ~ 'neither'
	))%>%
	ggplot(.,aes(x=occhange,color=TEchange_class))+geom_density(alpha=I(0.01))+theme_bw()+scale_x_continuous('Predicted Elongation Rate Change - Riboseq Occ')
dev.off()
normalizePath('plots/figures/figure2/te_change_vs_occscorechange.pdf')


occchange_vs_te_df%>%filter(!down)%>%{split(.$occhange,.$up)}%>%{t.test(.[[1]],.[[2]])}
occchange_vs_te_df%>%filter(!up)%>%{split(.$occhange,.$down)}%>%{t.test(.[[1]],.[[2]])}


pdf('plots/figures/figure2/te_change_vs_tRNAscorechange.pdf')
tRNAchange_vs_te_df%>%
	filter(!(up&down))%>%
	# filter(time!='P0')%>%
	mutate(TEchange_class = case_when(
		up==1 ~ 'up',
		down==1 ~ 'down',
		TRUE ~ 'neither'
	))%>%
	ggplot(.,aes(x=tRNA_score_change,color=TEchange_class))+geom_density(alpha=I(0.01))+theme_bw()+scale_x_continuous('Predicted Change in tRNA Abundance Score')
dev.off()
normalizePath('plots/figures/figure2/te_change_vs_tRNAscorechange.pdf')

tRNAchange_vs_te_df%>%filter(!down)%>%{split(.$tRNA_score_change,.$up)}%>%{t.test(.[[1]],.[[2]])}
tRNAchange_vs_te_df%>%filter(!up)%>%{split(.$tRNA_score_change,.$down)}%>%{t.test(.[[1]],.[[2]])}

###So is this the same genes???

tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%{txtplot(.$tRNA_score_change,.$occhange)}
tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%{cor.test(.$tRNA_score_change,.$occhange)}



#occ scores of things that change TE
pdf('tmp.pdf')
score_techange_df%>%
	filter(time=='E13')%>%
	filter(!(up&down))%>%
	mutate(TEchange_class = case_when(
		up==1 ~ 'up',
		down==1 ~ 'down',
		TRUE ~ 'neither'
	))%>%
	ggplot(.,aes(x=log(occ_score),color=TEchange_class))+geom_density(alpha=I(0.01))+theme_bw()
dev.off()
normalizePath('tmp.pdf')

allTEchangedf


###Do optimized looking CDS tend to have higher TE?
abste_opt_df<-bestmscountebayes$coef[,'TE']%>%enframe('uprotein_id','TE')%>%mutate(protein_id=uprotein_id%>%str_replace('_\\d+$',''))%>%
	left_join(pc_optcodons)


cor.test(abste_opt_df$TE,pc_optcodons$pc_opt)
txtplot(abste_opt_df$TE,pc_optcodons$pc_opt,width=100,pch='.')
#No.

#Do 


#Are occupancy and global frequency correlated?
tRNA_occ_df%>%
	left_join(enframe(overallcodonfreqs,'codon','freq'))%>%
	split(.,.$time)%>%
	map(~ cor.test(.$freq,.$occupancy))

#Are 
tRNA_occ_df%>%
	left_join(enframe(codon_is_optimal,'codon','optimality'))%>%
	mutate(AA=as.character(translate(DNAStringSet(codon))))%>%
	group_by(AA)%>%
	filter(n()>1)%>%
	mutate(least_occ=occupancy ==min(occupancy))%>%
	split(.,.$time)%>%
	map(~ identity(table(.$least_occ,.$optimality)))


################################################################################
########Now combine the tRNA and codon occupancy info
################################################################################ head(allcodsigmean_isomerge)


glutRNAabund<-allcodsigmean_isomerge%>%filter(codon%>%str_detect('CTC'))%>%.$signal
valtRNAabund<-allcodsigmean_isomerge%>%filter(codon%>%str_detect('AAC'))%>%.$signal
stopifnot(valtRNAabund<glutRNAabund)

allcodsigmean_isomerge<-allcodsigmean%>%
	filter(sample%>%str_detect('Total'))%>%
	# group_by(time,sample,codon)%>%
	# summarise(signal = -log2(sum(2^(-signal))))
	group_by(time,sample,codon,decoder)%>%
	nest%>%
	mutate(medsig = map_dbl(data,~median(.$signal)))%>%
	group_by(codon)%>%
	slice(which.max(medsig))%>%
	unnest

#get occupancies
codonoccs<-codonprofiles%>%inner_join(offsets)%>%
	filter(position <= -offset, position >= -offset-8)%>%
	# filter(position <= -offset, position >= -offset-5)%>%
	# filter(position <= -offset-3, position >= -offset-5)%>%
	filter(readlen%in%c('rl29'))%>%
	ungroup%>%
	mutate(time=str_extract(sample,'[^_]+'))%>%
	group_by(time,codon)%>%
	summarise(occupancy=mean(signal,na.rm=T))

#Now merge with the tRNA data
tRNA_occ_df<-allcodsigmean_isomerge%>%
	mutate(codon=str_replace(codon,'\\w+\\-',''))%>%
	left_join(codonoccs)

normtRNA_occ_df<-tRNA_occ_df%>%
	group_by(codon)%>%
	# mutate(signal = signal-mean(signal),occupancy=occupancy-mean(occupancy))
	mutate(signal = signal-mean(signal),occupancy=occupancy-occupancy[time=="E13"])

normtRNA_occ_df%>%
	# filter(! sample %>% str_detect('13')) %>% 
	# filter(! sample %>% str_detect('P0')) %>% 
	{cor.test(.$occupancy,.$signal,use='complete',method='pearson')}


pdf('tmp.pdf',w=24,h=14)
tRNA_occ_df%>%
	# filter(codon%>%str_detect(c('TTC|GTC|CAC|AAC|ATG')))%>%
	# filter(time%>%str_detect("E13"))%>%
	{
	#
	labelpos = group_by(.,codon)%>%summarise(occupancy=mean(occupancy),signal=mean(signal))
	#
	# filter(codon%>%str_detect("ATG"))%>%
	ggplot(.,aes(x=signal,y=occupancy,color=time,group=codon))+
	geom_point(size=8)+
	geom_line(color=I('grey'))+
	facet_wrap(nrow=3,time ~ . )+
	scale_color_manual(values=stagecols)+
	geom_text(show.legend=F,data=labelpos,aes(label=codon,color=NULL),size=9)+
	scale_x_continuous(name='Summed tRNA Expression')
}
dev.off()
normalizePath('tmp.pdf')

tRNA_occ_df %>% {cor.test(.$occupancy,.$signal,use='complete')}
tRNA_occ_df%>%filter(time=='P0') %>% {cor.test(.$occupancy,.$signal,use='complete')}

tRNA_occ_df%>%
	filter(! codon %>% str_detect('AAC')) %>% 
	{cor.test(.$occupancy,.$signal,use='complete')}

normtRNA_occ_df%>%filter(signal==0)



##Plot of change over time
pdf('tmp.pdf',w=24,h=14)
normtRNA_occ_df%>%
	
	# filter(! sample %>% str_detect('E13')) %>% 
	# filter(! sample %>% str_detect('P0')) %>% 
	# filter(codon%>%str_detect(c('TTC|GTC|CAC|AAC|ATG')))%>%
	# filter(time%>%str_detect("E13"))%>%
	{
	#
	labelpos = group_by(.,codon)%>%summarise(occupancy=mean(occupancy),signal=mean(signal))
	#
	# filter(codon%>%str_detect("ATG"))%>%
	ggplot(.,aes(x=signal,y=occupancy,color=time,group=codon))+
	geom_point(size=8)+
	# geom_line(color=I('grey'))+
	# facet_wrap(nrow=3,time ~ . )+
	scale_color_manual(values=stagecols)+
	# geom_text(show.legend=F,data=labelpos,aes(label=codon,color=NULL),size=9)+
	scale_x_continuous(name='Summed tRNA Expression - Median')+
	scale_y_continuous(name='Summed Occupancy - Median')+
	theme_bw()
}
dev.off()
normalizePath('tmp.pdf')

	{cor.test(.$occupancy,.$signal,use='complete')}

# codonstbl<-read.table('https://raw.githubusercontent.com/zhanxw/anno/master/codon.txt')
# codonstbl%>%write_tsv('ext_data/codons.txt')
codonstbl<-read_tsv('ext_data/codons.txt')
#segregate the 
codonnames <- paste0(codonstbl[[2]],'-',codonstbl[[1]])[match(names(codons2scan),codonstbl[[1]])]

#codoncovdf<-
profilemats%<>%map_df(.id='bam',.%>%map_df(.id='codon',.%>%enframe('position','signal')))

bamsamples<-bams%>%basename%>%str_remove('.bam')

sampleparams<-'src/sample_parameter.csv'%>%fread

profilemats%>%head
profilemats$sample = profilemats$bam%>%match(bams)%>%bamsamples[.]
profilemats$time<-sampleparams$time[match(profilemats$sample,sampleparams$sample_id)]



profilemats$codon%<>%as.numeric
profilemats$codon<-profilemats$codon%>%codonnames[.]

#summarise to get rid of the periodicity signal
codonprofmat<-profilemats%>%mutate(codonpos=ceiling(position/3))%>%group_by(codon,position,bam)

codonprofmat$codonpos%<>%subtract(FLANKCODS+1)
codonprofmat%<>%ungroup


require(plotly)

samplestage[codonprofmat$sample]%>%unique

highcod = 'AAC'
lowcod = 'GTC'

pdf('tmp.pdf')
codonprofmat%>%
	filter(codon%>%str_detect(c('GTC|AAC|ATG')))%>%
	# filter(codon%>%str_detect(c('Glu-TTC|GTC|Val-CAC|AAC|ATG')))%>%
	filter(sample%>%str_detect(c('ribo')))%>%
	{print(ggplot(.,aes(codonpos,signal,group=bam,
		color=samplestage[sample]))+
		scale_color_manual(values=displaystagecols)+
		facet_grid(codon~.)+
		geom_line(aes())+
		theme_bw())}
dev.off()
normalizePath('tmp.pdf')




save.image('data/codon_coverage.Rdata')


#' Convert from Rle to one column matrix
#'
Q
setAs("Rle", "Matrix", function(from) {
    rv <- runValue(from)
    nz <- rv != 0
    i <- which(as.vector(from !=0))
    x <- rep(rv[nz], runLength(from)[nz])

	length(i)    

    sparseMatrix(i= i, p=c(0L, length(x)), x=x,
                 dims=c(length(from), 1))
})

#' Convert from DataFrame of Rle to sparse Matrix
#'
setAs("DataFrame", "Matrix", function(from) {
  mat = do.call(cbind, lapply(from, as, "Matrix"))
  colnames(mat) <- colnames(from)
  rownames(mat) <- rownames(from)
  mat
})



