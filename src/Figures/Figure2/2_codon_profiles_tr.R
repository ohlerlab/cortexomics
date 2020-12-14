
################################################################################
########Do with deepshape etc instead
################################################################################
library(GenomicFeatures)
base::source(here::here('src/R/Rprofile.R'))
if(!exists('iso_tx_countdata')) load('data/1_integrate_countdata.R')
if(!exists('cdsgrl')) base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
# base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
trid2gid = cds%>%mcols%>%as.data.frame%>%select(transcript_id,gene_id)%>%{safe_hashmap(.[[1]],.[[2]])}


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

ribocovtrs = read_tsv('pipeline/scikitribotrs.txt',col_names=F)[[1]]

cds2use <- cdsgrl[ribocovtrs]

cdsseq <- cds2use%>%{GenomicFeatures::extractTranscriptSeqs(.,x=fafile)}
allcodons=getGeneticCode()

highcountcovtrs = ribocovtrs[trid2gid[[ribocovtrs]]%in%highcountgenes]
trlens = exonsgrl%>%width%>%sum

bamtbls = Sys.glob('pipeline/deepshapebamdata/*.bam.reformat')
names(bamtbls) = bamtbls%>%str_extract('[^/]*?(?=.bam.ref)')
bamtbl=bamtbls[1]

offsets <- read_tsv('ext_data/offsets_manual.tsv')

if(!file.exists(here('data/fpcovlist.rds'))){
	fpcovlist = bamtbls%>%mclapply(mc.cores=4,function(bamtbl){
		bamtbl%>%
			str_interp('grep -e $')%>%
			fread(select=c(2,6,7))%>%
			set_colnames(c('transcript_id','start','readlen'))%>%
			mutate(transcript_id=trimids(transcript_id))%>%
			filter(transcript_id%in%highcountcovtrs)%>%
			filter(between(readlen,25,31))%>%
			{GRanges(.$transcript_id,IRanges(.$start,w=1),readlen=.$readlen)}%>%
			{seqlevels(.) = highcountcovtrs ;.}%>%
			{seqlengths(.) = trlens[highcountcovtrs] ;.}%>%
			split(.,.$readlen)%>%
			lapply(coverage)
	})
	saveRDS(fpcovlist,here('data/fpcovlist.rds'))
}else{
	fpcovlist<-readRDS(here('data/fpcovlist.rds'))
}

#GRanges()%>%{seqlevels(.)=names(trlens[highcountcovtrs]);seqlengths(.)=trlens[highcountcovtrs];.}%>%coverage

cdsstarts = cdsgrl[highcountcovtrs]%>%sort_grl_st%>%resize_grl(1)%>%unlist%>%
	pmapToTranscripts(exonsgrl[names(.)]%>%sort_grl_st)%>%
	{setNames(start(.),as.character(seqnames(.)))}
cdsends = cdsgrl[highcountcovtrs]%>%sort_grl_st%>%resize_grl(1,'end')%>%unlist%>%
	pmapToTranscripts(exonsgrl[names(.)]%>%sort_grl_st)%>%
	{setNames(start(.),as.character(seqnames(.)))}
trcds = GRanges(names(cdsstarts),IRanges(cdsstarts,cdsends))

toptrs = iso_tx_countdata$abundance%>%
  as.data.frame%>%rownames_to_column('transcript_id')%>%
  pivot_longer(-transcript_id,names_to='dataset',values_to='TPM')%>%
  separate(dataset,c('time','assay','rep'))%>%
  group_by(transcript_id)%>%
  filter(assay=='ribo')%>%
  filter(transcript_id%>%is_in(ribocovtrs))%>%
  summarise(TPM=mean(TPM))%>%
  arrange(desc(TPM))
toptrs=toptrs%>%head(5e3)%>%.$transcript_id
exonseq = exonsgrl[toptrs]%>%extractTranscriptSeqs(x=fafile)

i=1
if(!file.exists(here('data/codmatchwindowlist.rds'))){
	codmatchwindowlist <- lapply(seq_along(allcodons)%>%setNames(names(allcodons)),function(i){
		#	
		codon=names(allcodons)[[i]]
		message(codon)
		codmatches<-vmatchPattern(pattern=codon,exonseq[toptrs])#exclude the start ccodon
		#
		nmatches = 	codmatches%>%elementNROWS 
		#
		matchgr<-codmatches%>%unlist%>%GRanges(names(.),.)
		matchgr$cdspos = start(matchgr) - cdsstarts[as.vector(seqnames(matchgr))]
		matchgr%<>%subset(cdspos %%3 == 0)
		seqlengths(matchgr) = exonsgrl%>%width%>%.[seqlevels(matchgr)]%>%sum
		codmatchwindows<-matchgr%>%resize(width(.)+(2*(3*FLANKCODS)),'center')
		codmatchwindows <- codmatchwindows[!is_out_of_bounds(codmatchwindows)]
		codmatchwindows%<>%subsetByOverlaps(trcds)
		codmatchwindows
	})
	saveRDS(codmatchwindowlist,here('data/codmatchwindowlist.rds'))
}else{
	codmatchwindowlist<-readRDS(here('data/codmatchwindowlist.rds'))
}

#
trseqinfo = Seqinfo(seqnames=names(trlens),seqlengths=trlens)
#

# startsigs = 


reduce=purrr::reduce
startproflist = 
	imap(fpcovlist['E13_ribo_1'],function(sampfpcov,sampname){
		trsums = sampfpcov%>%map(sum)%>%reduce(`+`)#sum over counts for that transcript
		sampfpcov%>%imap(function(rlfpcov,rl){
			rl=as.numeric(rl)
			rlfpcov = rlfpcov/(trsums)
			# rlfpcov = sampfpcov[[rl]]
			rloffset = offsets%>%filter(length==rl)%>%.$offset
			startwinds = cdsstarts%>%
				enframe('seqnames','start')%>%mutate(end=start)%>%
				GRanges%>%
				{seqinfo(.)=trseqinfo[seqlevels(.)];.}%>%
				{suppressWarnings({resize(.,3,'start')%>%resize(6,'end')%>%resize(45+3,'start')%>%shift(-rloffset)})}%>%
				.[!is_out_of_bounds(.)]
			#
			message('.')
			startwinds = startwinds%>%subset(seqnames%in%names(trsums)[trsums!=0])
			startwindsums = rlfpcov[startwinds]%>%as.matrix%>%colMeans
		})
	})

sampname = 'E13_ribo_1'
sampfpcov = fpcovlist[[1]]

endproflist = 
	imap(fpcovlist['E13_ribo_1'],function(sampfpcov,sampname){
		trsums = sampfpcov%>%map(sum)%>%reduce(`+`)#sum over counts for that transcript
		sampfpcov%>%imap(function(rlfpcov,rl){
			rl=as.numeric(rl)
			rlfpcov = rlfpcov/(trsums)
			# rlfpcov = sampfpcov[[rl]]
			rloffset = offsets%>%filter(length==rl)%>%.$offset
			startwinds = cdsends%>%
				enframe('seqnames','start')%>%mutate(end=start)%>%
				GRanges%>%
				{seqinfo(.)=trseqinfo[seqlevels(.)];.}%>%
				{suppressWarnings({resize(.,3,'end')%>%resize(6,'start')%>%resize(45+3,'end')%>%shift(-rloffset)})}%>%
				.[!is_out_of_bounds(.)]
			#
			message('.')
			startwinds = startwinds%>%subset(seqnames%in%names(trsums)[trsums!=0])
			startwindsums = rlfpcov[startwinds]%>%as.matrix%>%colMeans
		})
	})
stop()

startproflist[[1]][['28']][1:9]%>%txtplot
endproflist[[1]][['28']][(48-2-6):48]%>%txtplot

startproflist[[1]][['29']][1:9]%>%txtplot
endproflist[[1]][['29']][(48-2-6):48]%>%txtplot



startproflist[[1]]%>%map(.%>%matrix(byrow=TRUE,ncol=3)%>%{cbind(.,.-.)}%>%t%>%{txtplot(.,ylim=c(1,max(.)),height=20,width=100)})

#
startwinds = cdsstarts%>%
	enframe('seqnames','start')%>%mutate(end=start)%>%
	GRanges%>%
	{seqinfo(.)=trseqinfo[seqlevels(.)];.}%>%
	{suppressWarnings({resize(.,3,'start')%>%resize(30,'start')%>%shift(-rloffset)})}%>%
	.[!is_out_of_bounds(.)]
#

startproflist[[1]][['28']]%>%matrix(byrow=TRUE,ncol=3)%>%{cbind(.,.-.)}%>%t%>%{txtplot(.,ylim=c(min(.)*2,max(.)),height=20,width=100)}


startproflist[[1]][['27']][1:9]%>%txtplot
startproflist[[1]][['28']][1:9]%>%txtplot
startproflist[[1]][['29']][1:9]%>%txtplot
startproflist[[1]][['30']][1:9]%>%txtplot


startproflist[[1]]%>%reduce(`+`)%>%txtplot

endproflist[[1]]%>%reduce(`+`)%>%txtplot


allcodlist=codmatchwindowlist%>%GRangesList%>%unlist
cods = names(allcodlist)%>%str_split('\\.')%>%map_chr(1)
sampfpcov=fpcovlist[[1]]
# if(!file.exists(here('data/fpprofilelist.rds'))){
fpprofilelist = 
	imap(fpcovlist,function(sampfpcov,sampname){
		trsums = sampfpcov%>%map(sum)%>%reduce(`+`)#sum over counts for that transcript
		sampfpcov%>%lapply(function(rlfpcov){
			rlfpcov = rlfpcov/(trsums)
			allcodlistnz = allcodlist%>%subset(seqnames%in%names(trsums)[trsums!=0])
			message('.')
			rlfpcov[allcodlistnz]%>%split(cods)%>%lapply(as.matrix)%>%map(colMeans)
		})
	})


saveRDS(fpprofilelist,here('data/fpprofilelist.rds'))

fpprofilelist<-readRDS(here('data/fpprofilelist.rds'))


################################################################################
########testing the pca based a-site calls
################################################################################

codonprofiledat = fpprofilelist%>%map_depth(3,.%>%enframe('position','count'))%>%map_df(.id='sample',.%>%map_df(.id='readlen',.%>%bind_rows(.id='codon')))
#
windowoffsets = fpprofilelist[[1]][[1]][[1]]%>%length%>%seq(1,.)%>%subtract(1+(FLANKCODS*3))%>%multiply_by(-1)
windowpos = fpprofilelist[[1]][[1]][[1]]%>%length%>%seq(1,.)%>%subtract(1+(FLANKCODS*3))
# codonprofiledat$position = windowpos[codonprofiledat$position]
codonprofiledat%<>%mutate(position = position - 1 - (FLANKCODS*3))
codonprofiledat%<>%group_by(sample,readlen,codon)%>%mutate(count= count / median(count))
codonprofiledat%<>%filter(!codon %in% c('TAG','TAA','TGA'))
#
profvarpca = codonprofiledat%>%
	split(.,.$sample)%>%
	map_df(.id='sample',.%>%
		split(.,list(.$readlen))%>%
		# .[[1]]%>%
		map_df( .id='readlen',.%>%
			mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
			# filter(position> -numreadlen-1,position < 1)%>%
			filter(position> -numreadlen+6,position < -6)%>%
			# filter(sample=='E13_ribo_1')%>%
			ungroup%>%
			select(-numreadlen,-readlen,-sample)%>%
			# group_by(codon)%>%mutate(count/median(count))%>%
			spread(position,count)%>%
			{set_rownames(.[,-1],.$codon)}%>%
			princomp%>%{.$loadings[,1]}%>%{./.[which.max(abs(.))]}%>%enframe('position','pca1')
		)
	)
#
profvarpca%<>%select(sample,readlen,position,pca1)
profvarpca$readlen = paste0('rl',profvarpca$readlen)
library(rlang)
offsets <- read_tsv('ext_data/offsets_manual.tsv')
#

plotfile<-'plots/figures/figure2/trna_codons/fppos_vs_codon_pcascore.pdf'
offsets%<>%mutate(readlen=paste0('rl',length))
pdf(plotfile,w=12,h=12)
profvarpca%>%
	# slice_by(sample,c(1,2,3,4,5,6))%>%
	filter(sample%>%str_detect('ribo_1'))%>%
ggplot(data=.,aes(y=pca1,x=as.numeric(position)))+geom_point()+
	facet_grid(readlen~sample)+
		geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=offsets,aes(xintercept= -offset-5),color=I('green'),linetype=2)
dev.off()
normalizePath(plotfile)


{
offsets%<>%mutate(readlen=paste0(length))
pdf('plots/figures/figure2/trna_codons/fppos_vs_codon_variance.pdf',w=12,h=12)
#plotting variance amongst codons at each point.
codonprofiledat%>%
	ungroup%>%
	group_by(sample,readlen,position)%>%
	# filter(count==0)%>%.$codon%>%unique
	# filter(signal!=0)%>%
	# filter(position==1)%>%
	filter(sample%>%str_detect('ribo_1'))%>%
	filter(!is.nan(count))%>%
	summarise(sdsig=sd(count,na.rm=T)/median(count,na.rm=T))%>%
	separate(sample,c('time','assay','rep'))%>%
	group_by(readlen,time,assay,position)%>%
	summarise(sdsig=mean(sdsig))%>%
	mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
	filter(position> -numreadlen+6,position < -6)%>%
	filter(numreadlen>=25,numreadlen<=31)%>%
	arrange(position)%>%
	{
		qplot(data=.,x=position,y=sdsig)+
		theme_bw()+
		facet_grid(readlen~time)+
		scale_y_continuous('between codon variation (meannorm)')+
		scale_x_continuous('5 read position relative to codon ')+
		geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=offsets,aes(xintercept= -offset-5),color=I('green'),linetype=2)+
		ggtitle("variance of 5' read occurance vs position")
	}%>%print
dev.off()
normalizePath('plots/figures/figure2/trna_codons/fppos_vs_codon_variance.pdf')
}
pcaderrivedpos = profvarpca%>%
	# mutate(pca13wind = pca1+lead(pca1)+lead(pca1,2))%>%
	mutate(pca13wind = pca1)%>%
	group_by(readlen,position)%>%summarise(pca13wind=mean(pca13wind))%>%
	group_by(readlen)%>%slice(which.max(pca13wind))
# pcaderrivedpos=pcaderrivedpos%>%mutate(position = map(as.numeric(position),~ (.:(.+2))))%>%unnest

codonprofiles_pcawind = codonprofiles%>%semi_join(pcaderrivedpos%>%mutate(position=as.numeric(position)))
codonprofiles_pcawind = codonprofiles_pcawind%>%group_by(sample,codon,fraction)%>%
	summarise(occ_nonorm=sum(occ_nonorm),occupancy=sum(occupancy))%>%
	separate(sample,c('time','assay','rep'))%>%
	group_by(codon,time)%>%
	summarise(occ_nonorm=sum(occ_nonorm),occupancy=sum(occupancy))
	
codonprofiles_pcawind%>%inner_join(timecodon_tAb)%>%
	filter(.,time=='E13')%>%{
	quicktest(.$abundance,.$occ_nonorm)
}
