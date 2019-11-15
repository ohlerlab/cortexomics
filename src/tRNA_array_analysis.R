
#Notes - I want to pull out the 2 normalized values per probe from each 
readxl::excel_sheets(here('ext_data/tRNA_data/*Rep*/*')%>%Sys.glob%>%.[[1]])
alltRNAreps <- Sys.glob(here('ext_data/tRNA_data/*Rep*/*'))

allread <- alltRNAreps%>%
	str_subset('.xlsx$')%>%
	setNames(.,basename(.))%>%
	lapply(.%>%readxl::read_xlsx(.,sheet='Data pre-processing'))

testcolind = allread[[1]]%>%colnames%>%`==`('Test(normalization)')%>%which
ctlcolind = allread[[1]]%>%colnames%>%`==`('Control(normalization)')%>%which

startsamp<-allread['Total E145 VS Total E13.xlsx']%>%map_df(.id='file',.%>%tail(-1)%>%select(1,!!ctlcolind,!!ctlcolind+1))%>%
	mutate(sample=str_extract(file,regex('(?<= VS ).*(?=.xlsx)')))

othersamps<-allread%>%map_df(.id='file',.%>%tail(-1)%>%select(1,!!testcolind,!!testcolind+1))%>%
	mutate(sample=str_extract(file,regex('.*(?= VS )')))

allcodonsig <- list(startsamp,othersamps)%>%lapply(set_colnames,c('file','decoder','rep1','rep2','sample'))%>%bind_rows%>%
	gather(rep,signal,-file,-sample,-decoder)%>%group_by(file,sample,decoder,rep)%>%slice(1)

allcodonsig%<>%mutate(codon = decoder%>%str_replace('-\\d+$',''))
allcodonsig%<>%mutate(iscodon = str_detect(decoder,'\\w+-\\w+-\\d$'))
allcodonsig%<>%mutate(time = sample%>%str_extract('(?<= ).*$'))



allcodsigmean <- allcodonsig%>%
	filter(iscodon)%>%
	group_by(time,sample,codon,decoder)%>%
	summarise(signal=mean(as.numeric(signal)))

allcodonsig%>%filter(between(as.numeric(signal),11.8956,11.8957))
allcodonsig%>%filter(between(as.numeric(signal),11.8956,11.8957))%>%.$file%>%unique

alltRNAreps%>%str_subset('Poly P0 VS Poly E13.xlsx')%>%read_xlsx(sheet='Data pre-processing')

	filter(signal>12.0166,signal<12.0167)%>%as.data.frame%>%

allcodonsig%>%.$signal%>%table%>%sort%>%tail
	filter(signal>12.0166,signal<12.0167)%>%as.data.frame

allcodsigmean$signal%>%table%>%sort%>%tail

#how each decoder
plotfile<-'plots/figures/figure2/tRNA_sig_strip.pdf'
pdf(plotfile,w=9,h=16)
allcodsigmean%>%
	group_by(decoder)%>%
	mutate(dmean=mean(signal))%>%
	group_by(codon)%>%arrange(dmean)%>%
	mutate(decoder = as_factor(decoder))%>%
	filter(sample%>%str_detect('Total'))%>%
	group_by(codon)%>%
	# group_by(decoder,)
	ggplot(data=.,aes(x=decoder,color=stageconv[time],y=signal))+
	geom_point()+
	scale_color_manual(values=stagecols)+
	scale_y_continuous(breaks=seq(0,20,by=2.5))+
	coord_flip()+
	theme_bw()
	# theme(axis.text.x=element_text(angle=45,vjust=.5,size=6))
	# facet_grid(time ~ . )
dev.off()
normalizePath(plotfile)

allcodsigmean%>%.$codon%>%str_subset('Asp-')%>%unique

#group them on same row
plotfile<-'plots/figures/figure2/tRNA_sig_strip.pdf'
pdf(plotfile,w=9,h=16)
allcodsigmean%>%
	filter(sample%>%str_detect('Total'))%>%
	group_by(codon)%>%
	mutate(cmean=mean(signal))%>%
	ungroup%>%
	arrange(cmean)%>%
	mutate(codon = as_factor(codon))%>%
	# group_by(decoder,)
	ggplot(data=.,aes(x=codon,color=stageconv[time],y=signal))+
	geom_point()+
	scale_color_manual(values=stagecols)+
	scale_y_continuous(breaks=seq(0,20,by=2.5))+
	coord_flip()+
	theme_bw()
	# theme(axis.text.x=element_text(angle=45,vjust=.5,size=6))
	# facet_grid(time ~ . )
dev.off()
normalizePath(plotfile)



#group them on same row
plotfile<-'plots/figures/figure2/tRNA_sig_strip_poly.pdf'
pdf(plotfile,w=9,h=16)
allcodsigmean%>%
	filter(sample%>%str_detect('Poly'))%>%
	group_by(codon)%>%
	mutate(cmean=mean(signal))%>%
	ungroup%>%
	arrange(cmean)%>%
	mutate(codon = as_factor(codon))%>%
	# group_by(decoder,)
	ggplot(data=.,aes(x=codon,color=stageconv[time],y=signal))+
	geom_point()+
	scale_color_manual(values=stagecols)+
	scale_y_continuous(breaks=seq(0,20,by=2.5))+
	coord_flip()+
	theme_bw()
	# theme(axis.text.x=element_text(angle=45,vjust=.5,size=6))
	# facet_grid(time ~ . )
dev.off()
normalizePath(plotfile)

#df ranking codons by their mean difference from the mean signal at 
allcodsigmean%>%
	filter(sample=='Total E13')%>%
	separate(decoder,c('AA','codon','decod'))%>%
	group_by(time,sample,codon,decod)%>%
	group_by(sample)%>%
	mutate(allmean = mean(signal))%>%
	group_by(codon)%>%
	mutate(dev=min(abs(signal-allmean)))%>%
	arrange(desc(dev))

#This makes it seem like the codons tend to have a high and a low one 
#but no

pdf('tmp.pdf')
qplot()
dev.off()

#
cdstousetrna<-cds%>%split(.,.$protein_id)%>%.[highprotein_ids]%>%unlist

codons <- getSeq(cdstousetrna,x=FaFile(REF))%>%split(cdstousetrna$protein_id)%>%lapply(.%>%unlist%>%xscat)%>%DNAStringSet%>%vmatchPattern('TCT',.)%>%unlist%>%subset((start %% 3) == 1)

codonsgr <- GRanges(names(codons),codons)%>%mapFromTranscripts(cdstousetrna%>%split(.,.$protein_id))%>%subset(strand=='+')

codonsgr_exp <- codonsgr%>%resize(3+9+9,'center')


signalovercodons<-bams%>%str_subset('ribo')%>%mclapply(function(bam){
	riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=MAPQTHRESH,which=codonsgr_exp)
	reads <- readGAlignments(bam,param=riboparam)
	mcols(reads)$cdsshift <- get_cds_offsets(reads,psite_model$offsets,psite_model$compartments)
	#reads <- reads%>%subset(width==27)
	reads <- reads%>%subset(width %in% psite_model$offsets$length)
	mcols(reads)$length <- width(reads)
	reads%<>%subset(!is.na(cdsshift))
	psites <- apply_psite_offset(reads,c('cdsshift'))%>%as("GRanges")
	mcols(psites)$length <- mcols(reads)$length
	startwindpsites <- psites%>%mapToTranscripts(codonsgr_exp)
	start(startwindpsites)
})

codplotdf <- signalovercodons%>%
	lapply(table)%>%
	setNames(basename(bams%>%str_subset('ribo')))%>%
	map_df(.id='bam',enframe,'pos','count')%>%
	group_by(bam)%>%
	mutate(count=count/sum(count))





stagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
stageconv = names(stagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
bamcols <- bams%>%basename%>%str_extract('[^_]+')%>%stageconv[.]%>%stagecols[.]%>%setNames(basename(bams))


plotfile<-'plots/TCT_ribo_cov.pdf'%T>%pdf(h=4,w=8)
qplot(data=codplotdf,color=bam,x=as.numeric(pos)-10,y=count,geom='line')+
	scale_x_continuous(name='Position Relative to TCT')+
	# scale_color_discrete(guide=TRUE)+
	scale_color_manual(values=bamcols)+
	scale_y_continuous(name='Relative P Site Density')+
	theme_bw()
dev.off()
normalizePath(plotfile)%>%message

#TODO - make sure strand spec 

