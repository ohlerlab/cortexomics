library(readxl)
library(tidyverse)

################################################################################
########Reading the tRNA data
################################################################################
hkgenes2use<-c("5S rRNA",'18S rRNA')

alltRNAreps <- Sys.glob(here('ext_data/tRNA_data/*Rep*/*'))%>%str_subset('/[PT8][^/]+$')

#data is pulled from a combination of report and 'test/'ctrl' status
datasources <- tibble(file=c(alltRNAreps,alltRNAreps[1]),
	case=c(rep('Test',length(alltRNAreps)),'Control')
)

datasources%<>%mutate(sample = 
ifelse(case=='Test',
	str_extract(basename(file),regex('.*(?= VS )')),
	str_extract(file,regex('(?<= VS ).*(?=.xlsx)'))
)
)

#Now iterate over these
exprdf<-map2(datasources$file,datasources$case,.f=function(file,case){
	excel_sheets(file)
	hksheet <- readxl::read_xlsx(file,sheet='Choose Housekeeping Genes')
	hks2use_rows <- which(hksheet[[1]][1:4] %in% hkgenes2use)
	casecol <- if(case=='Test') 4 else 27


	hknormfacts <- list(rep1=mean(as.numeric(hksheet[[casecol]][hks2use_rows])),rep2=mean(as.numeric(hksheet[[casecol+1]][hks2use_rows])))
	# calibfactrow <- which(hksheet[[1]]=="Calibration factor")
	calibfactrow <- which(hksheet[[1]]=="PPC")
	ppcnormfacts <- list(rep1=mean(as.numeric(hksheet[[casecol]][calibfactrow])),rep2=mean(as.numeric(hksheet[[casecol+1]][calibfactrow])))

	exprsheet <- if(case=='Test') 'Test Sample Data' else 'Control Sample Data' 

	message(file)	
	suppressMessages({exprsheet<-read_xlsx(file,sheet=exprsheet)})

	exprdf = exprsheet[c(-1),1:4]%>%as.data.frame%>%head(-2)%>%set_colnames(c('decoder','well','rep1','rep2'))
	exprdf

})

allcodonsig<- exprdf%>%setNames(datasources$sample)%>%bind_rows(.id='sample')

allcodonsig%<>%mutate(codon = decoder%>%str_replace('-\\d+$',''))
allcodonsig%<>%mutate(iscodon = str_detect(decoder,'\\w+-\\w+-\\d$'))
allcodonsig%<>%mutate(time = sample%>%str_extract('(?<= ).*$'))
suppressMessages({allcodonsig$rep1%<>%as.numeric})
suppressMessages({allcodonsig$rep2%<>%as.numeric})
allcodonsig <- gather(allcodonsig,rep,signal,rep1:rep2)
#get normfacts
#normfacts <- allcodonsig%>%group_by(sample,rep)%>%summarise(normfact = mean(signal[decoder%in%hkgenes2use]))
#Now normalize

allcodonsig %<>% group_by(sample,rep) %>% mutate(signal = signal - mean(signal[decoder%in%hkgenes2use]))




#
 #=IF(ISERROR(AVERAGE(D9:D11));"";AVERAGE(D9:D11))

# #D17 this is a lookup - it gets the value in A17, (PPC), looks it up in the test sheet
# #Then 
# IF($A17="";"";
# 	IF(
# 		VLOOKUP($A17;'Test Sample Data'!$A$3:$V$194;D$16+2;FALSE)=0;
# 		"";
# 		VLOOKUP($A17;'Test Sample Data'!$A$3:$V$194;D$16+2;FALSE)
# 	)
# )

# #D18 Tis seems totally pointless - guess it's form other tables where it matters
# =IF(ISERROR(AVERAGE('Choose Housekeeping Genes'!D$17));"";AVERAGE('Choose Housekeeping Genes'!D$17))

# #D19 IPC overall averag, sums rows
# =IF(ISNUMBER(D18);
# 	AVERAGE(
# 		'Choose Housekeeping Genes'!$D$17:$W$17;
# 		'Choose Housekeeping Genes'!$AA$17:$AT$17);
# 	"")

# #forCalibration factor D20 - (a subtraction) use the correction factor, or 0 if it's not calculable
# IF(
# 	ISNUMBER('Choose Housekeeping Genes'!D$18-'Choose Housekeeping Genes'!D$19);
# 	'Choose Housekeeping Genes'!D$18-'Choose Housekeeping Genes'!D$19;
# 	0
# )

# #here's the excel formula dictating their numbers, refers to calibration factor on D/E20
# #This is pre normalization
# IF(
# 	SUM('Test Sample Data'!C$3:C$194)>10;
# 	IF(
# 		AND(
# 			ISNUMBER('Test Sample Data'!C3);
# 			'Test Sample Data'!C3<35;
# 			'Test Sample Data'!C3>0
# 		);
# 		'Test Sample Data'!C3-'Choose Housekeeping Genes'!D$20;35-'Choose Housekeeping Genes'!D$20
# 	);""
# )

# #the housekeeping norms have some tortured logic to deal with missing info
# =IF($A9="";"";IF(VLOOKUP($A9;'Choose Housekeeping Genes'!$A$3:$W$5;D$8+3;FALSE)=0;"";IF(ISNUMBER(VLOOKUP($A9;'Choose Housekeeping Genes'!$A$3:$W$5;D$8+3;FALSE));VLOOKUP($A9;'Choose Housekeeping Genes'!$A$3:$W$5;D$8+3;FALSE);"")))

# #Meanwhile the actual housekeeping gene average is here
# '=IF(ISERROR(AVERAGE(D9:D11));"";AVERAGE(D9:D11))'

# #which is used for The final number - after calibration using PPC, they subtrat 
# =IF(ISNUMBER(C3-'Choose Housekeeping Genes'!D$12);C3-'Choose Housekeeping Genes'!D$12;"")


# #Get the control values for each one
# #Now get the actual values for each one
# #Deal with the missing values
# #Produce output of the form
# set_colnames,c('file','decoder','rep1','rep2','sample')


# #Notes - I want to pull out the 2 normalized values per probe from each 
# tRNArepsheetnames <- c("Instructions", "Transcript table", "Test Sample Data", "Control Sample Data",
# "Choose Housekeeping Genes", "Data pre-processing", "All expressed genes",
# "Scatter Plot", "Volcano Plot", "up_Bar Graph", "down_Bar Graph"
# )

# readxl::excel_sheets(here('ext_data/tRNA_data/*Rep*/*')%>%Sys.glob%>%.[[1]])%>%dput

# allread <- alltRNAreps%>%
# 	str_subset('.xlsx$')%>%
# 	str_subset(neg=TRUE,'\\$')%>%
# 	setNames(.,basename(.))%>%
# 	lapply(.%>%readxl::read_xlsx(.,sheet='Data pre-processing'))

# testcolind = allread[[1]]%>%colnames%>%`==`('Test(normalization)')%>%which
# ctlcolind = allread[[1]]%>%colnames%>%`==`('Control(normalization)')%>%which

# startsamp<-allread['Total E145 VS Total E13.xlsx']%>%map_df(.id='file',.%>%tail(-1)%>%select(1,!!ctlcolind,!!ctlcolind+1))%>%
# 	mutate(sample=str_extract(file,regex('(?<= VS ).*(?=.xlsx)')))

# othersamps<-allread%>%map_df(.id='file',.%>%tail(-1)%>%select(1,!!testcolind,!!testcolind+1))%>%
# 	mutate(sample=str_extract(file,regex('.*(?= VS )')))

# allcodonsig <- list(startsamp,othersamps)%>%lapply(set_colnames,c('file','decoder','rep1','rep2','sample'))%>%bind_rows%>%
# 	gather(rep,signal,-file,-sample,-decoder)%>%group_by(file,sample,decoder,rep)%>%slice(1)

# allcodonsig%<>%mutate(codon = decoder%>%str_replace('-\\d+$',''))
# allcodonsig%<>%mutate(iscodon = str_detect(decoder,'\\w+-\\w+-\\d$'))
# allcodonsig%<>%mutate(time = sample%>%str_extract('(?<= ).*$'))



allcodsigmean <- allcodonsig%>%
	filter(iscodon)%>%
	group_by(time,sample,codon,decoder)%>%
	summarise(signal=mean(na.rm=T,as.numeric(signal)))

allcodsigmean%<>%mutate(fraction=str_extract(sample,'\\w+'))

allcodsigmean%>%filter(fraction=='Poly')

allcodsigmean%>%group_by(codon,time)%>%mutate(ispolyenriched=signal==max(signal[fraction=='Poly']))




allcodsigmean%>%filter(codon=='Val-AAC')

(2^(-8.33))
(2^(-35.33))

log2((2^(-8.33)) / 
(2^(-35.33)))


log2(2^(-8.33))
log2(2^(-35.33))

-log2(2^(-8.33))
-log2(2^(-35.33))

allcodsigmean_isomerge<-allcodsigmean%>%
	filter(sample%>%str_detect('Total'))%>%
	group_by(time,sample,codon)%>%
	# group_by(decoder)%>%group_slice(1)%>%
	summarise(signal = -log2(sum(2^(-signal),na.rm=T)))



allcodsigmean_isomerge%>%filter(codon=='Val-AAC')



stop()


# alltRNAreps%>%str_subset('Poly P0 VS Poly E13.xlsx')%>%read_xlsx(sheet='Data pre-processing')

	# filter(signal>12.0166,signal<12.0167)%>%as.data.frame%>%

# allcodonsig%>%.$signal%>%table%>%sort%>%tail
	# filter(signal>12.0166,signal<12.0167)%>%as.data.frame

# allcodsigmean$signal%>%table%>%sort%>%tail


################################################################################
########Now plot it
################################################################################
	
#how each decoder
plotfile<-'plots/figures/figure2/tRNA_sig_strip_allnorm.pdf'
pdf(plotfile,w=9,h=16)
allcodsigmean%>%
	group_by(decoder)%>%
	mutate(dmean=mean(signal))%>%
	group_by(codon)%>%arrange(dmean)%>%
	mutate(decoder = as_factor(decoder))%>%
	filter(sample%>%str_detect('Total'))%>%
	group_by(codon)%>%
	# group_by(decoder,)
	ggplot(data=.,aes(x=decoder,color=time,y=signal))+
	geom_point()+
	scale_color_manual(values=stagecols)+
	scale_y_continuous(breaks=seq(0,30,by=2.5))+
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


#group them on same row
plotfile<-'plots/figures/figure2/isomerge_tRNA_sig_strip_poly.pdf'
pdf(plotfile,w=9,h=16)
allcodsigmean_isomerge%>%
	filter(sample%>%str_detect('Total'))%>%
	group_by(codon)%>%
	mutate(cmean=mean(signal))%>%
	ungroup%>%
	arrange(cmean)%>%
	mutate(codon = as_factor(codon))%>%
	# group_by(decoder,)
	ggplot(data=.,aes(x=codon,color=time,y=signal))+
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

