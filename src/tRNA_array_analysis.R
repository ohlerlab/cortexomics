library(readxl)
library(here)
library(magrittr)
library(tidyverse)
library(conflicted)

################################################################################
########Reading the tRNA data
################################################################################
source(here('src/R/Rprofile.R'))

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
trnaexprdf<-map2(datasources$file,datasources$case,.f=function(file,case){
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

allcodonsig<- trnaexprdf%>%setNames(datasources$sample)%>%bind_rows(.id='sample')

codonfromdecoder = .%>%str_replace('-\\d+$','')%>%str_extract('')
allcodonsig%<>%mutate(anticodon = decoder%>%str_replace('-\\d+$',''))
allcodonsig%<>%mutate(iscodon = str_detect(decoder,'\\w+-\\w+-\\d$'))
allcodonsig%<>%mutate(time = sample%>%str_extract('(?<= ).*$'))
suppressMessages({allcodonsig$rep1%<>%as.numeric})
suppressMessages({allcodonsig$rep2%<>%as.numeric})
allcodonsig <- gather(allcodonsig,rep,signal,rep1:rep2)
#get normfacts
#normfacts <- allcodonsig%>%group_by(sample,rep)%>%summarise(normfact = mean(signal[decoder%in%hkgenes2use]))
#Now normalize

allcodonsig %<>% group_by(sample,rep) %>% mutate(signal = signal - mean(signal[decoder%in%hkgenes2use]))


allcodsigmean <- allcodonsig%>%
	filter(iscodon)%>%
	group_by(time,sample,anticodon,decoder)%>%
	summarise(signal=mean(na.rm=T,as.numeric(signal)))
#
allcodsigmean%<>%mutate(fraction=str_extract(sample,'\\w+'))
#
allcodsigmean_isomerge<-allcodsigmean%>%
	# filter(sample)%>%
	# filter(sample%>%str_detect('Poly'))%>%
	# mutate(AA = GENETIC_CODE[codon]%>%qs('S+'))%>%
	mutate(fraction = str_extract(sample,'\\w+'))%>%
	group_by(fraction,time,sample,anticodon)%>%
	summarise(signal = log2(sum(2^(-signal),na.rm=T)))
	# spread(fraction,signal)
#
allcodsigmean_isomerge%>%filter(anticodon=='Val-AAC')
glutRNAabund<-allcodsigmean_isomerge%>%filter(anticodon%>%str_detect('CTC'))%>%.$signal
valtRNAabund<-allcodsigmean_isomerge%>%filter(anticodon%>%str_detect('AAC'))%>%.$signal
#

#add codon info
allcodsigmean%<>%mutate(codon = str_extract(anticodon,'[^-]+$')%>%DNAStringSet%>%
		reverseComplement%>%
		as.character%>%qs('S+'))
allcodsigmean_isomerge%<>%mutate(codon = str_extract(anticodon,'[^-]+$')%>%DNAStringSet%>%
		reverseComplement%>%
		as.character%>%qs('S+'))
allcodsigmean%<>%mutate(AA = GENETIC_CODE[str_extract(codon,'[^\\-]+$')]%>%qs('S+'))
allcodsigmean_isomerge%<>%mutate(AA = GENETIC_CODE[str_extract(codon,'[^\\-]+$')]%>%qs('S+'))
#
allcodsigmean%>%mutate(as.character(translate(DNAStringSet(codon))))%>%head(1)%>%{stopifnot(.$anticodon=='Ala-AGC' & (.$codon=='GCT'))}
allcodsigmean_isomerge%>%mutate(as.character(translate(DNAStringSet(codon))))%>%head(1)%>%{stopifnot(.$anticodon=='Ala-AGC' & (.$codon=='GCT'))}
#
glutRNAabund<-allcodsigmean_isomerge%>%filter(codon%>%str_detect('CTC'))%>%.$signal
valtRNAabund<-allcodsigmean_isomerge%>%filter(codon%>%str_detect('AAC'))%>%.$signal
stopifnot(valtRNAabund<glutRNAabund)
#


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

################################################################################
########Now plot it
################################################################################
stagecolsdisplay <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
stageconv = names(stagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
stagecols <- stagecolsdisplay%>%setNames(c('E13','E145','E16','E175','P0'))
#how each decoder
plotfile<-'plots/figures/figure2/trna_codons/trna_sig_total_alldecod.pdf'
pdf(plotfile,w=9,h=16)
allcodsigmean%>%
	group_by(decoder)%>%
	mutate(dmean=mean(signal))%>%
	group_by(anticodon)%>%
	arrange(dmean)%>%
	mutate(decoder = as_factor(decoder))%>%
	filter(sample%>%str_detect('Total'))%>%
	group_by(anticodon)%>%
	# group_by(decoder,)
	ggplot(data=.,aes(x=decoder,color=time,y= - signal))+
	geom_point()+
	scale_color_manual(values=stagecols)+
	scale_y_continuous(name='- deltaCt',breaks=seq(0,-30,by= - 5))+
	coord_flip()+
	theme_bw()
	# theme(axis.text.x=element_text(angle=45,vjust=.5,size=6))
	# facet_grid(time ~ . )
dev.off()
normalizePath(plotfile)
#
#how each decoder
plotfile<-'plots/figures/figure2/trna_codons/trna_sig_poly_alldecod.pdf'
pdf(plotfile,w=9,h=16)
allcodsigmean%>%
	group_by(decoder)%>%
	mutate(dmean=mean(signal))%>%
	group_by(anticodon)%>%arrange(dmean)%>%
	mutate(decoder = as_factor(decoder))%>%
	filter(sample%>%str_detect('Poly'))%>%
	group_by(anticodon)%>%
	# group_by(decoder,)
	ggplot(data=.,aes(x=decoder,color=time,y= -signal))+
	geom_point()+
	scale_color_manual(values=stagecols)+
	scale_y_continuous(name='- deltaCt',breaks=seq(0,-30,by= - 5))+
	coord_flip()+
	theme_bw()
	# theme(axis.text.x=element_text(angle=45,vjust=.5,size=6))
	# facet_grid(time ~ . )
dev.off()
normalizePath(plotfile)


plotfile<-'plots/figures/figure2/trna_codons/trna_sig_total_mergedecod.pdf'
pdf(plotfile,w=9,h=16)
allcodsigmean_isomerge%>%
	filter(sample%>%str_detect('Total'))%>%
	group_by(anticodon)%>%
	mutate(cmean=mean(signal))%>%
	ungroup%>%
	arrange(cmean)%>%
	mutate(anticodon = as_factor(anticodon))%>%
	# group_by(decoder,)
	ggplot(data=.,aes(x=anticodon,color=stageconv[time],y=signal))+
	geom_point()+
	scale_color_manual(values=stagecols)+
	scale_y_continuous(name='- deltaCt',breaks=-seq(0,20,by=2.5))+
	coord_flip()+
	theme_bw()
	# theme(axis.text.x=element_text(angle=45,vjust=.5,size=6))
	# facet_grid(time ~ . )
dev.off()
normalizePath(plotfile)
#
plotfile<-'plots/figures/figure2/trna_codons/trna_sig_poly_mergedecod.pdf'
pdf(plotfile,w=9,h=16)
allcodsigmean_isomerge%>%
	filter(sample%>%str_detect('Poly'))%>%
	group_by(anticodon)%>%
	mutate(cmean=mean(signal))%>%
	ungroup%>%
	arrange(cmean)%>%
	mutate(anticodon = as_factor(anticodon))%>%
	# group_by(decoder,)
	ggplot(data=.,aes(x=anticodon,color=stageconv[time],y=signal))+
	geom_point()+
	scale_color_manual(values=stagecols)+
	scale_y_continuous(name='- deltaCt',breaks=-seq(0,20,by=2.5))+
	coord_flip()+
	theme_bw()
	# theme(axis.text.x=element_text(angle=45,vjust=.5,size=6))
	# facet_grid(time ~ . )
dev.off()
normalizePath(plotfile)


#group them on same row
plotfile<-'plots/figures/figure2/trna_codons/isomerge_tRNA_sig_strip_poly.pdf'
pdf(plotfile,w=9,h=16)
allcodsigmean_isomerge%>%
	filter(sample%>%str_detect('Total'))%>%
	group_by(anticodon)%>%
	mutate(cmean=mean(signal))%>%
	ungroup%>%
	arrange(cmean)%>%
	mutate(anticodon = as_factor(anticodon))%>%
	# group_by(decoder,)
	ggplot(data=.,aes(x=anticodon,color=time,y=signal))+
	geom_point()+
	scale_color_manual(values=stagecols)+
	scale_y_continuous(breaks=seq(0,20,by=2.5))+
	coord_flip()+
	theme_bw()
	# theme(axis.text.x=element_text(angle=45,vjust=.5,size=6))
	# facet_grid(time ~ . )
dev.off()
normalizePath(plotfile)








################################################################################
########Codon optimality scores based on codon frequencies, gene expression levels
################################################################################
cdswidths<-cds%>%split(.,.$protein_id)%>%width%>%sum
longestcdsids<-enframe(cdswidths,'protein_id','length')%>%left_join(cds%>%mcols%>%as_tibble%>%distinct(protein_id,gene_id))%>%group_by(gene_id)%>%slice(which.max(length))%>%.$protein_id
cds_noov<-cds%>%subset(protein_id%in%longestcdsids)
cds_noov%<>%split(.,.$protein_id)
# cds_noov<-cds_noov%>%subset(strand=='+')%>%split(.,.$protein_id)%>%head(10)%>%lapply(head,1)%>%GRangesList%>%unlist%>%resize(3)

times <- names(stagecols)%>%setNames(.,.)
library(Rsamtools)
codonfreqs<-mymemoise(oligonucleotideFrequency)(extractTranscriptSeqs(FaFile(REF),cds_noov),3,step=3)
rownames(codonfreqs)<-names(cds_noov)
overallcodonfreqs<-codonfreqs%>%colSums
bestcds<- cds%>%subset(protein_id%in%best_protein_ids)%>%split(.,.$protein_id)
usedcodonfreqs<-mymemoise(oligonucleotideFrequency)(extractTranscriptSeqs(FaFile(REF),bestcds),3,step=3)%>%set_rownames(names(bestcds))
codon_is_optimal <- overallcodonfreqs%>%enframe('codon','freq')%>%
	mutate(AA=as.character(translate(DNAStringSet(codon))))%>%
	group_by(AA)%>%mutate(optimal = freq==max(freq))%>%
	{setNames(.$optimal,.$codon)}

optimalcodons <- names(codon_is_optimal)[codon_is_optimal]
pc_optcodons <- ((usedcodonfreqs[,optimalcodons]%>%rowSums)/(usedcodonfreqs%>%rowSums))%>%setNames(rownames(usedcodonfreqs))%>%enframe('protein_id','pc_opt')




cdsexprvals <- (2^bestmscountvoom$E[,T])%>%{rownames(.)%<>%str_replace('_\\d+$','');.}

#length norm themj
cdsexprvals <- cdsexprvals / fData(countexprdata)$length[match(rownames(cdsexprvals),fData(countexprdata)$protein_id)]


cdsexprdf<- cdsexprvals[,1:20] %>% as.data.frame%>% rownames_to_column('protein_id')%>%gather(sample,signal,-protein_id)%>%
	separate(sample,c('time','assay','rep'))%>%group_by(protein_id,time,assay)%>%summarise(signal=mean(signal))%>%
	filter(assay=='ribo')%>%
	spread(time,signal)%>%
	arrange(match(protein_id,rownames(usedcodonfreqs)))

weighted_codon_usage <- lapply(times,function(itime)(usedcodonfreqs[,]*cdsexprdf[[itime]])%>%colSums%>%{./sum(.)}%>%enframe('codon','weightedusage'))%>%bind_rows(.id='time')

tRNA_abchange <- allcodsigmean_isomerge%>%
	group_by(fraction,codon)%>%
	filter(!is.na(signal))%>%
	filter(is.finite(signal),!is.na(signal))%>%
	nest%>%
	mutate(abundancechange = map_dbl(data,~{lm(data=.,signal ~ seq_along(time))$coef[2]}))

usage_v_abundance_df <- weighted_codon_usage%>%filter(time=='E13')%>%ungroup%>%select(-time)%>%left_join(tRNA_abchange%>%filter(is.finite(abundancechange))%>%select(fraction,codon,abundancechange))%>%
	filter(!is.na(abundancechange),!is.na(weightedusage))

wus_tab_cors<-usage_v_abundance_df%>%
	filter(!codon%in%c('TAG','TAA','TGA'))%>%
	group_by(fraction)%>%
	nest%>%mutate(
		cor = map_dbl(data,~ cor(.$abundancechange,.$weightedusage)),
		pval = map_dbl(data,~ cor.test(.$abundancechange,.$weightedusage)$p.value)
	)


fractions<-c('Poly','Total')
for(ifraction in fractions){
plotfile <- str_interp('plots/figures/figure2/trna_codons/trna_abchange_exprusage_${ifraction}.pdf')
pdf(plotfile)
usage_v_abundance_df%>%
	filter(!codon%in%c('TAG','TAA','TGA'))%>%
	filter(fraction == ifraction )%T>%{message(nrow(.))}%>%
	nest%>%
	mutate(labl = paste0(
		'r = ',map_dbl(data,~cor(.$abundancechange,.$weightedusage))%>%round(3),'\n',
		'p = ',map_dbl(data,~cor.test(.$abundancechange,.$weightedusage)$p.value%>%round(3))
	))%>%
	unnest%>%
	{
		print(
			qplot(data=.,x=.$abundancechange,y=.$weightedusage,label=codon,geom='blank')+
				scale_x_continuous('tRNA Abundance Signal Slope E13 - P0')+
				scale_y_continuous('RiboExpression Weighted E13 ')+
				geom_text(data=distinct(.,labl),hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))+
				geom_smooth(method='lm')+
				# facet_grid(tRNA_time~.)+
				geom_text()+
				geom_text(data=wus_tab_cors,aes(x=0,y=0,label=paste0('r = ',round(cor,2))))+
				theme_bw()
			)
	}
dev.off()
message(normalizePath(plotfile))
}





#are and global frequency correlated? - AA level?
pdf('plots/figures/figure2/trna_codons/AAusage_vs_summed_tRNAab.pdf')
allcodsigmean_isomerge%>%
	left_join(enframe(usedcodonfreqs%>%colSums,'codon','freq'))%>%
	group_by(fraction,time,AA)%>%
	# slice(which.max(signal))%>%
	# slice(which.max(freq))%>%
	summarise(signal = log2(sum(2^signal)),freq=sum(freq))%>%
	# summarise(signal = mean(signal),freq=sum(freq))%>%
	group_by(fraction,time)%>%
	nest%>%
	mutate(labl = paste0(
		'r = ',map_dbl(data,~cor(use='complete',.$freq,.$signal))%>%round(3),'\n',
		'p = ',map_dbl(data,~cor.test(use='complete',.$freq,.$signal)$p.value%>%round(3))
	))%>%
	unnest%>%
	{ggplot(.,aes(x=signal,y=freq))+geom_point()+facet_grid(time~fraction)+
		theme_bw()+geom_text(data=distinct(.,time,fraction,labl),hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))}+
	ggtitle('Amino acid usage frequency vs summed tRNA abundance')
dev.off()
normalizePath('plots/figures/figure2/trna_codons/AAusage_vs_summed_tRNAab.pdf')

#Are occupancy and global frequency correlated? - AA level??
pdf('plots/figures/figure2/trna_codons/AAweightedusage_vs_summed_tRNAab.pdf')
allcodsigmean_isomerge%>%
	left_join(weighted_codon_usage)%>%
	group_by(fraction,time,AA)%>%
	summarise(signal = log2(sum(2^signal)),weightedusage=sum(weightedusage))%>%
	filter(is.finite(signal))%>%
	group_by(fraction,time)%>%
	nest%>%
	mutate(labl = paste0(
		'r = ',map_dbl(data,~cor(use='complete',.$weightedusage,.$signal))%>%round(3),'\n',
		'p = ',map_dbl(data,~cor.test(use='complete',.$weightedusage,.$signal)$p.value%>%round(3))
	))%>%
	unnest%>%
	{ggplot(.,aes(x=signal,y=weightedusage))+geom_point()+facet_grid(time~fraction)+
		theme_bw()+geom_text(data=distinct(.,time,fraction,labl),hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))}+
	ggtitle('Amino acid expr weighted usage frequency vs summed tRNA abundance')
dev.off()
normalizePath('plots/figures/figure2/trna_codons/AAweightedusage_vs_summed_tRNAab.pdf')


#Are occupancy and global signal correlated?
pdf('plots/figures/figure2/trna_codons/codon_usage_vs_summed_tRNAab.pdf')
allcodsigmean_isomerge%>%
	mutate(signal = signal)%>%
	left_join(enframe(usedcodonfreqs%>%colSums,'codon','freq'))%>%
	filter(is.finite(signal))%>%
	group_by(fraction,time)%>%
	nest%>%
	mutate(labl = paste0(
		'r = ',map_dbl(data,~cor(use='complete',.$freq,.$signal))%>%round(3),'\n',
		'p = ',map_dbl(data,~cor.test(use='complete',.$freq,.$signal)$p.value%>%round(3))
	))%>%
	unnest%>%
	{ggplot(.,aes(x=signal,y=freq))+geom_point()+facet_grid(time~fraction)+
		theme_bw()+geom_text(data=distinct(.,time,fraction,labl),hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))}+
		ggtitle('Codon usage frequency vs summed tRNA abundance')
dev.off()
normalizePath('plots/figures/figure2/trna_codons/codon_usage_vs_summed_tRNAab.pdf')

#Are occupancy and global signal correlated?
pdf('plots/figures/figure2/trna_codons/codon_usage_vs_summed_tRNAab.pdf')
allcodsigmean_isomerge%>%
	mutate(signal = signal)%>%
	left_join(weighted_codon_usage)%>%
	filter(is.finite(signal))%>%
	group_by(fraction,time)%>%
	nest%>%
	mutate(labl = paste0(
		'r = ',map_dbl(data,~cor(use='complete',.$weightedusage,.$signal))%>%round(3),'\n',
		'p = ',map_dbl(data,~cor.test(use='complete',.$weightedusage,.$signal)$p.value%>%round(3))
	))%>%
	unnest%>%
	{ggplot(.,aes(x=signal,y=weightedusage))+geom_point()+facet_grid(time~fraction)+
		theme_bw()+geom_text(data=distinct(.,time,fraction,labl),hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))}+
		ggtitle('Codon usage frequency vs summed tRNA abundance')
dev.off()
normalizePath('plots/figures/figure2/trna_codons/codon_usage_vs_summed_tRNAab.pdf')

