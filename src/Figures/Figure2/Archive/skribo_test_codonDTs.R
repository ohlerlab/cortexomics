source("src/Figures/Figure2/8_tRNA_array_analysis.R")
allcodsigmean_isomerge%<>%dplyr::rename('abundance':='signal')
timecodon_tAb = allcodsigmean_isomerge%>%group_by(time,codon)%>%filter(fraction=='Total')%>%summarise(abundance=mean(abundance))
atimecodon_tAb = allcodsigmean_isomerge%>%group_by(codon)%>%filter(fraction=='Poly')%>%summarise(abundance=mean(abundance))

skribofiles = Sys.glob('pipeline/scikit_ribo/*/codons.csv')
skribofiles%<>%str_subset(neg=TRUE,'allstart')
names(skribofiles) = basename(dirname(skribofiles))

skribodts = skribofiles%>%map_df(.id='sample',fread)

timecodon_dt = skribodts%>%separate(sample,c('time','assay','rep'))%>%group_by(time,codon)%>%summarise(codon_dwell_time=mean(codon_dwell_time))
atimecodon_dt = skribodts%>%separate(sample,c('time','assay','rep'))%>%group_by(codon)%>%summarise(codon_dwell_time=mean(codon_dwell_time))


timecodon_tAb%>%left_join(timecodon_dt)%>%group_slice(1)%>%{
	quicktest(.$abundance,.$codon_dwell_time)
}

fread('pipeline/testnewsk/codons.csv')%>%
	mutate(sample='E13_ribo_1')%>%
	separate(sample,c('time','assay','rep'))%>%
	left_join(atimecodon_tAb)%>%{
		quicktest(.$abundance,.$codon_dwell_time)
	}

atimecodon_tAb%>%left_join(atimecodon_dt)%>%{
	quicktest(.$abundance,.$codon_dwell_time)
}



timecodon_tAb%>%left_join(timecodon_dt)%>%group_slice(1)%>%{
	quicktest(.$abundance,.$codon_dwell_time)
}



codonprofiles_pcawind%>%inner_join(timecodon_tAb)%>%
	filter(.,time=='E175')%>%{
	quicktest(.$abundance,.$occ_nonorm)
}


codonprofiles_pcawind%>%inner_join(timecodon_dt)%>%
	filter(.,time=='E13')%>%{
	quicktest(.$codon_dwell_time,.$occ_nonorm)
}


codonprofiles_pcawind%>%inner_join(timecodon_dt)%>%
	filter(.,time=='E13')%>%{
	quicktest(.$codon_dwell_time,.$occupancy)
}

codonoccs%>%left_join(timecodon_tAb)

trskribo = fread('pipeline/trtestcodondts.csv')%>%mutate(sample='E13_ribo_1')%>%separate(sample,c('time','assay','rep'))

timecodon_tAb%>%left_join(trskribo)%>%{
	quicktest(.$abundance,.$codon_dwell_time)
}

timecodon_dt%>%left_join(trskribo%>%select(tr_dt=codon_dwell_time,everything()))%>%{
	quicktest(.$codon_dwell_time,.$tr_dt)
}
