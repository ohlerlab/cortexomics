bamtbl=mainbamtbls[1]
if(!file.exists(here('data/sh_fpcovlist.rds'))){


shortoffsets <- tibble(
	offset=c(8),
	compartment='nucl',
	length=21,
)%>%mutate(readlen=paste0('rl',length))

trna_ab_df_samp = allcodsigmean_isomerge[c("fraction", "time", "sample", "anticodon", "abundance", "codon",
"weightedusage", "availability","rep")]

codon_data <- trna_ab_df_samp%>%
	filter(fraction=='Total')%>%select(-fraction,-sample)%>%
    select(time,rep,codon,abundance,availability)%>%
    group_by(time,codon)%>%
    summarise_at(vars(one_of(c('abundance','availability'))),list(mean))

sh_codprof%>%
	filter(position> -as.numeric(readlen)+6,position < -6)%>%
	separate(sample,c('time','assay','rep'))%>%
	left_join(codon_data)%>%
	group_by(position,time,rep)%>%
	nest%>%
	mutate(trnacor=map(data,~tidy(cor.test(.$abundance,.$count))))%>%
	unnest(trnacor)%>%
	arrange(p.value)%>%
	arrange(position)%>%.$estimate%>%txtplot




