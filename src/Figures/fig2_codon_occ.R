
library(conflicted)
library(rlang)
conflict_prefer("last", "dplyr")
conflict_prefer('setdiff','BiocGenerics')
conflict_prefer('intersect','BiocGenerics')
aatable <- read_tsv('ext_data/aa_properties.tsv')

################################################################################
########Plot occupancies at codon level
################################################################################



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
			spread(position,signal)%>%
			{set_rownames(.[,-1],.$codon)}%>%
			princomp%>%{.$loadings[,1]}%>%{./.[which.max(abs(.))]}%>%enframe('position','pca1')
		)
	)
profvarpca%<>%select(sample,readlen,position,pca1)

#
plotfile<-'plots/figures/figure2/trna_codons/fppos_vs_codon_pcascore.pdf'
offsets%<>%mutate(readlen=paste0('rl',length))
pdf(plotfile,w=12,h=12)
profvarpca%>%slice_by(sample,c(1,2,3,4,5,6))%>%
ggplot(data=.,aes(y=pca1,x=as.numeric(position)))+geom_point()+
	facet_grid(readlen~sample)+
		geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=offsets,aes(xintercept= -offset-5),color=I('green'),linetype=2)
dev.off()
normalizePath(plotfile)


# offsets%<>%mutate(readlen=paste0('rl',length))
# pdf('plots/figures/figure2/trna_codons/fppos_vs_codon_variance.pdf',w=12,h=12)
# #plotting variance amongst codons at each point.
# codonprofiles%>%
# 	ungroup%>%
# 	group_by(sample,readlen,position)%>%
# 	# filter(signal==0)%>%.$codon%>%unique
# 	# filter(signal!=0)%>%
# 	# filter(position==1)%>%
# 	filter(!is.nan(signal))%>%
# 	mutate(signal=occ_nonorm)%>%
# 	summarise(sdsig=sd(signal,na.rm=T)/mean(signal,na.rm=T))%>%
# 	separate(sample,c('time','assay','rep'))%>%
# 	group_by(readlen,time,assay,position)%>%
# 	summarise(sdsig=mean(sdsig))%>%
# 	mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
# 	filter(position> -numreadlen+6,position < -6)%>%
# 	filter(numreadlen>=25,numreadlen<=31)%>%
# 	arrange(position)%>%
# 	{
# 		qplot(data=.,x=position,y=sdsig)+
# 		theme_bw()+
# 		facet_grid(readlen~time)+
# 		scale_y_continuous('between codon variation (meannorm)')+
# 		scale_x_continuous('5 read position relative to codon ')+
# 		geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
# 		geom_vline(data=offsets,aes(xintercept= -offset-5),color=I('green'),linetype=2)+
# 		ggtitle("variance of 5' read occurance vs position")
# 	}%>%print
# dev.off()
# normalizePath('plots/figures/figure2/trna_codons/fppos_vs_codon_variance.pdf')

samplestage <- unique(codonprofiles$sample)%>%{setNames(stageconv[str_replace(.,'_.*?_.*?$','')],.)}
pdf('plots/figures/figure2/trna_codons/codonprof_rlindiv.pdf',w=24,h=14)
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
		geom_vline(data=offsets,aes(xintercept= -offset-5),color=I('green'),linetype=2)+
		# geom_rect(data=offsets,color=I('black'),alpha = I(0.1),aes(x=NULL,y=NULL,xmin= -8-offset, xmax = -offset, ymin = 0 ,ymax = Inf))+
		theme_bw())}
dev.off()
normalizePath('plots/figures/figure2/trna_codons/codonprof_rlindiv.pdf')

pdf('plots/figures/figure2/trna_codons/offsetwindow_rlmerge_codprofs.pdf',w=24,h=14)
codonprofiles%>%
	inner_join(offsets)%>%
	mutate(position = position+offset)%>%
	group_by(sample,codon,position)%>%summarise(signal=sum(signal))%>%
		filter(codon%>%str_detect(c('Glu-TTC|GTC|Val-CAC|AAC|ATG')))%>%
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
normalizePath('plots/figures/figure2/trna_codons/offsetwindow_rlmerge_codprofs.pdf')

# profilemats_unproc%>%saveRDS('data/fprofilemats_unproc.rds')

# profilemats <- profilemats_unproc%>%setNames(.,bams[1:length(.)])

codonoccs<-codonprofiles%>%inner_join(offsets)%>%
	filter(sample%>%str_detect('ribo'))%>%
	filter(position <= -offset-3, position >= -offset-5)%>%
	# filter(position <= -offset, position >= -offset-5)%>%
	# filter(position <= -offset-3, position >= -offset-5)%>%
	# filter(readlen%in%c('rl29'))%>%
	ungroup%>%
	mutate(time=str_extract(sample,'[^_]+'))%>%
	group_by(time,codon)%>%
	# group_slice(2)
	summarise(occupancy=mean(signal,na.rm=T))

codonoccs$time%>%unique%>%qs('S5')

################################################################################
########Now combine the tRNA and codon occupancy info
################################################################################ head(allcodsigmean_isomerge)

#Now merge with the tRNA data
tRNA_occ_df<-allcodsigmean_isomerge%>%
	ungroup%>%
	left_join(codonoccs)
#
tRNA_occ_df%<>%filter(is.finite(signal))
##Get the amino acid for each one
tRNA_occ_df%<>%mutate(AA = GENETIC_CODE[codon]%>%qs('S+'))
aatrna_occ_df<-tRNA_occ_df%>%
	group_by(time,AA)%>%
	summarise(occupancy = mean(occupancy),signal_aa = log2(sum(2^(-signal),na.rm=T)))
#
tRNA_occ_df%<>%	left_join(usedcodonfreqs%>%colSums%>%enframe('codon','freq'))%>%
	ungroup%>%
	mutate(common = freq>median(freq))

aaanova <- anova(
	tRNA_occ_df%>%lm(data=.,occupancy ~ 0+ time + AA)
)
library(txtplot)

anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + signal))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + codon))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA+ codon))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon+signal))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Poly'), occupancy ~ 1 + time + AA + codon+signal))

anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon+signal))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon+shared_signal))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon+availablity))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon+availablity_noshare))


anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon + signal))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + signal + codon))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + shared_signal+codon))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + availablity_noshare+codon))


anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Poly'), occupancy ~ 1 + time + AA + codon+signal))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Poly'), occupancy ~ 1 + time + AA + codon+availablity))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Poly'), occupancy ~ 1 + time + AA + availablity+codon))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Poly'), occupancy ~ 1 + time + AA + signal+codon))

anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Poly'), occupancy ~ 1 + time + AA + availablity + codon))



anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Poly'), occupancy ~ 1 + time + AA + codon + shared_signal))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Poly'), occupancy ~ 1 + time + AA + codon+availablity_noshare))


anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Poly'), occupancy ~ 1 + time + AA + codon))


lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1+time+AA+codon)

lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1+time+AA+signal)%>%anova
lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1+time+AA+shared_signal)%>%anova
lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1+time+AA+availablity_noshare)%>%anova
lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1+time+AA+availablity)%>%anova

anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + availablity+codon))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + signal+codon))

anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon+availablity))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon+signal))

anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon+signal+availablity))

anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + signal+availablity+codon))

anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon + availablity))

#if this small effect of availiabily is due to colinearity between codon and avail, then we should get 
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon + availablity))


anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon + signal))
anova(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + signal + codon))

summary(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + codon + availablity))
summary(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total'), occupancy ~ 1 + time + AA + availablity + codon))





tottRNA_occ_df<-tRNA_occ_df%>%safe_filter(fraction=='Total')
tottRNA_occ_df$timeaactlocc<-lm(data=tottRNA_occ_df, occupancy ~ 1 + time + AA)$residuals

tottRNA_occ_df%>%group_by(codon)%>%summarise(timeaactlocc=mean(timeaactlocc),signal=mean(signal))%>%
	{cor.test(.$timeaactlocc,.$signal)}

#
model_aa_time<-
	tRNA_occ_df%>%
	lm(data=.,occupancy ~ 1 + AA )
#
model_aacod_time<-
	tRNA_occ_df%>%
	lm(data=.,occupancy ~ 1  + time + AA + codon)
sigmodel_aacod_time<-
	tRNA_occ_df%>%
	lm(data=.,signal ~ 1  + time + AA + codon)


cor.test(use='complete',model_aacod_time$coef%>%.[names(.)%>%str_detect('codon')],
sigmodel_aacod_time$coef%>%.[names(.)%>%str_detect('codon')])


#do those tRNAs with more average occ, tend to have more tRNAsignal
tRNA_occ_df%>%group_by(fraction,time)%>%summarise(mean(signal))

#
tRNA_occ_df$codresid <- model_aa_time$residuals%>%
	{r=rep(0,as.numeric(last(names(.))));r[as.numeric(names(.))]<-.;r}%>%
	qs('R+')
#
tmeantRNA_occ_df <- tRNA_occ_df%>%
	filter(fraction=='Total')%>%
	group_by(codon)%>%summarise(
		signal = mean(signal),
		occupancy=mean(occupancy),
		availablity=mean(availablity),
		codresid=mean(codresid),
		shared_signal=mean(shared_signal),
		availablity_noshare=mean(availablity_noshare)
	)

model_aacod_time$coef%>%.[names(.)%>%str_detect('codon')]%>%enframe('codon','codoncoef')%>%
	mutate(codon=str_replace(codon,'codon',''))%>%
	left_join(tmeantRNA_occ_df)%>%
	{cor.test(.$codresid,.$availablity,use='complete')}

#
tRNA_occ_df$timerelabund <- tRNA_occ_df%>%lm(data=.,signal ~ time)%>%.$residuals%>%
	{r=rep(0,as.numeric(last(names(.))));r[as.numeric(names(.))]<-.;r}%>%
	qs('R+')

#
normtRNA_occ_df<-tRNA_occ_df%>%
	group_by(codon)%>%
	# filter(!is.na(occupancy))%>%
	# head(2)%>%
	# mutate(signal = signal-mean(signal),occupancy=occupancy-mean(occupancy))	
	# summarise(sum(time=='E13'))
	mutate(
		signal = signal-mean(signal),
		occupancy=occupancy-mean(occupancy),
		codresid=codresid-mean(codresid),
		timerelabund = timerelabund - mean(timerelabund)
	)%>%
	identity

#

#do occupancy and signal in general correalate
tmeantRNA_occ_df%>%{cor.test(.$codresid,.$availablity_noshare,use='complete',method='pearson')}
tmeantRNA_occ_df%>%{cor.test(.$codresid,.$signal,use='complete',method='pearson')}
tmeantRNA_occ_df%>%{cor.test(.$codresid,.$shared_signal,use='complete',method='pearson')}
tmeantRNA_occ_df%>%{cor.test(.$codresid,.$availablity,use='complete',method='pearson')}

#do occupancy and signal in general correalate
tmeantRNA_occ_df%>%{cor.test(.$occupancy,.$availablity_noshare,use='complete',method='pearson')}
tmeantRNA_occ_df%>%{cor.test(.$occupancy,.$signal,use='complete',method='pearson')}
tmeantRNA_occ_df%>%{cor.test(.$occupancy,.$shared_signal,use='complete',method='pearson')}
tmeantRNA_occ_df%>%{cor.test(.$occupancy,.$availablity,use='complete',method='pearson')}



#within time points
tRNA_occ_df%>%
	# filter(codon%>%str_detect(c('TTC|GTC|CAC|AAC|ATG')))%>%
	# filter(time%>%str_detect("E13"))%>%
	filter(abs(codresid)>0.0001)%>%
	group_by(time)%>%
	summarise(cor = cor(signal,codresid),pval=cor.test(signal,codresid)$p.value)

normtRNA_occ_df%>%
	# filter(codon%>%str_detect(c('TTC|GTC|CAC|AAC|ATG')))%>%
	# filter(time%>%str_detect("E13"))%>%
	filter(abs(codresid)>0.0001)%>%
	ungroup%>%
	summarise(cor = cor(signal,codresid),pval=cor.test(signal,codresid)$p.value)

normtRNA_occ_df%>%
	# filter(codon%>%str_detect(c('TTC|GTC|CAC|AAC|ATG')))%>%
	# filter(time%>%str_detect("E13"))%>%
	filter(abs(codresid)>0.0001)%>%
	ungroup%>%
	summarise(cor = cor(timerelabund,codresid),pval=cor.test(timerelabund,codresid)$p.value)


#plot with the AAs as colors
pdf('plots/figures/figure2/trna_codons/percodon_tRNAab_Riboocc.pdf',w=12,h=12)
tRNA_occ_df%>%
	# filter(codon%>%str_detect(c('TTC|GTC|CAC|AAC|ATG')))%>%
	# filter(sample%>%str_detect("Poly"))%>%
	# filter(abs(codresid)>0.0001)%>%
	mutate(timefrac = paste0(time,'_',fraction))%>%
	group_by(fraction,time)%>%
	nest%>%
	mutate(labl = paste0(
		'r = ',map_dbl(data,~cor(method='spearman',use='complete',.$occupancy,.$availablity_noshare))%>%round(3),'\n',
		'p = ',map_dbl(data,~cor.test(method='spearman',use='complete',.$occupancy,.$availablity_noshare)$p.value%>%round(3))
	))%>%
	unnest%>%
	{
	#
	labelpos = group_by(.,codon)%>%summarise(occupancy=mean(occupancy),availablity_noshare=mean(availablity_noshare))%>%
		mutate(AA = GENETIC_CODE[codon]%>%qs('S+'))
	#
	# filter(codon%>%str_detect("ATG"))%>%
	ggplot(filter(.,(codon%>%table(.)[.]) > 1),aes(x=availablity_noshare,y=occupancy,color=time))+
	geom_point(size=1)+
	# geom_line(color=I('grey'))+
	facet_wrap(nrow=3,timefrac ~ . )+
	scale_color_manual(values=stagecols)+
	geom_text(show.legend=F,data=labelpos,aes(label=codon,color=NULL),size=4)+
	geom_text(data=distinct(.,timefrac,time,fraction,labl),hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))+
	# geom_text(show.legend=F,data=labelpos,aes(label=AA,color=NULL),size=9)+
	# geom_line(aes(group=qs(GENETIC_CODE[codon],'S+'),color=NULL),size=1,linetype=2)+
	geom_smooth(method='lm')+
	theme_bw()+
	scale_x_continuous(name='Summed tRNA Expression')
}
dev.off()
normalizePath('plots/figures/figure2/trna_codons/percodon_tRNAab_Riboocc.pdf')


tRNA_occ_df%>%filter(AA=='D')

allcodsigmean%>%filter(AA=='D')
codonoccs$time%>%table

#plot with the AAs as colors
pdf('plots/figures/figure2/trna_codons/stripplot_aa_codon.pdf',w=12,h=5)
codonoccs%>%
	mutate(AA = GENETIC_CODE[codon])%>%
	# filter(fraction=='Total')%>%
	mutate(codonname = paste0(AA,'-',codon))%>%
	group_by(codonname)%>%mutate(codmean=mean(occupancy))%>%group_by(AA)%>%mutate(aamean=mean(occupancy))%>%
	ungroup()%>%arrange(aamean,codmean)%>%
	mutate(codonname=as_factor(codonname),AA=as_factor(AA),codonname=as_factor(codonname))%>%
	{
	ggplot(filter(.,(codon%>%table(.)[.]) > 1),aes(y=occupancy,x=codonname,color=time))+
	geom_point(size=1)+
	# geom_line(color=I('grey'))+
	scale_color_manual(values=stagecols)+
	facet_grid(.~AA,scale='free_x')+
	# geom_rect(xmax='R-AGG',xmin='R-CGT',ymin=-Inf,ymax=Inf,alpha=0.2)+
	theme_bw()+
	theme(axis.text.x=element_text(angle=45,vjust=.5))+
	scale_y_continuous(name='Normalized Occupancy Over Codon A site')
}
dev.off()
normalizePath('plots/figures/figure2/trna_codons/stripplot_aa_codon.pdf')


codonoccs$time%>%table


##Plot of change over time
pdf('plots/figures/figure2/trna_codons/percodon_AAreg_tRNAab_Riboocc.pdf',w=7,h=7)
normtRNA_occ_df%>%
	# safe_filter(sample%>%str_detect('Poly'))%>%
	left_join(usedcodonfreqs%>%colSums%>%enframe('codon','freq'))%>%
	ungroup%>%
	mutate(common = ifelse(freq>median(freq),'common','rare'))%>%
	{
	#
	labelpos = group_by(.,codon,fraction)%>%summarise(codresid=mean(codresid),signal=mean(signal))
	#
	# filter(codon%>%str_detect("ATG"))%>%
	ggplot(.,aes(x=signal,y=codresid,size=freq,color=time,group=codon))+
	geom_point(aes(size=freq))+
	# geom_line(color=I('grey'))+
	# facet_wrap(nrow=3,time ~ . )+
	scale_color_manual(values=stagecols)+
	facet_grid(common~fraction)+
	# geom_text(show.legend=F,data=labelpos,aes(label=codon,color=NULL),size=9)+
	scale_x_continuous(name='Temporal Change in Summed tRNA Expression')+
	scale_y_continuous(name='Temporal change in occupancy - AA effect')+
	theme_bw()
}
dev.off()
normalizePath('plots/figures/figure2/trna_codons/percodon_AAreg_tRNAab_Riboocc.pdf')


(lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + time + AA + signal + codon ))%T>%{anova(.)%>%print}



lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + time + AA + signal + codon )%T>%{anova(.)%>%print}

lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + time + AA  )%>%summary

lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + time + AA  + codon+availablity )%>%summary

lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + time + AA  + codon+signal )%>%summary

lm(data=
	tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), 
	occupancy ~ 1 + time + AA  + codon+signal )%>%confint

lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + time + AA  + signal+codon )%>%anova

lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total',time==unique(time)[1:5])%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + signal )%>%anova
lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total',time==unique(time)[4])%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + availablity )%>%anova
lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total',time==unique(time)[1])%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + availablity )%>%anova

lm(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + time + codon+  availablity )%>%anova

aov(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + time + signal+codon )%>%summary
aov(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + time + availablity+codon )%>%summary
aov(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + time + codon +availablity )%>%summary
aov(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + time + codon +availablity )%>%summary

aov(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + time +codon + signal )%>%confint

aov(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + signal + codon  )%>%summary
aov(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  +  codon + signal  )%>%summary

aov(data=tRNA_occ_df%>%safe_filter(fraction=='Total',time==unique(time)[5])%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + signal   )%>%summary


aov(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + time + codon +availablity )%>%summary

aov(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + time + availablity + codon )%>%summary

aov(data=tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal)), occupancy ~ 1 + AA  + time + availablity  )%>%summary

library(lme4)


tRNA_occ_df%>%safe_filter(fraction=='Total')%>%group_by(codon)%>%mutate(msig=mean(signal))%>%
{aov(lm(data=., signal ~ 1 + AA  + time)$residuals ~ .$codon)}


tot_tRNA_occ_df <- tRNA_occ_df %>%safe_filter(fraction=='Poly')%>%group_by(codon)%>%mutate(msig=mean(signal))
#So this is significant - signal matters before we control for codon.
aov(data=tot_tRNA_occ_df, occupancy ~ 1 + AA  + time + signal)%>%summary
#I would have thought this means that the mean levels for each codon were significantly
tot_tRNA_occ_df$prctlresids <- aov(data=tot_tRNA_occ_df, occupancy ~ 1 + AA  + time )$residuals
#Why is this not true.
tot_tRNA_occ_df%>%group_by(codon)%>%summarise(mr = mean(prctlresids),signal=mean(signal))%>%{cor.test(use='complete',.$mr,.$signal)}
#What if we look only at AAs with multiple codons?
tot_tRNA_occ_df%>%group_by(AA)%>%filter(n_distinct(codon)>1)%>%group_by(codon)%>%summarise(mr = mean(prctlresids),signal=mean(signal))%>%{cor.test(use='complete',.$mr,.$signal)}

tot_tRNA_occ_df%>%group_by(codon)%>%summarise(mr = mean(prctlresids),signal=mean(signal))%>%{txtplot(.$mr,.$signal)}
tot_tRNA_occ_df%>%group_by(AA)%>%filter(n_distinct(codon)>1)%>%group_by(codon)%>%summarise(mr = mean(prctlresids),signal=mean(signal))%>%{txtplot(.$mr,.$signal)}

#What if we look only at AAs with multiple codons?

tot_tRNA_occ_df%>%
	group_by(AA)%>%filter(n_distinct(codon)>1)%>%
	group_by(AA,codon)%>%summarise(mr = mean(occupancy),availablity=mean(availablity))%>%
	group_by(AA)%>%mutate(availablity = availablity - mean(availablity))%>%
	group_by(AA)%>%mutate(mr = mr - mean(mr))%>%
	{cor.test(use='complete',.$mr,.$availablity)}

tot_tRNA_occ_df%>%
	group_by(AA)%>%filter(n_distinct(codon)>1)%>%
	group_by(AA,codon)%>%summarise(mr = mean(occupancy),availablity=mean(availablity))%>%
	group_by(AA)%>%mutate(availablity = availablity - mean(availablity))%>%
	group_by(AA)%>%mutate(mr = mr - mean(mr))%>%
	{txtplot(.$mr,.$availablity)}



#total occ, row signal
for(fractioni in unique(tRNA_occ_df$fraction)){
	for(sigcol in syms(c('availablity','shared_signal','availablity_noshare','signal','weightedusage'))){
		tRNA_occ_df%>%
		filter(fraction==fractioni)%>%
		group_by(AA)%>%filter(n_distinct(codon)>1)%>%
		group_by(AA,codon)%>%summarise(mr = mean(occupancy),signal=mean(!!sigcol))%>%
		# group_by(AA)%>%mutate(signal = signal - mean(signal))%>%
		group_by(AA)%>%mutate(mr = mr - mean(mr))%>%
		identity%>%
	# {txtplot(.$mr,.$availablity)}
		{message(sigcol);message(fractioni);.}%>%
		{cor.test(use='complete',.$mr,.$signal)%>%print}%>%
		.$p_value
	}
}





# # codonstbl<-read.table('https://raw.githubusercontent.com/zhanxw/anno/master/codon.txt')
# # codonstbl%>%write_tsv('ext_data/codons.txt')
# codonstbl<-read_tsv('ext_data/codons.txt')
# #segregate the 
# codonnames <- paste0(codonstbl[[2]],'-',codonstbl[[1]])[match(names(codons2scan),codonstbl[[1]])]

# profilemats%<>%map_df(.id='bam',.%>%map_df(.id='codon',.%>%enframe('position','signal')))

# bamsamples<-bams%>%basename%>%str_remove('.bam')

# sampleparams<-'src/sample_parameter.csv'%>%fread

# profilemats%>%head
# profilemats$sample = profilemats$bam%>%match(bams)%>%bamsamples[.]
# profilemats$time<-sampleparams$time[match(profilemats$sample,sampleparams$sample_id)]



# profilemats$codon%<>%as.numeric
# profilemats$codon<-profilemats$codon%>%codonnames[.]

# #summarise to get rid of the periodicity signal
# codonprofmat<-profilemats%>%mutate(codonpos=ceiling(position/3))%>%group_by(codon,position,bam)

# codonprofmat$codonpos%<>%subtract(FLANKCODS+1)
# codonprofmat%<>%ungroup


# require(plotly)

# highcod = 'AAC'
# lowcod = 'GTC'

# pdf('tmp.pdf')
# codonprofmat%>%
# 	filter(codon%>%str_detect(c('GTC|AAC|ATG')))%>%
# 	# filter(codon%>%str_detect(c('Glu-TTC|GTC|Val-CAC|AAC|ATG')))%>%
# 	filter(sample%>%str_detect(c('ribo')))%>%
# 	{print(ggplot(.,aes(codonpos,signal,group=bam,
# 		color=samplestage[sample]))+
# 		scale_color_manual(values=displaystagecols)+
# 		facet_grid(codon~.)+
# 		geom_line(aes())+
# 		theme_bw())}
# dev.off()
# normalizePath('tmp.pdf')




# save.image('data/codon_coverage.Rdata')
# load('data/codon_coverage.Rdata')


# #' Convert from Rle to one column matrix
# #'
# Q
# setAs("Rle", "Matrix", function(from) {
#     rv <- runValue(from)
#     nz <- rv != 0
#     i <- which(as.vector(from !=0))
#     x <- rep(rv[nz], runLength(from)[nz])

# 	length(i)    

#     sparseMatrix(i= i, p=c(0L, length(x)), x=x,
#                  dims=c(length(from), 1))
# })

# #' Convert from DataFrame of Rle to sparse Matrix
# #'
# setAs("DataFrame", "Matrix", function(from) {
#   mat = do.call(cbind, lapply(from, as, "Matrix"))
#   colnames(mat) <- colnames(from)
#   rownames(mat) <- rownames(from)
#   mat
# })







pdf('plots/figures/figure2/trna_codons/tRNA_abchange_genomeusage.pdf')
usage_v_abundance_df%>%
	left_join(usedcodonfreqs%>%colSums%>%{./sum(.)}%>%enframe('codon','freq'))%>%
	filter(!codon%in%c('TAG','TAA','TGA'))%>%
	{
		print(
			qplot(data=.,x=.$abundancechange,y=.$freq,label=codon,geom='blank')+
				scale_x_continuous('tRNA Abundance Signal Slope E13 - P0')+
				scale_y_continuous('Genome Usage')+
				geom_smooth(method='lm')+
				# facet_grid(tRNA_time~.)+
				geom_text()+
				geom_text(data=wus_tab_cors,aes(x=0,y=0,label=paste0('r = ',round(cor,2))))+
				theme_bw()
			)
	}
dev.off()
normalizePath('plots/figures/figure2/trna_codons/tRNA_abchange_genomeusage.pdf')




codons4occ <- tRNA_occ_df$codon%>%unique%>%setdiff(c('TAG','TAA','TGA'))



##Score for cds at particular time points
times<-tRNA_occ_df$time%>%unique%>%setNames(.,.)
itime <- times[1]

allTEchangedf<-read_tsv('tables/manuscript/go_all_highcount_updown.tsv')

timeoccscores <- map_df(.id='time',times,function(itime){
	# codons4occ='CTC'
	timeoccvect <- tRNA_occ_df%>%filter(time==itime)%>%mutate(occupancy=occupancy-mean(na.rm=T,occupancy))%>%{setNames(.$occupancy,.$codon)[codons4occ]}
	(t(usedcodonfreqs[,codons4occ])*(timeoccvect[codons4occ]))%>%colMeans(na.rm=T)%>%
	enframe('protein_id','elongscore')
})

timeSigscores <- map_df(.id='time',times,function(itime){
	# codons4occ='CTC'
	timetrrnavect <- tRNA_occ_df%>%filter(time==itime)%>%mutate(signal=signal-mean(na.rm=T,signal))%>%{setNames(.$signal,.$codon)[codons4occ]}
	(t(usedcodonfreqs[,codons4occ])*(timetrrnavect[codons4occ]))%>%colMeans(na.rm=T)%>%
	enframe('protein_id','tRNA_expr_score')
})

score_techange_df<-timeoccscores%>%left_join(timeSigscores)%>%left_join(ms_id2protein_id%>%distinct(protein_id,gene_name))%>%
	left_join(allTEchangedf)

occchange_vs_te_df <- score_techange_df%>%
	group_by(up,down,protein_id)%>%
	# filter(time!='P0')%>%
	# group_slice(1:2)%>%
	nest%>%
	mutate(elongchange = map_dbl(data,~{lm(data=.,elongscore ~ seq_along(time))$coef[2]}))%>%
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

pdf('plots/figures/figure2/trna_codons/te_change_vs_occscorechange.pdf')
occchange_vs_te_df%>%
	filter(!(up&down))%>%
	mutate(TEchange_class = case_when(
		up==1 ~ 'up',
		down==1 ~ 'down',
		TRUE ~ 'neither'
	))%>%
	ggplot(.,aes(x=elongchange,color=TEchange_class))+geom_density(alpha=I(0.01))+theme_bw()+
	scale_x_continuous('Predicted Elongation Rate Change - Riboseq Occ')
dev.off()
normalizePath('plots/figures/figure2/trna_codons/te_change_vs_occscorechange.pdf')

occchange_vs_te_df%>%filter(!down)%>%{split(.$elongchange,.$up)}%>%{t.test(.[[1]],.[[2]])}
occchange_vs_te_df%>%filter(!up)%>%{split(.$elongchange,.$down)}%>%{t.test(.[[1]],.[[2]])}


pdf('plots/figures/figure2/trna_codons/te_change_vs_tRNAscorechange.pdf')
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
normalizePath('plots/figures/figure2/trna_codons/te_change_vs_tRNAscorechange.pdf')

tRNAchange_vs_te_df%>%filter(!down)%>%{split(.$tRNA_score_change,.$up)}%>%{t.test(.[[1]],.[[2]])}
tRNAchange_vs_te_df%>%filter(!up)%>%{split(.$tRNA_score_change,.$down)}%>%{t.test(.[[1]],.[[2]])}

###So is this the same genes???
tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%{txtplot(.$tRNA_score_change,.$elongchange)}
tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%{cor.test(.$tRNA_score_change,.$elongchange)}


#occ scores of things that change TE
pdf('plots/figures/figure2/trna_codons/abs_elongscore_vs_tRNAscorechange.pdf')
score_techange_df%>%
	# filter(time=='E13')%>%
	filter(!(up&down))%>%
	mutate(TEchange_class = case_when(
		up==1 ~ 'up',
		down==1 ~ 'down',
		TRUE ~ 'neither'
	))%>%
	ggplot(.,aes(x=log(elongscore),color=TEchange_class,fill=TEchange_class))+
	geom_density(alpha=I(0.01))+
	theme_bw()+
	scale_x_continuous(limits=c(-8,2))+
	facet_grid(time~.)
dev.off()
normalizePath('plots/figures/figure2/trna_codons/abs_elongscore_vs_tRNAscorechange.pdf')


###Do optimized looking CDS tend to have higher TE?
abste_opt_df<-bestmscountebayes$coef[,'TE']%>%enframe('uprotein_id','TE')%>%mutate(protein_id=uprotein_id%>%str_replace('_\\d+$',''))%>%
	left_join(pc_optcodons)


cor.test(abste_opt_df$TE,pc_optcodons$pc_opt)
txtplot(abste_opt_df$TE,pc_optcodons$pc_opt,width=100,pch='.')

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



#Do 
#Are occupancy and global frequency correlated?
tRNA_occ_df%>%
	left_join(enframe(overallcodonfreqs,'codon','freq'))%>%
	split(.,.$time)%>%
	map(~ cor.test(.$freq,.$signal))

#Do 
#Are occupancy and global frequency correlated?
tRNA_occ_df%>%
	left_join(enframe(usedcodonfreqs%>%colSums,'codon','freq'))%>%
	split(.,.$time)%>%
	map(~ cor.test(.$freq,.$signal))


#Are occupanc
################################################################################
########Rank order plots of codons occupancy
################################################################################

#Are occupancy and global signal correlated?
pdf('plots/figures/figure2/trna_codons/occupancy_rankplot_bystage_poly.pdf',h=12,w=12)
tRNA_occ_df%>%
	filter(sample%>%str_detect('Poly'))%>%
	identity%T>%checkDataFrame%>%
	left_join(enframe(usedcodonfreqs%>%colSums,'codon','freq'))%>%
	left_join(weighted_codon_usage)%>%
	gather(freqtype,freq,freq,weightedusage)%>%
	group_by(time,freq)%>%
	inner_join(aatable)%>%
	filter(AA!='*')%>%
	split(.,list(.$time,.$freqtype))%>%
	lapply(function(x){
		x%>%arrange(occupancy/signal)%>%
		mutate(codon=as_factor(codon))%>%
		{ggplot(.,aes(x=codon,color=paste0(hydro,'_',pol),y=occupancy/signal))+
		geom_point()+
		scale_x_discrete(limits = .$codon,labels=paste0(.$codon,'-',.$AA))+
		# facet_grid(freqtype~time,scales='free')+
		scale_color_manual(values=c('mod_+'='lightblue','mod_N'='lightgrey','phil_+'='blue','phil_N'='grey','phob_N'='black','phil_-'='red'))+
		theme(axis.text.x=element_text(angle=45))+coord_flip()+
		scale_y_continuous(name=paste0('occupancy/signal'))+
		ggtitle(unique(.$sample))+theme_bw()
	}
	})%>%{ggpubr::ggarrange(plotlist=.,ncol=3,nrow=2)}
dev.off()
normalizePath('plots/figures/figure2/trna_codons/occupancy_rankplot_bystage_poly.pdf')
#
#Are occupancy and global signal correlated?
pdf('plots/figures/figure2/trna_codons/occupancy_rankplot_bystage_total.pdf',h=12,w=18)
tRNA_occ_df%>%
	safe_filter(sample%>%str_detect('Total'))%>%
	left_join(enframe(usedcodonfreqs%>%colSums,'codon','freq'))%>%
	left_join(weighted_codon_usage)%>%
	gather(freqtype,freq,freq,weightedusage)%>%
	group_by(time,freq)%>%
	inner_join(aatable)%>%
	filter(AA!='*')%>%
	split(.,list(.$time,.$freqtype))%>%
	lapply(function(x){
		x%>%arrange(occupancy/signal)%>%
		mutate(codon=as_factor(codon))%>%
		{ggplot(.,aes(x=codon,color=paste0(hydro,'_',pol),y=occupancy/signal))+
		geom_point()+
		scale_x_discrete(limits = .$codon,labels=paste0(.$codon,'-',.$AA))+
		# facet_grid(freqtype~time,scales='free')+
		scale_color_manual(values=c('mod_+'='lightblue','mod_N'='lightgrey','phil_+'='blue','phil_N'='grey','phob_N'='black','phil_-'='red'))+
		theme(axis.text.x=element_text(angle=45))+coord_flip()+
		scale_y_continuous(name=paste0('occupancy/signal'))+
		ggtitle(unique(.$sample))+theme_bw()
	}
	})%>%{ggpubr::ggarrange(plotlist=.,ncol=5,nrow=2)}
dev.off()
normalizePath('plots/figures/figure2/trna_codons/occupancy_rankplot_bystage_total.pdf')

aatable%>%filter(pol=='-')%>%
	left_join(allcodsigmean%>% mutate(AA = GENETIC_CODE[str_extract(codon,'[^\\-]+$')]%>%qs('S+')))


################################################################################
########Test if the protein linear effects back up the hypothesis
########That the elongation rate driven 
################################################################################

elonguppids<-occchange_vs_te_df%>%filter(elongchange > 0.2,up)%>%.$protein_id
noelongupids<-occchange_vs_te_df%>%filter(abs(elongchange) < 0.5,up)%>%.$protein_id

timechange <- contrasts.fit(bestmscountebayes,cbind(alltimeeff[,-1],timeTEeffect[,-1]))$coef%>%{rownames(.)%<>%str_replace('_\\d+','');.}
MStimechange <- contrasts.fit(bestmscountebayes,cbind(timeMSeffect[,]))$coef%>%{rownames(.)%<>%str_replace('_\\d+','');.}
bestrownames <- MStimechange%>%rownames
elonguppids%<>%intersect(bestrownames)
noelongupids%<>%intersect(bestrownames)

library(matrixStats)
protein_dev_classes <- MStimechange%>%{setNames(rowVars(.)%>%sqrt,rownames(.))}%>%enframe('protein_id','MS_timevar')%>%
	mutate(
		set=case_when(
			protein_id%in%elonguppids ~ 'TE Up, Elongation Score Change > 0.2\n(Getting Slower)',
			protein_id%in%noelongupids ~ 'TE Up, Elongation Score Change < 0.05 ',
			protein_id%in%(occchange_vs_te_df%>%filter(up==1)%>%.$protein_id) ~ 'TE Up, other',
			TRUE ~ 'NA'
		)
)

pdf('plots/figures/figure2/trna_codons/protein_dev_teclasssse_density.pdf')
protein_dev_classes%>%
	filter(set!='NA',set!='TE Up, other')%>%
	{qplot(data=.,x=MS_timevar,fill=set,geom='blank')+
	# scale_x_log10()+
	geom_density()+
	geom_histogram()+
	theme_bw()+facet_grid(set ~ . ,scale='free')}%>%print
dev.off()
normalizePath('plots/figures/figure2/trna_codons/protein_dev_teclasssse_density.pdf')


pdf('plots/figures/figure2/trna_codons/protein_ddev_scatter_TEclassees.pdfs')
MStimechange[,]%>%{setNames(rowVars(.)%>%sqrt,rownames(.))}%>%enframe('protein_id','MS_timevar')%>%left_join(occchange_vs_te_df)%>%
	mutate(
		set=case_when(
			protein_id%in%elonguppids ~ 'TE Up, Elongation Score Change > 0.2\n(Getting Slower)',
			protein_id%in%noelongupids ~ 'TE Up, Elongation Score Change < 0.05 ',
			protein_id%in%(occchange_vs_te_df%>%filter(up==1)%>%.$protein_id) ~ 'TE Up, other',
			TRUE ~ 'NA'
		))%>%
	filter(up==1)%>%
	qplot(data=.,x=log10(elongchange),color=set,y=log10(MS_timevar),geom='point',alpha=I(0.5))+theme_bw()
dev.off()
normalizePath('plots/figures/figure2/trna_codons/protein_ddev_scatter_TEclassees.pdfs')

protein_dev_classes%>%{split(.$MS_timevar%>%log,.$set)}%>%{wilcox.test(.[["TE Up, Elongation Score Change < 0.05 "]],.[['TE Up, Elongation Score Change > 0.2\n(Getting Slower)']])}

MStimechange

