library(here)
if(!exists('allcodsig_isomerge')) base::source(here('src/Figures/Figure3/3_tRNA_array_analysis.R'))
if(!exists('kl_df')) base::source(here('src/subset_dwell_times.R'))



# #does amino acid variability determine the asite effects?
# codonprofiledat%>%left_join(offsets)%>%filter(position== -offset)%>%
# 	mutate(AA=GENETIC_CODE[codon])%>%
# 	lm(data=.,count~AA)%>%summary
# codonprofiledat%>%left_join(offsets)%>%
# 	filter(between(position, -offset-5,-offset-3))%>%
# 	group_by(sample,codon)%>%summarise(count=sum(count))%>%
# 	mutate(AA=GENETIC_CODE[codon])%>%
# 	lm(data=.,count~AA)%>%summary

{
allcodsig_isomerge%<>%addcodon
trna_ab_df_samp = allcodsig_isomerge[c("fraction", "time", "sample", "anticodon", "abundance", "codon",
"weightedusage", "availability","rep")]
#get data on the trnas
trna_dat <- trna_ab_df_samp%>%
	select(-sample)%>%
    select(fraction,time,rep,codon,abundance,availability)%>%
    group_by(fraction,time,codon)%>%
    summarise_at(vars(one_of(c('abundance','availability'))),list(mean),na.rm=T)
# varoffsets2use<-codonvaroffsets
varoffsets2use<-kl_offsets2plot%>%
	# group_by(readlen,time)%>%slice(1)%>%
	mutate(phase=position%%3)%>%
	# select(offset,readlen,phase,sample)
	select(-position)
	# rloffsets<-offsets%>%select(readlen,length,offset)%>%mutate(readlen=str_replace(length,'^(\\d)','rl\\1'))
# lexp = 3+3 #include positions for positions corresponding to bigger offsets
# rexp = 3#include positions for positions corresponding to smaller offsets
# codonprofiles%>%
precodondata <- 
	# frustprofilelist%>%
	frustprofilelist%>%
    mutate_at('sample',clean_fr_sampnames)%>%
    group_by(sample)
codondata_notrna<-precodondata%>%
	mutate(phase=position%%3)%>%
	# inner_join(varoffsets2use,by=c('sample','readlen','phase'),allow_missing=T)%>%
	inner_join(varoffsets2use,by=c('readlen','phase'),allow_missing=T)%>%
    separate(sample,c('time','assay','rep'))%>%
    # filter(between(position, -offset-5,-offset))%>%
    # group_by(assay,time,codon,rep)%>%
    # filter(readlen%in%c('rl26','rl27','rl28','rl29','rl30'))%>%
    filter(readlen%in%c('rl29'))%>%
    group_by(time,assay,rep,readlen)%>%
    select(time,assay,rep,readlen,phase,codon,count,position,offset)
#
codondata_notrna%>%filter(is.na(offset))
#
#seperate a and p site occ
codondata_notrna%<>%    
	group_by(time,assay,rep,codon)%>%
	summarise(
    	# p_site_occ = sum(count[position== -offset]),
    	# a_site_occ = sum(count[position== -offset-3]),
    	# p_site_occ = sum(count[position== -11]),
    	# a_site_occ = sum(count[position== -11-3]),
    	# p_site_occ = sum(count[position== -offset]),
    	# p_site_occ = sum(count[position%>%between(-offset-2,-offset)]),
    	# a_site_occ = sum(count[position== -offset-3]),
    	# a_site_occ = sum(count[position%>%between(-offset-5,-offset-3)])
    	p_site_occ = sum(count[position%>%between(-offset,-offset)]),
    	a_site_occ = sum(count[position%>%between(-offset-3,-offset-3)]),
    	# a_site_occ = sum(count[position%>%between(-offset-3-3,-offset+1-3-3)]),
    )%>%
    arrange(assay!='ribo')
stopifnot(!is.na(codondata_notrna$a_site_occ))
stopifnot(!is.na(codondata_notrna$p_site_occ))
}
{
#codondatasave<-codondata
#seperate poly and total tRNA info
total_trnadat <- trna_dat%>%ungroup%>%filter(fraction=='Total')%>%select(time,codon,abundance,availability)
poly_trnadat <- trna_dat%>%ungroup%>%filter(fraction=='Poly')%>%
	select(time,codon,poly_abundance=abundance,poly_availability=availability)
#add tRNA info, 
#AA and aa corrected dwell time
codondata<-codondata_notrna%>%
    left_join(total_trnadat,by=c('time','codon'))%>%
    left_join(poly_trnadat,by=c('time','codon'))%>%
    mutate(AA=GENETIC_CODE[codon])%>%
    group_by(assay,time,rep,AA)%>%
    arrange(assay!='ribo')%>%
    mutate(aacor_p_site_occ=ifelse(n()==1,NA,p_site_occ-mean(p_site_occ)))
#summarise the replicates
repsumcodondata<-codondata%>%
	group_by(assay,time,codon,AA)%>%
	summarise_at(vars(abundance,availability,poly_abundance,poly_availability,a_site_occ,p_site_occ,aacor_p_site_occ),mean)%>%
    arrange(assay!='ribo')
#check all codons there
stopifnot(codondata%>%group_by(assay,time,rep)%>%tally%>%.$n%>%`==`(61)%>%all)
stopifnot(repsumcodondata%>%group_by(assay,time)%>%tally%>%.$n%>%`==`(61)%>%all)

# codondata%>%filter(assay=='ribo')%>%group_by(assay,time)%>%group_slice(1)%>%{quicktest(.$p_site_occ,.$availability)}
repsumcodondata%>%filter(assay=='ribo',time=='E13')%>%group_by(assay,time)%>%group_slice(1)%>%{quicktest(.$a_site_occ,.$availability)}
repsumcodondata%>%filter(assay=='ribo',time=='E13')%>%group_by(assay,time)%>%group_slice(1)%>%{quicktest(.$p_site_occ,.$availability)}
}
# codondata%>%filter(assay=='ribo')%>%group_by(assay,time)%>%group_slice(2)%>%{quicktest(.$p_site_occ,.$availability)}
# repsumcodondata%>%filter(assay=='ribo')%>%group_by(assay,time)%>%group_slice(2)%>%{quicktest(.$p_site_occ,.$availability)}
# repsumcodondata%>%filter(assay=='Polyribo')%>%group_by(assay,time)%>%group_slice(1)%>%{quicktest(.$p_site_occ,.$poly_abundance)}
# codondata %>% saveRDS(here('data/codondata.rds'))
#r epsumcodondata %>% saveRDS(here('data/repsumcodondata.rds'))

totrepsumcodondata <- repsumcodondata%>%filter(assay=='ribo')
polyrepsumcodondata <- repsumcodondata%>%filter(assay=='Polyribo')

#now plot
plotfile<- here(paste0('plots/','a_occ_stage_dtdist','.pdf'))
pdf(plotfile)
p1=totrepsumcodondata%>%
	ggplot(.,aes(x=a_site_occ,fill=time))+
	geom_density(alpha=I(0.5))+
	scale_fill_manual(name='stage',values=stagecols)+
	scale_x_continuous(paste0('a_site_occ_'))+
	ggtitle(paste0(' dt_distribution'))+
	theme_bw()
p2=totrepsumcodondata%>%
	ggplot(.,aes(x=a_site_occ,y=time,color=time))+
	geom_point(alpha=I(0.5),position='jitter')+
	scale_color_manual(name='stage',values=stagecols)+
	scale_x_continuous(paste0('a_site_occ_'))+
	ggtitle(paste0(' dt_distribution'))+
	theme_bw()
ggarrange(plotlist=list(p1,p2),nrow=2)
dev.off()
message(normalizePath(plotfile))
#now plot
plotfile<- here(paste0('plots/','p_occ_stage_dtdist','.pdf'))
pdf(plotfile)
p1=totrepsumcodondata%>%
	ggplot(.,aes(x=p_site_occ,fill=time))+
	geom_density(alpha=I(0.5))+
	scale_fill_manual(name='stage',values=stagecols)+
	scale_x_continuous(paste0('p_site_occ_'))+
	ggtitle(paste0(' dt_distribution'))+
	theme_bw()
p2=totrepsumcodondata%>%
	ggplot(.,aes(x=p_site_occ,y=time,color=time))+
	geom_point(alpha=I(0.5),position='jitter')+
	scale_color_manual(name='stage',values=stagecols)+
	scale_x_continuous(paste0('p_site_occ_'))+
	ggtitle(paste0(' dt_distribution'))+
	theme_bw()
ggarrange(plotlist=list(p1,p2),nrow=2)
dev.off()
message(normalizePath(plotfile))

# asiteoccs = subfpprofilelist[['allhigh']]%>%filter(sample%>%str_detect('E13_ribo'),readlen=='rl29')%>%filter(position== -11-3)%>%
	# group_by(sample,codon)%>%summarise_at(vars(count),mean)%>%
	# group_by(codon)%>%summarise_at(vars(count),mean)
# asiteoccs%>%rename('a_site_occ'=count)%>%
	# left_join(psiteoccs%>%rename('p_site_occ'=count),by='codon')%>%
	# mutate(foo=1)%>%
	# make_quantcompplot_fac(a_site_occ,p_site_occ,foo,fname='plots/atimes_adjecency_vs_psite_occ.pdf')

################################################################################
########
################################################################################
	
make_quantcompplot_fac <- function(compdf, col1, col2, facetvar=NULL, fname){
	require(LSD)
	base::source(here('Applications/LSD/R/LSD.heatscatter.R'))
	require(broom)
	col1<-enquo(col1)
	col2<-enquo(col2)
	facetvar<-enquo(facetvar)
	compdf%<>%dplyr::rename(facet=!!facetvar)
	corlabel = compdf%>%
		group_by(facet)%>%
		filter(is.finite(!!col1),is.finite(!!col2))%>%
		summarise(tidy(cor.test(!!col1, !!col2)))
	corlabel = corlabel%>%
		mutate(
			pformat=format(p.value,format='e',digits=4),
			pvalstring = ifelse(p.value > 0.001,round(p.value,4),pformat),
			labl=paste0('rho = ',round(estimate,3),'\n','pval = ',pvalstring))
	#
	nlabel= corlabel%>%
		mutate(labl=paste0('N=',parameter))
	facetnum <- n_distinct(compdf$facet)
	# compdf%<>%mutate(facet=!!facetvar)
	# pdf(fname,h=3,w=5*facetnum)
	# gplot = heatscatter(ggplot=TRUE,
			# compdf[[quo_name(col1)]],compdf[[quo_name(col2)]])+
	gplot=compdf%>%ggplot(aes(x=!!col1,y=!!col2))+
		geom_point()+
		scale_x_continuous(quo_name(col1))+
		scale_y_continuous(quo_name(col2))+
		facet_grid(.~facet)+
		geom_smooth(method='lm')+
		ggtitle(basename(fname))+
		geom_text(show.legend=F,data=corlabel,
			hjust=1,vjust=1,x= Inf,y=Inf,aes(label=labl))+
		geom_text(show.legend=F,data=nlabel,
			hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))+
		theme_bw()
	# dev.off()
	pdf(fname,h=5,w=5*facetnum)
	print(gplot)
	dev.off()
	message(normalizePath(fname))
}

totrepsumcodondata <- repsumcodondata%>%filter(assay=='ribo')
polyrepsumcodondata <- repsumcodondata%>%filter(assay=='Polyribo')
# dir.create('plots/p_a_redux/')
fname= here(paste0('plots/p_a_redux/','a_occ_ab_vs_dt','.pdf'))
totrepsumcodondata%>%make_quantcompplot_fac(abundance,a_site_occ,time,fname)
fname= here(paste0('plots/p_a_redux/','a_occ_av_vs_dt','.pdf'))
totrepsumcodondata%>%make_quantcompplot_fac(availability,a_site_occ,time,fname)
#
fname= here(paste0('plots/p_a_redux/','p_occ_ab_vs_dt','.pdf'))
totrepsumcodondata%>%make_quantcompplot_fac(abundance,p_site_occ,time,fname)
fname= here(paste0('plots/p_a_redux/','p_occ_av_vs_dt','.pdf'))
totrepsumcodondata%>%make_quantcompplot_fac(availability,p_site_occ,time,fname)
#
# fname= here(paste0('plots/p_a_redux/','p_site_occ_vs_poly_av_vs_dt','.pdf'))
# polyrepsumcodondata%>%
	# left_join(poly_trnadat,by=c('time','codon'))%>%
	# make_quantcompplot_fac(poly_availability,p_site_occ,time,fname)
#
fname= here(paste0('plots/p_a_redux/','abundance_freq','.pdf'))
trna_dat%>%left_join(enframe(overallcodonfreqs,'codon','freq'),by='codon')%>%
	filter(fraction=='Total')%>%
	group_by(time,fraction)%>%
	# mutate(abundance = replace_na(abundance,min(abundance%>%keep(is.finite))))%>%
	make_quantcompplot_fac(abundance,freq,time,fname)
#
fname= here(paste0('plots/p_a_redux/','ab_v_wusage','.pdf'))
trna_dat%>%
	left_join(weighted_codon_usage)%>%
	filter(fraction=='Total')%>%
	group_by(time,fraction)%>%
	make_quantcompplot_fac(abundance,weightedusage,time,fname)

stop()
################################################################################
########Process these
################################################################################

# if(all(fprofilemats_unproc$sample %in% 1:100)) fprofilemats_unproc$sample %<>% {names(fpsitelist)[as.numeric(.)]}

# codonprofiles<-fprofilemats_unproc
# codonprofiles%<>%filter(!codon %in% c('TAG','TAA','TGA'))
# codonprofiles%<>%mutate(position = position - 1 - (FLANKCODS*3))
# codonprofiles%<>%group_by(readlen)%>%filter(any(signal!=0))
# codonprofiles%<>%group_by(readlen,codon,sample)%>%
# 	mutate(occ_nonorm=signal)%>%
# 	mutate(occupancy = signal / median(signal))
# codonprofiles%<>%select(-signal)
# stopifnot(codonprofiles$occupancy%>%is.finite%>%all)
# codonprofiles %>%saveRDS('data/codonprofiles.rds')
# codonprofiles <- readRDS('data/codonprofiles.rds')

# if(!exists('codonprofiledat')) base::source('src/Figures/Figure1/3_1a_metacodons_tr.R')
# stopifnot((codonprofiledat$sample%>%n_distinct)==22)

# if(!file.exists(here('data/codonprofiledat.rds'))){

# }else{
	# codonprofiledat<-readRDS(here('data/codonprofiledat.rds'))
# }

################################################################################
########Adjacency predicts p-site occupancy?
################################################################################
sixmermat <- oligonucleotideFrequency(cdsseq,step=6,width=6)

adjacency_freq_df <- sixmermat %>% colSums%>%
	enframe('sixmer','freq')%>%
	mutate(cod1=str_extract(sixmer,'^\\w{3}'),cod2=str_extract(sixmer,'\\w{3}$'))%>%
	select(-sixmer)%>%
	identity

adjacency_freq_df%<>%group_by(cod1)%>%mutate(freq = freq/sum(freq))

psiteoccs = subfpprofilelist[['allhigh']]%>%filter(sample%>%str_detect('E13_ribo'),readlen=='rl29')%>%filter(position== -11)%>%
	group_by(sample,codon)%>%summarise_at(vars(count),mean)%>%
	group_by(codon)%>%summarise_at(vars(count),mean)
asiteoccs = subfpprofilelist[['allhigh']]%>%filter(sample%>%str_detect('E13_ribo'),readlen=='rl29')%>%filter(position== -11-3)%>%
	group_by(sample,codon)%>%summarise_at(vars(count),mean)%>%
	group_by(codon)%>%summarise_at(vars(count),mean)
asiteoccs%>%rename('a_site_occ'=count)%>%
	left_join(psiteoccs%>%rename('p_site_occ'=count),by='codon')%>%
	mutate(foo=1)%>%
	make_quantcompplot_fac(a_site_occ,p_site_occ,foo,fname='plots/atimes_adjecency_vs_psite_occ.pdf')

# adjacencycompplot <- 
# adjacency_freq_df%>%
# 		left_join(asiteoccs,by=c('cod2'='codon'))%>%
# 		group_by(cod1)%>%
# 		summarise(wsum_asiteocc=weighted.mean(count,freq,na.rm=T))%>%
# 		left_join(psiteoccs,by=c('cod1'='codon'))%>%
#         # {quicktest(.$wsum_asiteocc,.$count)}
#         identity

# adjacencycompplot%>%mutate(foo=1)%>%make_quantcompplot_fac(wsum_asiteocc,count,foo,fname='plots/atimes_adjecency_vs_psite_occ.pdf')


psiteoccs = subfpprofilelist[['allhigh']]%>%filter(sample%>%str_detect('P0'),readlen=='rl29')%>%filter(position== -11)%>%
	group_by(sample,codon)%>%summarise_at(vars(count),mean)%>%
	group_by(codon)%>%summarise_at(vars(count),mean)
asiteoccs = subfpprofilelist[['allhigh']]%>%filter(sample%>%str_detect('P0'),readlen=='rl29')%>%filter(position== -11)%>%
	group_by(sample,codon)%>%summarise_at(vars(count),mean)%>%
	group_by(codon)%>%summarise_at(vars(count),mean)

adjacencycompplot <- 
adjacency_freq_df%>%
		left_join(asiteoccs,by=c('cod2'='codon'))%>%
		group_by(cod1)%>%
		summarise(wsum_asiteocc=weighted.mean(count,freq,na.rm=T))%>%
		left_join(psiteoccs,by=c('cod1'='codon'))%>%
        # {quicktest(.$wsum_asiteocc,.$count)}
        identity


asiteoccs%>%left_join(psiteoccs,by='codon')%>%mutate(foo=1)%>%make_quantcompplot_fac(wsum_asiteocc,count,foo,fname='plots/atimes_adjecency_vs_psite_occ.pdf')


################################################################################
########Gingold codon types
################################################################################
gingold_types = 'ext_data/gingold_etal_2014_trna_types.tsv'%>%read_tsv
gingold_types=gingold_types%>%mutate(codon = DNAStringSet(anticodon)%>%reverseComplement%>%as.character)
if(!'type'%in%colnames(repsumcodondata)){
	repsumcodondata%<>%left_join(gingold_types%>%select(codon,type),by='codon')
}
#now plot
plotfile<- here(paste0('plots/','gingold_type_dist','.pdf'))
pdf(plotfile)
	repsumcodondata%>%
	filter(assay=='ribo',!is.na(type))%>%
	ggplot(.,aes(fill=type,x=p_site_occ))+
	geom_density(alpha=I(0.5))+
	facet_wrap(time~.)+
    # scale_fill_manual(values = stagecols) +
	scale_x_continuous(paste0('P_site_occ'))+
	ggtitle(paste0('DT vs gingold types'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))

lapply(unique(repsumcodondata$time),function(itp){
	repsumcodondata%>%
		filter(assay=='ribo',!is.na(type))%>%
		filter(time==itp,type%in%c('prol','diff'))%>%
		{split(.$p_site_occ,.$type)}%>%{t.test(.[[1]],.[[2]])}
})%>%map_df(.id='time',tidy)%T>%write_tsv('tables/gingold_ttest_vals.tsv')
message(normalizePath('tables/gingold_ttest_vals.tsv'))



repsumcodondata

################################################################################
########Okay so do poly DT and total DT even agree??
################################################################################
# left_join(totrepsumcodondata,polyrepsumcodondata,by=c('time','codon'),suffix=c('_tot','_poly'))%>%
	# make_quantcompplot_fac(a_site_occ_tot,a_site_occ_poly,time,
		# fname=paste0('plots/bDT_dtpoly_vsdt.pdf'))


#do the tRNAs which drop out have unusual dwell times?
repsumcodondata%>%mutate(abundancedropout=is.na(abundance))%>%
	group_by(assay,time)%>%nest%>%
	mutate(test=map(data,~tidy(split(.$dwell_time,.$abundancedropout)%>%{t.test(.[[1]],.[[2]])})))%>%
	unnest(test)

#now plot
plotfile<- here(paste0('plots/','ab_dropout_vs_dt','.pdf'))
pdf(plotfile)
repsumcodondata%>%
	mutate(abundancedropout=is.na(abundance))%>%
	ggplot(.,aes(x=dwell_time,fill=abundancedropout))+
	# geom_density(alpha=I(0.5))+
	geom_histogram(alpha=I(0.5))+
	# scale_color_discrete(name='colorname',colorvals)+
	facet_grid(assay~time)+
	scale_x_continuous(paste0('Dwell_time'))+
	# scale_y_continuous(paste0('yname'))+
	ggtitle(paste0('Dwell_time vs Abundance Dropout'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))




################################################################################
########Linear modeling
################################################################################
	
totrepsumcodondata%>%group_by(time)%>%nest%>%
	summarise(map(data,~tidy(lm(data=.,dwell_time~AA+abundance))))


totrepsumcodondata%>%split(.,.$time)%>%.[[1]]%>%lm(data=.,dwell_time~abundance+AA)%>%anova
totrepsumcodondata%>%split(.,.$time)%>%.[[1]]%>%lm(data=.,dwell_time~AA+abundance)%>%anova

totrepsumcodondata%>%split(.,.$time)%>%.[[1]]%>%lm(data=.,dwell_time~availability+AA)%>%anova
totrepsumcodondata%>%split(.,.$time)%>%.[[1]]%>%lm(data=.,dwell_time~AA+availability)%>%anova


# LRT(
	# totrepsumcodondata%>%split(.,.$time)%>%.[[1]]%>%lm(data=.,abundance~AA)%>%anova
	# totrepsumcodondata%>%split(.,.$time)%>%.[[1]]%>%lm(data=.,aacor_dwell_time~abundance+AA)%>%anova
# )

gingold_types = 'ext_data/gingold_etal_2014_trna_types.tsv'%>%read_tsv
gingold_types=gingold_types%>%mutate(codon = DNAStringSet(anticodon)%>%reverseComplement%>%as.character)


#now plot
plotfile<- here(paste0('plots/','stage_dtdist','.pdf'))
pdf(plotfile)
p1=totrepsumcodondata%>%
	ggplot(.,aes(x=dwell_time,fill=time))+
	geom_density(alpha=I(0.5))+
	scale_fill_manual(name='stage',values=stagecols)+
	scale_x_continuous(paste0('dwell_time_broad'))+
	ggtitle(paste0('broad dt_distribution'))+
	theme_bw()
p2=totrepsumcodondata%>%
	ggplot(.,aes(x=dwell_time,y=time,color=time))+
	geom_point(alpha=I(0.5),position='jitter')+
	scale_color_manual(name='stage',values=stagecols)+
	scale_x_continuous(paste0('dwell_time_broad'))+
	ggtitle(paste0('broad dt_distribution'))+
	theme_bw()
ggarrange(plotlist=list(p1,p2),nrow=2)
dev.off()
message(normalizePath(plotfile))



#now plot
plotfile<- here(paste0('plots/','broad_dt_aa_stripplot','.pdf'))
pdf(plotfile,w=15,h=5)
totrepsumcodondata%>%
	group_by(AA)%>%
		mutate(AA = GENETIC_CODE[codon])%>%
		# filter(fraction=='Total')%>%
		mutate(codonname = paste0(AA,'-',codon))%>%
		group_by(codonname)%>%
		mutate(codmean=mean(dwell_time))%>%group_by(AA)%>%mutate(aamean=mean(dwell_time))%>%
		ungroup()%>%arrange(aamean,codmean)%>%
		mutate(codonname=as_factor(codonname),AA=as_factor(AA),codonname=as_factor(codonname))%>%
	arrange(mean(dwell_time))%>%
	mutate(AA=as_factor(AA))%>%
	ggplot(.,aes(y=dwell_time,x=codon,color=time))+
	geom_point(alpha=I(0.5))+
	facet_grid(.~AA,scale='free')+
	scale_color_manual(name='stage',values=stagecols)+
	scale_y_continuous(paste0('dwell_time_broad'))+
	ggtitle(paste0('broad dt_aa_stripplot'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))





