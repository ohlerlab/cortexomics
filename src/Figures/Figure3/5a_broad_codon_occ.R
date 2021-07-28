if(!exists('allcodsig_isomerge')) base::source(here('src/Figures/Figure3/3_tRNA_array_analysis.R'))
if(!exists('codonprofiledat')) base::source(here('src/subset_dwell_times.R'))

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

psiteoccs = codonprofiledat%>%filter(sample%>%str_detect('E13_ribo'),readlen=='rl29')%>%filter(position== -11)%>%
	group_by(sample,codon)%>%summarise_at(vars(count),mean)%>%
	group_by(codon)%>%summarise_at(vars(count),mean)
asiteoccs = codonprofiledat%>%filter(sample%>%str_detect('E13_ribo'),readlen=='rl29')%>%filter(position== -11)%>%
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

		
adjacencycompplot%>%make_quantcompplot(wsum_asiteocc,count,'plots/atimes_adjecency_vs_psite_occ.pdf')
################################################################################
########Adjacency predicts p-site occupancy?
################################################################################
	
adjacency_freq_df <- sixmermat %>% colSums%>%
	enframe('sixmer','freq')%>%
	mutate(cod1=str_extract(sixmer,'^\\w{3}'),cod2=str_extract(sixmer,'\\w{3}$'))%>%
	select(-sixmer)%>%
	identity

adjacency_freq_df%<>%group_by(cod1)%>%mutate(freq = freq/sum(freq))

psiteoccs = codonprofiledat%>%filter(sample%>%str_detect('P0'),readlen=='rl29')%>%filter(position== -11)%>%
	group_by(sample,codon)%>%summarise_at(vars(count),mean)%>%
	group_by(codon)%>%summarise_at(vars(count),mean)
asiteoccs = codonprofiledat%>%filter(sample%>%str_detect('P0'),readlen=='rl29')%>%filter(position== -11)%>%
	group_by(sample,codon)%>%summarise_at(vars(count),mean)%>%
	group_by(codon)%>%summarise_at(vars(count),mean)


adjacency_freq_df%>%
		left_join(asiteoccs,by=c('cod2'='codon'))%>%
		group_by(cod1)%>%
		summarise(wsum_asiteocc=weighted.mean(count,freq,na.rm=T))%>%
		left_join(psiteoccs,by=c('cod1'='codon'))%>%slice(which.max(count))
		{quicktest(.$wsum_asiteocc,.$count)}

clean_fr_sampnames<-function(x) x%>%str_replace('.*_(Poly|80S)(.*)_()','\\2_\\1ribo_\\3')


codonvarprofiles <-  codonprofiledat%>%
# rustprofiledat%>%
	ungroup%>%
	group_by(sample,readlen,position)%>%
	mutate(count = count / mean(count))%>%
	# filter(count==0)%>%.$codon%>%unique
	# filter(signal!=0)%>%
	# filter(position==1)%>%
	# filter(sample%>%str_detect('ribo_1'))%>%
	filter(sample%>%is_in(mainsamps))%>%
	filter(!is.nan(count))%>%
	# summarise(sdsig=sd(count,na.rm=T)/median(count,na.rm=T))%>%
	separate(sample,c('time','assay','rep'))%>%
	group_by(rep,time,readlen,position)%>%
	mutate(count = count / mean(count))%>%
	summarise(sdsig=sd(count,na.rm=T))%>%
	group_by(readlen,time,position)%>%
	summarise(sdsig=mean(sdsig))%>%
	mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
	# filter(position> -6,position < 6)%>%
	filter(position> -numreadlen+6,position< -6)%>%
	filter(numreadlen>=25,numreadlen<=31)%>%
	arrange(position)

codonvaroffsets = codonvarprofiles%>%group_by(readlen,time)%>%slice(which.max(sdsig))%>%
	mutate(offset=-(position+3))%>%
	as.data.frame

codonvaroffsets = codonprofiledat%>%
	ungroup%>%
	distinct(sample)%>%
	separate(sample,c('time','assay','rep'),remove=F)%>%
	left_join(codonvaroffsets)%>%
	select(sample,readlen,numreadlen,offset)

{
offsets%<>%mutate(readlen=paste0(length))
plotfile='plots/sh_fppos_vs_codon_variance.pdf'
pdf(plotfile,w=12,h=3*5)
#plotting variance amongst codons at each point.
# sh_codprof%>%
codonvarprofiles%>%
	{
		qplot(data=.,x=position,y=sdsig)+
		theme_bw()+
		facet_grid(readlen~time)+
		scale_y_continuous('between codon variation (meannorm)')+
		scale_x_continuous('5 read position relative to codon ')+
		geom_vline(data=filter(offsets,readlen%in%.$readlen),aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=filter(offsets,readlen%in%.$readlen),aes(xintercept= -offset),color=I('green'),linetype=2)+
		# geom_vline(xintercept= 0,color=I('blue'),linetype=2)+
		# geom_vline(xintercept= -5,color=I('green'),linetype=2)+
		ggtitle("variance of 5' read occurance vs position")
	}%>%print
dev.off()
normalizePath(plotfile)
}



#does amino acid variability determine the asite effects?
codonprofiledat%>%left_join(offsets)%>%filter(position== -offset)%>%
	mutate(AA=GENETIC_CODE[codon])%>%
	lm(data=.,count~AA)%>%summary

codonprofiledat%>%left_join(offsets)%>%
	filter(between(position, -offset-5,-offset-3))%>%
	group_by(sample,codon)%>%summarise(count=sum(count))%>%
	mutate(AA=GENETIC_CODE[codon])%>%
	lm(data=.,count~AA)%>%summary

codonprofiledat


allcodsig_isomerge%<>%addcodon
trna_ab_df_samp = allcodsig_isomerge[c("fraction", "time", "sample", "anticodon", "abundance", "codon",
"weightedusage", "availability","rep")]
#get data on the trnas
trna_dat <- trna_ab_df_samp%>%
	select(-sample)%>%
    select(fraction,time,rep,codon,abundance,availability)%>%
    group_by(fraction,time,codon)%>%
    summarise_at(vars(one_of(c('abundance','availability'))),list(mean),na.rm=T)
#

rloffsets<-offsets%>%select(readlen,length,offset)%>%mutate(readlen=str_replace(length,'^(\\d)','rl\\1'))
lexp = 3+3 #include positions for positions corresponding to bigger offsets
rexp = 3#include positions for positions corresponding to smaller offsets
# codonprofiles%>%
codondata <- 
	# rustprofiledat%>%
	codonprofiledat%>%
	# filter(subgroup=='nochangehighe')%>%
	# filter(subgroup=='all')%>%
# codondata <- rustprofiledat%>%
	# filter(sample%>%str_detect('ribo'))%>%
	mutate_at(vars(sample),clean_fr_sampnames)%>%
    group_by(sample)%>%
    # mutate(readlen=paste0('rl',length))%>%
    # safe_left_join(rloffsets)%>%
    safe_left_join(codonvaroffsets)%>%
    separate(sample,c('time','assay','rep'))%>%
    # filter(readlen=='rl29')%>%
    # filter(readlen%in%c('rl29','rl30','rl28','rl27'))%>%
    filter(readlen%in%c('rl29','rl31','rl30','rl27','rl26'))%>%
    # mutate(offset=0)%>%#for the shifted data
    # filter(position <= -(offset-rexp))%>%
    # filter(position >= -(offset+lexp))%>%
    # filter(position == -offset)%>%
    filter(between(position, -offset-5,-offset))%>%#psite
    # filter(between(position, -offset-5,-offset-3))%>%
    group_by(assay,time,codon,rep)%>%
    summarise(
    	p_site_occ = sum(count[position%>%between(-offset-2,-offset)]),
    	a_site_occ = sum(count[position%>%between(-offset-5,-offset-3)])
    )%>%
    arrange(assay!='ribo')
#seperate poly and total tRNA info
total_trnadat <- trna_dat%>%filter(fraction=='Total')%>%select(time,codon,abundance,availability)
poly_trnadat <- trna_dat%>%filter(fraction=='Poly')%>%
	select(time,codon,poly_abundance=abundance,poly_availability=availability)
#add tRNA info, 
#AA and aa corrected dwell time
codondata%<>%
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
repsumcodondata%>%filter(assay=='ribo',time=='E13')%>%group_by(assay,time)%>%group_slice(1)%>%{quicktest(.$p_site_occ,.$availability)}
# codondata%>%filter(assay=='ribo')%>%group_by(assay,time)%>%group_slice(2)%>%{quicktest(.$p_site_occ,.$availability)}
# repsumcodondata%>%filter(assay=='ribo')%>%group_by(assay,time)%>%group_slice(2)%>%{quicktest(.$p_site_occ,.$availability)}
# repsumcodondata%>%filter(assay=='Polyribo')%>%group_by(assay,time)%>%group_slice(1)%>%{quicktest(.$p_site_occ,.$poly_abundance)}
codondata %>% saveRDS(here('data/codondata.rds'))
repsumcodondata %>% saveRDS(here('data/repsumcodondata.rds'))

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
########
################################################################################
	
repsumcodondata%>%left_join(codon_usage%>%select(codon,freq,cAI))%>%group_by(time)%>%group_slice(1)%>%{quicktest(.$abundance,.$freq)}

make_quantcompplot <- function(compdf, col1, col2, facetvar=NULL, fname){
	require(LSD)
	base::source(here('Applications/LSD/R/LSD.heatscatter.R'))
	require(broom)
	col1<-enquo(col1)
	col2<-enquo(col2)
	facetvar<-enquo(facetvar)
	compdf%<>%rename(facet=!!facetvar)
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
fname= here(paste0('plots/','bDT_ab_vs_dt','.pdf'))
totrepsumcodondata%>%make_quantcompplot(abundance,dwell_time,time,fname)

fname= here(paste0('plots/','bDT_av_vs_dt','.pdf'))
totrepsumcodondata%>%make_quantcompplot(availability,dwell_time,time,fname)

fname= here(paste0('plots/','bDT_aacordt_vs_dt','.pdf'))
totrepsumcodondata%>%make_quantcompplot(availability,aacor_dwell_time,time,fname)

fname= here(paste0('plots/','bDT_cAI_vs_dt','.pdf'))
totrepsumcodondata%>%
	left_join(codon_usage%>%select(codon,cAI))%>%
	make_quantcompplot(availability,cAI,time,fname)
fname= here(paste0('plots/','bDT_wcAI_vs_dt','.pdf'))
totrepsumcodondata%>%
	left_join(weighted_codon_usage%>%select(time,codon,w_cAI))%>%
	make_quantcompplot(availability,w_cAI,time,fname)

fname= here(paste0('plots/','bDT_poly_av_vs_dt','.pdf'))
totrepsumcodondata%>%
	left_join(poly_trnadat,by=c('time','codon'))%>%
	make_quantcompplot(poly_availability,dwell_time,time,fname)

fname= here(paste0('plots/','bDT_inpoly_av_vs_dt','.pdf'))
totrepsumcodondata%>%
	left_join(poly_trnadat,by=c('time','codon'))%>%
	filter(!is.na(poly_availability))%>%
	make_quantcompplot(availability,dwell_time,time,fname)


fname= here(paste0('plots/','bDT_poly_ab_vs_dt','.pdf'))
totrepsumcodondata%>%
	left_join(poly_trnadat,by=c('time','codon'))%>%
	make_quantcompplot(poly_abundance,dwell_time,time,fname)


#now do all that with the polysomal
fname= here(paste0('plots/','bDT_ab_vs_dtpoly','.pdf'))
totrepsumcodondata%>%make_quantcompplot(abundance,dwell_time,time,fname)
fname= here(paste0('plots/','bDT_av_vs_dtpoly','.pdf'))
totrepsumcodondata%>%make_quantcompplot(availability,dwell_time,time,fname)
fname= here(paste0('plots/','bDT_aacordt_vs_dtpoly','.pdf'))
totrepsumcodondata%>%make_quantcompplot(availability,aacor_dwell_time,time,fname)
fname= here(paste0('plots/','bDT_cAI_vs_dtpoly','.pdf'))
totrepsumcodondata%>%
	left_join(codon_usage%>%select(codon,cAI))%>%
	make_quantcompplot(availability,cAI,time,fname)
fname= here(paste0('plots/','bDT_wcAI_vs_dtpoly','.pdf'))
totrepsumcodondata%>%
	left_join(weighted_codon_usage%>%select(time,codon,w_cAI))%>%
	make_quantcompplot(availability,w_cAI,time,fname)
fname= here(paste0('plots/','bDT_poly_av_vs_dtpoly','.pdf'))
totrepsumcodondata%>%
	make_quantcompplot(poly_availability,dwell_time,time,fname)
fname= here(paste0('plots/','bDT_inpoly_av_vs_dtpoly','.pdf'))
totrepsumcodondata%>%
	filter(!is.na(poly_availability))%>%
	make_quantcompplot(availability,dwell_time,time,fname)
fname= here(paste0('plots/','bDT_poly_ab_vs_dtpoly','.pdf'))
totrepsumcodondata%>%
	make_quantcompplot(poly_abundance,dwell_time,time,fname)


################################################################################
########Okay so do poly DT and total DT even agree??
################################################################################
left_join(totrepsumcodondata,polyrepsumcodondata,by=c('time','codon'),suffix=c('_tot','_poly'))%>%
	make_quantcompplot(dwell_time_tot,dwell_time_poly,time,
		fname=paste0('plots/bDT_dtpoly_vsdt.pdf'))


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


LRT(
	totrepsumcodondata%>%split(.,.$time)%>%.[[1]]%>%lm(data=.,abundance~AA)%>%anova
	totrepsumcodondata%>%split(.,.$time)%>%.[[1]]%>%lm(data=.,aacor_dwell_time~abundance+AA)%>%anova
)

gingold_types = 'ext_data/gingold_etal_2014_trna_types.tsv'%>%read_tsv
gingold_types=gingold_types%>%mutate(codon = DNAStringSet(anticodon)%>%reverseComplement%>%as.character)


totrepsumcodondata%>%left_join()

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





