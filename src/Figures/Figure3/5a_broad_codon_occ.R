
################################################################################
########Process these
################################################################################

if(all(fprofilemats_unproc$sample %in% 1:100)) fprofilemats_unproc$sample %<>% {names(fpsitelist)[as.numeric(.)]}

codonprofiles<-fprofilemats_unproc
codonprofiles%<>%filter(!codon %in% c('TAG','TAA','TGA'))
codonprofiles%<>%mutate(position = position - 1 - (FLANKCODS*3))
codonprofiles%<>%group_by(readlen)%>%filter(any(signal!=0))
codonprofiles%<>%group_by(readlen,codon,sample)%>%
	mutate(occ_nonorm=signal)%>%
	mutate(occupancy = signal / median(signal))
codonprofiles%<>%select(-signal)
stopifnot(codonprofiles$occupancy%>%is.finite%>%all)
codonprofiles %>%saveRDS('data/codonprofiles.rds')
codonprofiles <- readRDS('data/codonprofiles.rds')



trna_ab_df_samp = allcodsigmean_isomerge[c("fraction", "time", "sample", "anticodon", "abundance", "codon",
"weightedusage", "availability","rep")]

codon_data <- trna_ab_df_samp%>%
	filter(fraction=='Total')%>%select(-fraction,-sample)%>%
    select(time,rep,codon,abundance,availability)%>%
    
    group_by(time,AA)%>%
    mutate(aacor_dwell_time =)
    group_by(time,codon)%>%
    summarise_at(vars(one_of(c('abundance','availability'))),list(mean))

repsumcodondata<-codondata%>%
	group_by(time,codon)%>%summarise_at(vars(abundance,availability,dwell_time,aacor_dwell_time),mean)

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
	nlabel= compdf%>%
		group_by(facet)%>%
		tally%>%
		summarise(labl=paste0('N=',n))
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


fname= here(paste0('plots/','bDT_ab_vs_dt','.pdf'))
repsumcodondata%>%make_quantcompplot(abundance,dwell_time,time,fname)
fname= here(paste0('plots/','bDT_av_vs_dt','.pdf'))
repsumcodondata%>%make_quantcompplot(availability,dwell_time,time,fname)
fname= here(paste0('plots/','bDT_aacordt_vs_dt','.pdf'))
repsumcodondata%>%make_quantcompplot(availability,dwell_time,time,fname)




pdf(h=15,w=5,plotfile)
col1=sym('aacor_dwell_time')
col2=sym('availability')
corlabel = repsumcodondata%>%filter(is.finite(!!col1),is.finite(!!col2))%>%
		group_by(time)%>%
		summarise(tidy(cor.test(!!col1, !!col2)))
corlabel = corlabel%>%
	mutate(
		pformat=format(p.value,format='e',digits=4),
		pvalstring = ifelse(p.value > 0.001,round(p.value,4),pformat),
		labl=paste0('rho = ',round(estimate,3),'\n','pval = ',pvalstring))
nlabel=repsumcodondata%>%group_by(time)%>%summarise(labl=paste0('N=',n()))
repsumcodondata%>%
	ggplot(.,aes(y=availability,x=aacor_dwell_time))+
	geom_point(alpha=I(0.5))+
	geom_smooth(method='lm')+
	facet_grid(time~.)+
	scale_fill_manual(name='stage',values=stagecols)+
	scale_x_continuous(paste0('Broad Dwell Time'))+
		geom_text(show.legend=F,data=corlabel,
			hjust=1,vjust=1,x= Inf,y=Inf,aes(color=NULL,label=labl))+
	geom_text(show.legend=F,data=nlabel,
			hjust=0,vjust=1,x= -Inf,y=Inf,aes(color=NULL,label=labl))
	# scale_y_continuous(paste0(''))+
	ggtitle(paste0('Dwell Time Dist'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))
