#first run
#src/R/Figures/Figure2/codon_coverage.R
if(!exists('here')) base::source(here::here('src/R/Rprofile.R'))
 conflict_prefer("intersect", "BiocGenerics")
if(!exists('codonprofiles')) load('data/codon_coverage.Rdata')
if(!exists('allcodsigmean_isomerge')) source(here('src/tRNA_array_analysis.R'))


stagecolsdisplay <- c(E12.5 = "#214098", E14 = "#2AA9DF", E15.5 = "#F17E22", E17 = "#D14E28", P0 = "#ED3124")
stagecols <- stagecolsdisplay %>% setNames(c("E13", "E145", "E16", "E175", "P0"))
stageconv <- names(stagecols) %>% setNames(c("E13", "E145", "E16", "E175", "P0"))


codonprofiles <- readRDS('data/codonprofiles.rds')

codonprofiles$fraction = case_when(
	codonprofiles$sample%>%str_detect('ribo') ~ 'total',
	codonprofiles$sample%>%str_detect('Poly') ~ 'poly',
	TRUE ~ 'NA'
)
stopifnot(!any(codonprofiles$fraction=='NA'))


offsets <- read_tsv('ext_data/offsets_manual.tsv')
allTEchangedf<-read_tsv('tables/manuscript/go_all_highcount_updown.tsv')

GENETIC_CODE<-Biostrings::GENETIC_CODE

fig2 <- function(){

	offsets %<>% select(offset, compartment, length, readlen)
	codonprofiles %<>% select(sample,codon,readlen,position,occ_nonorm,occupancy,fraction)
	allcodsigmean_isomerge%<>%.[c("fraction", "time", "sample", "anticodon", "abundance", "codon",
"AA",  "weightedusage", "balance")]

	offsets%<>%mutate(readlen=paste0('rl',length))

	codonoccs<-codonprofiles%>%
		inner_join(offsets,by='readlen')%>%
		filter(sample%>%str_detect('ribo'))%>%
		filter(position <= -offset-3, position >= -offset-5)%>%
		# filter(position <= -offset, position >= -offset-5)%>%
		# filter(position <= -offset-3, position >= -offset-5)%>%
		# filter(readlen%in%c('rl29'))%>%
		ungroup%>%
		mutate(time=str_extract(sample,'[^_]+'))%>%
		group_by(time,codon)%>%
		# group_slice(2)
		summarise(occupancy=mean(occupancy ,na.rm=T))

	##Occupancy over the span of the read

	pdf('plots/figures/figure2/trna_codons/fppos_vs_codon_variance.pdf',w=12,h=12)
	#plotting variance amongst codons at each point.
	codonprofiles%>%
		filter(fraction=='total')%>%
		ungroup%>%
		group_by(sample,readlen,position)%>%
		# filter(abundance ==0)%>%.$codon%>%unique
		# filter(abundance !=0)%>%
		# filter(position==1)%>%
		filter(!is.nan(occupancy ))%>%
		mutate(occupancy =occ_nonorm)%>%
		summarise(sdsig=sd(occupancy ,na.rm=T)/mean(occupancy ,na.rm=T))%>%
		separate(sample,c('time','assay','rep'))%>%
		group_by(readlen,time,assay,position)%>%
		summarise(sdsig=mean(sdsig))%>%
		mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
		filter(position> -numreadlen+6,position < -6)%>%
		filter(numreadlen>=25,numreadlen<=31)%>%
		arrange(position)%>%
		# filter(numreadlen==29)%>%
		{
			ggplot(data=.,aes(,x=position,y=sdsig))+
			geom_point()+
			theme_bw()+
			facet_grid(readlen~time)+
			scale_y_continuous('between codon variation (meannorm)')+
			scale_x_continuous('5 read position relative to codon',breaks=seq(-30,0,3))+
			geom_rect(data=offsets,aes(y=NULL,x=NULL,xmin=0.5-3-offset,xmax= 0.5-offset),
					x=NULL,ymin=-Inf,ymax=Inf,alpha=I(0.5),fill=I('darkgreen'))+
			geom_rect(data=offsets,aes(y=NULL,x=NULL,xmin=0.5-6-offset,xmax= 0.5-3-offset),
					x=NULL,ymin=-Inf,ymax=Inf,alpha=I(0.5),fill=I('purple'))+
			geom_rect(data=offsets,aes(y=NULL,x=NULL,xmin=0.5-3-offset+3,xmax= 0.5-offset+3),
					x=NULL,ymin=-Inf,ymax=Inf,alpha=I(0.5),fill=I('orange'))+
			ggtitle("variance of 5' read occurance vs position")
		}%>%print
	dev.off()
	normalizePath('plots/figures/figure2/trna_codons/fppos_vs_codon_variance.pdf')

	#occpancy vs AA

	#plot with the AAs as colors
	pdf('plots/figures/figure2/trna_codons/stripplot_aa_codon.pdf',w=12,h=5)
	codonoccs%>%
		mutate(AA = GENETIC_CODE[codon])%>%
		# filter(fraction=='Total')%>%
		mutate(codonname = paste0(AA,'-',codon))%>%
		group_by(codonname)%>%
		mutate(codmean=mean(occupancy))%>%group_by(AA)%>%mutate(aamean=mean(occupancy))%>%
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

	################################################################################
	########Now tRNA only figures
	################################################################################


	plotfile <- "plots/figures/figure2/trna_codons/trna_sig_allfrac_mergedecod.pdf"
	pdf(plotfile, w = 9*2, h = 16)
	allcodsigmean_isomerge %>%
	    # filter(sample %>% str_detect("Total")) %>%
	    group_by(anticodon) %>%
	    mutate(cmean = mean(abundance)) %>%
	    ungroup() %>%
	    arrange(cmean) %>%
	    mutate(anticodon = as_factor(anticodon)) %>%
	    # group_by(decoder,)
	    mutate(fraction=fraction%>%factor(.,c('Total','Poly','80S')))%>%
	    ggplot(data = ., aes(x = anticodon, color = stageconv[time], y = abundance)) +
	    facet_grid( . ~ fraction,scale='free_x')+
	    geom_point() +
	    scale_color_manual(values = stagecols) +
	    scale_y_continuous(name = "- deltaCt", breaks = -seq(0, 20, by = 2.5)) +
	    coord_flip() +
	    theme_bw()
	# theme(axis.text.x=element_text(angle=45,vjust=.5,size=6))
	# facet_grid(time ~ . )
	dev.off()
	normalizePath(plotfile)

	fracpairs = list(c('Poly','Total'),c('80S','Total'),c('80S','Poly'))	
	pdf("plots/figures/figure2/trna_codons/tRNA_poly_mono_cor.pdf",w=21)
	lapply(fracpairs,function(fracpair){allcodsigmean_isomerge%>%ungroup%>%
		filter(fraction%in%c(fracpair[1],fracpair[2]))%>%
		mutate(fraction = fraction%>%{ifelse(.!=fracpair[2],'alt','base')})%>%
		select(time,fraction,codon,abundance)%>%
		spread(fraction,abundance)%>%
		filter(!is.na(alt))%>%
		filter(is.finite(base),is.finite(alt))%>%
		mutate(acor = list(tidy(cor.test(base,alt))))%>%
		# mutate(a = list(tidy(cor.test(Total,`80S`))))%>%
		unnest(col=c(acor),names_sep='_')%>%
		mutate(fraction=paste0(fracpair[1],' vs ',fracpair[2]))%>%
		mutate(acor_label=paste0('rho = ',round(acor_estimate,3),'\n','pval = ',round(acor_p.value,3)))%>%
		# gather(`Poly`,`80S`,a_frac,a_val)
		{
			ggplot(data=.,aes(x=base,y=alt,label=acor_label))+
			scale_x_continuous(name=fracpair[2])+
			scale_y_continuous(name=fracpair[1])+
			# geom_point()+
			ggtitle(unique(.$fraction))+
			geom_text(aes(label=codon))+
			geom_text(data=distinct(.,fraction,acor_label),
				show.legend=F,hjust=0,vjust=1,aes(label=acor_label,x= -Inf,y=Inf),size=4)+
			theme_bw()
		}
	})%>%ggarrange(plotlist=.,ncol=3)
		dev.off()
		normalizePath('plots/figures/figure2/trna_codons/tRNA_poly_mono_cor.pdf')

	plotfile <- "plots/figures/figure2/trna_codons/trna_sig_allfrac_mergedecod.pdf"
	pdf(plotfile, w = 9*2, h = 16)
	allcodsigmean_isomerge %>%
	    # filter(sample %>% str_detect("Total")) %>%
	    group_by(anticodon) %>%
	    mutate(cmean = mean(abundance)) %>%
	    ungroup() %>%
	    arrange(cmean) %>%
	    mutate(anticodon = as_factor(anticodon)) %>%
	    # group_by(decoder,)
	    mutate(fraction=fraction%>%factor(.,c('Total','Poly','80S')))%>%
	    ggplot(data = ., aes(x = anticodon, color = stageconv[time], y = abundance)) +
	    facet_grid( . ~ fraction,scale='free_x')+
	    geom_point() +
	    scale_color_manual(values = stagecols) +
	    scale_y_continuous(name = "- deltaCt", breaks = -seq(0, 20, by = 2.5)) +
	    coord_flip() +
	    theme_bw()
	# theme(axis.text.x=element_text(angle=45,vjust=.5,size=6))
	# facet_grid(time ~ . )
	dev.off()
	normalizePath(plotfile)

	################################################################################
	########now both
	################################################################################
			
	codonoccs%<>%rename('dwell_time':=occupancy)
	#Now merge with the tRNA data
	tRNA_occ_df<-allcodsigmean_isomerge%>%
		ungroup%>%
		left_join(codonoccs)
	tRNA_occ_df%<>%filter(is.finite(abundance))
	##Get the amino acid for each one
	tRNA_occ_df%<>%mutate(AA = GENETIC_CODE[codon]%>%qs('S+'))
	aatrna_occ_df<-tRNA_occ_df%>%
		group_by(time,AA)%>%
		summarise(dwell_time = mean(dwell_time),abundance_aa = log2(sum(2^(-abundance),na.rm=T)))
	#
	tRNA_occ_df%<>%	left_join(usedcodonfreqs%>%colSums%>%enframe('codon','freq'))%>%
		ungroup%>%
		mutate(common = freq>median(freq))
    tRNA_occ_df_en <- tRNA_occ_df%>%left_join(tRNAenrichdf)
    #
 	tRNA_occ_df_en%<>%group_by(time,fraction,AA)%>%mutate(aacor_dwell_time = dwell_time - mean(dwell_time))
	codons4occ <- tRNA_occ_df$codon%>%unique%>%setdiff(c('TAG','TAA','TGA'))

	aa_corrected_occ_av_df <- tRNA_occ_df%>%
		# safe_filter(fraction=='Poly')%>%
		# filter(time=='E13')%>%
		# group_by(AA)%>%filter(n_distinct(codon)>1)%>%
		group_by(fraction,time,AA,codon)%>%summarise(dwell_time = mean(dwell_time),balance=mean(balance))%>%
		# group_by(AA)%>%mutate(balance = balance - mean(balance))%>%
		# group_by(AA)%>%mutate(balance = balance - mean(balance))%>%
		group_by(fraction,time,AA)%>%mutate(dwell_time_aacor = dwell_time - mean(dwell_time))


	tmp1 <- tRNA_occ_df%>%
		safe_filter(fraction=='Poly')%>%
		filter(time=='E13')%>%
		# group_by(AA)%>%filter(n_distinct(codon)>1)%>%
		group_by(AA,codon)%>%summarise(mr = mean(dwell_time),balance=mean(balance))%>%
		# group_by(AA)%>%mutate(balance = balance - mean(balance))%>%
		# group_by(AA)%>%mutate(balance = balance - mean(balance))%>%
		group_by(AA)%>%mutate(mr = mr - mean(mr))

	pdf('tmp.pdf')
	aa_corrected_occ_av_df%>%
		ggplot(data=.,aes(x=dwell_time_aacor,y=balance))+
		geom_point(size=1)+
		facet_grid(fraction ~ time)+
		geom_text(show.legend=F,aes(label=codon,color=NULL),size=4)
	dev.off()
	normalizePath('tmp.pdf')

	aa_corrected_occ_av_df%>%filter(time=='E13',fraction=='80S')%>%filter(balance< -1500)

	sigcols <- c('balance','abundance','weightedusage','abundance_enrich','balance_enrich')
    sigcols <- syms(sigcols)%>%setNames(sigcols)
    sigcol=sigcols[1]
    fractioni='Poly'


	#total occ, row abundance
	testvals <- map_df(.id='fraction',unique(tRNA_occ_df_en$fraction)%>%setNames(.,.),function(fractioni){
		map_df(.id='sigcol',sigcols,function(sigcol){
			message(fractioni)
			message(sigcol)
			selvals<-tRNA_occ_df_en%>%
			filter(fraction==fractioni)%>%
				filter(time=='E13')%>%
				select(AA,codon,dwell_time,!!sigcol)%>%
				group_by(AA)%>%
				filter(n_distinct(codon)>1)%>%
				group_by(AA)%>%mutate(aacor_dwell_time = dwell_time - mean(dwell_time))%>%
				select(aacor_dwell_time,!!sigcol)%>%
				# group_by(AA)%>%summarise(dwell_time = mean(dwell_time),abundance=mean(!!sigcol))%>%
				# group_by(AA)%>%mutate(abundance = abundance - mean(abundance))%>%
				identity%>%
				ungroup
			if(all(is.na(selvals[[quo_text(sigcol)]]))) return(NULL)
			selvals%>%{tidy(cor.test(use='complete',.$aacor_dwell_time,.[[quo_text(sigcol)]]))}
		})
	})
 	


	#make the grid of plots for the different codons
	for(ifraction in c('Poly','Total')){
		for(timeset in list(c('E13','P0'),tps)) {
		tRNA_occ_df_en$weightedusage =tRNA_occ_df_en$weightedusage
		sigcols <- c('abundance','weightedusage','balance','aacor_dwell_time','dwell_time')
	    sigcols <- syms(sigcols)%>%setNames(sigcols)
	    sigcol=sigcols[1]
	    fractioni='Poly'
	    sigcolx <- sigcols[[3]]
		sigcoly <- sigcols[[3]]
		plots<-lapply(sigcols,function(sigcoly){
			lapply(sigcols,function(sigcolx){
				if(sigcolx==sigcoly){	
				    plot <- tRNA_occ_df_en%>%filter(time%in%timeset,fraction==ifraction)%>%
				    	{qplot(data=.,x=!!sigcolx,fill=time,geom='blank')+
				    	# geom_histogram(aes(y=..density..))+
				    	geom_density(aes(y=..density..),alpha=I(0.5))+
				    	scale_fill_manual(values=stagecols)}
				}else{
				    plot <- tRNA_occ_df_en%>%filter(time%in%timeset,fraction==ifraction)%>%
				    	{qplot(data=.,x=!!sigcolx,color=time,y=!!sigcoly,geom='point',position='')+scale_color_manual(values=stagecols)} 
				}
				plot <- plot+theme_bw()
			})
		})
		plots <- ggarrange(plotlist=plots%>%purrr::flatten(.),ncol=length(sigcols),nrow=length(sigcols))
		plots <- plots+theme_bw()
		plotfile = paste0('plots/figures/figure2/trna_codons/codon_stat_grid_',ifraction,'_',paste0(timeset,collapse='-'),'.pdf')
	    pdf(plotfile,h=20,w=20)
	    print(plots)
	    dev.off()
		normalizePath(plotfile)%>%message
		}
	}

	tRNA_occ_df_en%>%filter(time=='E13',fraction=='Poly')%>%lm(data=.,dwell_time ~ AA+abundance)	%>%anova
	tRNA_occ_df_en%>%filter(fraction=='Poly')%>%lm(data=.,dwell_time ~ AA+time+codon+abundance)	%>%anova



	Sys.glob('pipeline/star/reports/*/*idx*')%>%str_subset('total')%>%str_subset('.bam.',neg=T)%>%setNames(.,.)%>%map_dbl(.%>%fread%>%filter(V1%>%str_detect('chr'))%>%group_by(V1=='chrM')%>%summarise(n=sum(V3))%>%{.$n[2]/sum(.$n)})%>%
		txtplot

	Sys.glob('pipeline/star/reports/*/*idx*')%>%str_subset('ribo')%>%str_subset('.bam.',neg=T)%>%setNames(.,.)%>%map_dbl(.%>%fread%>%filter(V1%>%str_detect('chr'))%>%group_by(V1=='chrM')%>%summarise(n=sum(V3))%>%{.$n[2]/sum(.$n)})%>%
		txtplot


	################################################################################
	########Now look at gene-level changes in TE vs 'elong'
	################################################################################
	


	get_tidy_pval <- function(x) x %>%tidy%>%.$p.value%>%round(4)%>%{ifelse(. == 0,' p < 0.0001 ',paste0('p = ',.))}	
	wtest <- function(df,contcol,factcol){
		contcol=quo_text(contcol)
		factcol=quo_text(factcol)
	    stopifnot(n_distinct(df[[factcol]])==2)
	    df%>%{split(.[[contcol]],.[[factcol]])}%>%{wilcox.test(.[[1]],.[[2]])}
	}
	get_double_pvalstring <- function(df, contcol,fct1,fct2){
		contcol = enquo(contcol)
		stopifnot(quo_text(contcol)%in%colnames(df))
		fct1 = enquo(fct1)
		fct2 = enquo(fct2)
		paste(sep='\n',
			paste0('Up: wilcox test', get_tidy_pval(wtest(df%>%filter(!(!!(fct1))), contcol, fct2 ) )),
			paste0('down: wilcox test', get_tidy_pval(wtest(df%>%filter(!(!!fct2)), contcol, fct1) ) )
		)
	}
	make_density_threefact_fig <- function(df,classcolname,contcol,fct1,fct2,xname,plotfile){
		fct1 = enquo(fct1)
		fct2 = enquo(fct2)
		classcolname = enquo(classcolname)
		contcol = enquo(contcol)
		fct1text = quo_text(fct1)
		fct2text = quo_text(fct2)
		stopifnot(quo_text(contcol)%in%colnames(df))
		#
		# xlabs = xbreaks[-c(1,length(xbreaks))]
		#
		lclip <- df[[quo_text(contcol)]]%>%quantile(0.01)
		rclip <- df[[quo_text(contcol)]]%>%quantile(1-0.01)
		df%<>%filter(dplyr::between(!!contcol,lclip,rclip))
		#
		xlims <- df[[quo_text(contcol)]]%>%range
		xspan <- xlims%>%dist
		tickwidth <- 10^(round(log10(xspan))-1) * 2.5
		xlimlow <- floor(xlims[1]/(tickwidth))*tickwidth
		xlimhigh <- ceiling(xlims[2]/(tickwidth))*tickwidth
		xlims <- c(xlimlow,xlimhigh)
		xbreaks <- seq(xlimlow,xlimhigh,tickwidth)
		#
		pdf(plotfile)
		{df%>%
			filter(!( (!!fct1) & (!!fct2) ))%>%
			# filter(time!='P0')%>%
			mutate(!!classcolname := case_when(
				!!fct1==1 ~ fct1text,
				!!fct2==1 ~ fct2text,
				TRUE ~ 'neither'
			))%>%
			ggplot(.,aes(x=!!contcol,color=!!classcolname))+geom_density(alpha=I(0.01))+theme_bw()+
				scale_x_continuous(name=xname,limits=xlims,labels=xbreaks,breaks=xbreaks )+
				geom_text(data=data.frame(label=get_double_pvalstring(df,!!contcol,!!fct1,!!fct2)),aes(label=label),color=I('black'),x=Inf,y=Inf,vjust=1,hjust=1)
			}%>%print
		dev.off()
		normalizePath(plotfile)
	}

	timeoccscores <- map_df(.id='time',times,function(itime){
		timeoccvect <- tRNA_occ_df%>%filter(time==itime)%>%mutate(dwell_time=dwell_time-mean(na.rm=T,dwell_time))%>%{setNames(.$dwell_time,.$codon)[codons4occ]}
		(t(usedcodonfreqs[,codons4occ])*(timeoccvect[codons4occ]))%>%colMeans(na.rm=T)%>%
		enframe('protein_id','elongscore')
	})

	fractions <- unique(tRNA_occ_df$fraction) 
	
	timeSigscores <- map_df(.id='time',times,function(itime){
		map_df(.id='fraction',fractions%>%setNames(.,.),function(fractioni){
			timetrrnavect <- tRNA_occ_df%>%
				filter(time==itime,fraction==fractioni)%>%
				 mutate(abundance=ifelse(abundance %in% -Inf,min(abundance[is.finite(abundance)],na.rm=T)-1,abundance))%>%
				mutate(abundance=abundance-mean(na.rm=T,abundance))%>%				
				{setNames(.$abundance,.$codon)[codons4occ]}
			(t(usedcodonfreqs[,codons4occ])*(timetrrnavect[codons4occ]))%>%colMeans(na.rm=T)%>%
			enframe('protein_id','abundance')
		})
	})

	timeaAvailscores <- map_df(.id='time',times,function(itime){
		map_df(.id='fraction',fractions%>%setNames(.,.),function(fractioni){
			timetrrnavect <- tRNA_occ_df%>%
				filter(time==itime,fraction==fractioni)%>%
				 mutate(balance=ifelse(balance %in% -Inf,min(balance[is.finite(balance)],na.rm=T)-1,balance))%>%
				 mutate(balance=balance-mean(na.rm=T,balance))%>%				
				{setNames(.$balance,.$codon)[codons4occ]}
			(t(usedcodonfreqs[,codons4occ])*(timetrrnavect[codons4occ]))%>%colMeans(na.rm=T)%>%
			enframe('protein_id','balance')
		})
	})

	#
	score_techange_df<-timeoccscores%>%
		left_join(timeSigscores,by=c('time','protein_id'))%>%
		left_join(timeaAvailscores,by=c('time','fraction','protein_id'))%>%
		left_join(ms_id2protein_id%>%distinct(protein_id,gene_name))%>%
		left_join(allTEchangedf)

	#Data frame showing elongation change over time for each gene
	occchange_vs_te_df <- score_techange_df%>%
		group_by(up,down,protein_id)%>%
		filter(fraction=='Total')%>%
		# filter(time!='P0')%>%
		# group_slice(1:2)%>%
		nest%>%
		mutate(elongchange = map_dbl(data,~{lm(data=.,elongscore ~ seq_along(time))$coef[2]}))%>%
		mutate(tps = map_dbl(data,nrow))%>%
		select(-data)

	tRNAchange_vs_te_df <- score_techange_df%>%
	filter(fraction=='Total')%>%
		group_by(fraction,up,down,protein_id)%>%
		# filter(time!='P0')%>%
		# group_slice(1:2)%>%
		nest%>%
		mutate(tRNA_score_change = map_dbl(data,~{lm(data=.,abundance ~ seq_along(time))$coef[2]}))%>%
		mutate(tps = map_dbl(data,nrow))%>%
		select(-data)

	avail_change_vs_te_df <- score_techange_df%>%
		filter(fraction=='Total')%>%
		group_by(fraction,up,down,protein_id)%>%
		# filter(time!='P0')%>%
		# group_slice(1:2)%>%
		mutate(ntime=seq_along(time))%>%
		# select(ntime,balance)%>%
		filter(is.finite(balance))%>%
		nest%>%
		mutate(avail_change = map_dbl(data,~lm.fit(matrix(.$balance),matrix(c(rep(1,5),.$ntime),nrow=5,ncol=2))$coef[2]))%>%
		mutate(tps = map_dbl(data,nrow))%>%
		select(-data)



	make_density_threefact_fig(df=occchange_vs_te_df,
		classcolname = TEchange_class,
		contcol = elongchange,
		fct1=up,fct2=down,
		plotfile='plots/figures/figure2/trna_codons/te_change_vs_occscorechange.pdf',
		xname='Predicted Elongation Rate Change - Riboseq Occ'
	)


	make_density_threefact_fig(df=tRNAchange_vs_te_df,
		classcolname = TEchange_class,
		contcol = tRNA_score_change,
		fct1=up,fct2=down,
		plotfile='plots/figures/figure2/trna_codons/te_change_vs_tRNAscorechange.pdf',
		xname='Predicted Change in tRNA Abundance Score'
	)

	make_density_threefact_fig(df=avail_change_vs_te_df,
		classcolname = TEchange_class,
		contcol = avail_change,
		fct1=up,fct2=down,
		plotfile='plots/figures/figure2/trna_codons/te_change_vs_availchange.pdf',
		xname='Predicted Change in tRNA balance Score'
	)

	###So is this the same genes???
	library(txtplot)
	tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%{txtplot(.$tRNA_score_change,.$elongchange)}
	avail_change_vs_te_df%>%left_join(occchange_vs_te_df)%>%{cor.test(.$avail_change,.$elongchange)}

	tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%{txtplot(.$tRNA_score_change,.$elongchange)}
	avail_change_vs_te_df%>%left_join(occchange_vs_te_df)%>%{cor.test(.$tRNA_score_change,.$elongchange)}

	tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%
		filter(elongchange>2,tRNA_score_change<4)



	testvals%>%filter(p.value<0.05)

	library(broom)


	################################################################################
	########Now let's do the same, but with transcriptional change rather than TE change.
	################################################################################
	

	foldchangecatdf <- readRDS(here('data/foldchangecatdf.rds'))

	timeoccscores%<>%select(-matches('gene_id'))
	timeoccscores%<>%safe_left_join(mcols(cds)%>%as.data.frame%>%distinct(protein_id,gene_id))

	allTXNchangedf <- foldchangecatdf%>%group_by(gene_id)%>%
		mutate(up = any((transcriptional_logFC>0)&(transcriptional_adj.P.Val<0.05) ))%>%
		mutate(down = any((transcriptional_logFC<0)&(transcriptional_adj.P.Val<0.05) ))%>%
		summarise(up=any(up),down=any(down))%>%mutate_all(list(~replace_na(.,FALSE)))%>%
		filter(!(up&down))


	score_txnchange_df<-timeoccscores%>%
		left_join(timeSigscores,by=c('time','protein_id'))%>%
		left_join(timeaAvailscores,by=c('time','fraction','protein_id'))%>%
		left_join(ms_id2protein_id%>%distinct(protein_id,gene_name))%>%
		left_join(allTXNchangedf)

	#Data frame showing elongation change over time for each gene
	occchange_vs_txn_df <- score_txnchange_df%>%
		group_by(up,down,protein_id)%>%
		filter(fraction=='Total')%>%
		# filter(time!='P0')%>%
		# group_slice(1:2)%>%
		nest%>%
		mutate(elongchange = map_dbl(data,~{lm(data=.,elongscore ~ seq_along(time))$coef[2]}))%>%
		mutate(tps = map_dbl(data,nrow))%>%
		select(-data)

	tRNAchange_vs_txn_df <- score_txnchange_df%>%
		filter(fraction=='Total')%>%
		group_by(fraction,up,down,protein_id)%>%
		# filter(time!='P0')%>%
		# group_slice(1:2)%>%
		nest%>%
		mutate(tRNA_score_change = map_dbl(data,~{lm(data=.,abundance ~ seq_along(time))$coef[2]}))%>%
		mutate(tps = map_dbl(data,nrow))%>%
		select(-data)

	avail_change_vs_txn_df <- score_txnchange_df%>%
		filter(fraction=='Total')%>%
		group_by(fraction,up,down,protein_id)%>%
		# filter(time!='P0')%>%
		# group_slice(1:2)%>%
		mutate(ntime=seq_along(time))%>%
		select(ntime,balance)%>%
		filter(is.finite(balance))%>%
		arrange(protein_id,fraction)%>%
		nest%>%
		# .$data%>%.[1]
		# %>%select(balance,ntime)
		mutate(avail_change = map_dbl(data,~lm.fit(matrix(.$balance),matrix(c(rep(1,length(.$ntime)),.$ntime),nrow=length(.$ntime),ncol=2))$coef[2]))%>%
		mutate(tps = map_dbl(data,nrow))%>%
		select(-data)



	make_density_threefact_fig(df=occchange_vs_txn_df,
		classcolname = TXN_change,
		contcol = elongchange,
		fct1=up,fct2=down,
		plotfile='plots/figures/figure2/trna_codons/txn_change_vs_occscorechange.pdf',
		xname='Predicted Elongation Rate Change - Riboseq Occ'
	)
	make_density_threefact_fig(df=tRNAchange_vs_txn_df,
		classcolname = TXN_change,
		contcol = tRNA_score_change,
		fct1=up,fct2=down,
		plotfile='plots/figures/figure2/trna_codons/txn_change_vs_tRNAscorechange.pdf',
		xname='Predicted Change in tRNA Abundance Score'
	)
	make_density_threefact_fig(df=avail_change_vs_txn_df,
		classcolname = TXN_change,
		contcol = avail_change,
		fct1=up,fct2=down,
		plotfile='plots/figures/figure2/trna_codons/txn_change_vs_availchange.pdf',
		xname='Predicted Change in tRNA balance Score'
	)

	pdf()
	tRNAchange_vs_txn_df%>%
		filter(!(up&down))%>%
		# filter(time!='P0')%>%
		mutate(Txn_change_class = case_when(
			up==1 ~ 'up',
			down==1 ~ 'down',
			TRUE ~ 'neither'
		))%>%
		ggplot(.,aes(x=tRNA_score_change,color=Txn_change_class))+geom_density(alpha=I(0.01))+theme_bw()+
			scale_x_continuous()
	dev.off()
	normalizePath('plots/figures/figure2/trna_codons/txn_change_vs_tRNAscorechange.pdf')


	#abundance  vs TE change
	#
}

options(error=dump.frames)
purely(fig2,throw_error=F)()

save.image('data/fig2_codon_occ_2.Rdata')