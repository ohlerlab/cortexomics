
# ## Gene level Changes

	isoformab = iso_tx_countdata$abundance%>%
		as.data.frame%>%
		rownames_to_column('tr_id')%>%
		gather(contrast,TPM,-tr_id)%>%
		separate(contrast,into=c('time','assay','rep'))%>%
		group_by(time,assay,tr_id)%>%
		summarise(TPM=mean(TPM))

	trcodusage = codonfreqs%>%
		as.data.frame%>%
		rownames_to_column('tr_id')%>%
		gather(codon,count,-tr_id)

		if(!file.exists(here('data/gcodusagedf.rds'))){
			gcodusagedf = isoformab%>%mutate(g_id=trid2gid[[tr_id]])%>%filter(assay=='total')%>%
				left_join(trcodusage)%>%
				group_by(g_id,time,codon)%>%
				summarise(usage=sum(count*TPM))
			gcodusagedf%<>%group_by(g_id,time)%>%mutate(usage = usage/sum(usage))
			saveRDS(gcodusagedf,here('data/gcodusagedf.rds'))
		}else{
			gcodusagedf<-readRDS(here('data/gcodusagedf.rds'))
		}


		slowcodons = tRNA_occ_df_en%>%filter(time=='E13',fraction=='total')%>%filter(dwell_time>1.5)%>%.$codon
		trcodusage$AA=NULL
		discr_codons = tRNA_occ_df_en%>%distinct(codon,AA)%>%mutate(slow=codon%in%slowcodons)%>%group_by(AA)%>%filter(!all(slow),!all(!slow))
		all_codons = tRNA_occ_df_en%>%distinct(codon,AA)%>%mutate(slow=codon%in%slowcodons)

		slowfracdf = trcodusage%>%
			inner_join(discr_codons,by='codon')%>%
			group_by(tr_id,AA)%>%
			summarise(slowfrac=mean(slow*count)/sum(count))

		slowusagemat = spread(slowfracdf,AA,slowfrac)%>%ungroup%>%select(-tr_id)%>%as.matrix
		slowusagemat = slowusagemat[!apply(slowusagemat,1,function(x)any(is.na(x))),]
		#looks like controlling for amino acid usage doesn't work well.
		slowusagemat%>%princomp%>%.$loadings

		allslowfracdf = trcodusage%>%
			group_by(tr_id)%>%
			mutate(frac=count/sum(count))
		allslowfracdf%>%select(tr_id,codon,frac)%>%spread(codon,frac)
			inner_join(all_codons)%>%
			summarise(slowfrac=mean(slow*count)/sum(count))

		ncodonfreqmat = codonfreqs%>%sweep(1,FUN='/',STAT=rowSums(.))
		
		slowfrac = log2(ncodonfreqmat[,all_codons%>%filter(slow)%>%.$codon]%>%rowSums)
		nonslowfrac = log2(ncodonfreqmat[,all_codons%>%filter(!slow)%>%.$codon]%>%rowSums)


		longtrids = cdsgrl%>%width%>%sum%>%enframe('tr_id','length')%>%mutate(gid=trid2gid[[tr_id]])%>%group_by(gid)%>%slice(which.max(length))%>%.$tr_id

 		geneslowfrac = (slowfrac[longtrids]-nonslowfrac[longtrids])%>%.[is.finite(.)]

 		genemoreslow = geneslowfrac%>%{. > median(.,na.rm=T)}

 		onts = c('BP','MF','CC')
 		
 		source('src/R/Functions/go_term_funcs.R')
		stopifnot(exists("GTOGO"))
		if(!'gene_id'%in%colnames(GTOGO))GTOGO%<>%mutate(gene_id=ensembl_gene_id)
 		get_cluster_gos = function(clustvect){
		    out=map_df(.id='ontology',onts%>%setNames(.,.),function(ont){
		      lapply(unique(clustvect)%>%setNames(.,.),(function(clustval){
		            gids <- names(clustvect)
		            stopifnot(mean(names(clustvect)%in%GTOGO$gene_id)>.9)
		            filtGTOGO <- GTOGO %>%filter(gene_id %in%names(clustvect))
		             projmemoise(rungo)(
		              gids[clustvect==clustval],
		              filtGTOGO,
		              ont,
		              algo='classic'
		            )
		        }))%>%map('result')%>%bind_rows(.id='cluster')
		    })
		    out
		}

 		genemoreslow%>%{names(.)%<>%trid2gid[[.]];.}%>%get_cluster_gos

		codpca$loadings[,1]


		trcodusage%>%group_by(tr_id)%>%group_slice(1:10)%>%group_by(tr_id,AA)

		
		countpred_df<-readRDS('data/countpred_df.rds')
		gidtedf = countpred_df%>%select(g_id=gene_id,contrast,TE=logFC)%>%filter(str_detect(contrast,'TE'))%>%separate(contrast,c('time','tmp'))
		gcodusagedf%>%left_join(gidtedf)




		teform = as.formula(paste0('TE ~ ',paste0(collapse='+',codons4occ)))
		gcodus_te_model = gcodusagedf%>%filter(time=='E13')%>%spread(codon,usage)%>%left_join(gidtedf)%>%lm(data=.,teform)
		temodcoefdf = gcodus_te_model%>%.$coefficients%>%.[codons4occ]%>%enframe('codon','temodel_coef')%>%mutate(time=itime)

		temodcoefdf%>%left_join(codonoccs)%>%{quicktest(log2(.$temodel_coef),.$dwell_time)}
		temodcoefdf%>%left_join(trna_ab_df%>%filter(fraction=='total'))%>%{quicktest(log2(.$temodel_coef),.$dwell_time)}
		stop()
		gcodavaildf = gcodusagedf%>%
			# filter(time=='E13')%>%
			safe_left_join(allow_missing=TRUE,tRNA_occ_df%>%filter(fraction=='total')%>%select(time,codon,availability))%>%
			filter(!is.na(availability))%>%
			group_by(g_id,time)%>%summarise(totalavail = logSumExp(log2(usage)+availability))%>%
			left_join(gidtedf)


		teavail_glm4 = gcodavaildf%>%group_by(g_id)%>%filter(!is.na(TE))%>%group_slice(1:10000)%>%filter(n()==5)%>%lm(data=.,TE ~ g_id+totalavail)
		anova(teavail_glm4)

		glmfit = glm4(count ~ 0 + gene + codon+a_codon+phase, data=sitedf,family=negative.binomial(mtheta),MXITER=400,doFit=T, sparse=T, verbose=T)
		gcodavaildf%>%lm(data=.,TE ~ g_id+totalavail)%>%
			confint


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

	usedcodonfreqsnorm <- usedcodonfreqs[,]%>%{./colSums(t(.))}
	timeoccscores <- map_df(.id='time',times,function(itime){
		timeoccvect <- tRNA_occ_df%>%
			filter(time==itime)%>%
			mutate(dwell_time=dwell_time-mean(na.rm=T,dwell_time))%>%
			{setNames(.$dwell_time,.$codon)[codons4occ]}
		usedcodonfreqsnorm[,names(timeoccvect)]%*%matrix(timeoccvect,ncol=1)%>%.[,1]%>%
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
			timetrrnavect <- timetrrnavect[!is.na(timetrrnavect)]
			
			# (t(usedcodonfreqsnorm[,names(timetrrnavect)])%*%(timetrrnavect))%>%
			# 	colMeans(na.rm=T)%>%
			# 	enframe('protein_id','abundance')

		usedcodonfreqsnorm[,names(timetrrnavect)]%*%matrix(2^timetrrnavect,ncol=1)%>%.[,1]%>%
		log2%>%
		enframe('protein_id','abundance')

		})
	})

	timeaAvailscores <- map_df(.id='time',times,function(itime){
		map_df(.id='fraction',fractions%>%setNames(.,.),function(fractioni){
			timetrrnavect <- tRNA_occ_df%>%
				filter(time==itime,fraction==fractioni)%>%
				 mutate(availability=ifelse(availability %in% -Inf,min(availability[is.finite(availability)],na.rm=T)-1,availability))%>%
				 mutate(availability=availability-mean(na.rm=T,availability))%>%				
				{setNames(.$availability,.$codon)[codons4occ]}
timetrrnavect <- timetrrnavect[!is.na(timetrrnavect)]

				usedcodonfreqsnorm[,names(timetrrnavect)]%*%matrix(2^timetrrnavect,ncol=1)%>%.[,1]%>%
				log2%>%
				enframe('protein_id','availability')

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
		mutate(elongchange = map_dbl(data,~{lm(data=.,elongscore/mean(elongscore) ~ seq_along(time))$coef[2]}))%>%
		mutate(tps = map_dbl(data,nrow))%>%
		select(-data)

	tRNAchange_vs_te_df <- score_techange_df%>%
	filter(fraction=='Total')%>%
		group_by(fraction,up,down,protein_id)%>%
		# filter(time!='P0')%>%
		# group_slice(1:2)%>%
		nest%>%
		mutate(tRNA_score_change = map_dbl(data,~{lm(data=.,(abundance) ~ seq_along(time))$coef[2]}))%>%
		mutate(tps = map_dbl(data,nrow))%>%
		select(-data)



	avail_change_vs_te_df <- score_techange_df%>%
		filter(fraction=='Total')%>%
		group_by(fraction,up,down,protein_id)%>%
		# filter(time!='P0')%>%
		# group_slice(1:2)%>%
		mutate(ntime=seq_along(time))%>%
		# select(ntime,availability)%>%
		filter(is.finite(availability))%>%
		nest%>%
		mutate(avail_change = map_dbl(data,~lm.fit(matrix(.$availability),matrix(c(rep(1,5),.$ntime),nrow=5,ncol=2))$coef[2]))%>%
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
		xname='Predicted Change in tRNA availability Score'
	)

	###So is this the same genes???
	library(txtplot)
	# tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%{txtplot(.$tRNA_score_change,.$elongchange)}
	# avail_change_vs_te_df%>%left_join(occchange_vs_te_df)%>%{cor.test(.$avail_change,.$elongchange)}

	# tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%{txtplot(.$tRNA_score_change,.$elongchange)}
	# avail_change_vs_te_df%>%left_join(occchange_vs_te_df)%>%{cor.test(.$tRNA_score_change,.$elongchange)}

	# tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%
	# 	filter(elongchange>2,tRNA_score_change<4)



	# testvals%>%filter(p.value<0.05)

	# library(broom)


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
		select(ntime,availability)%>%
		filter(is.finite(availability))%>%
		arrange(protein_id,fraction)%>%
		nest%>%
		# .$data%>%.[1]
		# %>%select(availability,ntime)
		mutate(avail_change = map_dbl(data,~lm.fit(matrix(.$availability),matrix(c(rep(1,length(.$ntime)),.$ntime),nrow=length(.$ntime),ncol=2))$coef[2]))%>%
		mutate(tps = map_dbl(data,nrow))%>%
		select(-data)

#

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
		xname='Predicted Change in tRNA availability Score'
	)

	plotfile='plots/figures/figure2/trna_codons/txn_change_vs_tRNAscorechange.pdf'
	pdf(plotfile)
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
	normalizePath(plotfile)
