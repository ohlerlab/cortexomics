get_predictions <- function(bestmscountebayes,mscountvoomdesign){
	#df of predictions from limma model
	datagroup_names <- mscountvoomdesign%>%.$dataset%>%str_replace('_\\d+$','')%>%unique%>%setNames(.,.)
	sample_contrasts<-bestmscountebayes$design%>%unique%>%set_rownames(datagroup_names)%>%t
	datagroup<-datagroup_names[1]
	prediction_df<-	lapply(datagroup_names,function(datagroup){
		message(datagroup)
		topTable(contrasts.fit(bestmscountebayes,sample_contrasts[,datagroup,drop=F]),coef=1,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('uprotein_id')
	})%>%bind_rows(.id='datagroup')
	prediction_df%<>%as_tibble
	prediction_df%<>%separate(datagroup,c('time','assay'))
}


# maxwidth <- max(trajectoryplots_range_widths)
# ymaxs_exp <- trajectoryplots_range_centers + (maxwidth/2) 
# ymins_exp <- trajectoryplots_range_centers - (maxwidth/2) 
# ymins_exp_breaks <- ymins_exp%>%divide_by(breakint)%>%floor%>%multiply_by(breakint)
# ymaxs_exp_breaks <- ymaxs_exp%>%divide_by(breakint)%>%floor%>%multiply_by(breakint)
 

####Now we select the 




uprotein_ms_diffs <- contrasts.fit(bestmscountebayes,contrasts = timeMSeffect[,2:5])%>%topTable(coef=which(tps=='P0')-1,number=1e9,confint=0.95)%>%as.data.frame%>%
	rownames_to_column('uprotein_id')%>%safe_left_join(ms_id2protein_id%>%distinct(ms_id,gene_name,gene_id,uprotein_id))%>%group_by(gene_id)%>%
	mutate(bestprotmatch = abs(logFC) == min(abs(logFC)))


testtp='E175'


ms_genes_w_sig_TE <- lapply(tps[-1],function(testtp){
	eBayes(contrasts.fit(lmFit(mscountvoom[best_uprotein_ids,]),contrasts = timeTEeffect[,2:5]))%>%topTable(coef=which(tps==testtp)-1,number=1e9,confint=0.95)%>%rownames_to_column('uprotein_id')%>%
	safe_left_join(ms_id2protein_id%>%distinct(ms_id,gene_name,gene_id,uprotein_id))%>%
	filter(uprotein_id%in%best_uprotein_ids)%>%
	filter(adj.P.Val < 0.05)
})%>%bind_rows

#Fix the to top table here
count_te_fit <- eBayes(contrasts.fit(lmFit(countvoom[,]),contrasts = head(timeTEeffect[-3,2:5],-4)))
count_te_coefs <- lapply(tps[-1]%>%setNames(.,.),function(testtp){count_te_fit%>% topTable(coef=which(tps==testtp)-1,number=1e9,confint=0.95)%>%rownames_to_column('protein_id')})%>%bind_rows(.,.id='time')
count_te_coefs%<>%safe_left_join(mcols(cds)%>%as.data.frame%>%distinct(gene_name,gene_id,protein_id),by='protein_id')

genes_w_sig_TE_df <- count_te_coefs %>%
		group_by(gene_id)%>%
		slice(which.min(adj.P.Val))%>%
		filter(adj.P.Val < 0.05)%>%as.data.frame



allmscontr<-(alltimeeff+timeTEeffect+timeMSeffect)[,-1]

genes_w_sig_ms_changedf <- lapply(tps[-1],function(testtp){
	tmp<- 
			contrasts.fit(bestmscountebayes,contrasts = allmscontr)%>%
			topTable(coef=which(tps==testtp)-1,number=1e9,confint=0.95)%>%
			rownames_to_column('uprotein_id')
	tmp %>% safe_left_join(ms_id2protein_id%>%distinct(gene_name,gene_id,uprotein_id),by='uprotein_id')%>%
		group_by(gene_id)%>%
		slice(which.min(adj.P.Val))%>%
		filter(adj.P.Val < 0.05)
})%>%bind_rows


genes_w_sig_ms_changedf <- lapply(tps[-1],function(testtp){
	tmp<- 
			contrasts.fit(bestmscountebayes,contrasts = allmscontr)%>%
			topTable(coef=which(tps==testtp)-1,number=1e9,confint=0.95)%>%
			rownames_to_column('uprotein_id')
	tmp %>% safe_left_join(ms_id2protein_id%>%distinct(gene_name,gene_id,uprotein_id),by='uprotein_id')%>%
		group_by(gene_id)%>%
		slice(which.min(adj.P.Val))%>%
		filter(adj.P.Val < 0.05)
})%>%bind_rows


genes_w_sig_ms_dev <- lapply(tps[-1],function(testtp){
	tmp<- 
		mscountebayescontr%>%
			topTable(.,coef=colnames(timeMSeffect)%>%str_subset(testtp),number=1e9,confint=0.95)%>%
			rownames_to_column('uprotein_id')
	tmp %>% safe_left_join(ms_id2protein_id%>%distinct(gene_name,gene_id,uprotein_id),by='uprotein_id')%>%
		group_by(gene_id)%>%
		slice(which.min(adj.P.Val))%>%
		filter(adj.P.Val < 0.05)
})%>%bind_rows

test_that("Satb2 isoform worries are over?",{
	genes_w_sig_ms_changedf$gene_id%>%n_distinct
	genes_w_sig_ms_dev$gene_id%>%n_distinct
	ms_genes_w_sig_TE$gene_id%>%n_distinct

	'Satb2' %in% genes_w_sig_ms_changedf$gene_name
	genes_w_sig_ms_changedf%>%filter(gene_name=='Satb2')%>%as.data.frame
	genes_w_sig_TE_df%>%filter(gene_name=='Satb2')%>%as.data.frame
	count_te_coefs%>%filter(gene_name=='Satb2')%>%as.data.frame%>%arrange(protein_id)%>%filter(time=='P0')

	ms_id2protein_id%>%filter(gene_name=='Satb2')%>%distinct(gene_id,uprotein_id,gene_name)
	satb2cds$gene_id%>%unique
	satb2cds$protein_id%>%unique
	cds%>%mcols%>%as.data.frame%>%subset(gene_name=='Satb2')%>%distinct(protein_id)
	cds%>%mcols%>%as.data.frame%>%subset(gene_name=='Satb2')%>%distinct(transcript_id)


	satb2cds<-cds%>%subset(gene_name=='Satb2')
		satb2cds%>%split(.,.$protein_id)%>%width%>%sum
	djsatb2cds<-satb2cds%>%disjoin(with=TRUE)
	djsatb2cds$revmap%>%max
	# djsatb2cds$stack(djsatb2cds$revmap)%>%{split(unlist(as(satb2cds$protein_id,"CharacterList")[.$value]),.$name)}

	ms_id2protein_id%>%filter(protein_id=='ENSMUSP00000135163')


	library(rtracklayer)
	library(Rsamtools)

	satb2cdsaaseq <- satb2cds%>%setNames(.,.$protein_id)%>%getSeq(x=FaFile('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'),.)%>%
		split(satb2cds$protein_id)%>%lapply(do.call,what=xscat)%>%DNAStringSet%>%translate

	satb2cdsaaseq%>%writeXStringSet('pipeline/satb2cds.fa')
		'ENSMUSP00000110057'
		'ENSMUSP00000135391'
	satb2cdsaaseq%>%str_extract('HSSAAQA.*?VERVERE')


	##More of the variance is explained 
	bestpids <- best_uprotein_ids%>%str_replace('_\\d+$','')
	itimecoefs <- itime_mscountebayes$coef%>%colnames%>%str_subset('time')
	itimevarpca <- itime_mscountebayes$coef[,itimecoefs]%>%princomp
	ieigs <- itimevarpca$sdev^2
	ivarexplained <- ieigs / sum(ieigs)

	expect_true('Satb2' %in%(genes_w_sig_ms_changedf%>%filter(adj.P.Val<0.05)%>%.$gene_name))

})


test_that("I remember how linear combinations of normal dists work...",{
	a=3;sa=2
# b=2;sb=.1
# pgrid<-expand.grid(a=c(1,10,100),b=-2:2,sa=c(2,4,9),sb=c(.1,1,3))
# sim<-function(a,b,sa,sb){
# 	diffs<-rnorm(10e3,mean=a,sd=sa) - rnorm(10e3,mean=b,sd=sb)
# 	actualdiff<-mean(diffs)
# 	actualsd<-sd(diffs)
# 	tdiff <- a - b
# 	tsd <- sqrt( ((sa^2)) + ((1)*(sb^2)) )
# 	return(list(a,b,sa,sb,actualdiff,tdiff,actualsd,tsd))
# }
# lapply(1:nrow(pgrid),function(ii) do.call(sim,as.list(pgrid[ii,])))%>%simplify2array%>%t


})




################################################################################
########Trying to figure Satb2 diff out
################################################################################
#Maaaayyyyybe the number of multiple tests I'm carrying out here is screwing with the power?
# changing the weights so taht the MS doesn't effect teh count weights has the count
# and the MS TE change estimates very similiar to one another.
# We might need to estimate dispersion trends seperately
# We might also get better results when we redo the splines
# 
test_that("I've identified the cause of this Satb2 diff",{
	#so, is the old score the result of not limiting to cds?
	#oldest - kind of low P0
	satb2gid=nonredgnames%>%filter(gene_name=='Satb2')%>%.$gene_id%>%unique

	expect_true(c(7582) == ('pipeline/feature_counts_old/data/P0_total_2/P0_total_2.feature_counts'%>%fread%>%filter(Geneid==satb2gid)%>%.[[7]]))

	#do we need to look at variance inference for limma - can we see the effect with say xtail
	#so for xtail satb2 is definitely sig - and the effect resembles that in our limma model
	expect_true('pipeline/xtail/xtail_P0.txt'%>%fread%>%filter(feature_id==satb2gid)%>%.$p_value %>%`<`(0.05))

		#Comparing countss	
	maintotsamps<-colnames(allcountmat)%>%str_subset('total')
	'pipeline/exprdata/transformed_data.txt'%>%fread%>%filter(gene_name=='Satb2')
	satb2gid <- cds%>%subset(gene_name=='Satb2')%>%.$gene_id%>%head(1)
	'pipeline/feature_counts/all_feature_counts'%>%fread%>%filter(feature_id==satb2gid)%>%select(mainribosamps)%>%.[mainribosamps]
	allsegcounts%>%filter(protein_id==satb2ids[2])%>%filter(sample%in%mainribosamps)%>%select(sample,total)%>%spread(sample,total)%>%.[mainribosamps]
	allsegcounts%>%filter(protein_id==satb2ids[2])%>%filter(sample%in%mainribosamps)%>%select(sample,!!RIBOSIGCOL)%>%spread(sample,total)%>%.[mainribosamps]

	'pipeline/feature_counts/all_feature_counts'%>%fread%>%filter(feature_id==satb2gid)%>%select(maintotsamps)%>%.[maintotsamps]
	allsegcounts%>%filter(protein_id==satb2ids[1])%>%filter(sample%in%maintotsamps)%>%select(sample,total)%>%spread(sample,total)%>%.[maintotsamps]
		

	#what doe sthe final TP effect look like
	#somehow significant even though confidence intervals cross zero
	# contrasts.fit(mscountebayes,contrasts = timeTEeffect[,2:5])%>%topTable(coef=which(tps=='P0')-1,number=1e9,confint=0.95)%>%.[unique(satb2_uids),]
	
	# mscountebayescontr%>%topTable(coef=which(tps=='P0')-1,number=1e9,confint=0.95)%>%as.data.frame%>%
	# 	rownames_to_column('uprotein_id')%>%safe_left_join(ms_id2protein_id)%>%group_by(gene_name)

	# lowconflimsatb2_te <- mscountebayescontr%>%topTable(coef='TE_P0',number=1e9,confint=0.95)%>%.[unique(satb2_uids),]%>%.$CI.L%>%.[1]
	# expect_gt(lowconflimsatb2_te,0.5)

	#and the individual splines?




	test_that("Our use of the contrast function really is correct - Satb2 looks right",{})

})


satb2sig <- countsmatnorm[satb2_uids[1],]%>%.[c(1:4,17:20)]
satb2sig%>%{mean(.[1:2])-mean(.[3:4])}
satb2sig%>%{mean(.[5:6])-mean(.[7:8])}


{

prediction_df <- get_predictions(bestmscountebayes,mscountvoomdesign)
# countpred_df <- get_predictions(countebayes,mscountvoomdesign%>%filter(assay!='MS'))

#process our mass spec data for plotting
postprecdf<-postprecmat%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,precision,-uprotein_id)%>%mutate(rep=NA)
postmeandf <- postmeanmat%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,means,-uprotein_id)%>%mutate(rep=NA)
procmsdf <- postprecdf%>%left_join(postmeandf)
msdf <- matched_ms%>%left_join(ms_id2protein_id%>%distinct(ms_id,uprotein_id,protein_id))
msdf%<>%select(uprotein_id,protein_id,dataset,signal,time,rep=replicate)%>%mutate(assay='MS')
stopifnot(!ms_id2protein_id$uprotein_id%>%anyDuplicated)
msdf$signal%<>%log2

countsmatnorm <- allcountmat %>% {sweep(.,2,STATS = DESeq2::estimateSizeFactorsForMatrix(.),FUN='/')}
countsmatnorm <- mscountvoom$E[best_uprotein_ids,]
exprdf <- countsmatnorm%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,signal,-uprotein_id)%>%as_tibble
exprdf%<>%separate(dataset,c('time','assay','rep'))%>%left_join(ms_id2protein_id%>%distinct(uprotein_id,protein_id))%>%
	mutate(rep = as.numeric(rep))



#Extract confidence inttervals for effects
time_eff_contrasts <- contrastmatall
effects<-colnames(time_eff_contrasts)
istimete<-colnames(time_eff_contrasts)%>%str_detect('TE_')
time_eff_contrasts[,istimete] %<>%add( time_eff_contrasts[,'TE'])
istimeMSde<-colnames(time_eff_contrasts)%>%str_detect('MS_dev.')
time_eff_contrasts[,istimeMSde] %<>% add(time_eff_contrasts[,'TE'])
time_eff_contrasts[,istimeMSde] %<>% add(time_eff_contrasts[,'TE'])
time_eff_contrasts[,istimeMSde] %<>% add(time_eff_contrasts[,'MS_dev'])
time_eff_contrasts%<>%.[,istimete | istimeMSde]
timeeffnames<-colnames(time_eff_contrasts)%>%setNames(.,.)

time_eff_contrastsdf<-	lapply(timeeffnames,function(effect){
	message(effect)
	topTable(contrasts.fit(bestmscountebayes,time_eff_contrasts[,effect]),coef=1,number=Inf,confint=.95)%>%
	as.data.frame%>%rownames_to_column('uprotein_id')
})%>%bind_rows(.id='effect')
time_eff_contrastsdf%<>%separate(effect,c('assay','time'))
time_eff_contrastsdf%>%head

#Extract confidence intervals for time points
datagroup_names <- mscountvoomdesign%>%.$dataset%>%str_replace('_\\d+$','')%>%unique%>%setNames(.,.)
sample_contrasts<-bestmscountebayes$design%>%unique%>%set_rownames(datagroup_names)%>%t
#add in TE to these
datagroup_names = c(datagroup_names,paste0(tps,'_','TE'))%>%setNames(.,.)
sample_contrasts%<>%cbind(timeTEeffect%>%{.['TE',]<-1;.})
colnames(sample_contrasts)<-datagroup_names
datagroup<-datagroup_names[1]
prediction_df<-	lapply(datagroup_names,function(datagroup){
	message(datagroup)
	# prediction_ob$coef
	topTable(contrasts.fit(bestmscountebayes,sample_contrasts[,datagroup,drop=F]),coef=1,number=Inf,confint=.95)%>%
	as.data.frame%>%rownames_to_column('uprotein_id')
})%>%bind_rows(.id='datagroup')
prediction_df%<>%as_tibble
prediction_df%<>%separate(datagroup,c('time','assay'))

stopifnot(!any(prediction_df$assay%>%unique%>%is_in(tps)))

msrescale2lfq<-(postmeanmat%>%colMedians(na.rm=T)%>%median) - (mscountvoom$E[,21:25]%>%colMedians%>%median)
# msrescale2lfq<-0
}

test_that("predictions look sane!",{
	testorcuid<-'ENSMUSP00000048319_5829'
	orc3mspredisgdf <- prediction_df%>%filter(uprotein_id==testorcuid)
	orc3mspredisgdf%>%left_join(exprdf%>%filter(uprotein_id==testorcuid)%>%select(uprotein_id,signal))%>%
		mutate(within = (CI.L < signal) & (signal < CI.R))%>%
		select(assay,signal,logFC,CI.L,CI.R,within)%>%
		filter(assay=='MS')

	testi <- which(mscountvoom$E%>%rownames%>%`==`(testorcuid))

	exprdf%>%filter(uprotein_id==testorcuid)
	postmeanmat[testorcuid,]-msrescale2lfq
	exprdf%>%filter(assay=='MS',uprotein_id==testorcuid)

	#our data is going into the 
	expect_true(mscountvoom$E[testorcuid,]%>%tail(5)%>%{(.[1] - .[5]) > 1.5})
})

assaynames = c('MS'='Mass Spec','total'='RNAseq','ribo'='Ribo-Seq','TE'='TE')
tpnames = c('E12.5','E14','E15.5','E17','P0')%>%setNames(tps)
breakint <- 0.5
commonyspanplots  <- function(trajectoryplots,breakint=0.25,minrange=NULL){
	trajectoryplots_ranges <- map(trajectoryplots,~ggplot_build(.)$layout$panel_scales_y[[1]]$range$range)
	trajectoryplots_range_centers <- trajectoryplots_ranges%>%map_dbl(mean)
	centeredranges <- map2(trajectoryplots_ranges,trajectoryplots_range_centers,`-`)
	centbreakrangemin<-centeredranges%>%map_dbl(1)%>%min%>%divide_by(breakint)%>%floor%>%multiply_by(breakint)
	centbreakrangemax<-centeredranges%>%map_dbl(2)%>%max%>%divide_by(breakint)%>%ceiling%>%multiply_by(breakint)
	trajrangeminsnap <- (trajectoryplots_range_centers + centbreakrangemin ) %>% divide_by(breakint)%>%floor%>%multiply_by(breakint)
	trajrangemaxsnap <- (trajectoryplots_range_centers + centbreakrangemax ) %>% divide_by(breakint)%>%ceiling%>%multiply_by(breakint)
	i=1
	if(!is.null(minrange)){
		for(i in seq_along(trajrangeminsnap)){
			ydist = trajrangemaxsnap[i] - trajrangeminsnap[i]
			if(ydist < minrange){
				yshift = ceiling(abs(((minrange - ydist) / 2)) / breakint)
				trajrangeminsnap = trajrangeminsnap - (breakint*yshift)
				trajrangemaxsnap = trajrangemaxsnap + (breakint*yshift)
			}
		}
	}
	for(i in seq_along(trajectoryplots)){
		trajectoryplots[[i]] = 	trajectoryplots[[i]] + 
			# scale_y_continuous(name ='Log2 LFQ / Log2 Normalized Counts',breaks = seq(-10,10,by=breakint),labels = round(seq(-10,10,by=breakint), 3))+
			scale_y_continuous(name =ifelse(i==1,'Log Fold Change Vs Median',''),breaks = seq(-10,10,by=1),labels = round(2^seq(-10,10,by=1), 3))+
			coord_cartesian(ylim=c(trajrangeminsnap[[i]],trajrangemaxsnap[[i]]))
	}
	trajectoryplots
}
plotlistlist = list()
genenamelist = list()
genes2plot = c('Nes','Flna','Bcl11b','Satb2')
for(testname in genes2plot){
{
assert_that(testname %in% nonredgnames$gene_name)
testpids <- cds%>%subset(gene_name==testname)%>%.$protein_id%>%unique
test_uids<-unique(exprdf$uprotein_id)%>%str_subset(testpids%>%paste0(collapse='|'))
#get data for that test gene
ggdf <- exprdf%>%filter(uprotein_id%in%test_uids)
ggdf%<>%bind_rows(ggdf%>%filter(assay=='ribo')%>%mutate(assay='TE',signal = signal - (ggdf%>%filter(assay=='total')%>%.$signal)))
postmeanmat_scaled <- postmeanmat - msrescale2lfq
ggdf_msconf <- 
	safe_left_join(
		((postmeanmat_scaled[test_uids,,drop=F])-(1.96*postprecsdmat[test_uids,,drop=F]))%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,CI.L,-uprotein_id)%>%separate(dataset,c('time','assay')),
		((postmeanmat_scaled[test_uids,,drop=F])+(1.96*postprecsdmat[test_uids,,drop=F]))%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,CI.R,-uprotein_id)%>%separate(dataset,c('time','assay'))
	)%>%
	safe_left_join(
		((postmeanmat_scaled[test_uids,,drop=F]))%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,signal,-uprotein_id)%>%separate(dataset,c('time','assay')),
		)
# ggdf_msconf$time%<>%factor%>%time
#get ms-protein id pairs
test_uids<-ggdf$uprotein_id%>%unique
scaledata = function(x,assay2plot){
	if(assay2plot=='MS'){
		out = x%>%mutate_at(vars(one_of(c('signal','logFC','CI.R','CI.L'))),list(function(x)x+msrescale2lfq))
	}else{
		out = x	
	}
	rescale = if('signal' %in% colnames(x)){ median(x$signal[x$time=='E13']) } else {median(x$logFC[x$time=='E13'])}
	out = x %>%mutate_at(vars(one_of(c('signal','logFC','CI.R','CI.L'))),list(function(x)x-rescale))
	out
}
# plotdf%>%filter(assay==assay2plot)%>%scaledata(assay2plot='MS')
# scaledata = identity
#
trajfile = './plots/tmp.pdf'
for(testuid in test_uids){
	testuid %in% best_uprotein_ids
	uidfilt<-.%>%filter(uprotein_id==testuid)
	ms_id2protein_id%>%filter(uprotein_id==testuid)
	plotdf<-ggdf%>%uidfilt
	assays2plot <- unique(ggdf$assay)%>%sort%>%rev%>%.[order(.=='TE')]
	trajectoryplots<-lapply(assays2plot,function(assay2plot){
		if(assay2plot!='MS'){
			points = geom_point()
			linerange=NULL
		}else{
			points = geom_point(data=msdf%>%semi_join(ms_id2protein_id%>%filter(gene_name==testname))%>%scaledata('dont unscale'))
			linerange=geom_linerange(color=I('blue'),data=ggdf_msconf%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),aes(y=signal,ymin=CI.L,ymax=CI.R))
		}
		ggplot(
		data = plotdf%>%filter(assay==assay2plot)%>%scaledata(assay2plot),
		aes(
			x=as.numeric(as_factor(time)),
			y=signal
		))+
		points+linerange+
		geom_ribbon(data=prediction_df%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),aes(x=as.numeric(as_factor(time)),y=logFC,ymin=CI.L,ymax=CI.R),fill='darkgreen',alpha=I(0.5))+	
		geom_line(data=prediction_df%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),linetype=2,aes(x=as.numeric(as_factor(time)),y=logFC))+	
		# geom_line(data=prediction_df%>%uidfilt,aes(x=time,y=logFC))+	
		scale_x_continuous(name='Stage',labels=tpnames)+
		theme_bw()+
		ggtitle(assaynames[assay2plot])
		# facet_wrap( ~ assay,scales='free')+
	})
	trajectoryplots<-commonyspanplots(trajectoryplots,breakint=0.25,minrange=8)
	plotlistlist = append(plotlistlist,list(trajectoryplots))
	genenamelist = append(genenamelist,testname)
}
}
}
genenamelist%<>%unlist
unflatinds = rep(seq_along(plotlistlist),map_dbl(plotlistlist,length))
plotlistlist = plotlistlist%>%flatten%>%commonyspanplots%>%split(unflatinds)%>%setNames(genenamelist)

for(testname in genenamelist){
	trajectoryplot<-ggarrange(plotlist=plotlistlist[[testname]],ncol=4)
	trajectoryplot<-annotate_figure(trajectoryplot,top  = str_interp('Data vs Linear Model - ${testname}'))
	pdf(trajfile,w=12,h=4)
	print(trajectoryplot)
	dev.off()
	system(str_interp('cp ${trajfile} plots/figures/figure2/traject_${testname}.pdf'))
	message(normalizePath(trajfile))
	message(normalizePath(str_interp('plots/figures/figure2/traject_${testname}.pdf')))
}

'plots/trajectorys_splinelimma/'%>%dir.create
#
#
 




test_that('the linear modeling with the MS confidence intervals seems to have worked',{
	stop() # not convinced this is the case, e.g. Spr, they seem way to wide
	# I should plot the count only confidence interavls to see....
})

test_that("The confidence intervals for the MS specific effect really are properly dependent on the precision",{
	#TODO plot the confidence intervals of everything our linear model fits.
})




limmafit <- mscountlm
get_limmafit_predvals <- function(limmafit,designmatrix){
  (limmafit$coef %*% t(limmafit$design))%>%
	  set_colnames(designmatrix$dataset)%>%
	  as.data.frame%>%
	  rownames_to_column('gene_name')%>%
	  gather(dataset,signal,-gene_name)%>%
	  left_join(designmatrix)%>%
	  distinct(gene_name,time,assay,.keep_all = TRUE)%>%
	  select(-dplyr::matches('rep'))%>%
	  mutate(dataset=paste0(as.character(time),'_',assay))
}
get_limmafit_predvals(mscountlm,mscountvoomdesign)

ntp=n_distinct(allms$time)


cds%>%mcols%>%data.frame%>%distinct(gene_id,gene_name,transcript_id,protein_id)%>%saveRDS('pipeline/allids.txt')
allids <- readRDS('pipeline/allids.txt')
ms_id2protein_id%<>%safe_left_join(allids%>%distinct(protein_id,gene_name))
satb2ids <- ms_id2protein_id%>%filter(gene_name=='Satb2')%>%.$protein_id

allmscountmat[satb2ids,]
timeeffect[satb2ids[1],]




test_that("we can get variance estimate from our linear model",{
	timeeffect[1,]%>%select(dplyr::matches('^ns.as.num'))
})


test_that("Our linear modeling is a decent substitute for xtail etc.",{
	timeeffect[1,]%>%select(dplyr::matches('^ns.as.num'))
})






#protein ids in annotaiton
#protein ids counts
#protein IDs with no RNAseq
#protein IDs with no Riboseq
#protein IDs with no matching MS


#add the precision weights

#add the 


# limma_pred <- get_limmafit_predvals(
# 	limma::lmFit(,
# 	countdesign)




# ################################################################################
# ########Also try correlation of MS with periodicity
# ################################################################################
# specdata <- allsegcounts%>%select(protein_id,dataset=sample,signal=spec_coef)%>%filter(dataset%in%mainribosamps)
# c(specmat,specdesign)%<-% get_matrix_plus_design(specdata,protein_id,transform=identity,sigcol=signal)
# specmat<-specmat[!apply(specmat,1,function(x) any(is.nan(x))),]
# limma_pred <- get_limmafit_predvals(
# 	limma::lmFit(limma::voom(specmat,design=model.matrix(~ns(as.numeric(time),3), specdesign))),
# 	specdesign)
# specpredmat <- limma_pred%>%get_matrix_plus_design(gene_name)%>%.[[1]]
# specpreddsnested <- specpredmat%>%as.data.frame%>%rownames_to_column('protein_id')%>%group_by(protein_id)%>%nest
# specms_ribo_corddf <- ms_id2protein_id%>%inner_join(mspreddsnested)%>%inner_join(specpreddsnested,by='protein_id')%>%
# 	mutate( ms_cor = map2_dbl(data.x,data.y, ~ possibly(cor,NA)(unlist(.x),unlist(.y),use='complete') ))

# spec_count_comp_table <- specms_ribo_corddf%>%select(ms_id,protein_id,ms_cor)%>%left_join(ms_ribo_corddf%>%select(protein_id,ms_cor),by='protein_id')%>%
# 	group_by(ms_id)%>%arrange(desc(ms_cor.y))%>%slice(1)%>%mutate(spec_diff = ms_cor.x - ms_cor.y,specbetter = ms_cor.x > ms_cor.y)


# #######

# allsegcounts%>%filter(sample==	 'E13_ribo_2')%>%filter(protein_id=='ENSMUSP00000000001')
# allsegcounts%>%filter(sample==	 'E13_ribo_1')%>%filter(protein_id=='ENSMUSP00000000001')


# #Now collect the final matched dataset


# #Now we need to work out which of the protein ids is more correlated with the ms ID.
# myexprdata<-matchedms_mat
# designmatrix<-matched_ms_design
# {
# 	spline_n=3

# 	splinetimefit <- limma::lmFit(myexprmat,design=model.matrix(~ns(as.numeric(time),spline_n), designmatrix))
# 	isincomplete <- splinetimefit$coef%>%apply(1,function(x)sum(is.na(x)))!=0
# 	splinetimefit <- limma::lmFit(myexprmat[!incomplete_ids,],design=model.matrix(~ns(as.numeric(time),spline_n), designmatrix))


# 	timeeffect<-limma::topTable(bayesfit,number=sum(!isincomplete),coef=c(2,1+spline_n),confint=0.95)

# 	timechange_ms_ids <- rownames(timeeffect)[timeeffect$adj.P.Val<0.05]

# }

# countsmatrix <- countsmatrix %>% {sweep(.,2,STATS = sizefactors[colnames(countsmatrix)],FUN='/')}

# countsmatrix_snorm %>%{cbind(gene_name=rownames(.),as_data_frame(.))} %>% write_tsv(normcountstable)


