#!/usr/bin/env Rscript
	# facet_grid(scale='free',ifelse(model%in%names(assay2model),'Seq Data','MS') ~ . )+
message('loading libraries')
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(DESeq2))
suppressMessages(library(assertthat))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(here))
suppressMessages(library(biomaRt))
suppressMessages(library(testthat))
library(zeallot)
library(splines)
library(limma)

get_matrix_plus_design <- function(exprdata,idcol,transform=log2,sigcol=signal){
	idcol<-enquo(idcol)
	sigcol<-enquo(sigcol)

	designmatrix<-exprdata%>%ungroup%>%distinct(dataset)%>%separate(dataset,c('time','assay','rep'),remove=FALSE)
	designmatrix$time <- factor(designmatrix$time)


	myexprmat <- exprdata%>%ungroup%>%select(!!idcol,dataset,!!sigcol)%>%spread(dataset,signal)%>%{set_rownames((as.matrix(.[,-1])),.[[1]])}
	
	myexprmat <- transform(myexprmat)

	list(myexprmat,designmatrix)

}	

sizefactnorm<-function(countmat){
	sizefactors<-DESeq2::estimateSizeFactorsForMatrix(countmat)
	countmat <- countmat %>% {sweep(.,2,STATS = sizefactors[colnames(countmat)],FUN='-')}
}

message('...done' )


LOWCOUNTLIM <- 10

source(here('src/R/Rprofile.R'))

args <- c(
	segment_counts_list = paste0(collapse=',',Sys.glob(here('pipeline/riboseq_quant/data/*/segment_counts_df.tsv'))),
	msfile=here('pipeline',file.path('ms_tables/ms_LFQ_total_ms_tall.tsv')),
	gtf = here('pipeline/my_gencode.vM12.annotation.gtf'),
  transformdexprfile=here('pipeline',file.path('exprdata/transformed_data.txt')),
  designmatrixfile=here('pipeline',file.path('exprdata/designmatrix.txt')),
  normcountstable=here('pipeline','exprdata/allcounts_snorm.tsv')
)
#
for(i in names(args)) assign(i,args[i])
#
if(length(base::commandArgs(trailingOnly=TRUE))) args[] <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
args <- args[!is.na(args)]
message(capture.output(dput(args)))
for(i in names(args)) assign(i,args[i])
#get all the se
segment_counts_list%<>%str_split(string=.,pattern=',')%>%.[[1]]
allsegcounts <- segment_counts_list%>%setNames(.,basename(dirname(.)))%>%map_df(fread,.id='sample')

#get exons
if(!exists('gtf_gr')) gtf_gr<-rtracklayer::import(con=gtf,format='gtf')
exons <- gtf_gr%>%subset(type=='exon')
cds <- gtf_gr%>%subset(type=='CDS')
  


#add gene id
allsegcounts_nz <- allsegcounts%>%group_by(protein_id)%>%filter(any(centercount>10))
allsegcounts_nz %<>% safe_left_join(mcols(gtf_gr)[,c('gene_id','protein_id')]%>%as.data.frame%>%distinct(gene_id,protein_id))

cds_nz <- cds%>%subset(protein_id %in% unique(allsegcounts_nz$protein_id))

unique(allsegcounts_nz$protein_id)


mainribosamps <- '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/sample_parameter.csv'%>%
	fread%>%
	filter(is.na(fraction))%>%
	filter(assay=='ribo')%>%
	filter(sample_id%>%
	str_detect(neg=T,'test'))%>%
	.$sample_id




test_that("It looks like the counting worked!",{
	allsegcounts_nz%>%group_by(protein_id)%>%group_slice(1)%>%as.data.frame
	expect_true('E13_ribo_1' %in% mainribosamps)
	expect_true(all(mainribosamps %in% allsegcounts_nz$sample))
})




################################################################################
########Load the matching between protein IDS and mass spec ids
################################################################################

allms=data.table::fread(msfile)
#some formatting differences
allms%<>%select(ms_id=Protein_IDs,everything())
allms$time%<>%str_replace('p5','5')
allms$dataset%<>%str_replace('p5','5')
allms$dataset%<>%str_replace('_rep','_')
allms$dataset%<>%str_replace('^[^_]+_','')#no need to annotate what signal type it is
allms$dataset%<>%str_replace('total','MS')#no need to annotate what signal type it is


if(!exists('ms_id2protein_id')){
	ms_id2protein_id <- with(new.env(),{

		if(file.exists('data/ms_id2protein_id.rds')){
			ms_id2protein_id <- readRDS('data/ms_id2protein_id.rds')
			# ms_id2protein_id %>%saveRDS('data/ms_id2protein_id.rds')

		}else{
			source(here('src/R/Load_data/get_ms_gencode_tablematch.R'));	
		}
		ms_id2protein_id
	})
}


#keep the mass spec we've matched to protein IDs
matched_ms <- allms%>%semi_join(ms_id2protein_id%>%distinct(ms_id,protein_id))
matched_ms%<>%group_by(ms_id)%>%filter(!all(is.na(signal)))

#get our matrix and design df
c(matchedms_mat,matched_ms_design)%<-% get_matrix_plus_design(matched_ms,ms_id)

#get the satb2 ids we want to look at
allids <- readRDS('pipeline/allids.txt')
ms_id2protein_id%<>%safe_left_join(allids%>%distinct(protein_id,gene_name))
satb2ids <- ms_id2protein_id%>%filter(gene_name=='Satb2')%>%.$protein_id

test_that("We are matching stuff with dashes in the uniprot ids",{
  expect_true(ms_id2protein_id%>%filter(protein_id%in%satb2ids)%>%.$ms_id%>%n_distinct%>%`>`(1))
})
satb2ms_ids<-ms_id2protein_id%>%filter(protein_id%in%satb2ids)%>%.$ms_id%>%unique

test_that("the protein id mapping looks the way we expect",{
	expect_true(!any(ms_id2protein_id%>%distinct(ms_id,protein_id)%>%duplicated))
	expect_true(!any(ms_id2protein_id$uprotein_id%>%duplicated))
})

################################################################################
########Now, use proDD to get our proteomic values as posterios
################################################################################
c(posteriorsum,posteriors) %<-% with(new.env(),{
	if(!exists('posteriors')){
		if(file.exists('data/proDD.data')){
			load(file='data/proDD.data')
		}else{
			source(here('src/R/Load_data/protein_impute.R'))
		}
	}

	stopifnot(setequal(posteriorsum$ms_id, matched_ms$ms_id%>%unique))
	# matched_ms%>%group_by(ms_id,time)%>%summarise(mean(na.omit(signal)))
	list(posteriorsum,posteriors)
})
matched_ms%>%filter(ms_id %in% satb2ms_ids)%>%.$ms_id%>%unique
posteriorsum%>%filter(ms_id %in% satb2ms_ids)%>%.$ms_id%>%unique

posteriorsum%<>%mutate(assay='MS')

stop("need to deal with the non unique protein ids")
#maybe the count values get tagged, or the 
ms_id2protein_id%<>%ungroup%>%mutate(uprotein_id = paste0(protein_id,'_',as.numeric(factor(ms_id))))
assert_that(!ms_id2protein_id$uprotein_id%>%anyDuplicated)

#add the protein mean estimates
postmeanmat <- posteriorsum%>%
	left_join(ms_id2protein_id,allow_missing=TRUE,allow_dups=TRUE)%>%
	filter(!is.na(protein_id))%>%
	filter(protein_id %in% rownames(allcountmat))%>%
	select(uprotein_id,mean,time)%>%
	spread(time,mean)%>%
	{structure(as.matrix(.[,-1]),.Dimnames=list(.[,1],colnames(.)[-1]))}%>%
	{colnames(.)%<>%paste0('_MS');.}

postprecmat <- posteriorsum%>%
	left_join(ms_id2protein_id,allow_missing=TRUE,allow_dups=TRUE)%>%
	filter(!is.na(protein_id))%>%
	filter(protein_id %in% rownames(allcountmat))%>%
	select(uprotein_id,precision,time)%>%
	spread(time,precision)%>%
	{structure(as.matrix(.[,-1]),.Dimnames=list(.[,1],colnames(.)[-1]))}%>%
	{colnames(.)%<>%paste0('_MS');.}



################################################################################
########Perform linear modeling to select the best protein IDs
################################################################################

SPLINE_N <- 4

####
ribocounts <- allsegcounts_nz%>%select(protein_id,dataset=sample,signal=centercount)%>%filter(dataset%in%mainribosamps)%>%mutate(signal=replace_na(signal,0))
c(ribocountmat,ribodesign)%<-% get_matrix_plus_design(ribocounts,protein_id,transform=identity,sigcol=signal)

hasnoribo <- ribocountmat%>%rowSums%>%`==`(0)
noriboprotids <- names(hasnoribo)[hasnoribo]
noribobutms <- noriboprotids%>%intersect(ms_id2protein_id$protein_id)
ribocountmat <- ribocountmat[!hasnoribo,]
countvoom <- limma::voom(ribocountmat,design=model.matrix(~ns(as.numeric(time),2), ribodesign))

#vooom object for all counts
allcountcounts <- allsegcounts_nz%>%select(protein_id,dataset=sample,signal=centercount)%>%filter(dataset%>%str_detect('ribo|total'))%>%mutate(signal=replace_na(signal,0))
c(allcountmat,allcountdesign)%<-% get_matrix_plus_design(allcountcounts,protein_id,transform=identity,sigcol=signal)
allcountmat %<>% replace_na(0)
hasnoallcount <- allcountmat%>%rowSums%>%`==`(0)
noallcountprotids <- names(hasnoallcount)[hasnoallcount]
noallcountbutms <- noallcountprotids%>%intersect(ms_id2protein_id$protein_id)
allcountmat <- allcountmat[!hasnoallcount,]
allcountdesign$time%<>%as.numeric
countmodel = model.matrix(~ ribo*ns(time,SPLINE_N),data=allcountdesign%>%mutate(ribo=assay=='ribo'))
countmodel <- countmodel%>%set_colnames(countmodel%>%colnames%>%str_replace('riboTRUE','TE')%>%str_replace('MSTRUE','MS_dev'))


uprot2protinds <- rownames(postmeanmat)%>%str_replace('_\\d+$','')%>%match(rownames(allcountmat))
allmscountmat <- allcountmat%>%.[uprot2protinds,]

#create voom 
# mscountvoom <- limma::voom(allmscountmat,design=countmodel)
countvoom <- limma::voom(allcountmat,design=countmodel)
itimecountvoom <- limma::voom(allcountmat,design=model.matrix(~ 1+ribo*time,data=allcountdesign%>%mutate(time=tps[time],ribo=assay=='ribo')))







satb2cds <- cds%>%subset(protein_id %in% satb2ids)

satb2cdsdj <- satb2cds%>%disjoin(with=TRUE)

mergeByOverlaps(satb2cds,satb2cdsdj)%>%as.data.frame%>%group_by(protein_id)%>%summarise(width=sum(satb2cdsdj.width))

################################################################################
########construct voom object 
################################################################################
	
# #add our MS vals into the voom object
# voomnames <- rownames(mscountvoom$E)
# #add expression values
# mscountvoom$E%<>%cbind(postmeanmat[voomnames,])
# #add precision weights
# mscountvoom$weights%<>%cbind(postprecmat[voomnames,])
# #add the design matrix
# mscountvoomdesign <- mscountvoom$E%>%
# 	colnames%>%str_split('_')%>%
# 	map(head,2)%>%
# 	simplify2array%>%t%>%as.data.frame%>%
# 	set_colnames(c('time','assay'))%>%
# 	mutate(dataset=colnames(mscountvoom$E),time=factor(time),ribo = assay %in% c('ribo','MS'),MS = assay=='MS')
# mscountvoomdesign$time%<>%as.numeric
# modelmat=model.matrix(~ 1 + ribo+MS+ns(time,3)+ns(time,3):(ribo+MS),data=mscountvoomdesign )
# mscountvoom$design <- modelmat
# #add the lib sizes to the object
# mscountvoom$targets%<>%rbind(colSums(postmeanmat[voomnames,])%>%
# 	stack%>%
# 	{set_rownames(data.frame(lib.size=.$values),.$ind)})

# dim(mscountvoom$E)
# dim(mscountvoom$design)
# dim(mscountvoom$targets)
# dim(mscount)

#add the design matrix

dsets <- colnames(allcountmat)%>%c(colnames(postmeanmat))
mscountvoomdesign <- dsets%>%
	str_split('_')%>%
	map(head,2)%>%
	simplify2array%>%t%>%as.data.frame%>%
	set_colnames(c('time','assay'))%>%
	mutate(dataset=dsets,time=factor(time),ribo = assay %in% c('ribo','MS'),MS = assay=='MS')
mscountvoomdesign$time%<>%as.numeric%>%{.-median(.)}
mscountmodelmat=model.matrix(~ 1 + ribo+MS+ns(time,SPLINE_N)+ns(time,SPLINE_N):(ribo+MS),data=mscountvoomdesign )
mscountmodelmat <- mscountmodelmat%>%set_colnames(mscountmodelmat%>%colnames%>%str_replace('riboTRUE','TE')%>%str_replace('MSTRUE','MS_dev'))
uprot2protinds <- rownames(postmeanmat)%>%str_replace('_\\d+$','')%>%match(rownames(allcountmat))
mscountonlyvoom <- limma::voom(cbind(allcountmat[uprot2protinds,]%>%set_rownames(rownames(postmeanmat))))
countonlyweights <- mscountonlyvoom$weights
mscountvoom <- limma::voom(cbind(allcountmat[uprot2protinds,]%>%set_rownames(rownames(postmeanmat)),postmeanmat),design=mscountmodelmat)
#(countonlyweights - mscountvoom$weights[,1:20] )%>%txtdensity

voomeffects <- mscountvoom$design%>%colnames

assert_that(identical(c("(Intercept)", "TE", "MS_dev", "ns(time, SPLINE_N)1", "ns(time, SPLINE_N)2",
"ns(time, SPLINE_N)3", "ns(time, SPLINE_N)4", "TE:ns(time, SPLINE_N)1",
"TE:ns(time, SPLINE_N)2", "TE:ns(time, SPLINE_N)3", "TE:ns(time, SPLINE_N)4",
"MS_dev:ns(time, SPLINE_N)1", "MS_dev:ns(time, SPLINE_N)2", "MS_dev:ns(time, SPLINE_N)3",
"MS_dev:ns(time, SPLINE_N)4"),voomeffects))



################################################################################
########Manipulate the effects from spline space to linear space 
################################################################################
	
mscoln <- ncol(postprecmat)
totcoln <- ncol(mscountvoom$weights)

mscountvoom$weights%>%colnames
mscountvoom$weights[,1:(totcoln-mscoln)] <- countonlyweights
mscountvoom$weights[,(totcoln-mscoln+1):totcoln] <- postprecmat[,]

itime_modelmat<-model.matrix(~ 1 + ribo+MS+time+time:(ribo+MS),data=mscountvoomdesign%>%mutate(time=factor(tps[add(time,3)])))
itime_modelmat%<>%set_colnames(itime_modelmat%>%colnames%>%str_replace('riboTRUE','TE')%>%str_replace('MSTRUE','MS_dev'))

#for transforming back to time space from spline space
splinemat <- t(ns(1:5,SPLINE_N))
splinezeros <- matrix(0,ncol=5,nrow=SPLINE_N)

effect_is_time <- voomeffects%>%str_detect(negate=TRUE,'ns\\(')
stopifnot(length(rle(effect_is_time))==2)#all main effects, then time
nosplinen <- voomeffects%>%str_detect(negate=TRUE,'ns\\(')%>%which%>%tail(1)
maineffzeros <- matrix(0,ncol=5,nrow=nosplinen)

#contrast matrices transform back to time space
alltimeeff <- rbind(
	maineffzeros,
	t(ns(1:5,SPLINE_N)),
	splinezeros,
	splinezeros
)%>%set_rownames(voomeffects)%>%set_colnames(paste0('all_',tp))
#
timeTEeffect <- rbind(
	matrix(0,ncol=5,nrow=3),
	splinezeros,
	t(ns(-2:2,SPLINE_N)),
	splinezeros
)%>%set_rownames(voomeffects)%>%set_colnames(paste0('TE_',tps))
#
timeMSeffect <- rbind(
	matrix(0,ncol=5,nrow=3),
	splinezeros,
	splinezeros,
	t(ns(1:5,SPLINE_N))
)%>%set_rownames(voomeffects)%>%set_colnames(paste0('MS_dev',tps))

contrastmat <- diag(length(voomeffects))[,1:nosplinen]%>%set_colnames(voomeffects[1:nosplinen])%>%cbind(alltimeeff,timeTEeffect,timeMSeffect)
contrastmat%<>% . [,colSums(.!=0)!=0]#remove the zero cols at first time point

counteffs <- colnames(countvoom$design)
countnosplinen <- counteffs%>%str_detect(negate=TRUE,'ns\\(')%>%which%>%tail(1)
contrastmatcounts <- diag(length(counteffs))[,1:countnosplinen]%>%set_colnames(counteffs[1:countnosplinen])%>%
	cbind(head(alltimeeff[-3,],-SPLINE_N),head(timeTEeffect[-3,],-SPLINE_N),head(timeMSeffect[-3,],-SPLINE_N))
contrastmat%<>% . [,colSums(.!=0)!=0]#remove the zero cols at first time point










message('fitting linear model to counts and inferred MS CIs, using splines')
mscountlm <- lmFit(mscountvoom,design = mscountvoom$design)
mscountebayes <- eBayes(mscountlm)
mscountebayescontr <- eBayes(contrasts.fit(mscountlm,contrastmat))

message('fitting linear model to counts and inferred MS CIs, indep time points')
itime_mscountlm <- lmFit(mscountvoom,design = itime_modelmat)
itime_mscountebayes<-eBayes(itime_mscountlm)

message('fitting linear model to counts, with splines')
#create voom objects
countebayes <- eBayes(lmFit(countvoom))
message('fitting linear model to counts, indep timepoints')
itimecountebayes <- eBayes(lmFit(itimecountvoom)) 


timeeffect<-limma::topTable(itime_mscountebayes,number=1e9,confint=0.95)



#now we can calculate the minimum MS specific variance per gene.
ms_id2protein_id%<>%left_join(allids%>%distinct(gene_id,protein_id))


TODO - big contrast matrix with everything we want - so reparametrizing the splines.

uprotein_ms_diffs <- contrasts.fit(mscountebayes,contrasts = timeMSeffect[,2:5])%>%topTable(coef=which(tps=='P0')-1,number=1e9,confint=0.95)%>%as.data.frame%>%
	rownames_to_column('uprotein_id')%>%safe_left_join(ms_id2protein_id%>%distinct(ms_id,gene_name,gene_id,uprotein_id))%>%group_by(gene_id)%>%
	mutate(bestprotmatch = abs(logFC) == min(abs(logFC)))

bestuprotids <- uprotein_ms_diffs%>%filter(bestprotmatch)%>%.$uprotein_id

testtp='E175'

ms_genes_w_sig_TE <- lapply(tps[-1],function(testtp){
	eBayes(contrasts.fit(lmFit(mscountvoom[bestuprotids,]),contrasts = timeTEeffect[,2:5]))%>%topTable(coef=which(tps==testtp)-1,number=1e9,confint=0.95)%>%rownames_to_column('uprotein_id')%>%
	safe_left_join(ms_id2protein_id%>%distinct(ms_id,gene_name,gene_id,uprotein_id))%>%
	filter(uprotein_id%in%bestuprotids)%>%
	filter(adj.P.Val < 0.05)
})%>%bind_rows%>%.$gene_name%>%n_distinct

genes_w_sig_TE <- lapply(tps[-1],function(testtp){
	tmp<-eBayes(contrasts.fit(lmFit(countvoom[,]),contrasts = head(timeTEeffect[-3,2:5],-4)))%>%topTable(coef=which(tps==testtp)-1,number=1e9,confint=0.95)%>%rownames_to_column('protein_id')
	tmp %>% safe_left_join(mcols(cds)%>%as.data.frame%>%distinct(gene_name,gene_id,protein_id))%>%
		group_by(gene_id)%>%
		slice(which.min(adj.P.Val))%>%
		filter(adj.P.Val < 0.05)
})%>%bind_rows%>%.$gene_name%>%n_distinct



timecoefs = -c(1:3)
timevarpca <- mscountebayes$coef[bestuprotids,timecoefs]%>%princomp
summary(timevarpca)

timevarpca$loadings


##More of the variance is explained 
bestpids <- bestuprotids%>%str_replace('_\\d+$','')
itimecoefs <- itime_mscountebayes$coef%>%colnames%>%str_subset('time')
itimevarpca <- itime_mscountebayes$coef[,itimecoefs]%>%princomp
ieigs <- itimevarpca$sdev^2
ivarexplained <- ieigs / sum(ieigs)


mscountebayes%>%topTable


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
	expect_true(c(7582) == ('pipeline/feature_counts_old/data/P0_total_2/P0_total_2.feature_counts'%>%fread%>%filter(Geneid==satb2gid)%>%.[[7]]))

	#do we need to look at variance inference for limma - can we see the effect with say xtail
	#so for xtail satb2 is definitely sig - and the effect resembles that in our limma model
	'pipeline/xtail/xtail_P0.txt'%>%fread%>%filter(feature_id==satb2gid)

		#Comparing countss	
	maintotsamps<-colnames(allcountmat)%>%str_subset('total')
	'pipeline/exprdata/transformed_data.txt'%>%fread%>%filter(gene_name=='Satb2')
	satb2gid <- cds%>%subset(gene_name=='Satb2')%>%.$gene_id%>%head(1)
	'pipeline/feature_counts/all_feature_counts'%>%fread%>%filter(feature_id==satb2gid)%>%select(mainribosamps)%>%.[mainribosamps]
	allsegcounts%>%filter(protein_id==satb2ids[2])%>%filter(sample%in%mainribosamps)%>%select(sample,total)%>%spread(sample,total)%>%.[mainribosamps]
	allsegcounts%>%filter(protein_id==satb2ids[2])%>%filter(sample%in%mainribosamps)%>%select(sample,centercount)%>%spread(sample,centercount)%>%.[mainribosamps]

	'pipeline/feature_counts/all_feature_counts'%>%fread%>%filter(feature_id==satb2gid)%>%select(maintotsamps)%>%.[maintotsamps]
	allsegcounts%>%filter(protein_id==satb2ids[1])%>%filter(sample%in%maintotsamps)%>%select(sample,total)%>%spread(sample,total)%>%.[maintotsamps]
		

	#what doe sthe final TP effect look like
	#somehow significant even though confidence intervals cross zero
	contrasts.fit(mscountebayes,contrasts = timeTEeffect[,2:5])%>%topTable(coef=which(tps=='P0')-1,number=1e9,confint=0.95)%>%.[unique(satb2_uids),]
	
	mscountebayescontr%>%topTable(coef=which(tps=='P0')-1,number=1e9,confint=0.95)%>%as.data.frame%>%
		rownames_to_column('uprotein_id')%>%safe_left_join(ms_id2protein_id)%>%group_by(gene_name)

	lowconflimsatb2_te <- mscountebayescontr%>%topTable(coef='TE_P0',number=1e9,confint=0.95)%>%.[unique(satb2_uids),]%>%.$CI.L%>%.[1]
	expect_gt(lowconflimsatb2_te,0.5)

	#and the individual splines?


	BiocManager::install('variancePartition')



	test_that("Our use of the contrast function really is correct - Satb2 looks right",{})



	#linear effects
	(mscountebayes$coef %*% t(mscountebayes$design)) %>% set_colnames(mscountvoomdesign$dataset)

	(mscountebayes$coef %*% t(mscountebayes$design))


})

#object list









#process our mass spec data for plotting
postprecdf<-postprecmat%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,precision,-uprotein_id)%>%mutate(rep=NA)
postmeandf <- postmeanmat%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,means,-uprotein_id)%>%mutate(rep=NA)
procmsdf <- postprecdf%>%left_join(postmeandf)
msdf <- matched_ms%>%left_join(ms_id2protein_id%>%distinct(ms_id,uprotein_id,protein_id))
msdf%<>%select(uprotein_id,protein_id,dataset,signal,time,rep=replicate)%>%mutate(assay='MS')
stopifnot(!ms_id2protein_id$uprotein_id%>%anyDuplicated)
msdf$signal%<>%log2


countsmatnorm <- allcountmat %>% {sweep(.,2,STATS = DESeq2::estimateSizeFactorsForMatrix(.),FUN='/')}
normcountdf <- countsmatnorm%>%add(0.5)%>%log2%>%as.data.frame%>%rownames_to_column('protein_id')%>%gather(dataset,signal,-protein_id)%>%as_tibble

exprdf <- normcountdf%>%
	separate(dataset,c('time','assay','rep'))%>%
	left_join(ms_id2protein_id%>%distinct(uprotein_id,protein_id))%>%
	mutate(rep = as.numeric(rep))%>%
	bind_rows(msdf)%>%
	as_tibble

satb2sig <- countsmatnorm[satb2ids[1],]%>%.[c(1:4,17:20)]%>%log2
satb2sig%>%{mean(.[1:2])-mean(.[3:4])}
satb2sig%>%{mean(.[5:6])-mean(.[7:8])}


allcountcounts <- allsegcounts_nz%>%select(protein_id,dataset=sample,signal=centercount)%>%filter(dataset%>%str_detect('ribo|total'))%>%mutate(signal=replace_na(signal,0))
c(allcountmat,allcountdesign)%<-% get_matrix_plus_design(allcountcounts,protein_id,transform=identity,sigcol=signal)





testname <- 'Flna'
assert_that(testname %in% cds$gene_name)
testpids <- cds%>%subset(gene_name==testname)%>%.$protein_id%>%unique
test_uids<-unique(exprdf$uprotein_id)%>%str_subset(testpids%>%paste0(collapse='|'))
#get data for that test gene
ggdf <- exprdf%>%filter(uprotein_id%in%test_uids)
#get ms-protein id pairs
test_uids<-ggdf$uprotein_id%>%unique
#
trajfile = './plots/tmp.pdf'
#
pdf(trajfile)
for(testuid in test_uids){
	trajectoryplot<-ggplot(
		data = ggdf%>%filter(uprotein_id==testuid),
		aes(
			x=as.numeric(as_factor(time)),
			y=signal
		))+
		stat_summary(geom='line',fun.y=mean)+
		geom_point()+

		# geom_line(data=procmsdf,aes(y=log2(means),ymin=log2(means)-(1.96*precision),ymin=log2(means)+(1.96*precision)))+
		scale_x_continuous(name='Stage',labels=tps)+
		scale_y_continuous(name='Log2 LFQ / Log2 Normalized Counts')+
		theme_bw()+
		facet_wrap( ~ assay,scales='free')+
		ggtitle(label = str_interp('Data vs Linear Model - ${testname}')
	)
	print(trajectoryplot)
}
dev.off()
message(normalizePath(trajfile))


timeeffect
	ggplot(data=.,)+
	geom_linerange(
		data=posteriorsum%>%filter(ms_id==names(testrowneg)),
		aes(ymin=ymin,ymax=ymax,x=as.numeric(as.factor(time))))+
	theme_bw()+
	geom_point(data=.,aes(x=as.numeric(as.factor(time)),y=signal))+
	# scale_y_log10()+
	ggtitle(matched_ms%>%filter(ms_id==names(testrowneg))%>%.$gene_name%>%.[1])

test_that('the linear modeling with the MS confidence intervals seems to have worked',{

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


timeeffect<-limma::topTable(mscountebayes,number=sum(!isincomplete),coef=c(2,1+spline_n),confint=0.95)




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

countpredmat <- limma_pred%>%get_matrix_plus_design(gene_name)%>%.[[1]]

####
message("here I should check that the ones with zero coefficients also are just low count things")

countpreddsnested <- countpredmat%>%as.data.frame%>%rownames_to_column('protein_id')%>%group_by(protein_id)%>%nest
mspreddsnested <- mspredmat%>%as.data.frame%>%rownames_to_column('ms_id')%>%group_by(ms_id)%>%nest

ms_ribo_corddf <- ms_id2protein_id%>%inner_join(mspreddsnested)%>%inner_join(countpreddsnested,by='protein_id')%>%
	mutate( ms_cor = map2_dbl(data.x,data.y, ~ possibly(cor,NA)(unlist(.x),unlist(.y),use='complete') ))

top_protein_ids <- ms_ribo_corddf  %>%group_by(ms_id)%>%arrange(desc(ms_cor))%>%slice(1)%>%select(-data.x,-data.y)

exprdata <- rbind(
	top_protein_ids%>%select(ms_id,protein_id)%>%left_join(allsegcounts%>%select(protein_id,dataset=sample,signal=centercount)),
	top_protein_ids%>%select(ms_id,protein_id)%>%left_join(allms%>%select(ms_id,dataset,signal))
	)

c(exprmat,designmat) %<-% (exprdata%>%get_matrix_plus_design(protein_id,identity,signal))

exprmat[1:10,1]%>%{DESeq2::vst(.)}


library(limma)
splinemodel <- limma::lmFit(
	limma::voom(exprmat,design=model.matrix(~ns(as.numeric(time),3), designmat))
)





################################################################################
########Also try correlation of MS with periodicity
################################################################################
specdata <- allsegcounts%>%select(protein_id,dataset=sample,signal=spec_coef)%>%filter(dataset%in%mainribosamps)
c(specmat,specdesign)%<-% get_matrix_plus_design(specdata,protein_id,transform=identity,sigcol=signal)
specmat<-specmat[!apply(specmat,1,function(x) any(is.nan(x))),]
limma_pred <- get_limmafit_predvals(
	limma::lmFit(limma::voom(specmat,design=model.matrix(~ns(as.numeric(time),3), specdesign))),
	specdesign)
specpredmat <- limma_pred%>%get_matrix_plus_design(gene_name)%>%.[[1]]
specpreddsnested <- specpredmat%>%as.data.frame%>%rownames_to_column('protein_id')%>%group_by(protein_id)%>%nest
specms_ribo_corddf <- ms_id2protein_id%>%inner_join(mspreddsnested)%>%inner_join(specpreddsnested,by='protein_id')%>%
	mutate( ms_cor = map2_dbl(data.x,data.y, ~ possibly(cor,NA)(unlist(.x),unlist(.y),use='complete') ))

spec_count_comp_table <- specms_ribo_corddf%>%select(ms_id,protein_id,ms_cor)%>%left_join(ms_ribo_corddf%>%select(protein_id,ms_cor),by='protein_id')%>%
	group_by(ms_id)%>%arrange(desc(ms_cor.y))%>%slice(1)%>%mutate(spec_diff = ms_cor.x - ms_cor.y,specbetter = ms_cor.x > ms_cor.y)


#######

allsegcounts%>%filter(sample==	 'E13_ribo_2')%>%filter(protein_id=='ENSMUSP00000000001')
allsegcounts%>%filter(sample==	 'E13_ribo_1')%>%filter(protein_id=='ENSMUSP00000000001')


#Now collect the final matched dataset


#Now we need to work out which of the protein ids is more correlated with the ms ID.
myexprdata<-matchedms_mat
designmatrix<-matched_ms_design
{
	spline_n=3

	splinetimefit <- limma::lmFit(myexprmat,design=model.matrix(~ns(as.numeric(time),spline_n), designmatrix))
	isincomplete <- splinetimefit$coef%>%apply(1,function(x)sum(is.na(x)))!=0
	splinetimefit <- limma::lmFit(myexprmat[!incomplete_ids,],design=model.matrix(~ns(as.numeric(time),spline_n), designmatrix))


	timeeffect<-limma::topTable(bayesfit,number=sum(!isincomplete),coef=c(2,1+spline_n),confint=0.95)

	timechange_ms_ids <- rownames(timeeffect)[timeeffect$adj.P.Val<0.05]

}

countsmatrix <- countsmatrix %>% {sweep(.,2,STATS = sizefactors[colnames(countsmatrix)],FUN='/')}

countsmatrix_snorm %>%{cbind(gene_name=rownames(.),as_data_frame(.))} %>% write_tsv(normcountstable)




################################################################################
########Test that all this has worked
################################################################################
isUnique <- function(vector){return(!any(duplicated(vector))) ; !any()}
test_that("Most of our high signal MS aspec guys should have a matched protein ID,There should be a 1:1 correspondence",{
	expect_true(all(exprdata$ms_id %in% allms$ms_id))
	
	highsigms_ids <- allms%>%group_by(gene_name,ms_id)%>%summarise(signal=median(na.omit(signal)))%>%slice(which.max(signal))%>%.$ms_id
	expect_true(mean(highsigms_ids %in% ms_id2protein_id$ms_id) > 0.94)
	
	expect_true(!any(is.na(exprdata$ms_id)))
	expect_true(!any(is.na(exprdata$protein_id)))
	expect_true(isUnique(exprdata$ms_id))
	expect_true(isUnique(exprdata$protein_id))
	expect_true(isUnique(exprdata$gene_name))
})

test_that("each of our MS aspec guys should have a matched protein ID,There should be a 1:1 correspondence",{})


test_that("Other protein_id - ms_id links should also be accessible",{
	expect_true(all(exprdata$protein_id %in% countstable$protein_id))
	expect_true(countstable$ms_id %in% allms$ms_id)
	expect_true(setequal(
		countstable%>%group_by(ms_id)%>%slice(which.max(ms_cor))%>%pluck('protein_id'),
		exprdata%>%pluck('protein_id')
	))
})

#There should in general be a good correlation between RNAseq and Riboseq
#In terms of overall levels
test_that("Correlations should look right",{
	expect_true(all(between(countstable$ms_cor,-1,1)))#correlations are there and look right
	expect_true(mean(between(countstable$ms_cor,0,1)) > 0.9) #most of them are positive
})

#Satb2 should show increasing TE
isUnique <- function(vector){return(!any(duplicated(vector))) ; !any()}
test_that("Satb2 should have increasing TE",{
	
	

	cds%>%mcols%>%data.frame%>%distinct(gene_id,gene_name,transcript_id,protein_id)%>%saveRDS('pipeline/allids.txt')
	allids <- readRDS('pipeline/allids.txt')
	ms_id2protein_id%<>%safe_left_join(allids%>%distinct(protein_id,gene_name))
	satb2ids <- ms_id2protein_id%>%filter(gene_name=='Satb2')%>%.$protein_id

	expect_true(allmscountmat[satb2ids,]%>%colMeans%>%last%>%`!=`(0))	
	
	
	satb2timeeffects <- as.matrix(timeeffect[satb2ids,3:5]) %*% t(ns(1:ntp,3))%>%colMeans
	expect_true(satb2timeeffects%>%first < satb2timeeffects%>%last)

	'Satb2' %in% exprdata$gene_name
	satb2ddata <- (exprddata%>%filter(gene_name=='Satb2'))
	'Q8VI24' %in% satb2ddata$ms_id
	
	#the total expression levels should still agree with previous findidngs
	#we'll need a way to compare quantification moddes in our different genes
	te_table%>%filter(time=='P0')%$%{log2fc > 0 & padj < 0.05}
})

test_that("mostly just a few funny cases where we have ribo but no MS",{
	expect_true(length(noribobutms)==38)
	expect_true(all(has_low_ms))
})
#Test that our linear modeling worked



test_that("Mappability masks used",{

})

test_that("Xtail's TE should broadly agree with limma",{})
test_that("Xtail's TE should broadly agree with DESeq",{})

test_that("Spline based TE change implemented",{})

#There should in general be a decent correaltion between the spectral coef and the MS

#' But - shoul we be using correlation smoothed, or just normal corelation???	