#!/usr/bin/env Rscript
	# facet_grid(scale='free',ifelse(model%in%names(assay2model),'Seq Data','MS') ~ . )+
message('loading libraries')
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))
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
library(conflicted)
library(zeallot)
library(splines)
library(limma)

conflict_prefer('rowMedians','Biobase')
conflict_prefer('setequal','S4Vectors')
conflict_prefer("between", "dplyr")
conflict_prefer("matches", "dplyr")
matches <- dplyr::matches

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


LOWCOUNTLIM <- 100

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
allsegcounts%<>%separate(sample,into=c('time','assay','rep'))%>%group_by(protein_id,assay,time)%>%mutate(hascounts = sum(total)>LOWCOUNTLIM)
allsegcounts%<>%unite(sample,time,assay,rep)
allsegcounts_nz <- allsegcounts%>%group_by(protein_id)%>%filter(any(hascounts))
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
	assert_that(all(mainribosamps %in% allsegcounts_nz$sample))
})





################################################################################
########Load the matching between protein IDS and mass spec ids
################################################################################
{
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

nonredgnames <- mcols(cds)%>%as.data.frame%>%distinct(gene_name,gene_id)%>%group_by(gene_name)%>%
	mutate(new_gene_name=paste0(gene_name[1],c('',paste0('_',2:n())))[1:n()])%>%
	ungroup%>%
	select(gene_name=new_gene_name,gene_id)

ms_id2protein_id%<>%select(-matches('gene_name'))%>%safe_left_join(nonredgnames)

#keep the mass spec we've matched to protein IDs
matched_ms <- allms%>%semi_join(ms_id2protein_id%>%distinct(ms_id,protein_id))
matched_ms%<>%group_by(ms_id)%>%filter(!all(is.na(signal)))

#get our matrix and design df
c(matchedms_mat,matched_ms_design)%<-% get_matrix_plus_design(matched_ms,ms_id)

#get the satb2 ids we want to look at
allids <- readRDS('pipeline/allids.txt')
# allids2<-allids %>%bind_rows(.,mutate(.,gene_name=paste0(gene_name,'_2')))

# ms_id2protein_id%<>%safe_left_join(allids%>%distinct(protein_id,gene_name))

ms_id2protein_id%<>%select(-matches('gene_name'))%>%safe_left_join(mcols(cds)%>%as.data.frame%>%distinct(transcript_id,gene_id))%>% safe_left_join(nonredgnames,by='gene_id')

}

# ms_id2protein_id%>%mutate(id=seq_len(n()))%>%left_join(allids2%>%distinct(protein_id,gene_name))%>%group_by(id)%>%filter(n()>1)
# #now we can calculate the minimum MS specific variance per gene.
# ms_id2protein_id%>%left_join(allids%>%distinct(gene_id,protein_id))
# ms_id2protein_id%>%mutate(id=seq_len(n()))%>%left_join(allids2%>%distinct(protein_id,gene_name))%>%group_by(id)%>%filter(n()>1)%>%head(2)%>%as.data.frame
# ms_id2protein_id%>%as.data.frame%>%head

#deal with the few cases in which the same gene name is linked to multiple gene IDs


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
assert_that(posteriorsum%>%filter(ms_id %in% satb2ms_ids)%>%.$ms_id%>%unique%>%n_distinct%>%`>`(1))
assert_that(posteriorsum%>%filter(ms_id %in%'B7FAU9;Q8BTM8;B7FAV1')%>%.$precision%>%head(1)%>%`>`(100))



#Pick most precise MS
pickms_ids<-posteriorsum%>%
	group_by(ms_id)%>%summarise(mean_prec=mean(precision))%>%
	left_join(ms_id2protein_id%>%distinct(gene_id,ms_id))%>%
	group_by(gene_id)%>%
	mutate(mostprec_ms_id=mean_prec==max(mean_prec))

#always a most precise ms_id for a gene
assert_that((pickms_ids%>%filter(sum(mostprec_ms_id)>1)%>%nrow) ==0)
pickms_ids <- pickms_ids%>%filter(mostprec_ms_id)
posteriorsum%<>%mutate(assay='MS')

#maybe the count values get tagged, or the 
ms_id2protein_id%<>%ungroup%>%mutate(uprotein_id = paste0(protein_id,'_',as.numeric(factor(ms_id))))
assert_that(!ms_id2protein_id$uprotein_id%>%anyDuplicated)

msposteriorsumfilt <- posteriorsum%>%
	filter(ms_id%in%pickms_ids$ms_id)%>%
	left_join(ms_id2protein_id,allow_missing=TRUE,allow_dups=TRUE)%>%
	filter(!is.na(protein_id))%>%
	filter(protein_id %in% allsegcounts_nz$protein_id)

#add the protein mean estimates
postmeanmat <- msposteriorsumfilt%>%
	select(uprotein_id,mean,time)%>%
	spread(time,mean)%>%
	{structure(as.matrix(.[,-1]),.Dimnames=list(.[,1],colnames(.)[-1]))}%>%
	{colnames(.)%<>%paste0('_MS');.}

postprecmat <- msposteriorsumfilt%>%
	select(uprotein_id,precision,time)%>%
	spread(time,precision)%>%
	{structure(as.matrix(.[,-1]),.Dimnames=list(.[,1],colnames(.)[-1]))}%>%
	{colnames(.)%<>%paste0('_MS');.}

flnauid <- ms_id2protein_id%>%filter(ms_id=='B7FAU9;Q8BTM8;B7FAV1')%>%.$uprotein_id%>%head(1)
assert_that(postprecmat[flnauid,]%>%head(1)%>%`>`(100))

################################################################################
########Perform linear modeling to select the best protein IDs
################################################################################
tps <- mainribosamps%>%str_extract(.,'[^_]+')%>%unique
SPLINE_N <- 4
# RIBOSIGCOL <- 'centercount'
RIBOSIGCOL <- 'total'
allsegcounts_nz%>%colnames
best_satb2_uid <- 'ENSMUSP00000110057_4528'
best_satb2_pid <- 'ENSMUSP00000110057'
satb2cds <- cds%>%subset(protein_id %in% satb2ids)
satb2cdsdj <- satb2cds%>%disjoin(with=TRUE)

{
####
ribocounts <- allsegcounts_nz%>%select(protein_id,dataset=sample,signal=!!RIBOSIGCOL)%>%filter(dataset%in%mainribosamps)%>%mutate(signal=replace_na(signal,0))
c(ribocountmat,ribodesign)%<-% get_matrix_plus_design(ribocounts,protein_id,transform=identity,sigcol=signal)

hasnoribo <- ribocountmat%>%rowSums%>%`==`(0)
noriboprotids <- names(hasnoribo)[hasnoribo]
noribobutms <- noriboprotids%>%intersect(ms_id2protein_id$protein_id)
ribocountmat <- ribocountmat[!hasnoribo,]
countvoom <- limma::voom(ribocountmat,design=model.matrix(~ns(as.numeric(time),2), ribodesign))

#vooom object for all counts
allcountcounts <- allsegcounts_nz%>%select(protein_id,dataset=sample,signal=!!RIBOSIGCOL)%>%filter(dataset%>%str_detect('ribo|total'))%>%mutate(signal=replace_na(signal,0))
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

# oldsatb2counts<-allcountmat[best_satb2_pid,]
allcountmat[best_satb2_pid,]

}

{

featuredata = rownames(allcountmat)%>%data.frame(protein_id=.)%>%
	safe_left_join(as.data.frame(mcols(cds))%>%distinct(protein_id,gene_id,transcript_id))%>%
	safe_left_join(nonredgnames)%>%
	set_rownames(.$protein_id)

featuredata$is_gid_highest <- data.frame(
    gid = featuredata$gene_id,
    pid = featuredata$protein_id,
    rowmed = allcountmat%>%rowMedians
  )%>%  
  mutate( origorder = seq_len(n()))%>%
  group_by(gid)%>%
  sample_frac(1)%>%
  mutate(is_gid_highest = seq_len(n())==which.max(rowmed))%>%
  ungroup%>%arrange(origorder)%>%
  .$is_gid_highest

featuredata$length = width(cds)[match(featuredata$protein_id,cds$protein_id)]

library(txtplot)

allcountmat%>%as.data.frame%>%select(dplyr::matches('total'))%>%as.matrix%>%rowMedians%>%add(1)%>%log10%>%txtdensity


countexprdata <- ExpressionSet(
	allcountmat,
	AnnotatedDataFrame(allcountdesign%>%as.data.frame%>%set_rownames(colnames(allcountmat))),
	AnnotatedDataFrame(featuredata)
)
countexprdata%>%saveRDS(here('pipeline/exprdata/countexprset.rds'))

}


all(countvoom$E%>%rownames%in% fData(countexprdata)$protein_id)

################################################################################
########construct voom object 
################################################################################
	
{
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


mscountvoom <- limma::voom(
	cbind(
		allcountmat[uprot2protinds,]%>%set_rownames(rownames(postmeanmat)),
		2^postmeanmat
	),
	design=mscountmodelmat
)

mscoln <- ncol(postprecmat)
totcoln <- ncol(mscountvoom$weights)


ms_lib_sizes<-tail(mscountvoom$targets$lib.size/1e6,5)
postprecsdmat<-sqrt(1/postprecmat)
scaled_ms_sds <- sweep(postprecsdmat,2,ms_lib_sizes,`/`)

scaled_ms_vars <- scaled_ms_sds^2
scaled_ms_precisions <- 1 / (scaled_ms_vars)
# posteriorsum%>%left_join(ms_id2protein_id)%>%filter(uprotein_id=='ENSMUSP00000000049_2085')

mscountvoom$weights%>%colnames
mscountvoom$weights[,1:(totcoln-mscoln)] <- countonlyweights
# mscountvoom$weights[,(totcoln-mscoln+1):totcoln] <- scaled_ms_precisions[,]
mscountvoom$weights[,(totcoln-mscoln+1):totcoln] <- postprecmat[,]
}

test_that("the library size numbers in a limma object work as I think they do",{

	rnaunscaled <- mscountvoom$E[,1]%>%head%>%{2^.}%>%multiply_by(mscountvoom$targets$lib.size[1]/1000000)
	rnaactual <- allcountmat[uprot2protinds,][,1]%>%head

	msunscaled <- mscountvoom$E[,21]%>%tail%>%{2^.}%>%multiply_by(mscountvoom$targets$lib.size[21]/1000000)
	msactual <- postmeanmat[,][,1]%>%{2^.}%>%tail

	expect_true(between(mean(rnaactual-rnaunscaled),-1,1))
	expect_true(between(mean(abs((msactual-msunscaled)/msunscaled)),-0.001,0.001))

})

voomeffects <- mscountvoom$design%>%colnames

assert_that(identical(c("(Intercept)", "TE", "MS_dev", "ns(time, SPLINE_N)1", "ns(time, SPLINE_N)2",
"ns(time, SPLINE_N)3", "ns(time, SPLINE_N)4", "TE:ns(time, SPLINE_N)1",
"TE:ns(time, SPLINE_N)2", "TE:ns(time, SPLINE_N)3", "TE:ns(time, SPLINE_N)4",
"MS_dev:ns(time, SPLINE_N)1", "MS_dev:ns(time, SPLINE_N)2", "MS_dev:ns(time, SPLINE_N)3",
"MS_dev:ns(time, SPLINE_N)4"),voomeffects))



################################################################################
########Manipulate the effects from spline space to linear space 
################################################################################
	
{
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
)%>%set_rownames(voomeffects)%>%set_colnames(paste0('all_',tps))
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

ntps<-length(tps)
sntps <- seq_along(tps)
i_n <- 1

itimecontrasts=alltimeeff

stepwise_contrasts <- lapply(list(all=alltimeeff,TE=timeTEeffect,MS_dev=timeMSeffect),function(itimecontrasts){
	lapply(1:(ntps-1),function(i_n){
		cname <- colnames(itimecontrasts)[i_n+1]%>%str_replace('_','from_')
		(itimecontrasts[,i_n+1,drop=FALSE] - (itimecontrasts[,i_n,drop=FALSE])) %>%
		set_colnames(cname)
	})%>%do.call(cbind,.)
})%>%do.call(cbind,.)





contrastmatall <- diag(length(voomeffects))[,1:nosplinen]%>%set_colnames(voomeffects[1:nosplinen])%>%cbind(alltimeeff,timeTEeffect,timeMSeffect)
contrastmat <- contrastmatall%>% . [,colSums(.!=0)!=0]#remove the zero cols at first time point

counteffs <- colnames(countvoom$design)
countnosplinen <- counteffs%>%str_detect(negate=TRUE,'ns\\(')%>%which%>%tail(1)
contrastmatcounts <- diag(length(counteffs))[,1:countnosplinen]%>%set_colnames(counteffs[1:countnosplinen])%>%
	cbind(head(alltimeeff[-3,],-SPLINE_N),head(timeTEeffect[-3,],-SPLINE_N),head(timeMSeffect[-3,],-SPLINE_N))
contrastmatcounts%<>% . [,colSums(.!=0)!=0]#remove the zero cols at first time point

itimeeffs <- contrastmat%>%colnames

}

################################################################################
########model fitting
################################################################################	

{
#fit 

message('fitting linear model to counts and inferred MS CIs, using splines')
mscountlm <- lmFit(mscountvoom,design = mscountvoom$design)
mscountebayes <- eBayes(mscountlm)
# mscountebayescontr <- eBayes(contrasts.fit(mscountlm,contrastmat))

message('fitting linear model to counts and inferred MS CIs, indep time points')
itime_mscountlm <- lmFit(mscountvoom,design = itime_modelmat)
itime_mscountebayes<-eBayes(itime_mscountlm)

message('fitting linear model to counts, with splines')
#create voom objects
countebayes <- eBayes(lmFit(countvoom))
message('fitting linear model to counts, indep timepoints')
itimecountebayes <- eBayes(lmFit(itimecountvoom)) 

#Now test time dependence of MSDEV
best_uprotein_ids <- mscountebayes%>%topTable(coef=12:15,number=Inf)%>%rownames_to_column('uprotein_id')%>%
	safe_left_join(ms_id2protein_id%>%distinct(gene_id,uprotein_id))%>%
	group_by(gene_id)%>%slice(which.max(P.Value))%>%
	.$uprotein_id

#Now test time dependence of MSDEV
mscountebayes%>%topTable(coef=12:15,number=Inf)%>%rownames_to_column('uprotein_id')%>%
	safe_left_join(ms_id2protein_id%>%distinct(gene_id,uprotein_id))%>%
	group_by(gene_id)%>%slice(which.max(P.Value))%>%
	.$uprotein_id

#And redo the model with the best pairs
bestmscountebayes <- eBayes(lmFit(mscountvoom[best_uprotein_ids,]))


ebayes_stepwise <- contrasts.fit(bestmscountebayes,stepwise_contrasts)



stpallcoefdf <- ebayes_stepwise%>%topTable(number=Inf)%>%rownames_to_column('uprotein_id')%>%select(-AveExpr,-F,-P.Value,-adj.P.Val)


# totalsatb2steps<-stpallcoefdf%>%safe_left_join(ms_id2protein_id%>%distinct(uprotein_id,gene_name))%>%filter(gene_name=='Satb2')
stpallcoefdf%>%safe_left_join(ms_id2protein_id%>%distinct(uprotein_id,gene_name))%>%filter(gene_name=='Satb2')
ms_id2protein_id%>%filter(uprotein_id%in%best_uprotein_ids)%>%filter(gene_name=='Satb2')%>%as.data.frame

ebayes_stepwiseold<-ebayes_stepwise

}



stop()

####Now we select the 




uprotein_ms_diffs <- contrasts.fit(mscountebayes,contrasts = timeMSeffect[,2:5])%>%topTable(coef=which(tps=='P0')-1,number=1e9,confint=0.95)%>%as.data.frame%>%
	rownames_to_column('uprotein_id')%>%safe_left_join(ms_id2protein_id%>%distinct(ms_id,gene_name,gene_id,uprotein_id))%>%group_by(gene_id)%>%
	mutate(bestprotmatch = abs(logFC) == min(abs(logFC)))


testtp='E175'

ms_genes_w_sig_TE <- lapply(tps[-1],function(testtp){
	eBayes(contrasts.fit(lmFit(mscountvoom[bestuprotids,]),contrasts = timeTEeffect[,2:5]))%>%topTable(coef=which(tps==testtp)-1,number=1e9,confint=0.95)%>%rownames_to_column('uprotein_id')%>%
	safe_left_join(ms_id2protein_id%>%distinct(ms_id,gene_name,gene_id,uprotein_id))%>%
	filter(uprotein_id%in%bestuprotids)%>%
	filter(adj.P.Val < 0.05)
})%>%bind_rows%>%.$gene_name%>%n_distinct

#Fix the to top table here
count_te_fit <- eBayes(contrasts.fit(lmFit(countvoom[,]),contrasts = head(timeTEeffect[-3,2:5],-4)))
count_te_coefs <- lapply(tps[-1]%>%setNames(.,.),function(testtp){count_te_fit%>% topTable(coef=which(tps==testtp)-1,number=1e9,confint=0.95)%>%rownames_to_column('protein_id')})%>%bind_rows(.,.id='time')
count_te_coefs%<>%safe_left_join(mcols(cds)%>%as.data.frame%>%distinct(gene_name,gene_id,protein_id),by='protein_id')

genes_w_sig_TE_df <- count_te_coeffs %>%
		group_by(gene_id)%>%
		slice(which.min(adj.P.Val))%>%
		filter(adj.P.Val < 0.05)%>%as.data.frame



allmscontr<-(alltimeeff+timeTEeffect+timeMSeffect)[,-1]

genes_w_sig_ms_changedf <- lapply(tps[-1],function(testtp){
	tmp<- 
			contrasts.fit(mscountebayes,contrasts = allmscontr)%>%
			topTable(coef=which(tps==testtp)-1,number=1e9,confint=0.95)%>%
			rownames_to_column('uprotein_id')
	tmp %>% safe_left_join(ms_id2protein_id%>%distinct(gene_name,gene_id,uprotein_id),by='uprotein_id')%>%
		group_by(gene_id)%>%
		slice(which.min(adj.P.Val))%>%
		filter(adj.P.Val < 0.05)
})%>%bind_rows


genes_w_sig_ms_changedf <- lapply(tps[-1],function(testtp){
	tmp<- 
			contrasts.fit(mscountebayes,contrasts = allmscontr)%>%
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
	genes_w_sig_TE$gene_id%>%n_distinct

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
	djsatb2cds$ stack(djsatb2cds$revmap)%>%{split(unlist(as(satb2cds$protein_id,"CharacterList")[.$value]),.$name)}

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
	bestpids <- bestuprotids%>%str_replace('_\\d+$','')
	itimecoefs <- itime_mscountebayes$coef%>%colnames%>%str_subset('time')
	itimevarpca <- itime_mscountebayes$coef[,itimecoefs]%>%princomp
	ieigs <- itimevarpca$sdev^2
	ivarexplained <- ieigs / sum(ieigs)

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
	allsegcounts%>%filter(protein_id==satb2ids[2])%>%filter(sample%in%mainribosamps)%>%select(sample,!!RIBOSIGCOL)%>%spread(sample,centercount)%>%.[mainribosamps]

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





	test_that("Our use of the contrast function really is correct - Satb2 looks right",{})



	#linear effects
	(mscountebayes$coef %*% t(mscountebayes$design)) %>% set_colnames(mscountvoomdesign$dataset)

	(mscountebayes$coef %*% t(mscountebayes$design))

})


# exprdf <- normcountdf%>%
# 	separate(dataset,c('time','assay','rep'))%>%
# 	left_join(ms_id2protein_id%>%distinct(uprotein_id,protein_id))%>%
# 	mutate(rep = as.numeric(rep))%>%
# 	bind_rows(msdf)%>%
# 	as_tibble

satb2sig <- countsmatnorm[satb2_uids[1],]%>%.[c(1:4,17:20)]
satb2sig%>%{mean(.[1:2])-mean(.[3:4])}
satb2sig%>%{mean(.[5:6])-mean(.[7:8])}


# allcountcounts <- allsegcounts_nz%>%select(protein_id,dataset=sample,signal=centercount)%>%filter(dataset%>%str_detect('ribo|total'))%>%mutate(signal=replace_na(signal,0))
# c(allcountmat,allcountdesign)%<-% get_matrix_plus_design(allcountcounts,protein_id,transform=identity,sigcol=signal)



get_predictions <- function(mscountebayes,mscountvoomdesign){
	#df of predictions from limma model
	datagroup_names <- mscountvoomdesign%>%.$dataset%>%str_replace('_\\d+$','')%>%unique%>%setNames(.,.)
	sample_contrasts<-mscountebayes$design%>%unique%>%set_rownames(datagroup_names)%>%t
	datagroup<-datagroup_names[1]
	prediction_df<-	lapply(datagroup_names,function(datagroup){
		message(datagroup)
		topTable(contrasts.fit(mscountebayes,sample_contrasts[,datagroup,drop=F]),coef=1,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('uprotein_id')
	})%>%bind_rows(.id='datagroup')
	prediction_df%<>%as_tibble
	prediction_df%<>%separate(datagroup,c('time','assay'))
}

{

prediction_df <- get_predictions(mscountebayes,mscountvoomdesign)
countpred_df <- get_predictions(countebayes,mscountvoomdesign%>%filter(assay!='MS'))

#process our mass spec data for plotting
postprecdf<-postprecmat%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,precision,-uprotein_id)%>%mutate(rep=NA)
postmeandf <- postmeanmat%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,means,-uprotein_id)%>%mutate(rep=NA)
procmsdf <- postprecdf%>%left_join(postmeandf)
msdf <- matched_ms%>%left_join(ms_id2protein_id%>%distinct(ms_id,uprotein_id,protein_id))
msdf%<>%select(uprotein_id,protein_id,dataset,signal,time,rep=replicate)%>%mutate(assay='MS')
stopifnot(!ms_id2protein_id$uprotein_id%>%anyDuplicated)
msdf$signal%<>%log2

countsmatnorm <- allcountmat %>% {sweep(.,2,STATS = DESeq2::estimateSizeFactorsForMatrix(.),FUN='/')}
countsmatnorm <- mscountvoom$E[,]
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
	topTable(contrasts.fit(mscountebayes,time_eff_contrasts[,effect]),coef=1,number=Inf,confint=.95)%>%
	as.data.frame%>%rownames_to_column('uprotein_id')
})%>%bind_rows(.id='effect')
time_eff_contrastsdf%<>%separate(effect,c('assay','time'))
time_eff_contrastsdf$time%>%head

#Extract confidence intervals for time points
datagroup_names <- mscountvoomdesign%>%.$dataset%>%str_replace('_\\d+$','')%>%unique%>%setNames(.,.)
sample_contrasts<-mscountebayes$design%>%unique%>%set_rownames(datagroup_names)%>%t
datagroup<-datagroup_names[1]
prediction_df<-	lapply(datagroup_names,function(datagroup){
	message(datagroup)
	prediction_ob$coef
	topTable(contrasts.fit(mscountebayes,sample_contrasts[,datagroup,drop=F]),coef=1,number=Inf,confint=.95)%>%
	as.data.frame%>%rownames_to_column('uprotein_id')
})%>%bind_rows(.id='datagroup')
prediction_df%<>%as_tibble
prediction_df%<>%separate(datagroup,c('time','assay'))



msrescale2lfq<-(postmeanmat%>%colMedians(na.rm=T)%>%median) -(mscountvoom$E[,21:25]%>%colMedians%>%median)


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

{
testname <- 'Agap2'
assert_that(testname %in% cds$gene_name)
testpids <- cds%>%subset(gene_name==testname)%>%.$protein_id%>%unique
test_uids<-unique(exprdf$uprotein_id)%>%str_subset(testpids%>%paste0(collapse='|'))
#get data for that test gene
ggdf <- exprdf%>%filter(uprotein_id%in%test_uids)
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
#
trajfile = './plots/tmp.pdf'
#
pdf(trajfile)
for(testuid in test_uids){
	uidfilt<-.%>%filter(uprotein_id==testuid)
	plotdf<-ggdf%>%uidfilt
	assays2plot <- unique(ggdf$assay)%>%sort%>%rev
	trajectoryplots<-lapply(assays2plot,function(assay2plot){ggplot(
		data = plotdf%>%filter(assay==assay2plot),
		aes(
			x=as.numeric(as_factor(time)),
			y=signal
		))+
		geom_point()+
		geom_linerange(data=ggdf_msconf%>%uidfilt%>%filter(assay==assay2plot),aes(y=signal,ymin=CI.L,ymax=CI.R))+	
		geom_ribbon(data=prediction_df%>%uidfilt%>%filter(assay==assay2plot),aes(x=as.numeric(as_factor(time)),y=logFC,ymin=CI.L,ymax=CI.R),alpha=I(0.5))+	
		# geom_line(data=prediction_df%>%uidfilt,aes(x=time,y=logFC))+	
		scale_x_continuous(name='Stage',labels=tps)+
		scale_y_continuous(name='Log2 LFQ / Log2 Normalized Counts')+
		theme_bw()+ggtitle(assay2plot)
		# facet_wrap( ~ assay,scales='free')+
	})
	trajectoryplots<-commonyspanplots(trajectoryplots)
	trajectoryplot<-ggarrange(plotlist=trajectoryplots,ncol=3)
	trajectoryplot<-annotate_figure(trajectoryplot,top  = str_interp('Data vs Linear Model - ${testname}'))
	print(trajectoryplot)
}
dev.off()
message(normalizePath(trajfile))
'plots/trajectorys_splinelimma/'%>%dir.create

system(str_interp('cp ${trajfile} plots/trajectorys_splinelimma/${testname}.pdf'))
 }


test_that('the linear modeling with the MS confidence intervals seems to have worked',{
	stop() # not convinced this is the case, e.g. Spr, they seem way to wide
	# I should plot the count only confidence interavls to see....
})

test_that("The confidence intervals for the MS specific effect really are properly dependent on the precision",{
	#TODO plot the confidence intervals of everything our linear model fits.

})

breakint <- 0.5
commonyspanplots  <- function(trajectoryplots,breakint=0.5){
	trajectoryplots_ranges <- map(trajectoryplots,~ggplot_build(.)$layout$panel_scales_y[[1]]$range$range)
	trajectoryplots_range_centers <- trajectoryplots_ranges%>%map_dbl(mean)
	centeredranges <- map2(trajectoryplots_ranges,trajectoryplots_range_centers,`-`)
	centbreakrangemin<-centeredranges%>%map_dbl(1)%>%min%>%divide_by(breakint)%>%floor%>%multiply_by(breakint)
	centbreakrangemax<-centeredranges%>%map_dbl(2)%>%max%>%divide_by(breakint)%>%ceiling%>%multiply_by(breakint)
	trajrangeminsnap <- trajectoryplots_ranges%>%map_dbl(1)%>%map(divide_by,breakint)%>%map(floor)%>%map(multiply_by,breakint)
	trajrangemaxsnap <- trajrangeminsnap%>%map_dbl(add,centbreakrangemax-centbreakrangemin)
	for(i in seq_along(trajectoryplots)){
		trajectoryplots[[i]] = 	trajectoryplots[[i]] + coord_cartesian(ylim=c(trajrangeminsnap[[i]],trajrangemaxsnap[[i]]))
	}
	trajectoryplots
}

maxwidth <- max(trajectoryplots_range_widths)
ymaxs_exp <- trajectoryplots_range_centers + (maxwidth/2) 
ymins_exp <- trajectoryplots_range_centers - (maxwidth/2) 
ymins_exp_breaks <- ymins_exp%>%divide_by(breakint)%>%floor%>%multiply_by(breakint)
ymaxs_exp_breaks <- ymaxs_exp%>%divide_by(breakint)%>%floor%>%multiply_by(breakint)
 


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




################################################################################
########Test that all this has worked
################################################################################



test_that("mostly just a few funny cases where we have ribo but no MS",{
	expect_true(all(has_low_ms))
	gene_name_in_counts <- allms%>%select(ms_id,msgname=gene_name)%>%left_join(ms_id2protein_id%>%distinct(ms_id,protein_id),by='ms_id')%>%group_by(msgname)%>%summarise(any_incounts = any(!is.na(protein_id)))%>%
		{setNames(.$any_incounts,.$msgname)}
	allms%>%filter(gene_name %in% names(gene_name_in_counts[!gene_name_in_counts]))%>%head
	'Znf516'%in%cds$gene_name
})
#Test that our linear modeling worked


test_that("Mappability masks used",{

})

test_that("Xtail's TE should broadly agree with limma",{})
test_that("Xtail's TE should broadly agree with DESeq",{})

test_that("Spline based TE change implemented",{})


export_CIs <- countpred_df%>%
	select(protein_id=uprotein_id,time,assay,signal=logFC,CI.L,CI.R)%>%
	mutate(sd = (((CI.R - CI.L)/2)/1.96) )%>%select(-CI.R,-CI.L)%>%
	inner_join(ms_id2protein_id%>%select(uprotein_id,protein_id,gene_id,gene_name,ms_id)%>%filter(uprotein_id%in%best_uprotein_ids))%>%
	ungroup%>%
	select(-protein_id,-gene_id,-ms_id)%>%
	rbind(
		msposteriorsumfilt%>%
			filter(uprotein_id%in%best_uprotein_ids)%>%
			mutate(sd=sqrt(1/precision))%>%
			select(gene_name,time,assay,signal=mean,sd,uprotein_id)
	)

allgnameshavealldsets <- all(export_CIs%>%group_by(gene_name)%>%tally%>%.$n%>%`==`(length(datagroup_names)))
expect_true(allgnameshavealldsets)

pid_widths<-cds%>%split(.,.$protein_id)%>%width%>%sum%>%stack%>%as.data.frame
pid_widths%<>%set_colnames(c('width','protein_id'))

metainfo <- export_CIs%>%distinct(uprotein_id,gene_name)%>%
	mutate(protein_id = uprotein_id%>%str_replace('_\\d+$',''))%>%
	mutate(sig_TE = str_replace(uprotein_id,'_\\d+','')%in%genes_w_sig_TE$protein_id)%>%
	mutate(sig_MS_change = uprotein_id %in%genes_w_sig_ms_dev$uprotein_id)%>%
	left_join(pid_widths)

	
export_CIs%>%write_tsv(here('pipeline/exprdata/limma_proDD_CIs.tsv'))
metainfo%>%write_tsv(here('pipeline/exprdata/limma_genemetadata.tsv'))
