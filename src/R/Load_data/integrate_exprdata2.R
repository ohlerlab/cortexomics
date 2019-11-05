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
conflict_prefer("intersect", "BiocGenerics")
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

if(file.exists('data/ms_id2protein_id.rds')){
	ms_id2protein_id <- readRDS('data/ms_id2protein_id.rds')			# ms_id2protein_id %>%saveRDS('data/ms_id2protein_id.rds')
}else{
	source(here('src/R/Load_data/get_ms_gencode_tablematch.R'));	
}

nonredgnames <- mcols(cds)%>%as.data.frame%>%distinct(gene_name,gene_id)%>%group_by(gene_name)%>%
	mutate(new_gene_name=paste0(gene_name[1],c('',paste0('_',2:n())))[1:n()])%>%
	ungroup%>%
	select(gene_name=new_gene_name,gene_id)%>%
	left_join(cds%>%mcols%>%as.data.frame%>%distinct(transcript_id,gene_id))

ms_id2protein_id%<>%{select(.,-matches('gene_name'))}%>%safe_left_join(nonredgnames,by='transcript_id')

#keep the mass spec we've matched to protein IDs
matched_ms <- allms%>%semi_join(ms_id2protein_id%>%distinct(ms_id,protein_id))
matched_ms%<>%group_by(ms_id)%>%filter(!all(is.na(signal)))

#get our matrix and design df
c(matchedms_mat,matched_ms_design)%<-% get_matrix_plus_design(matched_ms,ms_id)

#get the satb2 ids we want to look at
allids <- readRDS('pipeline/allids.txt')
# allids2<-allids %>%bind_rows(.,mutate(.,gene_name=paste0(gene_name,'_2')))

# ms_id2protein_id%<>%safe_left_join(allids%>%distinct(protein_id,gene_name))

ms_id2protein_id%<>%select(-matches('gene_name'))
trgeneids <- mcols(cds)%>%as.data.frame%>%distinct(transcript_id,gene_id)
ms_id2protein_id%<>%safe_left_join(trgeneids)
ms_id2protein_id %<>% safe_left_join(nonredgnames%>%distinct(transcript_id,gene_name),by='transcript_id')

satb2ids <- ms_id2protein_id%>%filter(gene_name=='Satb2')%>%.$protein_id

}

# ms_id2protein_id%>%mutate(id=seq_len(n()))%>%left_join(allids2%>%distinct(protein_id,gene_name))%>%group_by(id)%>%filter(n()>1)
# #now we can calculate the minimum MS specific variance per gene.
# ms_id2protein_id%>%left_join(allids%>%distinct(gene_id,protein_id))
# ms_id2protein_id%>%mutate(id=seq_len(n()))%>%left_join(allids2%>%distinct(protein_id,gene_name))%>%group_by(id)%>%filter(n()>1)%>%head(2)%>%as.data.frame
# ms_id2protein_id%>%as.data.frame%>%head

#deal with the few cases in which the same gene name is linked to multiple gene IDs



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
mscountebayescontr <- eBayes(contrasts.fit(mscountlm,contrastmat))

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
best_protein_ids = best_uprotein_ids%>%str_replace('_\\d+$','')
# satb2ids %in% best_protein_ids
# satb2ids %>%intersect(best_protein_ids)


ebayes_stepwise <- contrasts.fit(bestmscountebayes,stepwise_contrasts)



stpallcoefdf <- ebayes_stepwise%>%topTable(number=Inf)%>%rownames_to_column('uprotein_id')%>%select(-AveExpr,-F,-P.Value,-adj.P.Val)


# totalsatb2steps<-stpallcoefdf%>%safe_left_join(ms_id2protein_id%>%distinct(uprotein_id,gene_name))%>%filter(gene_name=='Satb2')
stpallcoefdf%>%safe_left_join(ms_id2protein_id%>%distinct(uprotein_id,gene_name))%>%filter(gene_name=='Satb2')
ms_id2protein_id%>%filter(uprotein_id%in%best_uprotein_ids)%>%filter(gene_name=='Satb2')%>%as.data.frame

ebayes_stepwiseold<-ebayes_stepwise

}


################################################################################
########Test that all this has worked
################################################################################

test_that("mostly just a few funny cases where we have ribo but no MS",{
	# expect_true(all(has_low_ms))
	gene_name_in_counts <- allms%>%select(ms_id,msgname=gene_name)%>%left_join(ms_id2protein_id%>%distinct(ms_id,protein_id),by='ms_id')%>%group_by(msgname)%>%summarise(any_incounts = any(!is.na(protein_id)))%>%
		{setNames(.$any_incounts,.$msgname)}
	table(gene_name_in_counts)
	# allms%>%filter(gene_name %in% names(gene_name_in_counts[!gene_name_in_counts]))%>%.$gene_name%>%unique
	# 'Znf516'%in%cds$gene_name
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
