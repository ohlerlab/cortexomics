#!/usr/bin/env Rscript
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
library(zeallot)
library(splines)

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

#add gene id
allsegcounts %<>% safe_left_join(mcols(gtf_gr)[,c('gene_id','protein_id')]%>%as.data.frame%>%distinct(gene_id,protein_id))
allsegcounts_nz <- allsegcounts%>%group_by(protein_id)%>%filter(any(centercount>10))

unique(allsegcounts$protein_id)

#get exons
if(!exists('gtf_gr')) gtf_gr<-rtracklayer::import(con=gtf,format='gtf')
exons <- gtf_gr%>%subset(type=='exon')
cds <- gtf_gr%>%subset(type=='CDS')
  

cds_nz <- cds%>%subset(protein_id %in% unique(allsegcounts_nz$protein_id))

mainribosamps <- '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/sample_parameter.csv'%>%
	fread%>%
	filter(is.na(fraction))%>%
	filter(assay=='ribo')%>%
	filter(sample_id%>%
	str_detect(neg=T,'test'))%>%
	.$sample_id
stopifnot('E13_ribo_1' %in% mainribosamps)
mainribosamps %in% allsegcounts$sample

allsegcounts%>%filter(sample==	 'E13_ribo_2')%>%filter(protein_id=='ENSMUSP00000000001')
allsegcounts%>%filter(sample==	 'E13_ribo_1')%>%filter(protein_id=='ENSMUSP00000000001')



################################################################################
########Load the matching between protein IDS and mass spec ids
################################################################################

allms=data.table::fread(msfile)
#some formatting differences
allms%<>%select(ms_id=Protein_IDs,everything())
allms$dataset%<>%str_replace('p5','5')
allms$dataset%<>%str_replace('_rep','_')
allms$dataset%<>%str_replace('^[^_]+_','')#no need to annotate what signal type it is
allms$dataset%<>%str_replace('total','MS')#no need to annotate what signal type it is
	
ms_id2protein_id <- with(new.env(),{
	source(here('src/R/Load_data/get_ms_gencode_tablematch.R'));
	ms_id2protein_id
})
#keep the mass spec we've matched to protein IDs
matched_ms <- allms%>%semi_join(ms_id2protein_id%>%distinct(ms_id,protein_id))
matched_ms%<>%group_by(ms_id)%>%filter(!all(is.na(signal)))
ms_id2protein_id$ms_id%>%duplicated
ms_id2protein_id%<>%ungroup%>%filter(!duplicated(protein_id))
stopifnot(!any(ms_id2protein_id$protein_id%>%duplicated))


################################################################################
########Perform linear modeling to select the best protein IDs
################################################################################

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

get_limmafit_predvals <- function(limmafit,designmatrix){
  (limmafit$coef %*% t(limmafit$design))%>%
	  set_colnames(designmatrix$dataset)%>%
	  as.data.frame%>%
	  rownames_to_column('gene_name')%>%
	  gather(dataset,signal,-gene_name)%>%
	  left_join(designmatrix)%>%
	  distinct(gene_name,time,assay,.keep_all = TRUE)%>%
	  select(-rep)%>%
	  mutate(dataset=paste0(as.character(time),'_',assay))
}

library(limma)
stop()
#Get the spline smoothed estimates of abundance per timepoint for each ms id
c(matchedms_mat,matched_ms_design)%<-% get_matrix_plus_design(matched_ms,ms_id)
limma_pred <- get_limmafit_predvals(
	limma::lmFit(matchedms_mat,design=model.matrix(~ns(as.numeric(time),3), matched_ms_design)),
	matched_ms_design)

mspredmat <- limma_pred%>%get_matrix_plus_design(gene_name)%>%.[[1]]

#countdata<-allsegcounts%>%get_matrix_plus_design(countdata,protein_id)

#Get the spline smoothed estimates of abundance per timepoint for each gencode protein id
message('using the old data')


###Code for using the old data
# rawcountmat <- '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/feature_counts/all_feature_counts'%>%
# 	fread%>%
# 	select(feature_id,one_of(mainribosamps))%>%
# 	left_join(gtf_gr%>%mcols%>%as.data.frame%>%distinct(feature_id=gene_id,protein_id)%>%filter(!is.na(protein_id))%>%group_by(feature_id)%>%slice(1:n()))%>%
# 	ungroup%>%select(-feature_id,protein_id,one_of(mainribosamps))
# countdata<-rawcountmat%>%gather(dataset,signal,-protein_id)

countddata <- allsegcounts%>%select(protein_id,dataset=sample,signal=centercount)

c(countmat,countdesign)%<-% get_matrix_plus_design(countdata,protein_id,transform=identity,sigcol=signal)

limma_pred <- get_limmafit_predvals(
	limma::lmFit(limma::voom(countmat,design=model.matrix(~ns(as.numeric(time),3), countdesign))),
	countdesign)

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

stop()

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




# #and then transform the counts
# countsmatrix_snorm<-cbind(
#   DESeq2::vst(countsmatrix[,ribcols]),
#   DESeq2::vst(countsmatrix[,totcols])
# )





# top_protein_ids <- matched_ms%>%
# 	distinct(ms_id,protein_id)%>%
# 	left_join(allsegcounts)%>%
# 	group_by(ms_id,protein_id)%>%
# 	summarise(mcount = mean(total))%>%
# 	slice(which.max(mcount))

# top_protein_ids%>%
# 	select(-mcount)


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
	'Satb2' %in% exprdata$gene_name
	satb2ddata <- (exprddata%>%filter(gene_name=='Satb2'))
	'Q8VI24' %in% satb2ddata$ms_id
	
	#the total expression levels should still agree with previous findidngs
	#we'll need a way to compare quantification moddes in our different genes
	te_table%>%filter(time=='P0')%$%{log2fc > 0 & padj < 0.05}
})





test_that("Mappability masks used",{})

test_that("Xtail's TE should broadly agree with limma",{})
test_that("Xtail's TE should broadly agree with DESeq",{})

test_that("Spline based TE change implemented",{})

#There should in general be a decent correaltion between the spectral coef and the MS

#' But - shoul we be using correlation smoothed, or just normal corelation???	