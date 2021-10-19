
#!udsr/bin/env Rscript
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
#
conflict_prefer('setdiff','BiocGenerics')
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
	countmat
}

message('...done' )


LOWCOUNTLIM <- 32

source(here('src/Rprofile.R'))




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

 #make gene names non redundant
ids_nrgname <- mcols(cds)%>%as.data.frame%>%distinct(gene_name,gene_id)%>%group_by(gene_name)%>%
	mutate(new_gene_name=paste0(gene_name[1],c('',paste0('_',2:n())))[1:n()])%>%
	ungroup%>%
	select(gene_name=new_gene_name,gene_id)%>%
	left_join(cds%>%mcols%>%as.data.frame%>%distinct(transcript_id,gene_id,protein_id))

stopifnot(ids_nrgname%>%group_by(gene_id)%>%filter(n_distinct(gene_name)>1)%>%nrow%>%identical(0L))


{
#add gene id
allsegcounts%<>%separate(sample,into=c('time','assay','rep'))%>%group_by(protein_id,assay,time)%>%
	mutate(hascounts = sum(total[assay=='ribo'])>LOWCOUNTLIM)
allsegcounts%<>%unite(sample,time,assay,rep)
allsegcounts_nz <- allsegcounts%>%group_by(protein_id)%>%filter(any(hascounts))
allsegcounts_nz %<>% safe_left_join(mcols(gtf_gr)[,c('gene_id','protein_id')]%>%
		as.data.frame%>%distinct(gene_id,protein_id))
cds_nz <- cds%>%subset(protein_id %in% unique(allsegcounts_nz$protein_id))

allTEchangedf<-read_tsv('tables/manuscript/go_all_highcount_updown.tsv')

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

protid_gid_df <- cds%>%mcols%>%as.data.frame%>%distinct(protein_id,gene_id)

# allsegcounts_nz%>%safe_left_join(allids%>%select(protein_id,gene_name))%>%.$gene_id%>%n_distinct
allsegcounts%>%safe_left_join(protid_gid_df)%>%.$gene_id%>%n_distinct

highcountgenes<-allsegcounts_nz$gene_id%>%unique

setequal(highcountgenes,allTEchangedf$gene_id)

}
#true but takes ages
# stopifnot(all(allTEchangedf$gene_id%in%allsegcounts%>%safe_left_join(allids%>%distinct(gene_id,protein_id))%>%.$gene_id))
# stopifnot(all(allTEchangedf$gene_id%in%allsegcounts_nz%>%safe_left_join(allids%>%distinct(gene_id,protein_id))%>%.$gene_id))
# inclusiontable(allTEchangedf$gene_id,highcountgenes)
# inclusiontable(allTEchangedf$gene_name,ms_id2protein_id$gene_name)


# gapdh_protein_ids <- cds%>%subset(gene_name=='Gapdh')%>%.$protein_id%>%unique
# gapdh_protein_ids %in% allsegcounts$protein_id
# gapdh_protein_ids %in% allsegcounts_nz$protein_id




################################################################################
########Load the matching between protein IDS and mass spec ids
################################################################################


{

getms<-function(msfile){
	allms=data.table::fread(msfile)
	#some formatting differences
	allms%<>%select(ms_id=Protein_IDs,everything())
	allms$time%<>%str_replace('p5','5')
	allms$dataset%<>%str_replace('p5','5')
	allms$dataset%<>%str_replace('_rep','_')
	allms$dataset%<>%str_replace('^[^_]+_','')#no need to annotate what signal type it is
	allms$dataset%<>%str_replace('total','MS')#no need to annotate what signal type it is
	allms
}
allms <- getms(msfile)
allmsibaq <- getms(here('pipeline/ms_tables/ms_iBAQ_total_ms_tall.tsv'))
#allms <- allmsibaq
}

{
{
# allms%>%group_by(dataset)%>%group_slice(10)%>%filter(!is.na(ms_id))%>%safe_left_join(allmsibaq%>%select(dataset,ms_id,ibaq=signal))%>%{cor(log10(.$signal),log10(.$ibaq),use='complete')}


if(file.exists('data/ms_id2protein_id.rds')){
	ms_id2protein_id <- readRDS('data/ms_id2protein_id.rds')			# ms_id2protein_id %>%saveRDS('data/ms_id2protein_id.rds')
}else{
	source(here('src/R/Load_data/get_ms_gencode_tablematch.R'));	
}

ms_id2protein_id%<>%{select(.,-matches('gene_name'))}%>%
	safe_left_join(ids_nrgname%>%distinct(transcript_id,gene_name),by='transcript_id')

#get the satb2 ids we want to look at
allids <- readRDS('pipeline/allids.txt')
ms_id2protein_id%<>%select(-matches('gene_name'))
trgeneids <- mcols(cds)%>%as.data.frame%>%distinct(transcript_id,gene_id)
ms_id2protein_id%<>%safe_left_join(trgeneids)
ms_id2protein_id %<>% safe_left_join(ids_nrgname%>%distinct(transcript_id,gene_name),by='transcript_id')

satb2ids <- ms_id2protein_id%>%filter(gene_name=='Satb2')%>%.$protein_id

}


#keep the mass spec we've matched to protein IDs
matched_ms <- allms%>%semi_join(ms_id2protein_id%>%distinct(ms_id,protein_id))
matched_ms%<>%group_by(ms_id)%>%filter(!all(is.na(signal)))
#get our matrix and design df
c(matchedms_mat,matched_ms_design)%<-% get_matrix_plus_design(matched_ms,ms_id)

ibaqMatl2<-allmsibaq%>%
	semi_join(ms_id2protein_id%>%distinct(ms_id,protein_id))%>%
	group_by(ms_id)%>%filter(!all(is.na(signal)))%>%
	get_matrix_plus_design(.,ms_id)

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
{
conflict_prefer("setequal", "dplyr")
stopifnot(exists('matchedms_mat'))


#I was going for some kind of weird thing where I could use a source script
#as a function with caching. I didn't work out. I apologize to my future
#self

c(posteriorsum,posteriors) %<-% with(new.env(),{
	if(!exists('posteriors')){
		if(file.exists('data/proDD_ibaq.data')){
			load(file='data/proDD_ibaq.data')
		}else{
			source(here('src/R/Load_data/proDD.R'))
		}
	}

	stopifnot(setequal(posteriorsum$ms_id, matched_ms$ms_id%>%unique))
	# matched_ms%>%group_by(ms_id,time)%>%summarise(mean(na.omit(signal)))
	list(posteriorsum,posteriors)
})
#
stopifnot(posteriorsum$ms_id%>%n_distinct%>%identical(6267L))

assert_that(posteriorsum%>%filter(ms_id %in% satb2ms_ids)%>%.$ms_id%>%unique%>%n_distinct%>%`>`(1))
assert_that(posteriorsum%>%filter(ms_id %in%'B7FAU9;Q8BTM8;B7FAV1')%>%.$precision%>%head(1)%>%`>`(100))
assert_that(all(rownames(matchedms_mat)%>%is_in(posteriorsum$ms_id)))
assert_that(all(rownames(matchedms_mat)%>%is_in(posteriorsum$ms_id)))
assert_that(posteriorsum%>%filter(ms_id %in%'B7FAU9;Q8BTM8;B7FAV1')%>%.$precision%>%head(1)%>%`>`(100))


# inclusiontable(posteriorsum$ms_id, matched_ms$ms_id%>%unique)
	# matched_
(setdiff(posteriorsum$ms_id, matched_ms$ms_id%>%unique))
(setdiff(matched_ms$ms_id%>%unique, posteriorsum$ms_id))

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


genes_w_ms <- posteriorsum%>%safe_left_join(allow_dup=T,ms_id2protein_id%>%select(ms_id,gene_name))%>%.$gene_name

test_that("all of our ms posteriors are in the ms_id2protein_id object",{
	genes_w_ms <- posteriorsum%>%safe_left_join(allow_dup=T,ms_id2protein_id%>%select(ms_id,gene_name))%>%.$gene_name

	oldtegenes <- allTEchangedf%>%safe_left_join(allow_dup=T,ids_nrgname)%>%left_join(ms_id2protein_id)
	oldtegenes <- allTEchangedf$gene_name
	inclusiontable(genes_w_ms,oldtegenes)
	#there are currently 9 genes with ms that weren't in the oldte genes
	setdiff(genes_w_ms,oldtegenes)
	#is this the gene name redundancy?
	table(ms_id2protein_id$gene_name %in% ids_nrgname$gene_name)
	#so the above checks out, ms_id2protein_id gnames are all in ids_nrgname
})



}

################################################################################
########Perform linear modeling of counts only
################################################################################
tps <- mainribosamps%>%str_extract(.,'[^_]+')%>%unique
SPLINE_N <- 4
# RIBOSIGCOL <- 'centercount'
RIBOSIGCOL <- sym('total')
allsegcounts_nz%>%colnames
satb2_uids <- ms_id2protein_id%>%filter(gene_name=='Satb2')%>%.$uprotein_id%>%unique
best_satb2_uid <- 'ENSMUSP00000110057_4528'
best_satb2_pid <- 'ENSMUSP00000110057'
satb2cds <- cds%>%subset(protein_id %in% satb2ids)
satb2cdsdj <- satb2cds%>%disjoin(with=TRUE)

{

highcountpids <- allsegcounts_nz%>%group_by(gene_id,protein_id)%>%summarise(medcount=median(!!RIBOSIGCOL))%>%filter(medcount==max(medcount))%>%distinct(gene_id,protein_id)%>%mutate(is_gid_highest=T)

#Library size factors - these will be done with DESEQ, and need to be done only one with each cgene.
libsizenormfacts <- allsegcounts_nz%>%inner_join(highcountpids)%>%select(protein_id,sample,total)%>%spread(sample,total)%>%{DESeq2::estimateSizeFactorsForMatrix(as.matrix(.[,-1]))}


####
ribocounts <- allsegcounts_nz%>%select(protein_id,dataset=sample,signal=!!RIBOSIGCOL)%>%filter(dataset%in%mainribosamps)%>%mutate(signal=replace_na(signal,0))
c(ribocountmat,ribodesign)%<-% get_matrix_plus_design(ribocounts,protein_id,transform=identity,sigcol=signal)

hasnoribo <- ribocountmat%>%rowSums%>%`==`(0)
noriboprotids <- names(hasnoribo)[hasnoribo]
noribobutms <- noriboprotids%>%intersect(ms_id2protein_id$protein_id)
ribocountmat <- ribocountmat[!hasnoribo,]
countvoom <- limma::voom(ribocountmat,design=model.matrix(~ns(as.numeric(time),2), ribodesign))
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



allTEchangedf<-read_tsv('tables/manuscript/go_all_highcount_updown.tsv')


all(cds$gene_id %in% ids_nrgname$gene_id)
ids_nrgname$gene_name%>%anyDuplicated



ids_nrgname_pid<-ids_nrgname%>%distinct(transcript_id,gene_name)%>%safe_left_join(allids%>%select(transcript_id,protein_id))

countvoomgenenames<-data.frame(protein_id=rownames(countvoom))%>%
	safe_left_join(ids_nrgname_pid)%>%.$gene_name

#so all the old te changedf guys are in this count voom
stopifnot(all(allTEchangedf$gene_name %in% countvoomgenenames))


}




################################################################################
########construct voom object 
################################################################################
	
{
	conflict_prefer("between", "dplyr")
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

allTEchangedf<-read_tsv('tables/manuscript/go_all_highcount_updown.tsv')
mscountvoomgenenames<-data.frame(uprotein_id=rownames(mscountvoom))%>%
	safe_left_join(ms_id2protein_id%>%select(uprotein_id,gene_name))%>%.$gene_name

stopifnot(all(mscountvoomgenenames%in%allTEchangedf$gene_name))

}

################################################################################
########Fix the count expression data
################################################################################
	

{
conflict_prefer("rowMedians", "matrixStats")
featuredata = rownames(allcountmat)%>%data.frame(protein_id=.)%>%
	safe_left_join(as.data.frame(mcols(cds))%>%distinct(protein_id,gene_id,transcript_id))%>%
	safe_left_join(ids_nrgname)%>%
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

featuredata$length = sum(width(cds%>%split(.,.$protein_id)))[featuredata$protein_id]

library(txtplot)

countexprdata <- ExpressionSet(
	allcountmat,
	AnnotatedDataFrame(allcountdesign%>%as.data.frame%>%set_rownames(colnames(allcountmat))),
	AnnotatedDataFrame(featuredata)
)

countexprdata%>%saveRDS(here('pipeline/exprdata/countexprset.rds'))

all(featuredata$protein_id%in%allsegcounts_nz$protein_id)

}



################################################################################
########Manipulate the effects from spline space to linear space 
################################################################################
# design(allcountebayes)
{
itime_modelmat<-model.matrix(~ 1 + ribo+MS+time+time:(ribo+MS),data=mscountvoomdesign%>%mutate(time=factor(tps[add(time,3)])))
itime_modelmat%<>%set_colnames(itime_modelmat%>%colnames%>%str_replace('riboTRUE','TE')%>%str_replace('MSTRUE','MS_dev'))
#
#for transforming back to time space from spline space
splinemat <- t(ns(1:5,SPLINE_N))
splinezeros <- matrix(0,ncol=5,nrow=SPLINE_N)
#
effect_is_time <- voomeffects%>%str_detect(negate=TRUE,'ns\\(')
stopifnot(length(rle(effect_is_time))==2)#all main effects, then time
nosplinen <- voomeffects%>%str_detect(negate=TRUE,'ns\\(')%>%which%>%tail(1)
maineffzeros <- matrix(0,ncol=5,nrow=nosplinen)
#
#contrast matrices transform back to time space
alltimeeff <- rbind(
	maineffzeros,
	t(ns(1:5,SPLINE_N)),
	splinezeros,
	splinezeros
)%>%set_rownames(voomeffects)%>%set_colnames(paste0(tps,'_all'))
#
timeTEeffect <- rbind(
	matrix(0,ncol=5,nrow=3),
	splinezeros,
	t(ns(-2:2,SPLINE_N)),
	splinezeros
)%>%set_rownames(voomeffects)%>%set_colnames(paste0(tps,'_TE'))
#
timeMSeffect <- rbind(
	matrix(0,ncol=5,nrow=3),
	splinezeros,
	splinezeros,
	t(ns(1:5,SPLINE_N))
)%>%set_rownames(voomeffects)%>%set_colnames(paste0(tps,'_MSdev'))
#
ntps<-length(tps)
sntps <- seq_along(tps)
i_n <- 1
#
itimecontrasts=alltimeeff
#
stepwise_contrasts <- lapply(list(all=alltimeeff,TE=timeTEeffect,MS_dev=timeMSeffect),function(itimecontrasts){
	lapply(1:(ntps-1),function(i_n){
		cname <- colnames(itimecontrasts)[i_n+1]%>%str_replace('_','_step_')
		(itimecontrasts[,i_n+1,drop=FALSE] - (itimecontrasts[,i_n,drop=FALSE])) %>%
		set_colnames(cname)
	})%>%do.call(cbind,.)
})%>%do.call(cbind,.)
contrastmatall <- diag(length(voomeffects))[,1:nosplinen]%>%set_colnames(voomeffects[1:nosplinen])%>%cbind(alltimeeff,timeTEeffect,timeMSeffect)
contrastmat <- contrastmatall%>% . [,colSums(.!=0)!=0]#remove the zero cols at first time point
#
counteffs <- colnames(countvoom$design)
countnosplinen <- counteffs%>%str_detect(negate=TRUE,"ns\\(")%>%which%>%tail(1)
contrastmatcounts <- diag(length(counteffs))[,1:countnosplinen]%>%set_colnames(counteffs[1:countnosplinen])%>%
	cbind(head(alltimeeff[-3,],-SPLINE_N),head(timeTEeffect[-3,],-SPLINE_N),head(timeMSeffect[-3,],-SPLINE_N))
contrastmatcounts%<>% . [,colSums(.!=0)!=0]#remove the zero cols at first time point
#
itimeeffs <- contrastmat%>%colnames
#
}



################################################################################
########model fitting, in various spaces
################################################################################	

{
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
stop()
mscountebayes -> ibaqmscountebayes
best_uprotein_ids <- mscountebayes%>%topTable(coef=12:15,number=Inf)%>%rownames_to_column('uprotein_id')%>%
	safe_left_join(ms_id2protein_id%>%distinct(gene_id,uprotein_id))%>%
	group_by(gene_id)%>%slice(which.max(P.Value))%>%
	.$uprotein_id
}
#27.360168
}
bestpvals_lfq <- lfqmscountebayes%>%topTable(coef=12:15,number=Inf)%>%rownames_to_column('uprotein_id')%>%
	safe_left_join(ms_id2protein_id%>%distinct(gene_id,uprotein_id))%>%
	group_by(gene_id)%>%slice(which.max(P.Value))%>%select(F,P.Value)
bestpvals_ibaq <- ibaqmscountebayes%>%topTable(coef=12:15,number=Inf)%>%rownames_to_column('uprotein_id')%>%
	safe_left_join(ms_id2protein_id%>%distinct(gene_id,uprotein_id))%>%
	group_by(gene_id)%>%slice(which.max(P.Value))%>%select(F,P.Value)




#now plot


plotfile<- here(paste0('plots/','lfq_ibaq_comp','.pdf'))
pdf(plotfile)
	bestpvals_lfq%>%safe_left_join(bestpvals_ibaq%>%select(F_ibaq=F,pval_ibaq=P.Value))%>%
	ggplot(.,aes(x= - log10(P.Value),y = - log10(pval_ibaq)))+
	geom_point()
	# scale_color_discrete(name='colorname',colorvals)+
	# scale_x_continuous(paste0('xname'))+
	# scale_y_continuous(paste0('yname'))+
	ggtitle(paste0('LFQ vs. IBAQ'))+
	theme_bw()
dev.off()
normalizePath(plotfile)

best_satb2_pid%in% best_uprotein_ids


stop()
#And redo the model with the best pairs
#Both for MS and count
bestmscountvoom <- mscountvoom[best_uprotein_ids,]
bestmscountebayes <- eBayes(lmFit(mscountvoom[best_uprotein_ids,]))

#And for Counts only
best_protein_ids = best_uprotein_ids%>%str_replace('_\\d+$','')
bestonlycountebayes <- eBayes(lmFit(countvoom[best_protein_ids,]))


ms_id2protein_id%>%mutate(uprotein_id%in%best_uprotein_ids)

satb2ids %in% best_protein_ids
satb2ids %>%intersect(best_protein_ids)


ebayes_stepwise <- contrasts.fit(bestmscountebayes,stepwise_contrasts)



stpallcoefdf <- ebayes_stepwise%>%eBayes%>%topTable(number=Inf)%>%rownames_to_column('uprotein_id')%>%select(-AveExpr,-F,-P.Value,-adj.P.Val)


# totalsatb2steps<-stpallcoefdf%>%safe_left_join(ms_id2protein_id%>%distinct(uprotein_id,gene_name))%>%filter(gene_name=='Satb2')
stpallcoefdf%>%safe_left_join(ms_id2protein_id%>%distinct(uprotein_id,gene_name))%>%filter(gene_name=='Satb2')
ms_id2protein_id%>%filter(uprotein_id%in%best_uprotein_ids)%>%filter(gene_name=='Satb2')%>%as.data.frame

ebayes_stepwiseold<-ebayes_stepwise

# }

###############################################################################
#######Breakdown into confidence intervals over assays
###############################################################################

###
{
	#Extract confidence inttervals for effects
	time_eff_contrasts <- contrastmatall
	effects<-colnames(time_eff_contrasts)
	istimete<-colnames(time_eff_contrasts)%>%str_detect('_TE')
	time_eff_contrasts[,istimete] %<>%add( time_eff_contrasts[,'TE'])
	istimeMSde<-colnames(time_eff_contrasts)%>%str_detect('MSdev.')
	time_eff_contrasts[,istimeMSde] %<>% add(time_eff_contrasts[,'TE'])
	time_eff_contrasts[,istimeMSde] %<>% add(time_eff_contrasts[,'MS_dev'])
	time_eff_contrasts%<>%.[,istimete | istimeMSde]
	timeeffnames<-colnames(time_eff_contrasts)%>%setNames(.,.)
	effect = timeeffnames[1]
	time_eff_contrastsdf<-	lapply(timeeffnames,function(effect){
		message(effect)
		topTable(eBayes(contrasts.fit(bestmscountebayes,time_eff_contrasts[,effect])),coef=1,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('uprotein_id')
	})%>%bind_rows(.id='effect')
	#
	time_eff_contrastsdf%<>%separate(effect,c('assay','time'))
}
#
#
time_eff_contrastsdf%<>%select(assay, time, uprotein_id, logFC, CI.L, CI.R, AveExpr, t, P.Value, adj.P.Val, B)
#
#
get_contrast_cis <- function(ebayesob,contrastmat,rownamecol='uprotein_id',contrcol='datagroup'){
	datagroup_names <- colnames(contrastmat)
	#
	prediction_df<-	lapply(datagroup_names%>%setNames(.,.),function(datagroup){
		message(datagroup)
		# prediction_ob$coef
		topTable(eBayes(contrasts.fit(ebayesob,contrastmat[,datagroup,drop=F])),coef=1,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column(rownamecol)
	})%>%bind_rows(.id=contrcol)
	#
	#
	prediction_df%>%as_tibble%>%
	separate(datagroup,c('time','assay'))
}
#
pred_te_design <- (function(){
	#Extract confidence intervals for time points, including TE and MS
	datagroup_names <- mscountvoomdesign%>%.$dataset%>%str_replace('_\\d+$','')%>%unique%>%setNames(.,.)
	pred_te_design<-bestmscountebayes$design%>%unique%>%set_rownames(datagroup_names)%>%t
	#add in TE to these
	pred_te_design%<>%cbind(timeTEeffect%>%{.['TE',]<-1;.})
	# sample_contrasts%<>%cbind(.,timeMSeffectnm)
	pred_te_design
})()
#
#get the predictions - confidence intervals over actual datafor counts alone, and the whole shebang
countonly_pred_te_design<-
	pred_te_design%>%
	.[,colnames(.)%>%str_subset('total|ribo')]%>%
	.[rownames(.)%>%str_subset('MS_dev',negate=T),]%>%
	qs('M>5')
#
#
prediction_df_countonly <- get_contrast_cis(bestonlycountebayes,countonly_pred_te_design) %>% qs('D11')
prediction_df <- get_contrast_cis(bestmscountebayes,pred_te_design)
assert_that(all(prediction_df$assay%in%c('ribo','total','MS','TE','MSdev')))
#get the linear effects, differences
effcontrasts<-cbind(alltimeeff,timeTEeffect,
	timeMSeffect%>%set_colnames(colnames(.)%>%str_replace('_MSdev(\\w+)','\\1_MSdev'))
)

timeff_ciddf <- get_contrast_cis(bestmscountebayes,effcontrasts)

upids_w_sig_ms_dev <- timeff_ciddf%>%filter(adj.P.Val<0.05,assay%>%str_detect('MSdev'))%>%.$uprotein_id

stopifnot(timeff_ciddf$assay%>%unique%>%identical(c('all','TE','MSdev')))

###

# ################################################################################
# ########Test that all this has worked
# ################################################################################

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

#
allTEchangedf<-read_tsv('tables/manuscript/go_all_highcount_updown.tsv')
techange_prot_ids <- ms_id2protein_id%>%inner_join(allTEchangedf%>%filter(up|down))%>%.$protein_id
prediction_df%>%write_tsv('tables/prediction_df.tsv')

#Now, export linear CIs - a mix of the linear ones for our count data, and the
#proDA ones
export_CIs <- prediction_df_countonly%>%
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

#all genes with linear CIs are in the set tested for TE
stopifnot(all(export_CIs$gene_name%in%allTEchangedf$gene_name))

# allgnameshavealldsets <- all(export_CIs%>%group_by(gene_name)%>%tally%>%.$n%>%`==`(length(datagroup_names)))
# expect_true(allgnameshavealldsets)
#
#get the widths for our protein ids
pid_widths<-cds%>%split(.,.$protein_id)%>%width%>%sum%>%stack%>%as.data.frame
pid_widths%<>%set_colnames(c('width','protein_id'))
#collate some meta info - sig protein change, width

# metainfo <- export_CIs%>%distinct(uprotein_id,gene_name)%>%
# 	mutate(protein_id = uprotein_id%>%str_replace('_\\d+$',''))%>%
# 	mutate(protein_id %in% techange_prot_ids)%>%
# 	mutate(sig_MS_change = uprotein_id %in%genes_w_sig_ms_dev$uprotein_id)%>%
# 	safe_left_join(ms_id2protein_id%>%distinct(protein_id,gene_id,transcript_id))%>%
# 	left_join(pid_widths)

#there are unique tr and uprotein id pairs
stopifnot(n_distinct(ms_id2protein_id$uprotein_id)==n_distinct(ms_id2protein_id%>%{paste(.$uprotein_id,.$transcript_id)}))

metainfo$uprotein_id %in% best_uprotein_ids

allTEchangedf <- read_tsv('tables/manuscript/go_all_highcount_updown.tsv')
table(best_protein_ids%in%metainfo$protein_id)

#so meta info should have all cds ids, it should have logical columns annotating whether they appear in the
#segcounts, in the highcounts, in the ms, in the matched ms, counts, and in the chosen pairs of ms_counts
#whether they have sig TE change, sig MS change, 
#we also want info on dropouts - any droputs, whole missing stage

ms_id2protein_id$gene_name%>%setdiff(allids$gene_name)
ms_id2protein_id$gene_name%>%setdiff(ids_nrgname$gene_name)
allids$gene_name%>%setdiff(ids_nrgname$gene_name)
setdiff(ids_nrgname$gene_name,allids$gene_name)

n_stagemissingdf <- stagemissingdf%>%summarise(n_stagemissing=sum(missing))

dropoutuprotids <- matchedms_mat%>%apply(1,.%>%is.na%>%sum)%>%`!=`(0)%>%.[.]%>%names
stagemissingdf <- matchedms_mat%>%apply(1,.%>%is.na)%>%t%>%as.data.frame%>%rownames_to_column('ms_id')%>%gather(sample,isna,-ms_id)%>%
	separate(sample,c('time','assay','rep'))%>%group_by(ms_id,time)%>%summarise(missing=all(isna))


metainfo <- ids_nrgname%>%
	mutate(hascount = protein_id%in%allsegcounts$protein_id)%>%
	mutate(highcount = protein_id%in%allsegcounts_nz$protein_id)%>%	
	left_join(fData(countexprdata)%>%select(protein_id,is_gid_highest,width=length))%>%
	mutate(dTE = is_gid_highest & (gene_id%in%(allTEchangedf%>%filter(up|down)%>%.$gene_id)))%>%
	left_join(ms_id2protein_id%>%distinct(protein_id,ms_id,uprotein_id))%>%
	mutate(has_ms =!is.na(uprotein_id))%>%
	mutate(isbest = uprotein_id %in% best_uprotein_ids)%>%
	mutate(sig_MS_change = uprotein_id %in%upids_w_sig_ms_dev)%>%
	mutate(anydropout = ms_id%in%dropoutuprotids)%>%
	safe_left_join(allow_missing=TRUE,n_stagemissingdf%>%select(ms_id,n_stagemissing))

#each uprotein id is a unique transcript, exclude those lines that are just transcriopts without a uprotein id.
stopifnot(n_distinct(metainfo%>%filter(!is.na(uprotein_id))%>%.$uprotein_id)==n_distinct(metainfo%>%filter(!is.na(uprotein_id))%>%select(uprotein_id,transcript_id)))

metainfo%>%filter(!is.na(uprotein_id))%>%group_by(uprotein_id)%>%filter(n_distinct(transcript_id)>1)

#export for modeling
export_CIs%>%write_tsv(here('pipeline/exprdata/limma_proDD_CIs.tsv'))

metainfo%>%write_tsv(here('pipeline/exprdata/limma_genemetadata.tsv'))

#
exprdf <- bestmscountvoom$E%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,signal,-uprotein_id)%>%as_tibble
#
exprdf%<>%separate(dataset,c('time','assay','rep'))%>%
  left_join(ms_id2protein_id%>%
    distinct(uprotein_id,protein_id)
  )%>%
  mutate(rep = as.numeric(rep))
exprdf%>%filter(assay=='MS')


# countsamples<-allsegcounts_nz$sample%>%unique%>%str_subset('total|ribo')%>%str_subset('ribo')



# save.image('data/integrate_exprdata2.Rdata')
 load('data/integrate_exprdata2.Rdata')

