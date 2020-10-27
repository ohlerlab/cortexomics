
################################################################################
########This (hopefully final) version of the script uses deepshape prime
########to quantify Riboseq data
################################################################################
base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
}

#mcols(cds)[,c('transcript_id','gene_id')]%>%as.data.frame%>%
if(!exists('tx_countdata'))load('data/1_integrate_countdata.R')
gid2gnm = ids_nrgname%>%distinct(gene_id,gene_name)%>%{safe_hashmap(.[[2]],.[[1]])}

# msid2pid = metainfo%>%
# 	distinct(ms_id,protein_id)%>%
# 	{safe_hashmap(.[['ms_id']],.[['protein_id']])}
# msid2upid = metainfo%>%
# 	distinct(ms_id,uprotein_id)%>%
# 	{safe_hashmap(.[['ms_id']],.[['uprotein_id']])}
# upid2msid = metainfo%>%
# 	distinct(ms_id,uprotein_id)%>%
# 	{safe_hashmap(.[['uprotein_id']],.[['ms_id']])}



# # base::source(here::here('src/R/Rprofile.R'))
# # {

# # getms<-function(msfile){
# # 	allms=data.table::fread(msfile)
# # 	#some formatting differences
# # 	allms%<>%select(ms_id=Protein_IDs,everything())
# # 	allms$time%<>%str_replace('p5','5')
# # 	allms$dataset%<>%str_replace('p5','5')
# # 	allms$dataset%<>%str_replace('_rep','_')
# # 	allms$dataset%<>%str_replace('^[^_]+_','')#no need to annotate what signal type it is
# # 	allms$dataset%<>%str_replace('total','MS')#no need to annotate what signal type it is
# # 	allms
# # }
# # # allms <- getms(msfile)
# # allms <- getms(here('pipeline/ms_tables/ms_iBAQ_total_ms_tall.tsv'))
# # #allms <- allmsibaq
# # }


# # {
# # if(file.exists('data/ms_id2protein_id.rds')){
# # 	ms_id2protein_id <- readRDS('data/ms_id2protein_id.rds')			# ms_id2protein_id %>%saveRDS('data/ms_id2protein_id.rds')
# # }else{
# # 	source(here('src/R/Load_data/get_ms_gencode_tablematch.R'));	
# # }

# # ms_id2protein_id%<>%{select(.,-matches('gene_name'))}%>%
# # 	filter(transcript_id%in%alltrs)%>%
# # 	safe_left_join(ids_nrgname%>%distinct(transcript_id,gene_name),by='transcript_id')

# #get the satb2 ids we want to look at
# allids <- readRDS('pipeline/allids.txt')
# ms_id2protein_id%<>%select(-matches('gene_name'))
# trgeneids <- mcols(cds)%>%as.data.frame%>%distinct(transcript_id,gene_id)
# ms_id2protein_id%<>%safe_left_join(trgeneids)
# ms_id2protein_id %<>% safe_left_join(ids_nrgname%>%distinct(transcript_id,gene_name),by='transcript_id')

# satb2ids <- ms_id2protein_id%>%filter(gene_name=='Satb2')%>%.$protein_id

# }

# #get the ms which match one of our genes
# matched_ms <- allms%>%inner_join(ms_id2protein_id%>%distinct(ms_id,gene_id))
# #get the ms that isn't missing
# matched_ms%<>%group_by(ms_id)%>%filter(!all(is.na(signal)))
# #now get the ms_id which is the best for each gene
# #lowest number of missing, then highest median signal
# matched_ms = matched_ms%>%
# 	group_by(ms_id)%>%
# 	mutate(medsig=median(signal))%>%
# 	mutate(nmissing=sum(is.na(signal)))%>%
# 	group_by(gene_id)%>%
# 	arrange(nmissing,desc(medsig))%>%
# 	filter(ms_id==ms_id[1])%>%
# 	select(-medsig,-nmissing)
# #now name the mass spec according to our scheme
# matched_ms$gene_name = gid2gnm[[matched_ms$gene_id]]
# ms_mat = matched_ms%>%select(gene_id,signal,dataset)%>%pivot_wider(id_cols='gene_id',dataset,values_from=signal)%>%
# 	{set_rownames(as.matrix(.[,-1]),.[[1]])}
# stopifnot(all(rownames(ms_mat)%in%rownames(tx_countdata$abundance)))
# #get the ms ids we want
# goodmsids = matched_ms$ms_id%>%unique
# #get the transcripts for these
# ms_transcripts = ms_id2protein_id%>%filter(ms_id%in%goodmsids)%>%.$transcript_id
base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
  base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
}
# stop()



#
protgnmtrids<-readRDS('data/protgnmtrids.rds')

proteinmsdata<-readRDS('data/proteinmsdata.rds')
proteinmsdata$ms_id = proteinmsdata[[1]]
proteinmsdatamatch <- proteinmsdata
proteinmsdatamatch$cdsidx = match(proteinmsdatamatch$gene_name,fmcols(cdsgrl,gene_name))
proteinmsdatamatch <- proteinmsdatamatch%>%filter(!is.na(cdsidx))
ms_transcripts <- protgnmtrids$tr_id
#select iBAQ and input cols
sigcols <- proteinmsdatamatch%>%colnames%>%str_subset('iBAQ')%>%str_subset('input')
proteinmsdatamatch%<>%select(ms_id,gene_name,g_id,one_of(sigcols))
ms_mat <- proteinmsdatamatch[,-c(1:3)]%>%as.matrix
#
rownames(ms_mat) <- proteinmsdatamatch$ms_id
colnames(ms_mat) <- colnames(ms_mat)%>%
	str_replace('iBAQ\\.','')%>%
	str_replace('input','MS')%>%
	str_replace('rep','')%>%
	str_replace('p(\\d)','\\1')
ms_mat = proDA::median_normalization(ms_mat)
	
################################################################################
########Import it all to get tr length scaled counts
################################################################################
{	
library(tximport)
library(tidyverse)

rnasalmonfiles = Sys.glob(here('pipeline/salmon/data/*total*/quant.sf'))
dpoutfiles = Sys.glob(here('pipeline/deepshapeprime/fakesalmonfiles/*ribo*/*'))

allquantfiles = c(rnasalmonfiles,dpoutfiles)
names(allquantfiles) <- allquantfiles%>%dirname%>%basename

dptrs = dpoutfiles[[1]]%>%fread%>%.$Name
salmontrs = allquantfiles[[1]]%>%fread%>%.$Name%>%str_extract('ENSMUST\\w+')
# inclusiontable(dptrs,salmontrs)
trs = intersect(dptrs,salmontrs)


####KEY STEP
ms_trs = intersect(trs,ms_transcripts)

ms_tx_countdata = tximport(files=allquantfiles,
	ignoreTxVersion=TRUE,
	tx2gene=mcols(cds)[,c('transcript_id','gene_id')]%>%as.data.frame%>%distinct(transcript_id,gene_id),
	type='salmon',
	countsFromAbundance='scaledTPM',
	importer=function(file){
		read_tsv(file,col_types=cols())%>%
			mutate(Name=str_extract(Name,'ENSMUST\\w+'))%>%
			filter(Name%in%ms_trs)%>%arrange(match(Name,trs))
})

randomround = function(x)floor(x)+rbinom(length(x),1,x%%1)
ms_tx_countdata$counts%<>%as.data.frame%>%
	mutate_all(randomround)%>%
	as.matrix%>%
	set_rownames(rownames(ms_tx_countdata$counts))

ms_tx_countdata$counts%>%as.data.frame%>%
	rownames_to_column('gene_id')%>%
	mutate(gene_name = gid2gnm[[gene_id]])%>%
	select(-gene_id)%>%select(gene_name,everything())%>%
	mutate_at(vars(-gene_name),list(randomround))%>%
	write_tsv('data/ms_tx_scaled_countData.tsv')
}

{
sharedgenes = intersect(ms_tx_countdata$abundance%>%rownames,tx_countdata$abundance%>%rownames)

quicktest(
	log2(ms_tx_countdata$abundance[sharedgenes,1]),
	log2(tx_countdata$abundance[sharedgenes,1])
)
ms_matgidrows <- ms_mat%>%set_rownames(proteinmsdatamatch$g_id)
mssharedgenes <- intersect(rownames(ms_matgidrows),
	intersect(
		ms_tx_countdata$abundance[,]%>%rownames,
		tx_countdata$abundance[highcountgenes,]%>%rownames)
		)


quicktest(
	log2(ms_matgidrows[mssharedgenes,1]),
	log2(tx_countdata$abundance[mssharedgenes,'E13_total_1'])
)

quicktest(
	log2(ms_matgidrows[mssharedgenes,1]),
	log2(tx_countdata$abundance[mssharedgenes,'E13_ribo_1'])
)


quicktest(
	log2(ms_matgidrows[mssharedgenes,1]),
	log2(ms_tx_countdata$abundance[mssharedgenes,'E13_total_1'])
)

quicktest(
	log2(ms_matgidrows[mssharedgenes,1]),
	log2(ms_tx_countdata$abundance[mssharedgenes,'E13_ribo_1'])
)

#also test salmon tpms for ribo
salmonribo = tximport(files= Sys.glob(here('pipeline/salmon/data/*ribo*/quant.sf')),
	ignoreTxVersion=TRUE,
	tx2gene=mcols(cds)[,c('transcript_id','gene_id')]%>%as.data.frame%>%distinct(transcript_id,gene_id),
	type='salmon',
	countsFromAbundance='scaledTPM',
	importer=function(file){
		read_tsv(file,col_types=cols())%>%
			mutate(Name=str_extract(Name,'ENSMUST\\w+'))%>%
			filter(Name%in%ms_trs)%>%arrange(match(Name,trs))
})
salmonribogenes = salmonribo$abundance%>%rownames
srshgenes = intersect(mssharedgenes,salmonribogenes)

quicktest(
	(ms_tx_countdata$abundance[srshgenes,'E13_ribo_1']),
	(salmonribo$abundance[srshgenes,1])
)

quicktest(
	(ms_tx_countdata$counts[srshgenes,'E13_ribo_1']),
	(salmonribo$counts[srshgenes,1])
)
ms_tx_countdata$counts[srshgenes,'E13_ribo_1']%>%sum%>%divide_by(1e6)
salmonribo$counts[srshgenes,1]%>%sum%>%divide_by(1e6)
}






#get our matrix and design df
sample_info_df <- data.frame(name = colnames(ms_mat),
                             stringsAsFactors = FALSE)
sample_info_df$time <- sample_info_df$name%>%str_extract('[^_]+')
sample_info_df$assay <- sample_info_df$name%>%str_extract('(?<=_)[^_]+')

proDAfitms <- proDA::proDA(
			ms_mat, design = ~ time, 
            col_data = sample_info_df%>%filter(assay=='MS'),
            # reference_level = "E13"
            reference_level = NULL
)


#get predictions
predmat = proDA::predict(proDAfitms)
preddf = predmat%>%as.data.frame%>%
	rownames_to_column('ms_id')%>%
	gather(dataset,estimate,-ms_id)%>%
	separate(dataset,c('time','assay','rep'))

timevals = paste0('time',sample_info_df$time%>%unique%>%tail(-1))%>%setNames(.,.)

#now get the prediction ses
times = sample_info_df%>%filter(assay=='MS')%>%.$time%>%unique
# proda_pref_ses =imap(coefficient_variance_matrices,function(vcov,nm){
#by analogy to https://rdrr.io/bioc/limma/src/R/toptable.R
#(search for margin.error)
# proda_pref_ses =map(1:nrow(proDAfitms),function(i){
# 	vcov = proDAfitms$coefficient_variance_matrices[[i]]
# 	sqrt(proDAfitms$feature_parameters$s2[i]) *
# 	(diag(X%*%vcov%*%t(X))) *
# 	qt(alpha,proDAfitms$feature_parameters$df[i])
# })a
# proda_pref_ses = proda_pref_ses%>%map_df(.id='ms_id',~enframe(setNames(.,times),'time','se'))
#now create confidence intervals
preddf = preddf%>%safe_left_join(proda_pref_ses,by=c('ms_id','time'))%>%
	mutate(CI.L = estimate - (1.96*se))%>%
	mutate(CI.R = estimate + (1.96*se))

contrasts = lapply(1:5,function(.) as.numeric(1:5%in%c(1,.)))
contrasts%<>%setNames(times)
preddf<- map_df(.id='time',contrasts, ~proDA::test_diff(proDAfitms,contrast=.))
preddf%<>%rename('ms_id'='name')

preddf <- preddf %>% 	
	mutate(estimate=diff)%>%
	mutate(CI.L = diff - (se*1.96))%>%
	mutate(CI.R = diff + (se*1.96))
#also add the se se to the E13 ones, they don't have it for some reason
preddf%<>%mutate(se = ifelse(is.na(se),(CI.R-CI.L)/3.92,))

# prodalfcs <- prodalfcs%>%mutate(time=time%>%str_replace('time',''))%>%
# 	mutate(CI.L = diff - (1.96*se))%>%a
# 	mutate(CI.R = diff + (1.96*se))

# n_obs = sqrt(proDAfitms@elementMetadata@listData$n_approx)
# intstdev = proDAfitms$coefficient_variance_matrices%>%map_dbl(1)
# intse = intstdev / n_obs

# prodalfcs %<>% bind_rows(tibble(name = rownames(proDAfitms),time='E13',diff = 0)%>%
# 	mutate(CI.L = diff - (intse*1.96))%>%
# 	mutate(CI.R = diff + (intse*1.96)),.)


# preddf%>%group_by(gene_id)%>%group_slice(1)%>%{
# 	glm = lm(data=.,signal~time)
# 	predict(newdata=distinct(.,time),glm,interval='prediction')
# }


# X = model.matrix(~time ,sample_info_df%>%filter(assay=='MS'))%>%unique


# coefficient_variance_matrices = proDAfitms$coefficient_variance_matrices
# coefficient_variance_matrices%<>%setNames(rownames(proDAfitms))

# if(confint) {
# 		if(is.numeric(confint)) alpha <- (1+confint[1])/2 else alpha <- 0.975
# 		margin.error <- sqrt(eb$s2.post[top])*fit$stdev.unscaled[top,coef]*qt(alpha,df=eb$df.total[top])
# 		tab$CI.L <- M[top]-margin.error
# 		tab$CI.R <- M[top]+margin.error
# 	}



# prodalfcs %<>% rename('gene_name':='name')
# prodalfcs %<>% rename('gene_id':='gene_name')
# prodalfcs %<>% rename('gene_id':='gene_name')
#change of base formula
# prodalfcs%<>%mutate_at(vars(diff,CI.L,CI.R,se),list(~ . / log(2)))

# #make combined matrix
# proDAdatmat <- cbind(allcountmat[bestprotids,]%>%add(.1),ms_mat)%>%set_rownames(rownames(ms_mat))
proda_ms_pred <- proDAfitms$coefficients
proda_ms_pred[,2:5]=(proda_ms_pred[,2:5] + proda_ms_pred[,1])
colnames(proda_ms_pred)[1]='timeE13'
colnames(proda_ms_pred)%<>%str_replace('time(.*)','\\1_MS')

# re%>%slice(1)
# tgene = prodalfcs%>%slice(1)%>%.$gene_id%>%.[1]


ms_metadf <- tibble(
	gene_name = proteinmsdatamatch$gene_name,
	gene_id = proteinmsdatamatch$g_id,
	ms_id = proteinmsdatamatch$ms_id,
	n_missing = ms_mat%>%apply(2,is.na)%>%rowSums,
	median = ms_mat%>%matrixStats::rowMedians(.)
)

#add number of missing stages
ms_metadf <- ms_metadf%>%safe_left_join(by='gene_id',
	sel_ms_mat[,]%>%
    as.data.frame%>%
    rownames_to_column('gene_id')%>%
    # filter(!is.na(gene_id))%>%
    gather(dataset,signal,-gene_id)%>%
    separate(dataset,into=c('time','assay','replicate'))%>%
    group_by(gene_id,time)%>%summarise(missing=all(is.na(signal)))%>%
    summarise(n_stagemissing=sum(missing))
)

################################################################################
########Identify ribosomal proteins
################################################################################
  
#'get info on the ribosomal subu,nts from mats table
rids <- read_tsv(here('ext_data/riboprotids.tsv'))
lridssplit<-rids%>%filter(!`Mito-RP`,`RPL_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist
sridssplit<-rids%>%filter(!`Mito-RP`,`RPS_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist
ridssplit<-c(lridssplit,sridssplit)
# rids%>%filter(`Mito-RP`)%>%nrow

#get info on proteins so we can catagorize them

ms_metadf$is_rpl = sep_element_in(ms_metadf$ms_id,rids%>%filter(!`Mito-RP`,`RPL_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist)
ms_metadf$is_rps = sep_element_in(ms_metadf$ms_id,rids%>%filter(!`Mito-RP`,`RPS_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist)
ms_metadf$is_mito_rp = sep_element_in(ms_metadf$ms_id,rids%>%filter(`Mito-RP`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist)
# ms_metadf$is_mito_rps = sep_element_in(ms_metadf$ms_id,rids%>%filter(`Mito-RP`,`RPS_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist)

ms_metadf$is_rpl%>%table

preddf%<>%safe_left_join(distinct(ms_metadf,gene_name,gene_id,ms_id))

sel_ms_ids = ms_metadf%>%group_by(gene_name)%>%
	arrange(n_missing,desc(median))%>%slice(1)%>%.$ms_id


sel_ms_mat <- ms_mat
sel_ms_mat = sel_ms_mat[sel_ms_ids,]
rownames(sel_ms_mat) = proteinmsdatamatch$g_id[match(rownames(sel_ms_mat),proteinmsdatamatch$ms_id)]
sel_prodpreds <- preddf%>%filter(gene_id %in% rownames(sel_ms_mat))

proDAfitms%>%saveRDS('data/proDAfitms.rds')
sel_prodpreds%>%saveRDS('data/sel_prodpreds.rds')
sel_ms_mat%>%saveRDS('data/sel_ms_mat.rds')
ms_metadf%>%write_tsv('tables/ms_metadf.tsv')

if(F){
	testgname='Magohb'
	testgene = 'ENSMUSG00000030188'
	testms = proteinmsdatamatch$ms_id[match(testgene,proteinmsdatamatch$g_id)]
# > sel_ms_mat[testgene,]

}

################################################################################
########Do fold changes at the 
################################################################################
