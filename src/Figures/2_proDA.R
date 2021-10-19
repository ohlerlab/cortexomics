
################################################################################
########This (hopefully final) version of the script uses deepshape prime
########to quantify Riboseq data
################################################################################
if(!exists('fafileob')) {
	base::source("src/Figures/load_annotation.R")
}
#
if(!file.exists(here('data/protgnmtrids.rds'))){
	source('src/Figures/Figure0/1_ms_data_import.R')
}
protgnmtrids<-readRDS(here('data/protgnmtrids.rds'))
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
ms_mat_nonorm <- ms_mat 
ms_mat = proDA::median_normalization(ms_mat)

################################################################################
########Now run proDA
################################################################################

#get our matrix and design df
sample_info_df <- data.frame(name = colnames(ms_mat),
                             stringsAsFactors = FALSE)
sample_info_df$time <- sample_info_df$name%>%str_extract('[^_]+')
sample_info_df$assay <- sample_info_df$name%>%str_extract('(?<=_)[^_]+')


if(!file.exists(here('data/proDAfitms.rds'))){
	proDAfitms <- proDA::proDA(
			ms_mat, design = ~ time, 
            col_data = sample_info_df%>%filter(assay=='MS'),
            # reference_level = "E13"
            reference_level = NULL
)
	saveRDS(proDAfitms,here('data/proDAfitms.rds'))
}else{
	proDAfitms<-readRDS(here('data/proDAfitms.rds'))
}
tps = sample_info_df$time%>%unique
sample_info_df[[tps[[1]]]]=TRUE
for(i in seq_along(tps)){
	sample_info_df[[tps[i]]] = sample_info_df[['time']]%in%tps[length(tps):i]
}
stepdesign = as.formula(paste0('~',paste0(collapse='+',tps[-1])))

if(!file.exists(here('data/stepproDAfitms.rds'))){
		
	stepproDAfitms <- proDA::proDA(
				ms_mat, design = stepdesign, 
	            col_data = sample_info_df%>%filter(assay=='MS'),
	            # reference_level = "E13"
	            reference_level = NULL
	)
	saveRDS(stepproDAfitms,here('data/stepproDAfitms.rds'))
	
}else{
	stepproDAfitms<-readRDS(here('data/stepproDAfitms.rds'))
}

timevals = paste0('time',sample_info_df$time%>%unique%>%tail(-1))%>%setNames(.,.)

#now get the prediction ses
times = sample_info_df%>%filter(assay=='MS')%>%.$time%>%unique

contrasts = lapply(1:5,function(.) as.numeric(1:5%in%c(1,.)))
contrasts%<>%setNames(times)
preddf<- map_df(.id='time',contrasts, ~proDA::test_diff(proDAfitms,contrast=.))
preddf%<>%rename('ms_id'='name')

preddf <- preddf %>% 	
	mutate(estimate=diff)%>%
	mutate(CI.L = diff - (se*1.96))%>%
	mutate(CI.R = diff + (se*1.96))
#also add the se se to the E13 ones, they don't have it for some reason
preddf%<>%mutate(se = ifelse(is.na(se),(CI.R-CI.L)/3.92,se))


contrasts = lapply(1:5,function(.) as.numeric(1:5%in%c(.)))%>%setNames(tps)
contrdf<- map_df(.id='time',contrasts[-1], ~proDA::test_diff(proDAfitms,contrast=.))
contrdf%<>%rename('ms_id'='name')

contrdf <- contrdf %>% 	
	mutate(estimate=diff)%>%
	mutate(CI.L = diff - (se*1.96))%>%
	mutate(CI.R = diff + (se*1.96))
#also add the se se to the E13 ones, they don't have it for some reason
contrdf%<>%mutate(se = ifelse(is.na(se),(CI.R-CI.L)/3.92,se))

contrasts = paste0(tps[-1],'TRUE')%>%setNames(tps[-1])
stepstepprotcontrdf<- map_df(.id='time',contrasts, ~ proDA::test_diff(stepproDAfitms,contrast=.))
stepstepprotcontrdf$time%<>%str_replace('TRUE','')
stepstepprotcontrdf%<>%rename('ms_id'='name')

stepstepprotcontrdf <- stepstepprotcontrdf %>% 	
	mutate(estimate=diff)%>%
	mutate(CI.L = diff - (se*1.96))%>%
	mutate(CI.R = diff + (se*1.96))
#also add the se se to the E13 ones, they don't have it for some reason
stepstepprotcontrdf%<>%mutate(se = ifelse(is.na(se),(CI.R-CI.L)/3.92,se))

# #make combined matrix
proda_ms_pred <- proDAfitms$coefficients
proda_ms_pred[,2:5]=(proda_ms_pred[,2:5] + proda_ms_pred[,1])
colnames(proda_ms_pred)[1]='timeE13'
colnames(proda_ms_pred)%<>%str_replace('time(.*)','\\1_MS')

ms_metadf <- tibble(
	gene_name = proteinmsdatamatch$gene_name,
	gene_id = proteinmsdatamatch$g_id,
	ms_id = proteinmsdatamatch$ms_id,
	n_missing = ms_mat%>%apply(2,is.na)%>%rowSums,
	median = ms_mat%>%matrixStats::rowMedians(.)
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
rplids = rids%>%filter(!`Mito-RP`,`RPL_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist
rpsids = rids%>%filter(!`Mito-RP`,`RPS_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist
mrpids = rids%>%filter(`Mito-RP`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist
ms_metadf$is_rpl = sep_element_in(ms_metadf$ms_id,rplids)
ms_metadf$is_rps = sep_element_in(ms_metadf$ms_id,rpsids)
ms_metadf$is_mito_rp = sep_element_in(ms_metadf$ms_id,mrpids)

ms_metadf$is_rpl%>%table

preddf%<>%safe_left_join(distinct(ms_metadf,gene_name,gene_id,ms_id))

sel_ms_ids = ms_metadf%>%group_by(gene_name)%>%
	arrange(n_missing,desc(median))%>%slice(1)%>%.$ms_id
ms_metadf$best = ms_metadf$ms_id %in% sel_ms_ids

contrdf = contrdf%>%safe_left_join(ms_metadf%>%select(ms_id,gene_id,gene_name))%>%
	filter(ms_id %in% sel_ms_ids)

stepstepprotcontrdf <- stepstepprotcontrdf%>%
	safe_left_join(ms_metadf%>%select(ms_id,gene_id,gene_name))%>%
	filter(ms_id %in% sel_ms_ids)

sel_ms_mat <- ms_mat
sel_ms_mat = sel_ms_mat[sel_ms_ids,]
rownames(sel_ms_mat) = proteinmsdatamatch$g_id[match(rownames(sel_ms_mat),proteinmsdatamatch$ms_id)]
sel_prodpreds <- preddf%>%filter(ms_id %in% sel_ms_ids)

proDAfitms%>%saveRDS('data/proDAfitms.rds')
sel_prodpreds%>%saveRDS('data/sel_prodpreds.rds')
sel_ms_mat%>%saveRDS('data/sel_ms_mat.rds')
contrdf%>%saveRDS('data/contrdf.rds')
stepstepprotcontrdf%>%saveRDS('data/stepcontrdf.rds')


proDAfitms<-readRDS('data/proDAfitms.rds')
sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
sel_ms_mat<-readRDS('data/sel_ms_mat.rds')
contrdf<-readRDS('data/contrdf.rds')



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

#add number of missing stages
# ms_metadf <- ms_metadf%>%safe_left_join(by='gene_id',
anymissingdf = sel_ms_mat[,]%>%
    as.data.frame%>%
    rownames_to_column('gene_id')%>%
    # filter(!is.na(gene_id))%>%
    gather(dataset,signal,-gene_id)%>%
    separate(dataset,into=c('time','assay','replicate'))%>%
    group_by(gene_id,time)%>%summarise(missing=any(is.na(signal)))%>%
    mutate(time = paste0(time,'_anymissing'))%>%
    spread(time,missing)

ms_metadf <- ms_metadf%>%safe_left_join(anymissingdf,by='gene_id')

ms_metadf%>%saveRDS('data/ms_metadf.rds')

ms_metadf<-readRDS('data/ms_metadf.rds')

bestmeta <- ms_metadf%>%filter(best)
sel_ms_mat_nonorm <- ms_mat_nonorm[bestmeta$ms_id,]%>%
	set_rownames(bestmeta$gene_name)

sel_ms_mat_nonorm%>%saveRDS('data/sel_ms_mat_nonorm.rds')
sel_ms_mat%>%colMedians(na.rm=T)%>%txtplot
ms_mat%>%colMedians(na.rm=T)%>%txtplot
median_normalization(ms_mat)%>%colMedians(na.rm=T)%>%txtplot
