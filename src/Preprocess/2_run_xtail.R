#from : /Users/dharnet/bih_clustersrc/R/Load_data/integrate_countdata.R
{
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,here))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,stringr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tibble))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,magrittr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,assertthat))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,data.table))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tidyverse))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,DESeq2))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,xtail))

if(!exists('tx_countdata')) load('data/1_integrate_countdata.R')
# if(!exists('tx_countdata')) tx_countdata<-readRDS('data/tx_countdata.rds')
count_tbl <- tx_countdata$counts%>%as.data.frame%>%rownames_to_column('gene_id')
sampfile<-here('pipeline/sample_parameter.csv')
timepoints = fread(sampfile)$time%>%unique
samples = fread(sampfile)%>%
  filter(is.na(fraction))%>%
  filter(!str_detect(sample_id,'test'))%>%
  .$sample_id
gid2gnmv <-  ids_nrgname%>%distinct(gene_id,gene_name)%>%{setNames(.$gene_name,.$gene_id)}
count_tbl$gene_name <- gid2gnmv[count_tbl$gene_id]
count_tbl = count_tbl%>%select(gene_name,one_of(samples))

# uorfcounts <- fread(uORFcountfile)
# uorfcounts = uorfcounts%>%select(gene_name,one_of(samples))
# count_tbl <- count_tbl%>%rbind(uorfcounts)
xtailfiles = here('pipeline',paste0('xtail','/xtail_',timepoints[-1],'.txt'))%>%setNames(timepoints[-1])
timepoints = timepoints
my_contrasts = timepoints[-1] %>% setNames(.,.)
matches = dplyr::matches
# dir.create(xtailfiles[1]%>%dirname)
}


for(tp2 in (timepoints[-c(1)])){
    #samples to use
    allxtailcounts =   
      count_tbl%>%
      dplyr::select(matches(paste0('gene_name|',tp1,'|',tp2)))
    
    xtailcounts = allxtailcounts[!apply(allxtailcounts,1,.%>%is.na%>%any),]
    mrnatab <- xtailcounts%>%select(matches('total'))%>%as.data.frame%>%set_rownames(xtailcounts$gene_name)
    ribotab <- xtailcounts%>%select(matches('ribo'))%>%as.data.frame%>%set_rownames(xtailcounts$gene_name)
    
    conditionvect <-    ifelse(str_detect(colnames(xtailcounts%>%select(matches('total'))),tp1),tp1,tp2)   

    xtailres <- xtail::xtail(mrnatab,ribotab,conditionvect,threads=20)

    xtailtable <- xtailres$resultsTable%>%rownames_to_column%>%set_colnames(
      c("gene_name","mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1", "pvalue_v1", paste0(tp1,"_log2TE"), 
      paste0(tp2,"_log2TE"), "log2FC_TE_v2", "pvalue_v2", "log2fc", 
      "p_value", "adj_p_value")
  )
    write_tsv(xtailtable,xtailfiles[tp2])
}

isgenenames <- function(v) v%>%head%>%str_length%>%`<`(8)%>%sum%>%`>`(8)
# 
for(tp2 in (timepoints[-c(1)])){
  xtailtable <-  read_tsv(xtailfiles[tp2])
  base_means <- rowMeans(assayData(countexprdata)[['exprs']]%>%sweep(.,2,FUN='/',STAT=DESeq2::estimateSizeFactorsForMatrix(.)))
  base_means <- base_means[fData(countexprdata)$protein_id[fData(countexprdata)$is_gid_highest]]
  names(base_means) <- fData(countexprdata)$gene_name[fData(countexprdata)$is_gid_highest]
  xtailtable$base_mean <- base_means[xtailtable$gene_name]
  if(!isgenenames(xtailtable$gene_name)) {

    xtailtable$gene_name <- xtailtable$gene_name
    xtail_exprmatch <- match(xtailtable$gene_name,fData(countexprdata)$gene_name)
    xtailtable$gene_name <- fData(countexprdata)$gene_id[xtail_exprmatch]
  } 

  write_tsv(xtailtable,xtailfiles[[tp2]])

}



#now do things stepwise
for(i in seq_along(timepoints)[-1]){
    tp1 = timepoints[i-1]
    tp2 = timepoints[i]
    #samples to use
    allxtailcounts =   
      count_tbl%>%
      dplyr::select(matches(paste0('gene_name|',tp1,'|',tp2)))
    #
    xtailcounts = allxtailcounts[!apply(allxtailcounts,1,.%>%is.na%>%any),]
    mrnatab <- xtailcounts%>%select(matches('total'))%>%as.data.frame%>%set_rownames(xtailcounts$gene_name)
    ribotab <- xtailcounts%>%select(matches('ribo'))%>%as.data.frame%>%set_rownames(xtailcounts$gene_name)
    #
    conditionvect <-    ifelse(str_detect(colnames(xtailcounts%>%select(matches('total'))),tp1),tp1,tp2)   
    #
    nz_mrna1 <- mrnatab[,1:2]%>%rowSums%>%`!=`(0)
    nz_mrna2 <- mrnatab[,3:4]%>%rowSums%>%`!=`(0)
    # nz_ribo1 <- ribotab[,1:2]%>%rowSums%>%`!=`(0)
    # nz_ribo2 <- ribotab[,3:4]%>%rowSums%>%`!=`(0)
    nz_c1 <- cbind(
      ribotab[,conditionvect==conditionvect[1]],
      mrnatab[,conditionvect==conditionvect[1]]
    )%>%rowSums%>%`!=`(0)
    nz_c2 <- cbind(
      ribotab[,conditionvect==conditionvect[2]],
      mrnatab[,conditionvect==conditionvect[2]]
    )%>%rowSums%>%`!=`(0)

    nz_c1['Try5']
    nz_c2['Try5']

    # nz_ribo <- ribotab%>%rowSums%>%`!=`(0)
    # nz_ribo <- ribotab[,12]%>%rowSums%>%`!=`(0)
    nz_both = nz_c1 & nz_c2 & nz_mrna1 & nz_mrna2
    nz_bothgnms = intersect(names(nz_both[nz_both]),highcountgnms)
    # xtailres <- xtail(mrnatab[nz_both,],ribotab[nz_both,],conditionvect,threads=20)
    xtailres <- xtail(mrnatab[nz_bothgnms,],ribotab[nz_bothgnms,],conditionvect,threads=20)
    #
    xtailtable <- xtailres$resultsTable%>%as.data.frame%>%rownames_to_column%>%set_colnames(
      c("gene_name","mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1", "pvalue_v1", paste0(tp1,"_log2TE"), 
      paste0(tp2,"_log2TE"), "log2FC_TE_v2", "pvalue_v2", "log2fc", 
      "p_value", "adj_p_value")
    )
    #
    xtailfile = here('pipeline',paste0('xtail','/xtail_',tp2,'_v_',tp1,'.txt'))
    #
    write_tsv(xtailtable,xtailfile)
}




isgenenames <- function(v) v%>%head%>%str_length%>%`<`(8)%>%sum%>%`>`(8)

for(tp2 in (timepoints[-c(1)])){
  xtailtable <-  read_tsv(xtailfiles[tp2])
  base_means <- rowMeans(assayData(countexprdata)[['exprs']]%>%sweep(.,2,FUN='/',STAT=DESeq2::estimateSizeFactorsForMatrix(.)))
  base_means <- base_means[fData(countexprdata)$protein_id[fData(countexprdata)$is_gid_highest]]
  names(base_means) <- fData(countexprdata)$gene_name[fData(countexprdata)$is_gid_highest]
  xtailtable$base_mean <- base_means[xtailtable$gene_name]
  if(!isgenenames(xtailtable$gene_name)) {

    xtailtable$gene_name <- xtailtable$gene_name
    xtail_exprmatch <- match(xtailtable$gene_name,fData(countexprdata)$gene_name)
    xtailtable$gene_name <- fData(countexprdata)$gene_id[xtail_exprmatch]
  } 

  write_tsv(xtailtable,xtailfiles[[tp2]])

}

