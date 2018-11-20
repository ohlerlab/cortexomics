suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,stringr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tibble))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,magrittr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,assertthat))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,data.table))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tidyverse))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,DESeq2))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,xtail))

args <- c(
  countfile='feature_counts/all_feature_counts',
  uORFcountfile='SaTAnn/uORFs.feature_counts',
  outdir= 'xtail'
)

args[] = commandArgs(trailingOnly=TRUE)
for(i in names(args)) assign(i,args[i])

featurecountsagg <- data.table::fread(countfile)
uorfcounts <- fread(uORFcountfile)

timepoints = fread('sample_parameter.csv')$time%>%unique
samples = fread('sample_parameter.csv')%>%
  filter(is.na(fraction))%>%
  filter(!str_detect(sample_id,'test'))%>%
  .$sample_id
featurecountsagg = featurecountsagg%>%select(feature_id,one_of(samples))
uorfcounts = uorfcounts%>%select(feature_id,one_of(samples))

featurecountsagg <- featurecountsagg%>%rbind(uorfcounts)

xtailfiles = file.path(paste0(outdir,'/xtail_',timepoints[-1],'.txt'))%>%setNames(timepoints[-1])

timepoints = timepoints
my_contrasts = timepoints[-1] %>% setNames(.,.)

tp1 = timepoints[1]

for(tp2 in timepoints[-1]){
    #samples to use
    xtailcounts =   
      featurecountsagg%>%
      dplyr::select(matches(paste0('feature_id|',tp1,'|',tp2)))
    
    xtailcounts = xtailcounts[!apply(xtailcounts,1,.%>%is.na%>%any),]
    mrnatab <- xtailcounts%>%select(matches('total'))%>%as.data.frame%>%set_rownames(xtailcounts$feature_id)
    ribotab <- xtailcounts%>%select(matches('ribo'))%>%as.data.frame%>%set_rownames(xtailcounts$feature_id)
    conditionvect <-       ifelse(str_detect(colnames(xtailcounts%>%select(matches('total'))),tp2),tp1,tp2)   

    xtailres <- xtail::xtail(mrnatab,ribotab,conditionvect,threads=20)

    xtailtable <- xtailres$resultsTable%>%rownames_to_column%>%set_colnames(
      c("feature_id","mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1", "pvalue_v1", "E145_log2TE", 
      "E13_log2TE", "log2FC_TE_v2", "pvalue_v2", "log2fc", 
      "p_value", "adj_p_value")
    )

    write_tsv(xtailtable,xtailfiles[tp2])
  
}


# save.image(file.path('prepimage.R'))


#let's make some toy images. Genes



# colnames(xtailres)%<>%str_replace_all('[\\(\\s\\)]','_')
# 

# annotation_gr <- rtracklayer::import("~/bih_cluster/projects/cubit/current/static_data/annotation/GENCODE/M12/GRCm38/gencode.vM12.annotation.gff3")
# 
# ribodiffres<-left_join(
#   ribodiffres,
#   annotation_gr%>%mcols%>%as.data.frame%>%select(gene_id,gene_name)%>%distinct(gene_id,gene_name),
#   by=c('geneIDs'='gene_id')
# )%>%  select(-disper,-pval,-X8)%>%
#   select(gene_name,everything())
# 
# ribodiffressig = ribodiffres %>% filter(padj<0.05,!is.nan(padj))
# 
# ribodiffressig%>%arrange(log2FC_TE_P0_vs_E13_)%>%distinct(gene_name,.keep_all=TRUE)
# ribodiffressig%>%arrange(-log2FC_TE_P0_vs_E13_)%>%distinct(gene_name,.keep_all=TRUE)
# 
# ribodiffres %<>% mutate(l2fc_te = `log2FC_TE(P0 vs E13)`)
# ribodiffres %<>% mutate(log10_mean_te = log10(TEE13+TEP0) / 2)
# 
# ribodiffres$aveLogCPM <-
#   edgeR::aveLogCPM(ribodiffcounts%>%select(-Entry)%>%as.matrix%>%as.integer)%>%
#   setNames(ribodiffcounts$Entry)%>%
#   .[ribodiffres$geneIDs]
#   
# ribodiffres%>%{qplot(data=.,x =aveLogCPM,y=l2fc_te,color=padj<0.05)}
