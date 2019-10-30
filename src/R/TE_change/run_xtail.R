suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,here))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,stringr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tibble))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,magrittr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,assertthat))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,data.table))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tidyverse))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,DESeq2))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,xtail))

args <- c(
  countfile='pipeline/exprdata/countexprset.rds',
  uORFcountfile='SaTAnn/uORFs.feature_counts',
  outdir= 'xtail'
)

if(!interactive()) args[] = commandArgs(trailingOnly=TRUE)
for(i in names(args)) assign(i,args[i])

countexprdata <- readRDS(countfile)
# count_tbl <- data.table::fread(countfile)
count_tbl <- assayData(countexprdata)[['exprs']][countexprdata@featureData$is_gid_highest,]%>%
  as.data.table%>%mutate(feature_id=fData(countexprdata)$gene_name[fData(countexprdata)$is_gid_highest])
# oldcounttbl<-count_tbl
count_tbl$feature_id%>%table%>%.[.>1]
count_tbl%>%{.$feature_id=='Satb2'}%>%which

sampfile<-here('pipeline/sample_parameter.csv')
timepoints = fread(sampfile)$time%>%unique
samples = fread(sampfile)%>%
  filter(is.na(fraction))%>%
  filter(!str_detect(sample_id,'test'))%>%
  .$sample_id
count_tbl = count_tbl%>%select(feature_id,one_of(samples))

# uorfcounts <- fread(uORFcountfile)
# uorfcounts = uorfcounts%>%select(feature_id,one_of(samples))
# count_tbl <- count_tbl%>%rbind(uorfcounts)

xtailfiles = here('pipeline',paste0(outdir,'/xtail_',timepoints[-1],'.txt'))%>%setNames(timepoints[-1])

timepoints = timepoints
my_contrasts = timepoints[-1] %>% setNames(.,.)

tp1 = timepoints[1]

matches = dplyr::matches
# dir.create(xtailfiles[1]%>%dirname)



# for(tp2 in timepoints[2:]){
for(tp2 in (timepoints[-c(1)])){

    #samples to use
    allxtailcounts =   
      count_tbl%>%
      dplyr::select(matches(paste0('feature_id|',tp1,'|',tp2)))
    
    xtailcounts = xtailcounts[!apply(xtailcounts,1,.%>%is.na%>%any),]
    mrnatab <- xtailcounts%>%select(matches('total'))%>%as.data.frame%>%set_rownames(xtailcounts$feature_id)
    ribotab <- xtailcounts%>%select(matches('ribo'))%>%as.data.frame%>%set_rownames(xtailcounts$feature_id)
    
    conditionvect <-    ifelse(str_detect(colnames(xtailcounts%>%select(matches('total'))),tp1),tp1,tp2)   

    stop()



    xtailres <- xtail::xtail(mrnatab,ribotab,conditionvect,threads=20)

    xtailtable <- xtailres$resultsTable%>%rownames_to_column%>%set_colnames(
      c("feature_id","mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1", "pvalue_v1", paste0(tp1,"_log2TE"), 
      paste0(tp2,"_log2TE"), "log2FC_TE_v2", "pvalue_v2", "log2fc", 
      "p_value", "adj_p_value")
    )

    write_tsv(xtailtable,xtailfiles[tp2])
  
}
isgenenames <- function(v) v%>%head%>%str_length%>%`<`(8)%>%sum%>%`>`(8)

for(tp2 in (timepoints[-c(1)])){
  xtailtable <-  read_tsv(xtailfiles[tp2])
  base_means <- rowMeans(assayData(countexprdata)[['exprs']]%>%sweep(.,2,FUN='/',STAT=DESeq2::estimateSizeFactorsForMatrix(.)))
  base_means <- base_means[fData(countexprdata)$protein_id[fData(countexprdata)$is_gid_highest]]
  names(base_means) <- fData(countexprdata)$gene_name[fData(countexprdata)$is_gid_highest]
  xtailtable$base_mean <- base_means[xtailtable$feature_id]
  if(!isgenenames(xtailtable$feature_id)) {

    xtailtable$gene_name <- xtailtable$feature_id
    xtail_exprmatch <- match(xtailtable$feature_id,fData(countexprdata)$gene_name)
    xtailtable$feature_id <- fData(countexprdata)$gene_id[xtail_exprmatch]
  } 

  write_tsv(xtailtable,xtailfiles[[tp2]])

}


# count_tbl[countexprdata@featureData[countexprdata@featureData$is_gid_highest]$gene_name=='Satb2',]
# xtailtable[countexprdata@featureData[countexprdata@featureData$is_gid_highest]$gene_name=='Satb2',]

# fData(countexprdata)%>%data.frame%>%.[11552,]


# library(forcats)

# coldatadeseq <- phenoData(countexprdata)@data
# coldatadeseq$assay%<>%as_factor%>%fct_relevel(c('total','ribo'))
# coldatadeseq$time%<>%timepoints[.]%>%factor

# DESeqte<-DESeqDataSetFromMatrix(as.data.frame(count_tbl%>%.[,rownames(phenoData(countexprdata))]%>%as.matrix%>%set_rownames(count_tbl[[1]])),coldatadeseq,design= ~ time*assay)
# tedds <- DESeq(DESeqte)
# resultsNames(tedds)

# deseqte <- results(tedds,c(0,0,0,0,0,0,0,0,0,1))

# library(txtplot)

# deseqte%>%as.data.frame%>%rownames_to_column('gene_name')%>%left_join(xtailtable)%>%
#   {cor(y=.$log2FoldChange,x=.$log2fc,use='complete')}

# deseqte%>%as.data.frame%>%rownames_to_column('gene_name')%>%left_join(xtailtable)%>%
#   {cor(use='complete',y=-log10(.$padj),x=-log10(.$adj_p_value))}




# #####NOw with splines
# library(splines)
# deseqmat <- as.data.frame(count_tbl%>%.[,rownames(phenoData(countexprdata))]%>%as.matrix%>%set_rownames(count_tbl[[1]]))
# coldatadeseqspl <- coldatadeseq%>%mutate(time = as.numeric(time))
# deseqte_spl <- DESeqDataSetFromMatrix(deseqmat,coldatadeseqspl,design= ~ ns(time,3)*assay)
# tedds <- DESeq(deseqte_spl)

# resultsNames(tedds)

# deseqte <- results(tedds,c(0,0,0,0,0,ns(1:5,3)[5,]))

# library(txtplot)

# deseqte%>%as.data.frame%>%rownames_to_column('gene_name')%>%left_join(xtailtable)%>%
#   {cor(y=.$log2FoldChange,x=.$log2fc,use='complete')}

# deseqte%>%as.data.frame%>%rownames_to_column('gene_name')%>%left_join(xtailtable)%>%
#   {cor(use='complete',y=-log10(.$padj),x=-log10(.$adj_p_value))}

# deseqte$padj%>%`<`(0.05)%>%table

# xtailtable$adj_p_value%>%`<`(0.05)%>%table

# txtplot(-log10(deseqte$padj),-log10(xtailtable$adj_p_value))

# xtailtable%>%filter(gene_name=='Satb2')%>%as.data.frame
# deseqte%>%as.data.frame%>%rownames_to_column('gene_name')%>%left_join(xtailtable%>%select(feature_id,gene_name ))%>%filter(gene_name=='Satb2')%>%as.data.frame

# deseqte






stopifnot(xtailfiles%>%file.exists)

library(txtplot)

xtailtable%>%.$E145_log2TE%>%na.omit%>%txtdensity

satb2fdata <- fData(countexprdata)[fData(countexprdata)$gene_name=='Satb2',]
exprs(countexprdata)[fData(countexprdata)$gene_name=='Satb2',]
xtailtable[count_tbl$feature_id %in% satb2fdata$protein_id,]




xtailfiles%>%map_df(fread)%>%group_by(feature_id)%>%filter(any(adj_p_value<0.05))%>%slice(1)

techange <- xtailfiles%>%map_df(fread)%>%group_by(feature_id)%>%filter(any((adj_p_value<0.05) & (abs(log2fc) > 0.32)))%>%slice(1)

#I'd also like to try DESeq2 with splines
oldxtailfiles <- xtailfiles%>%str_replace('xtail/' ,'xtail_old/')


ovtable(
oldxtailfiles[[4]]%>%fread%>%filter(adj_p_value<0.05)%>%.$feature_id,
xtailfiles[[4]]%>%fread%>%mutate(feature_id = countexprdata@featureData$gene_id[as.numeric(as.character(.$feature_id))])%>%filter(adj_p_value<0.05)%>%.$feature_id
)

#Is satb2 delta TE?
xtailfiles[[4]]%>%fread%>%mutate(feature_id = countexprdata@featureData$gene_name[as.numeric(as.character(.$feature_id))])%>%filter(feature_id%>%str_detect('Satb2'))


ovtable<-function(x,y){
  all=union(x,y)
  table(all%in%unique(x),all%in%unique(y))
}
xtailfiles

fread(here('pipeline',oldxtailfiles[[1]]))






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
