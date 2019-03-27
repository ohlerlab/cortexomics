library(rtracklayer)
library(stringr)
library(magrittr)
library(assertthat)
library(tidyverse)
library(here)

select =dplyr::select
source('~/Dropbox/code_new/df_functions.R')

kallistodir = here("./data/kallisto/")
metacols = c('time','assay','replicate')

annotation_gr <- rtracklayer::import("~/bih_cluster/projects/cubit/current/static_data/annotation/GENCODE/M12/GRCm38/gencode.vM12.annotation.gff3")
annotation_gr$stripped_gene_ID <- annotation_gr$gene_id%>%str_replace('\\.\\d+','')
transcript_gr <- annotation_gr %>% subset( (.$type=='transcript'))


kallistocounts <- kallistodir %>%
  list.files(full=TRUE, patt="abundance.tsv",recurs=TRUE)%>%
  setNames(.,basename(dirname(dirname(.))))%>%
  map(read_tsv,col_types=cols())%>%
  map(select,target_id,tpm)%>%
  bind_rows(.id = 'dataset')

#we need to ensure each entry of our kallisto table has a unique gene ID
stopifnot(kallistocounts$target_id%>%str_count('ENSMUSG\\d+\\.\\d+')%>%equals(1)%>%all)
kallistocounts <- kallistocounts %>% mutate(ID = str_extract(target_id,'ENSMUST\\d+\\.\\d+'))


markergenestbl <- readr::read_tsv(here::here("./data/marker_genes.tsv"),col_names = TRUE)
markergenestbl$gencode%<>%str_replace('\\.\\d+','')
markergenes <- markergenestbl$gencode %>% str_split(',') %>% flatten_chr
markergene_gr <- annotation_gr%>% subset( (.$type=='transcript') & (stripped_gene_ID %in% markergenes) )
# stopifnot(all(markergenes %in% markergene_gr$stripped_gene_ID))
# markergenes[!markergenes %in% markergene_gr$stripped_gene_ID]

#get counts for specific genes
# markergene_counts <- markergene_gr%>%mcols%>%as.data.frame%>%select(ID,gene_name)%>%left_join(htseqcounts)
kallistocounts <- kallistocounts %>% mutate(ID = str_extract(target_id,'ENSMUST\\d+\\.\\d+'))
transcript_tpm <- transcript_gr%>%mcols%>%as.data.frame%>%select(ID,transcript_name,gene_name)%>%left_join(kallistocounts%>%select(-target_id),by='ID')

#metadata for our datasets
transcript_tpm <- 
  transcript_tpm%>%
  mutate(tmp = dataset%>%str_replace('_28_30','-28-30'))%>%
  separate(tmp,metacols,sep='_')%>%
  mutate(assay = str_replace(assay,'-28-30',''))

validvalue <- function(x) is.finite(x) & Negate(is.na)(x) & Negate(is.nan)(x) & `>`(x,0) 

# './data/ms_cortexomics_dec_2017/325_new_2_all_log2_LFQ_n7464.txt'%>%readLines(10)

"/fast/users/harnettd_c/cortexomics/gdrive/MS_Ribo+Total_Proteome/"

ms_data=
  # read_tsv(here('./data/ms_cortexomics_dec_2017/325_new_1_all_n7464.txt'),comment = '#')%>%
  read_tsv(here('./data/ms_cortexomics_dec_2017/325_new_2_all_log2_LFQ_n7464.txt'),comment = '#')%>%
  {colnames(.) = str_replace_all(colnames(.),' ','_');.}%>%
  select(gene_name=Gene_names,everything())

delete_japanese<-function(dblstring) str_extract(dblstring,'(NA)|(NaN)|(\\d+\\.?\\d*)')%>%as.numeric

ms_data<-ms_data%>%mutate_at(vars(matches('^(E\\d\\d|P0)')),funs(delete_japanese))

ms_data$gene_name%>%str_subset(regex('foxp2|CAGH44|SPCH1|TNRC10|2810043D05Rik|D0Kist7',ignore_case=TRUE))
ms_data$gene_name%>%str_subset(regex('fezf2|fezl|Zfp312',ignore_case=TRUE))
ms_data$gene_name%>%str_subset(regex('cux2|cutl2|cux-2',ignore_case=TRUE))


gene_has_term <- function(gene_has_term){
  hip_cor_res.t<-as.data.frame(hip_cor_res)
  hip_cor_res.t$genes<-rownames(hip_cor_res.t)
  diff.regul.genes<-rep(0,nrow(hip_cor_res.t))
  names(diff.regul.genes)<-hip_cor_res.t$genes
  diff.regul.genes[dplyr::filter(hip_cor_res.t, abs(log2FoldChange)>0.5, padj < 0.01)$genes] <- 1
  tgd.bp.diff <- new( "topGOdata", ontology = "BP", allGenes = factor(diff.regul.genes), 
                      nodeSize = 5, annot = annFUN.org, 
                      mapping = "org.Mm.eg.db", ID = "ensembl")
  ann.genes <- genesInTerm(tgd.bp.diff)
  n=0
  for (i in hip_cor_bp$table$elim$GO.ID[1:20]){
    n=n+1
    j <- filter(hip_cor_res.t, hip_cor_res.t$genes %in% ann.genes[[i]], hip_cor_res.t$padj<0.01,  abs(hip_cor_res.t$log2FoldChange)>0.5)
    # filename=paste(i, ".xls", sep = "")
    #  filename=paste("GO_tables/BP/", filename, sep="")
    
    filename2=paste(hip_cor_bp$table$elim$Term[n], ".xls", sep = "")
    filename2=paste("GO_tables/BP/", filename2, sep="")
    
}
# #i need to deal witht he m ultiple protein groups somehow
# ms_data%>%
#   group_by(gene_name)%>%
# mutate(is_group = str_detect(Protein_IDs,';'))

ms_tall <- 
  ms_data %>%
  # select(gene_name,matches('iBAQ_'))%>%
  select(gene_name,matches('^(E\\d\\d|P0)'))%>%
  filter(!is.na(gene_name))%>%
  gather(dataset,signal,-gene_name)%>%
  mutate(dataset=paste0('LFQ_',dataset))%>%
  mutate(signal = 2^signal)

ms_tall <- ms_tall %>%
  separate(dataset,c('tmp','time','replicate'),remove = FALSE)%>%
  select(-tmp)%>%
  mutate(signal = ifelse(signal==0,NA,signal))%>%
  group_by(time,replicate)%>%
  mutate(signal = (1e6*signal) / sum (na.omit(signal)) )

stopifnot(ms_data$Protein_IDs%>%str_split(';')%>%unlist%>%duplicated%>%any%>%not)
# ms_data
# ms_data%>%colnames%>%grep(value=TRUE,inver=TRUE,pattern='iBAQ|Oxidation|Deamination',x=.)%>%sort
# ms_data%>%colnames%>%grep(value=TRUE,pattern='iBAQ',x=.)

# t=assert_that(ms_tall %has_name% 'IMS')
# t=assert_that(ms_tall %has_name% 'iBAQ')
t=assert_that(ms_tall %has_name% 'signal')
t=assert_that(ms_tall %has_name% 'time')
t=assert_that(ms_tall %has_name% 'gene_name')

ms_data$Fasta_headers%>%sample(1)



#now join the transcript and 
gene_tpm<-
  transcript_tpm%>%
  group_by(gene_name,time,assay,replicate,dataset)%>%
  summarise(TPM=sum(tpm))

#overlap in gene names
list(g=tolower(gene_tpm$gene_name),m=tolower(ms_tall$gene_name))%>%
{table(unique(c(.$g,.$m)) %in% .$g,unique(c(.$g,.$m)) %in% .$m)}

commongname = intersect(gene_tpm$gene_name,ms_tall$gene_name)

common_gene_tpm <- gene_tpm %>% semi_join(data.frame(gene_name= commongname ))%>%rename(signal=TPM)


expression_table<-
  ms_tall%>%
  filter(gene_name %in% commongname)%>%
  select(signal,gene_name,dataset,replicate,time = time)%>%
  mutate(assay='MS')%>%
  bind_rows(.,common_gene_tpm[colnames(.)])

expression_table$assay<-expression_table$assay%>%recode_factor(.,'total'='TotalRNASeq','ribo'='RiboSeq','MS'='MassSpec')

#pag4?
expression_table%>%filter(gene_name=='Pa2g4')
#format for modelling
expression_table$log_signal = log2(expression_table$signal)
expression_table$ntime = expression_table$time%>%factor%>%as.integer

