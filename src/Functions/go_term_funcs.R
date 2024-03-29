#This will screen the results of our DEseeq analysis using GO terms.
library(biomaRt)
library(tidyverse)
library(magrittr)
library(assertthat)
select<-dplyr::select

#' GO Enrichment Plot
#'
#' @param go_table a data frame with GO enrichment results as produced by
#'   'run_go_enrich'
#' @param sort_var string, the variable in 'go_table' used to sort rows
#' @param title_str string, a title for the enrichment plot
#' @param num_top integer, number of top GO terms to be included
#'
#' @return a ggplot2 object
#' @export
plot_go_enrich <- function(go_table, sort_var, title_str, num_top = 10){
  assert_that(is.string(sort_var))
  assert_that(is.string(title_str))
  assert_that(is.count(num_top))
  assert_that(go_table %has_name% sort_var)
  assert_that(go_table %has_name% "Significant")
  assert_that(go_table %has_name% "Annotated")
  assert_that(go_table %has_name% "Term")
  #
  go_table$p_value <- go_table[[sort_var]]
  stopifnot(!all(is.na(go_table$p_value))) 
  #add plotting variables
  got_table <-
    go_table %>%
    dplyr::arrange(p_value) %>%
    dplyr::mutate(gene_ratio = Significant/Annotated) %>%
    head(num_top) %>%
    filter(Significant!=0)%>%
    dplyr::mutate(Term = factor(Term, levels = rev(Term)))
  
  #log cale breaks for the pvalues
  coltranslims = got_table$p_value %>%na.omit%>% log10%>%multiply_by(-1)%>%{c(floor(min(.)),ceiling(max(.)))}
  stopifnot(is.numeric(coltranslims))
  n_colbreaks = min(5,(coltranslims[2]-coltranslims[1])+1) # how many breaks int he pvalue col scale?
  col_label_format <- function(x) format(1/(10^x),digits = 2) 
  #log scale breaks for the number of genes
  sizetranslims = got_table$Significant %>%na.omit%>% log2%>%{c(floor(min(.)),ceiling(max(.)))}
  message(sizetranslims)
  n_sizebreaks = min(5,(sizetranslims[2]-sizetranslims[1])+1) # how many breaks int he pvalue col scale?
  size_label_format <- function(x) format((2^x),digits = 2) 
  

  got_table$Term %<>% as.character%>%linebreaklongsent%>%as.factor
  got_table$Term %<>% forcats::fct_reorder(got_table$gene_ratio)

  #now plot
  go_plot <-
    ggplot(got_table,aes(size = log2(Significant) , y = Term, color = -log10(p_value), x = gene_ratio)) + 
    geom_point() + 
    ggplot2::theme_bw() + 
    ggplot2::labs(size = "# significant  genes", y = "go term", color = "p-value", x = "gene-ratio", title=title_str) + 
    #apply transformed color scale
    ggplot2::scale_color_continuous(
      limits = coltranslims,
      breaks = round(seq(from=min(coltranslims),to=max(coltranslims),len=n_colbreaks)),
      labels = col_label_format,high = 'darkgreen',low='lightgreen'
    )+
    #apply transformed size scale
    ggplot2::scale_size_continuous(
      limits = sizetranslims,
      breaks = round(seq(from=min(sizetranslims),to=max(sizetranslims),len=n_sizebreaks)),
      labels = size_label_format
    )+
    theme(axis.text=element_text(size=7), axis.title=element_text(size=14,face="bold"))
  #and return
  return(go_plot)
}

linebreaklongsent <-function(sent,n=3){
  for(i in seq_along(sent)){
    sent[i]%<>%str_split(' ')%>%map(.%>%{suppressWarnings(split(.,floor((seq(0,length(.))/3))))}%>%map(paste0,collapse=' ')%>%paste0(collapse='\n'))
  }
  sent%>%unlist
}
linebreaklongsent(sent=c('fdasfdas fdsafdas fdsafdsa fdsafdsa fdsafdas','fdsafdas fdsafdsa'))
#For a predfined list of genes
#
if(!file.exists('data/GTOGO.rds')){
  require(biomaRt)
  mart<-useMart('ENSEMBL_MART_ENSEMBL')
  mart<-useDataset('mmusculus_gene_ensembl',mart)
  GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id", "go_id","external_gene_name"), mart = mart)
  GTOGO%<>%select(gene_name=external_gene_name,go_id,ensembl_gene_id)
  dir.create('data')
  saveRDS(GTOGO,'data/GTOGO.rds')
}else{
  GTOGO<-readRDS('data/GTOGO.rds')
}

GTOGO %<>% select(gene_name, go_id, ensembl_gene_id)

rungo <- function(int.genes,GTOGO,my_ontology,background=NULL,algo='elim'){
  if(is.logical(int.genes)){
    stopifnot(any(int.genes))
    background = names(int.genes)
    int.genes = names(int.genes)[int.genes]
  } 
  if(!is.null(background)) GTOGO <- GTOGO%>%filter(ensembl_gene_id%in%background)
  geneID2GO <- by(GTOGO$go_id, GTOGO$ensembl_gene_id, function(x) as.character(x))
  
  require(topGO)


  all.genes <- sort(unique(as.character(GTOGO$ensembl_gene_id)))
  int.genes <- factor(as.integer(all.genes %in% int.genes))
  stopifnot(any(int.genes==1))
  names(int.genes) = all.genes

  
  go.obj <- new("topGOdata", ontology=my_ontology
                , allGenes = int.genes
                , annot = annFUN.gene2GO
                , gene2GO = geneID2GO,
                nodeSize=5
  )
  
  #?runTest
  results     <- runTest(go.obj, algorithm = algo, statistic = "fisher")
  results.tab <- GenTable(object = go.obj, elimFisher = results,topNodes = 100)
  results.tab%<>%mutate(Enrichment = Significant / Expected )
  results.tab%<>%mutate(elimFisher = as.numeric(elimFisher) )
  results.tab
}

plot_go_tree <- function(toptermgotable,ont='BP',tit_postfix=''){
  par(mar=c(1, 4,1,20),cex=0.8)
  c('Term','gene_name') %in% toptermgotable
  terms <- toptermgotable%>%filter(ontology==ont)
  stopifnot(nrow(terms)>0)
  
  
  terms%>%
    distinct(Term,gene_name)%>%
    mutate(tmp=TRUE)%>%
    spread(gene_name,tmp)%>%
    map_df(replace_na,F)%>%
    {set_rownames(as.matrix(.[,-1]),.$Term)}%>%
    dist('manhattan')%>%
    hclust%>%
    as.dendrogram%>%
    plot( cex = 0.2, horiz=TRUE, main = paste0('GO Terms clustered by overlap in significant genes:\n ',ont,' ',tit_postfix))
}

GTOGO%<>%dplyr::mutate(gene_id=ensembl_gene_id)

onts = c('BP','MF','CC')

get_cluster_gos <- function(clustvect){
    names(clustvect) %<>% gnm2gidv[.]
    out=map_df(.id='ontology',onts%>%setNames(.,.),function(ont){
      lapply(unique(clustvect)%>%setNames(.,.),safely(function(clustval){
            gids <- names(clustvect)
            stopifnot(mean(names(clustvect)%in%GTOGO$gene_id)>.9)
            filtGTOGO <- GTOGO %>%filter(gene_id %in%names(clustvect))
             projmemoise(rungo)(
              gids[clustvect==clustval],
              filtGTOGO,
              ont,
              algo='classic'
            )
        }))%>%map('result')%>%bind_rows(.id='cluster')
    })
    out
}

go_comparison_plot <-function(x){x%>%map_df(.id='cluster',.%>%arrange(elimFisher)%>%head(10)%>%select(Term,elimFisher))%>%
  mutate(Term = as_factor(Term))%>%
  mutate(cluster = as_factor(cluster))%>%
  ggplot(.,aes(x=cluster,y=Term,color=-log10(elimFisher),size=-log10(elimFisher)))+
  geom_point()+
  theme_bw()
}

