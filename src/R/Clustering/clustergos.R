
#mark satb2
clustergos <- clustlist%>%map('cluster')%>%
  lapply(function(vct){names(vct) %<>% gnm2gid[[.]];vct})%>%
  lapply(function(clustvect){
    map_df(.id='ontology',onts%>%setNames(.,.),function(ont){
      mclapply(unique(clustvect)%>%setNames(.,.),safely(function(clustval){
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
})

source('src/R/Functions/go_term_funcs.R')
stopifnot(exists("GTOGO"))

GTOGO%<>%dplyr::mutate(gene_id=ensembl_gene_id)

gnm2gid = metainfo%>%filter(isbest)%>%distinct(gene_name,gene_id)%>%{safe_hashmap(.[[1]],.[[2]])}
onts = c('BP','MF','CC')


#clustergos%<>%.[names(clustlist)]
#memoise
go_comparison_plot <- projmemoise(go_comparison_plot)

