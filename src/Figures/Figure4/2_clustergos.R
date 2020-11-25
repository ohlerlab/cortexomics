

source('src/R/Functions/go_term_funcs.R')
stopifnot(exists("GTOGO"))

GTOGO%<>%dplyr::mutate(gene_id=ensembl_gene_id)

# gnm2gid = metainfo%>%filter(isbest)%>%distinct(gene_name,gene_id)%>%{safe_hashmap(.[[1]],.[[2]])}
onts = c('BP','MF','CC')

#mark satb2
# ont='BP'
# clustvect=hclustob$cluster
# clustval = unique(clustvect)%>%setNames(.,.)%>%head(1)
get_cluster_gos <- function(clustvect){
    names(clustvect) %<>% gnm2gid[[.]]
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

#clustergos%<>%.[names(clustlist)]
#memoise

go_comparison_plot <-function(x)x%>%map_df(.id='cluster',.%>%arrange(elimFisher)%>%head(10)%>%select(Term,elimFisher))%>%
  mutate(Term = as_factor(Term))%>%
  mutate(cluster = as_factor(cluster))%>%
  ggplot(.,aes(x=cluster,y=Term,color=-log10(elimFisher),size=-log10(elimFisher)))+
  geom_point()+
  # ggtitle(mapname)+
  theme_bw()


# go_comparison_plot <- projmemoise(go_comparison_plot)

# #now plot
# plotfile<- here(paste0('plots/','cluster_go_bp','.pdf'));cairo_pdf(h=21,w=21,plotfile)
# go_comparison_plot(clustergos[[1]]%>%filter(ontology=='BP')%>%{split(.,.$cluster)})+ggtitle('Biological Process')
# dev.off()
# message(normalizePath(plotfile))
# plotfile<- here(paste0('plots/','cluster_go_mf','.pdf'));cairo_pdf(h=21,w=21,plotfile)
# go_comparison_plot(clustergos[[1]]%>%filter(ontology=='MF')%>%{split(.,.$cluster)})+ggtitle('Molecular Function')
# dev.off()
# message(normalizePath(plotfile))
# plotfile<- here(paste0('plots/','cluster_go_cc','.pdf'));cairo_pdf(h=21,w=21,plotfile)
# go_comparison_plot(clustergos[[1]]%>%filter(ontology=='CC')%>%{split(.,.$cluster)})+ggtitle('Molecular Function')
# dev.off()
# message(normalizePath(plotfile))



