#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(assertthat))
message('...done')

defaultargs <- c(
  clusterdata='exprdata/limma_fold_changes.txt',
  outdir='clusters/limma_fold_changes/kmeans/',
  maxclustnum=10
)

args <- coalesce(
  commandArgs(trailingOnly=TRUE)%>%
      {
        aarg=grepl(.,pattern='^--args$')
        if(any(aarg)){
          .[-c(0:aarg)]
        }else{
          .
        }
    } %>%
    .[1:length(defaultargs)] %>%
    setNames(names(defaultargs)),
  defaultargs
)
for(i in names(args)) assign(i,args[i])

outdir%>%dir.create(rec=TRUE,showWarnings = FALSE)
#read in the cluster data
clusterdata %<>% read_tsv

clusterdata %<>% select(gene_name, everything())
assert_that(map_chr(clusterdata,class)[1] == 'character')
assert_that(all(map_chr(clusterdata,class)[-1] == 'numeric'))

#
notallmissing <- clusterdata[,-1]%>%as.matrix%>%apply(2,Negate(is.na))%>%apply(1,all)
clusterdata <- clusterdata[notallmissing,]

clustobs<-lapply(1:maxclustnum,kmeans,x=as.matrix(clusterdata[,-1]))

kmeans_ss <- tibble(ncluster = seq_along(clustobs),
                  ss = map_dbl(clustobs, function(x) x$betweenss / x$totss))
elbowplot<-kmeans_ss %>%
ggplot(aes(x = ncluster, y = ss)) +
geom_point() + geom_line() +
scale_x_continuous(breaks = unique(kmeans_ss$ncluster))+theme_bw()

ggsave(plot=elbowplot,file=paste0(outdir,'/elbowplot.pdf'))

#
gene_cluster_df<-map(clustobs,'cluster')%>%
	setNames(.,paste0('k',seq_along(.)))%>%
	as_data_frame%>%
	bind_cols(clusterdata[TRUE,'gene_name'],.)


gene_cluster_df%>%write_tsv(file.path(outdir,'perplex_cluster_table.txt'))


