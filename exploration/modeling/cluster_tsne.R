#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(limma))
suppressMessages(library(DESeq2))
message('...done')

defaultargs <- c(
  clusterdatafile='exprdata/limma_fold_changes_ribodfilt.txt',
  outdir='clusters/limma_fold_changes_ribodfilt/tsne/',
  maxclustnum=10
)

PERPLEX = 50
PERPLEXITIES = c(PERPLEX/10,PERPLEX,PERPLEX*10)

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

message(capture.output(print(args)))
#args <- 'Rscript ../exploration/modeling/cluster_tsne.R exprdata/cent_scaled_exprdata.txt clusters/cent_scaled_exprdata/tsne/'%>%str_split(' ')%>%.[[1]]%>%.[-1]

for(i in names(args)) assign(i,args[i])
PERPLEX%<>%as.numeric
maxclustnum%<>%as.numeric


dir.create(outdir,showWarn=F,rec=TRUE)
#read in the cluster datafg

clusterdata <- read_tsv(clusterdatafile)


clusterdata %<>% select(gene_name, everything())
assert_that(map_chr(clusterdata,class)[1] == 'character')
assert_that(all(map_chr(clusterdata,class)[-1] == 'numeric'))

# clusterdatamat <- clusterdata  { set_rownames(as.matrix(.[,-1]),.[[1]]) }

#and cluster with kmeans
notallmissing <- clusterdata[,-1]%>%as.matrix%>%apply(2,Negate(is.na))%>%apply(1,all)
clusterdata <- clusterdata[notallmissing,]
for(i in seq_along(PERPLEXITIES)){

  tsnedimfile <- paste0(outdir,'/tsne_dims_p',PERPLEXITIES[[i]],'.txt')
  
  if(!file.exists(tsnedimfile)){
    tsne_dims = tsne::tsne(clusterdata[,-1], perplexity=PERPLEXITIES[[i]],max_iter=300,epoch=50)
    tsne_dims<-cbind(clusterdata[[1]],tsne_dims)
    write_tsv(as.data.frame(tsne_dims),tsnedimfile)
  }else{
    tsne_dims <- read_tsv(tsnedimfile)
  }

  clustob<-lapply(1:maxclustnum,kmeans,x=tsne_dims[,-1])

  kmeans_ss <- tibble(ncluster = seq_along(clustob),
                      ss = map_dbl(clustob, function(x) x$betweenss / x$totss))
  elbowplot<-kmeans_ss %>%
    ggplot(aes(x = ncluster, y = ss)) +
    geom_point() + geom_line() +
    scale_x_continuous(breaks = unique(kmeans_ss$ncluster))+theme_bw()

  ggsave(plot=elbowplot,file=paste0(outdir,'/perplexity_',PERPLEXITIES[[i]],'_elbowplot.pdf'))

  #
  gene_cluster_df<-map(clustob,'cluster')%>%
   	setNames(.,paste0('k',seq_along(.)))%>%
   	as_data_frame%>%
   	bind_cols(clusterdata[TRUE,'gene_name'],.)

  for(k in seq_len(maxclustnum)){
  	TSNE_kmean_cluster <- as.factor(clustob[[k]]$cluster)
  	tsneplot <- qplot(x=tsne_dims[[2]],y=tsne_dims[[3]],color=TSNE_kmean_cluster)+geom_point()+
  		ggtitle(paste0('Tsne Clusters, k=: ',k,'\nInput Data:',clusterdatafile),sub=paste0('perplexity = ',PERPLEX))

  	pdf(paste0(outdir,'/','clusters_',k,'perp',PERPLEXITIES[[i]],'.pdf')%>%normalizePath%T>%message)
  	print(tsneplot);dev.off()
  	# ggsave(plot=tsneplot,file=paste0(outdir,'/','clusters_',k,'.pdf')%>%normalizePath%T>%message)
  }

  gene_cluster_df%>%write_tsv(file.path(outdir,paste0('perplex_',PERPLEXITIES[[i]],'cluster_table.txt')))
}