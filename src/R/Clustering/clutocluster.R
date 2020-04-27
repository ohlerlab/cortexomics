library(here)
source(here::here('src/R/Rprofile.R'))

# suppressMessages({library(svglite)})
suppressMessages({library(readr)})
#suppressMessages({library(psd)})
suppressMessages({library(txtplot)})
suppressMessages({library(rtracklayer)})
suppressMessages({library(stringr)})
suppressMessages({library(data.table)})
suppressMessages({library(assertthat)})
suppressMessages({library(parallel)})
suppressMessages({library(dplyr)})
suppressMessages({library(purrr)})
suppressMessages({library(here)})
suppressMessages({library(magrittr)})
suppressMessages({library(stringr)})
suppressMessages({library(tidyverse)})
# suppressMessages({library(bamsignals)})
suppressMessages({library(zeallot)})
suppressMessages({library(stringr)})
suppressMessages({library(here)})
suppressMessages({library(assertthat)})
suppressMessages({library(multitaper)})
suppressMessages({library(here)})
suppressMessages({library(ggpubr)})
suppressMessages({library(conflicted)})

# suppressMessages({map(lsf.str("package:BiocGenerics"),.f=conflict_prefer,'BiocGenerics')})
conflict_prefer("intersect", "BiocGenerics")
reduce <- GenomicRanges::reduce
#
MAPQTHRESH <- 200
USEPHASE <- FALSE
USERIBOSEQC <- FALSE


CLUTOBIN <- here('Applications/cluto-2.1.2/Linux-x86_64/vcluster')

file.copy('pipeline/clustering/cluto_clustering/tmp_clusterdata.txt','ext_data/old_clutodata.txt')
oldclutomat <- 'ext_data/old_clutodata.txt'%>%fread%>%as.matrix
oldclutoclusts<-readRDS('/fast/groups/ag_ohler/work/dharnet_m/cortexomics/data/oldclustercluto.rds')

#Now the new data
metainfo<-read_tsv(here('pipeline/exprdata/limma_genemetadata.tsv'))
#
prediction_df<-read_tsv('tables/prediction_df.tsv')%>%safe_left_join(metainfo%>%distinct(uprotein_id,gene_name))
##Now let's compare the full model to 3df spline - a 4 df spline is simply a perfect fit
limmapredf <- prediction_df%>% select(time, assay, uprotein_id,gene_name, logFC, CI.L, CI.R, adj.P.Val)%>%
  filter(assay%in%c('ribo','total','MS'))
gnames <- limmapredf$gene_name%>%unique
#
#get the predictions for a set of assays, subtract out avarage by gene/assay
getassaypreds <- function(limmapredf,selassays){
  clutomat<- limmapredf %>% select(time, assay, gene_name, logFC)%>%
      filter(!is.na(logFC))%>%
      filter(assay %in% selassays)%>%
      group_by(gene_name,assay)%>%
      mutate(signal = logFC - mean(na.rm=T,logFC))%>%
      ungroup%>%
      mutate(dataset=paste0(time,assay))%>%
      select(dataset,signal,gene_name)%>%
      spread(dataset,signal)%>%
      {set_rownames(as.matrix(.[,-1]),.$gene_name)}%>%
      as.data.frame
    clutomat
}

get_cluto_clusts<-function(clutomat,clustnum,
                           clutoinput= 'pipeline/clustering/cluto_clustering/tmp_clusterdata.txt',
                           clutodir = 'pipeline/clustering/cluto_clustering/',
                           clustfnpar = '',
                           agglofrom=50,
                           clutobin=CLUTOBIN
){
  #create folder for cluto
  qs(clutomat,'M+')
  clutoinput%>%dirname%>%dir.create(rec=T,showW=F)
  clutoinput<-clutoinput%>%{file.path(normalizePath(dirname(.),mustWork = F),basename(.))}
  #put our table with a 'dims' line before it
  cat(file=clutoinput,x=str_interp('${nrow(clutomat)} ${ncol(clutomat)}\n'))
  as.data.frame(clutomat)%>%write_tsv(append = T, clutoinput,col_names = F)
  #  
  system('cd ${clutodir} ;${clutobin} ${clutoinput} ${clustnum} -agglofrom=${agglofrom} -sim=cos ${clustfnpar} -clmethod=graph -showtree -fulltree -showfeatures -showsummaries=cliques -labeltree'%>%str_interp%T>%message)
  Sys.sleep(2)
  clusteringfile <- clutoinput%>%{paste0(.,'.clustering.',clustnum)}
  stopifnot(file.exists(clusteringfile))
  clutoclusts <- clusteringfile%>%scan%>%setNames(rownames(clutomat))
  list(cluster=clutoclusts)
}




################################################################################
########First let's test with old clusters
################################################################################
# load('data/cluto_clusts.Rdata')
# oldclutomat<-clutomat%>%as.matrix
clutmat_RNA <- oldclutomat%>%.[,1:5]

clutoclustlist_RNA <- lapply(4%>%setNames(.,.),get_cluto_clusts,clutomat=oldclutomat,clutoinput='pipeline/clustering/testcluster_clutmat_RNAonly_exp.txt')

clutoclustlist_RNA <- lapply(30%>%setNames(.,.),get_cluto_clusts,clutomat=oldclutomat,clutoinput='pipeline/clustering/testcluster_clutmat_RNAonly_exp.txt',agglofrom=30)

clutoclustlist_RNA[[1]]%>%table

clutoclustlist_RNA <- lapply(5%>%setNames(.,.),get_cluto_clusts,clutomat=oldclutomat,clutoinput='pipeline/clustering/testcluster_clutmat_RNAonly_exp.txt',clustfnpar = ' -cstype=best ')
clutoclustlist_RNA[[1]]%>%table

clutmat_RNA <- oldclutomat%>%.[,colnames(.)%>%str_subset('total')]
clutmat_RNA_Ribo <- oldclutomat%>%.[,colnames(.)%>%str_subset('total|ribo')]
clutmat_RNA_MS <- oldclutomat%>%.[,colnames(.)%>%str_subset('total|MS')]
clutmat_RNA_Ribo_MS <- oldclutomat%>%.[,]
clutmat_Ribo_MS <- oldclutomat%>%.[,colnames(.)%>%str_subset('ribo|MS')]
clutoclustlist_RNA <- lapply(13%>%setNames(.,.),get_cluto_clusts,clutomat=clutmat_RNA,clutoinput='pipeline/clustering/testcluster_clutmat_RNAonly_exp.txt')
stop()
clutoclustlist_RNA_Ribo <- lapply(13%>%setNames(.,.),get_cluto_clusts,clutomat=clutmat_RNA_Ribo,clutoinput='pipeline/clustering/testcluster_clutmat_RNA_Ribo_exp.txt')
clutoclustlist_RNA_MS <- lapply(13%>%setNames(.,.),get_cluto_clusts,clutomat=clutmat_RNA_MS,clutoinput='pipeline/clustering/testcluster_clutmat_RNA_Ribo_exp.txt')
clutoclustlist_Ribo_MS <- lapply(13%>%setNames(.,.),get_cluto_clusts,clutomat=clutmat_Ribo_MS,clutoinput='pipeline/clustering/testcluster_clutmat_Ribo_MS_exp.txt')

clutoclustlist_RNA_Ribo_MS <- lapply(13%>%setNames(.,.),get_cluto_clusts,clutomat=clutmat_RNA_Ribo_MS,clutoinput='pipeline/clustering/testcluster_clutmat_RNA_Ribo_MS_exp.txt')
#is this new clustering similiar to the old?
oldnames <- oldclutoclusts[[1]][[1]]%>%names

newcl<-clutoclustlist_RNA_Ribo_MS[['13']]$cluster%>%setNames(oldnames)
oldcl<-oldclutoclusts[['13']][[1]]%>%setNames(oldnames)

confusiontable<-function(newcl,oldcl){
conftab<-enframe(newcl,'gene','newcl')%>%left_join(oldcl%>%enframe('gene','oldcl'))%>%{table(.$newcl,.$oldcl)}%>%
  {.[order(rowSums(.)),]}%>%
  {.[,order(colSums(.))]}
conftab  
}

confusiontable(newcl,oldcl)

conftabold

conftab->conftabold
#Yes, acceptably so.


#Does it make a differeence, to take away the MS?
confusiontable(clutoclustlist_RNA_Ribo[['13']]$cluster,clutoclustlist_RNA_Ribo_MS[['13']]$cluster)
#Yes (no suprises)

#What if we then take away the Ribo?
confusiontable(clutoclustlist_RNA_Ribo[['13']]$cluster,clutoclustlist_RNA[['13']]$cluster)
#That does too.

#Perhaps not the best measure...
confusiontable(clutoclustlist_RNA[['13']]$cluster,clutoclustlist_Ribo[['13']]$cluster)

#This still doesn't really satisfy me - I'm not conviced these clusters are in any way stable....

#okay so what if I compare distance on the Ribo and mRNA levels
assayi='total'
cl1=1
cl2=13
assaydists<-map_df(.id='assay',c('total','ribo','MS')%>%setNames(.,.),function(assayi){
  map_df(.id='cl1',0:12%>%setNames(.,.),function(cl1){
    map_df(.id='cl2',0:12%>%setNames(.,.),function(cl2){
      colset <- colnames(clutomat)%>%str_subset(assayi)
      clustmean1 <- clutomat[clutoclustlist_RNA_Ribo_MS[['13']]$cluster==cl1,colset]%>%colMeans
      clustmean2 <- clutomat[clutoclustlist_RNA_Ribo_MS[['13']]$cluster==cl2,colset]%>%colMeans
      data_frame(dist=sum(abs(clustmean1 - clustmean2)^2))
    })
  })
})

assaydists%<>%spread(cl2,dist)

cor(
  assaydists%>%filter(assay=='total')%>%.[,-c(1:2)]%>%unlist,
  assaydists%>%filter(assay=='ribo')%>%.[,-c(1:2)]%>%unlist
)

cor(
  assaydists%>%filter(assay=='total')%>%.[,-c(1:2)]%>%unlist,
  assaydists%>%filter(assay=='MS')%>%.[,-c(1:2)]%>%unlist
)

cor(
  assaydists%>%filter(assay=='ribo')%>%.[,-c(1:2)]%>%unlist,
  assaydists%>%filter(assay=='MS')%>%.[,-c(1:2)]%>%unlist
)



txtplot(
  assaydists%>%filter(assay=='total')%>%.[,-c(1:2)]%>%unlist,
  assaydists%>%filter(assay=='ribo')%>%.[,-c(1:2)]%>%unlist
)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)

set.seed(123)

mat = cbind(rbind(matrix(rnorm(16, -1), 4), matrix(rnorm(32, 1), 8)),
            rbind(matrix(rnorm(24, 1), 4), matrix(rnorm(48, -1), 8)))

rownames(mat) = paste0("R", 1:12)
colnames(mat) = paste0("C", 1:10)

tothmp<-Heatmap(column_title='Distance Between Cluto Clusters - RNA',assaydists%>%filter(assay=='total')%>%.[,-c(1:2)]%>%as.matrix%>%{(.-sum(.))/sd(.)},row_order=0:12,column_order=0:12,cluster_rows=F,cluster_columns=F)
ribohmp<-Heatmap(column_title='Distance Between Cluto Clusters - Ribo',assaydists%>%filter(assay=='ribo')%>%.[,-c(1:2)]%>%as.matrix%>%{(.-sum(.))/sd(.)},row_order=0:12,column_order=0:12,cluster_rows=F,cluster_columns=F)
tothmp+ribohmp

Q

clutoclustlist_RNA_Ribo_MS[['13']]$cluster



#If I cluster iwth the new DATA, but on the same gene set, does it look the same,

#If I cluster with all the new data, does it still lookg the same



################################################################################
########Now try cluto on new data
################################################################################
  
n_clutmat_RNAonly <- getassaypreds(limmapredf,'total')
n_clutmat_RNA_Ribo <- getassaypreds(limmapredf,c('total','ribo'))
n_clutmat_RNA_Ribo_MS <- getassaypreds(limmapredf,c('total','ribo','MS'))
n_clutmat_Ribo_MS <- getassaypreds(limmapredf,c('ribo','MS'))


clutoclustlist <- lapply(1:2%>%setNames(.,.),get_cluto_clusts,clutomat=clutmat_RNAonly,clutoinput='pipeline/clustering/testcluster_clutmat_RNAonly_exp.txt')
clutoclustlist <- lapply(1:2%>%setNames(.,.),get_cluto_clusts,clutomat=clutmat_RNA_Ribo,clutoinput='pipeline/clustering/testcluster_clutmat_RNA_Ribo_exp.txt')
clutoclustlist <- lapply(1:2%>%setNames(.,.),get_cluto_clusts,clutomat=clutmat_RNA_Ribo_MS,clutoinput='pipeline/clustering/testcluster_clutmat_RNA_Ribo_MS_exp.txt')
clutoclustlist <- lapply(1:2%>%setNames(.,.),get_cluto_clusts,clutomat=clutmat_Ribo_MS,clutoinput='pipeline/clustering/testcluster_clutmat_Ribo_MS_exp.txt')

clutoclusts <- clutoclustlist[['20']]$cluster


# clutoclustlist%>%saveRDS('pipeline/clustering/cluto_clusters.rds')
clutoclustlist <- readRDS('pipeline/clustering/cluto_clusters.rds')
conflict_prefer('rowMedians','Biobase')

aggnesclust<-cluster::agnes(clutomat)

# 
# BiocManager::install('factoextra')
# BiocManager::install('NbClust')

#First let's excludefrom the matrix and our clusterlist, the -1 guys, who have no cluster
noclustnames<-clutoclustlist[[4]]$cluster%>%`<`(0)%>%which
clutomat_filt <- clutomat[setdiff(gnames,noclustnames),]
clutoclustlist_filt <- clutoclustlist%>%map(function(x){x$cluster=x$cluster+1;x$cluster%<>%keep(`>`,0) ;x})
# Elbow method
cosine_sim<-function(DF){
  Matrix <- as.matrix(DF)
  sim <- Matrix / sqrt(rowSums(Matrix * Matrix))
  sim <- sim %*% t(sim)
  sim
}

#' ### Assessing the Number of clusters
#' There are a bunhc
library(factoextra)

nbclustplot2<-factoextra::fviz_nbclust(as.data.frame(clutomat_filt), function(mat,n) clutoclustlist_filt[[n]], method = "wss",k.max = 20) +
  labs(subtitle = "Elbow method - processed data")
nbclustplot2

nbclustplot2_dist<-fviz_nbclust(cosine_sim(clutomat_filt), function(mat,n) clutoclustlist_filt[[n]], method = "wss",k.max = 20) +
  labs(subtitle = "Elbow method - cosine distances")
nbclustplot2_dist

nbclustplot_gap<-fviz_nbclust(as.data.frame(clutomat_filt), function(mat,n) clutoclustlist_filt[[n]],  method = "gap_stat", nboot = 200, k.max=20) +
  labs(subtitle = "Gap_stat Method")

nbclustplot_gap_dist<-fviz_nbclust(cosine_sim(clutomat_filt), function(mat,n) clutoclustlist_filt[[n]],  method = "gap_stat", nboot = 200, k.max=20) +
  labs(subtitle = "Gap_stat Method - cosine distances")
nbclustplot_gap_dist



####Now, decon style, nearest neghbor clustering
#on only differential TE 


make_cluster_trajplots<-function(clutoclusts,indata=modelpreddf){
  indata%>%left_join(data_frame(gene_name=names(clutoclusts),cluster=clutoclusts))%>%
    select(gene_name,time,assay,predicted_signal_full_step,cluster)%>%
    filter(!is.na(predicted_signal_full_step))%>%
    filter(!is.na(cluster))%>%
    group_by(cluster)%>%mutate(clustern=paste0('Cluster_',LETTERS[cluster],' n = ',n_distinct(gene_name)))%>%
    group_by(gene_name,assay)%>%
    mutate(predicted_signal_full_step = predicted_signal_full_step-mean(predicted_signal_full_step,na.rm=T))%>%
    ggplot(aes(x=as_factor(time),y=predicted_signal_full_step,group=gene_name,color=as.factor(clustern)))+
    #    geom_line(alpha=I(0.1))+
    geom_line(alpha=I(0.1))+
    facet_grid(assay~clustern)+
    coord_cartesian(ylim=c(-2,2))+
    scale_y_continuous('Centered log2-signal')+
    stat_summary(aes(x=as_factor(time),group=clustern,color=as.factor(clustern)),alpha=I(1),color=I('black'),fun.y=median, fun.ymin=median, fun.ymax = median, geom='line',linetype=2)+
    theme_bw()+
    ggtitle(paste0('Cluto Clusters'))
}


make_cluster_trajplots(clutoclustlist_filt[[5]]$cluster)
make_cluster_trajplots(clutoclustlist_filt[[6]]$cluster)
make_cluster_trajplots(clutoclustlist_filt[[9]]$cluster)
make_cluster_trajplots(clutoclustlist_filt[[10]]$cluster)
make_cluster_trajplots(clutoclustlist_filt[[11]]$cluster)
make_cluster_trajplots(clutoclustlist_filt[[12]]$cluster)

cairo_pdf(w=16,h=12,'./plots/clustering/trajplots_cluto_13.pdf'%T>%{normalizePath(.)%>%message})
make_cluster_trajplots(clutoclustlist_filt[[13]]$cluster)
dev.off()


make_cluster_trajplots(kmeans$cluster)
# make_cluster_trajplots((myhclusts,13))
source('src/R/Functions/go_term_funcs.R')

#' ### Now run a go enrichment on each of the clusters, with the set as a whole as the background.
goenrichmentsbp<-clutoclustlist_filt[[13]]$cluster%>%split(.,.)%>%
  map(~rungo((names(.)),GTOGO%>%mutate(ensembl_gene_id=gene_name)%>%filter(ensembl_gene_id %in%gnames),'BP'))
goenrichmentsmf<-clutoclustlist_filt[[13]]$cluster%>%split(.,.)%>%
  map(~rungo((names(.)),GTOGO%>%mutate(ensembl_gene_id=gene_name)%>%filter(ensembl_gene_id %in%gnames),'MF'))
goenrichmentscc<-clutoclustlist_filt[[13]]$cluster%>%split(.,.)%>%

  map(~rungo((names(.)),GTOGO%>%mutate(ensembl_gene_id=gene_name)%>%filter(ensembl_gene_id %in%gnames),'CC'))

#get all go terms
allterms <- goenrichmentsbp%>%map('Term')%>%unlist%>%unique

#now get the similiarity matrix in go terms

source("src/R/Functions/go_term_funcs.R")

gomat <- Matrix::Matrix(data=FALSE,sparse=TRUE,
  nrow=n_distinct(GTOGO$ensembl_gene_id),
  ncol=n_distinct(GTOGO$go_id)
)

rownames(gomat)<-unique(GTOGO$ensembl_gene_id)
colnames(gomat)<-unique(GTOGO$go_id)
indmat<-data.frame(GTOGO$ensembl_gene_id%>%as.factor%>%as.numeric,GTOGO$go_id%>%as.factor%>%as.numeric)%>%as.matrix
gomat[indmat]<-TRUE
#discard terms that are always zero in our data
gomat <- gomat[,Matrix::colSums(gomat[metainfo$gene_id,]) > 10]


#we weant to go from a gene x term matrix to a cluster x term matrix

str(gomat[match(clustgeneids[[2]],rownames(gomat)),TRUE])

#so we could get the entropy of each cluster's distribution over the terms, giving us n entropies
#or the entropy of each term's distribution over the clusters
cnum_termentrops<-lapply(oldclutoclusts,function(clustvect){
  clustgeneids <- clustvect$cluster%>%
    enframe('gene_name','cluster')%>%
    left_join(metainfo%>%
    distinct(gene_name,gene_id))%>%
    filter(!is.na(gene_id))%>%
    {split(.$gene_id,.$cluster)}

  term_clust_mat <- lapply(clustgeneids,function(gids){
    Matrix::colSums(gomat[match(gids,rownames(gomat)),TRUE])
  })%>%simplify2array

  term_entropies <- term_clust_mat%>%apply(1,function(x){y=x[x>0];y=y/sum(y);sum(y * log(y))})

  term_entropies

})

pdf('tmp.pdf')
cnum_termentrops%>%setNames(paste0('n_',seq_along(.)))%>%enframe('ncluster','term_entropy')%>%
  mutate(ncluster = factor(ncluster,unique(ncluster)))%>%
  unnest%>%
  qplot(data=.,x=term_entropy,color=ncluster,geom='density')+facet_grid(ncluster~.,scale='free')+theme_bw()
dev.off()
normalizePath('tmp.pdf')

nrow(gomat)


library(gridExtra)

#Crete image table with top3 GO terms
top3clustergos<-goenrichmentsbp%>%map(~arrange(.,elimFisher))%>%map(head,5)%>%map(~select(.,Term,elimFisher))

library(here)

ont='BP'
clustnum=13
gonum=3
method='CLUTO'
topgo3plot<-here('plots/clustering/topgoterms/',paste0(ont,clustnum,method,gonum,'.topgos.pdf'))

topgo3plot%>%dirname%>%dir.create(rec=T)


texttheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 0.5)),
    colhead = list(fg_params=list(cex = 0.5)),
    rowhead = list(fg_params=list(cex = 0.5)))
#
pdf(topgo3plot,w=16,h=12)
# --- Graph 1 : If you want ONLY the table in your image :
# First I create an empty graph with absolutely nothing :
plist<-lapply(1:clustnum, function(n){
  qplot(1:10, 1:10, geom = "blank",main=paste0(method,' cluster number: ',n)) + 
  theme_bw() + theme(line = element_blank(), text = element_blank()) +
# Then I add my table :
  annotation_custom(grob = tableGrob(top3clustergos[[n]],rows=NULL,theme=texttheme))
})
ggpubr::annotate_figure(ggpubr::ggarrange(plotlist=plist),top=paste0('Top GO ',ont,' catagories ',method,' clusters '))
dev.off()


##Now look at enrichment fo rTE genes
alltegenes <- Sys.glob(here('pipeline/xtail/xtail_*.txt'))%>%map(.%>%
  fread%>%filter(adj_p_value<0.05)%>%left_join(fread('pipeline/ids.txt'),by=c(feature_id='gene_id'))%>%.$gene_name
)%>%unlist
clusttenumplot<-here('plots/clustering/topgoterms/',paste0(ont,clustnum,method,gonum,'.TEprops.pdf'))
#
pdf(clusttenumplot,w=16,h=12)
# --- Graph 1 : If you want ONLY the table in your image :
# First I create an empty graph with absolutely nothing :
n=1
plist<-lapply(1:clustnum, function(n){
  clutoclustlist_filt[[13]]$cluster%>%keep(.==n)%>%names %>% is_in(alltegenes)%>%table%>%.[order(names(.))]%>%as.numeric%>%tibble(n=.)%>%mutate(set=c('static TE','changing TE'))%>%
    qplot(data=.,fill=set,y=n,x='')+geom_bar(stat='identity')+theme_bw()+coord_polar("y", start=0)+
    ggtitle('Numbers of TE changing genes')
})
plotlabels=LETTERS[(1:13)+1]
ggpubr::annotate_figure(ggpubr::ggarrange(plotlist=plist,labels=plotlabels),top=paste0('TE proportions ',method,' clusters '))
dev.off()
normalizePath(clusttenumplot)





#' #### Check - cluto with raw values
#' The above clustering seems to work well - I posit this is because the
#' mean centering and cosine similiarity procedure subtracts out the magnitude of change over time and 
#' allows the algorithm to cluster things which are under distinct regulatory influences.
#' As a check, I now try clustering with less processed data. First - raw FPKMs

clutomatraw<-2^clutomat

#now cluster.
clutoclustlistraw <- lapply(1:30%>%setNames(.,.),get_cluto_clusts,clutomat=clutomatraw)
clutoclustsraw <- clutoclustlistraw[['20']]$cluster

noclustnames<-clutoclustlistraw[[4]]$cluster%>%`<`(0)%>%which
clutomat_filtraw <- clutomat[setdiff(gnames,noclustnames),]
clutoclustlist_filtraw <- clutoclustlistraw%>%map(function(x){x$cluster=x$cluster+1;x$cluster%<>%keep(`>`,0) ;x})


#and do elbow plot
nbclustplot2_raw<-fviz_nbclust((clutomat_filtraw), function(mat,n) clutoclustlist_filtraw[[n]], method = "wss",k.max = 30) +
  labs(subtitle = "Elbow method - Raw data")

nbclustplot2_raw

nbclustplot2_distraw<-fviz_nbclust(cosine_sim(clutomat_filtraw), function(mat,n) clutoclustlist_filtraw[[n]], method = "wss",k.max = 20) +
  labs(subtitle = "Elbow method - cosine distances, Raw data")

install.packages("clustree")
library(clustree)

clutoclustlist_filt%>%map(.%>%.$cluster%>%enframe('gene_name','cluster'))%>%bind_rows(.id='K')%>%mutate(K=paste0('K',K))%>%spread(K,cluster)%>%
    clustree('K')

warnings()

#and display the clusters as trajectories.
make_cluster_trajplots(clutoclustlistraw[[13]]$cluster)

clutoclustlistraw[[13]]%>%table%>%sort
clutoclustlist[[13]]%>%table%>%sort

#' So it looks like clustering on more 'raw' values works too - this disproves my hypothesis that log-fold changes centering was doing the work.


#' #### Check - cluto with coefficients
#' The above clustering seems to work well - I posit this is because the

clutomatraw<-2^clutomat

#now cluster.
clutoclustlistraw <- lapply(1:20%>%setNames(.,.),get_cluto_clusts,clutomat=clutomatraw)
clutoclustsraw <- clutoclustlistraw[['20']]$cluster

noclustnames<-clutoclustlistraw[[4]]$cluster%>%`<`(0)%>%which
clutomat_filtraw <- clutomat[setdiff(gnames,noclustnames),]
clutoclustlist_filtraw <- clutoclustlistraw%>%map(function(x){x$cluster=x$cluster+1;x$cluster%<>%keep(`>`,0) ;x})


#and do elbow plot
nbclustplot2_raw<-fviz_nbclust((clutomat_filtraw), function(mat,n) clutoclustlist_filtraw[[n]], method = "wss",k.max = 20) +
  labs(subtitle = "Elbow method - Raw data")

nbclustplot2_raw

nbclustplot2_distraw<-fviz_nbclust(cosine_sim(clutomat_filtraw), function(mat,n) clutoclustlist_filtraw[[n]], method = "wss",k.max = 20) +
  labs(subtitle = "Elbow method - cosine distances, Raw data")

#
make_cluster_trajplots(clutoclustlistraw[[13]]$cluster)

clutoclustlistraw[[13]]%>%table%>%sort
clutoclustlist[[13]]%>%table%>%sort


load('data/cluto_clusts.Rdata')

