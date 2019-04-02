library(tidyverse)
library(magrittr)
library(data.table)
library(stringr)
library(magrittr)
library(splines)

get_limmafit_predvals <- function(limmafit,designmatrix){
  (limmafit$coefficients %*% t(limmafit$design))%>%
    set_colnames(designmatrix$dataset)%>%
    as.data.frame%>%
    rownames_to_column('gene_name')%>%
    gather(dataset,predicted_signal,-gene_name)%>%
    left_join(designmatrix)%>%
    distinct(gene_name,time,assay,.keep_all = TRUE)
}


get_limmafit_stdevs <- function(limmafit,designmatrix){
  (((limmafit$stdev.unscaled)^2) %*% t(limmafit$design))%>%
    set_colnames(designmatrix$dataset)%>%
    as.data.frame%>%
    rownames_to_column('gene_name')%>%
    gather(dataset,var_signal,-gene_name)%>%
    mutate(sd_signal = sqrt(var_signal))%>%
    left_join(designmatrix)%>%
    distinct(gene_name,time,assay,.keep_all = TRUE)
}

getlimmapredictions <- function(modfit,modname,designmatrix){
  predictedvals <- get_limmafit_predvals(modfit,designmatrix)
  predictedstdevs <- get_limmafit_stdevs(modfit,designmatrix)
  predictedvals%<>%left_join(predictedstdevs)
  predictedvals%<>%select(-matches('var_'))
  predictedvals%<>%filter(rep==1)
  predictedvals%<>%select(-rep)
  colnames(predictedvals) %<>%str_replace('predicted_signal',paste0('predicted_signal_',modname))
  colnames(predictedvals) %<>%str_replace('sd_signal',paste0('sd_signal_',modname))
  predictedvals
}



limmafits <- list()
designmatrix$time %<>% as_factor

#fit the full model
limmafits[['full']] = limma::lmFit(exprmatrix,
                                   design=model.matrix( ~ time*(ribo+MS), designmatrix)
)

#with 4 splines
limmafits[['spline_4']] = limma::lmFit(exprmatrix,
                                        design=model.matrix( ~ ns(as.numeric(time),4)*(ribo+MS), designmatrix)
)

#with fewer splines
limmafits[['spline_3']] = limma::lmFit(exprmatrix,
                                       design=model.matrix( ~ ns(as.numeric(time),3)*(ribo+MS), designmatrix)
)

#with fewer splines
limmafits[['spline_2']] = limma::lmFit(exprmatrix,
                                       design=model.matrix( ~ ns(as.numeric(time),2)*(ribo+MS), designmatrix)
)

###And the stepwise model
tps <- designmatrix$time%>%levels
designmatrix%<>%.[,colnames(.)%>%setdiff(tps)]
designmatrix %<>% cbind(as.numeric(designmatrix$time)%>%
  {sweep(replicate(length(.),1:length(tps)),STATS = .,MARGIN = 2,FUN='<=')}%>%
  t%>%set_colnames(designmatrix$time%>%levels))

limmafits[['full_step']] = limma::lmFit(exprmatrix,
                                        design=model.matrix( as.formula(paste0(' ~ (',paste0(tps[-1],collapse='+'),') * (ribo+MS)'))
, designmatrix)
)

##Now, let's do some model comparisons
###First let's confirm that the stepwise model gives similiar results to the jumps one
ribomsdesign<-designmatrix%>%filter(!assay=='total')

#get predictions of different models 
modelpreddf <- exprdata %>%select(gene_name,dataset,signal)%>%
  left_join(getlimmapredictions(limmafits[['full']],'full',designmatrix)) %>%
  left_join(getlimmapredictions(limmafits[['full_step']],'full_step',designmatrix)) %>%
  left_join(getlimmapredictions(limmafits[['spline_4']],'spline_4',designmatrix))%>%
  left_join(getlimmapredictions(limmafits[['spline_3']],'spline_3',designmatrix))%>%
  left_join(getlimmapredictions(limmafits[['spline_2']],'spline_2',designmatrix))

modelpreddf%>%colnames
modelpreddf$predicted_signal_full

##
modelpreddf%>%colnames

get_varexplained_df <- function(modelpreddf,fullpredname,redpredname){
  myformula = as.formula(paste0(fullpredname,' ~ ',redpredname))

    modelpreddf%>%group_by(gene_name)%>%
    nest%>%
    mutate(varexplained = map(.$data,~lm(data = ., myformula)%>%anova%>%.$"Sum Sq"%>%{.[1]/sum(.)} ))%>%
    identity%>%
    select(-data)%>%
    unnest
}
varexpldf<-get_varexplained_df(modelpreddf,fullpredname = 'predicted_signal_full',redpredname = 'predicted_signal_full_step')
##This confirms the stepwise model gives identical predictions to the non stepwise model

varexpldf$varexplained%>%table


##Now let's compare the full model to 3df spline - a 4 df spline is simply a perfect fit
varexpldf<- get_varexplained_df(modelpreddf,fullpredname = 'predicted_signal_full',redpredname = 'predicted_signal_spline_4')

varexpldf<- get_varexplained_df(modelpreddf,fullpredname = 'predicted_signal_full',redpredname = 'predicted_signal_spline_3')

##Now let's look at our PCA plots using the step wise

allcoefftbl<- limmafits[['full_step']]$coefficients%>%as.data.frame%>%
    rownames_to_column('gene_name')%>%
    select(-`(Intercept)`,-MSTRUE,-riboTRUE)
  

#perform pca

pcafit <- princomp(allcoefftbl%>%select(-gene_name))

#make plots of the pcas
svglite::svglite('../plots/pcafit_limmafcs.svg')
plot(pcafit,main='Fold Change Over Time - PCA')
plot(pcafit$scores[,1:2],ylim=pcafit$scores[,1:2]%>%range,xlim=pcafit$scores[,1:2]%>%range)
plot(pcafit$scores[,2:3],ylim=pcafit$scores[,2:3]%>%range,xlim=pcafit$scores[,2:3]%>%range)
dev.off()
getwd()%>%message

svg(h=5,w=8,'../plots/limmafc_pca_loadings.svg'%T>%{normalizePath(.)%>%message})
ggpubr::ggarrange(ncol=2,
                  pcafit$loading[,1]%>%enframe('dimension','loading')%>%
                    arrange(str_detect(dimension,'MS'),str_detect(dimension,'ribo'),str_detect(dimension,'P0'),str_detect(dimension,'E17'),str_detect(dimension,'E16'),str_detect(dimension,'E14'))%>%
                    mutate(dimension=factor(dimension,unique(dimension)))%>%
                    ggplot(aes(x=dimension,y=loading))+stat_identity(geom='bar')+
                    theme_minimal()+
                    theme(axis.text.x=element_text(angle=45,vjust=0.5))+
                    ggtitle('Developmental Fold Changes \n PCA1')
                  ,
                  pcafit$loading[,2]%>%enframe('dimension','loading')%>%
                    arrange(str_detect(dimension,'MS'),str_detect(dimension,'ribo'),str_detect(dimension,'P0'),str_detect(dimension,'E17'),str_detect(dimension,'E16'),str_detect(dimension,'E14'))%>%
                    mutate(dimension=factor(dimension,unique(dimension)))%>%
                    ggplot(aes(x=dimension,y=loading))+stat_identity(geom='bar')+
                    theme_minimal()+
                    theme(axis.text.x=element_text(angle=45,size=8,vjust=0.5))+
                    ggtitle('Developmental Fold Changes \n PCA2')
)

dev.off()


hmapfile <- file.path(paste0('heatmaps/',basename(transformdexprfile),'_limma_cs.pdf'))
hmapfile%>%dirname%>%dir.create

#
hmapdata <- as.matrix(allcoefftbl%>%select(-gene_name))%>%pmax(.,-3)%>%pmin(.,3)
#
cairo_pdf(normalizePath(hmapfile) %T>% message,w = 10, h = 10)
#
cbreaks <-
  cut_number(hmapdata,100,ordered=TRUE)%>%levels%>%str_split(',')%>%map(str_replace_all,'[^0-9\\-\\.e]','')%>%{c(.[[1]][[1]],map_chr(.,2))}%>%as.numeric
# 
par(lend = 1,rend=1)
par(mar=c(2, 2, 2, 2))
gplots::heatmap.2(hmapdata,
                  trace='none',
                  scale='none',
                  # breaks=cbreaks,
                  # scale="row",
                  # dendrogram="none",
                  # Colv = FALSE,
                  # Rowv = FALSE,
                  labRow=FALSE,
                  cexRow=1,
                  cexCol=0.75,srtCol=45,
                  main=paste0("Limma Fold Change Matrix \n "),
                  # RowSideColors=protcolors,
                  # ColSideColors=protcolors,
                  # key.xlab="Log2(norm_iBAQ ratio)"
)
# legend(x=0,y=0.85, legend=catcolors$pcat,fill=catcolors$color,cex=0.7)
dev.off()


library(dbscan)

coefhdb <- dbscan::hdbscan(allcoefftbl[,-1],10)
diffTEnames<-'/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/xtail/xtail_P0.txt'%>%read_tsv%>%filter(adj_p_value < 0.05)%>%
  select(gene_id=feature_id)%>%
  left_join(by='gene_id',read.table(header=T,'/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/ids.txt'))%>%
  .$gene_name

difftecoeffs <- allcoefftbl[,-1]%>%set_rownames(allcoefftbl[,1])%>%.[intersect(rownames(.),diffTEnames),]

pcafit <- princomp(difftecoeffs)

difftecoefhdb<-dbscan::hdbscan(difftecoeffs,10)


#make plots of the pcas
svglite::svglite('../plots/pcafit_limmafcs.svg')
plot(pcafit,main='Fold Change Over Time - PCA')
plot(pcafit$scores[,1:2],ylim=pcafit$scores[,1:2]%>%range,xlim=pcafit$scores[,1:2]%>%range)
plot(pcafit$scores[,2:3],ylim=pcafit$scores[,2:3]%>%range,xlim=pcafit$scores[,2:3]%>%range)
dev.off()
getwd()%>%message

#on only differential TE 
modelpreddf%>%left_join(data_frame(gene_name=rownames(difftecoeffs),cluster=difftecoefhdb$cluster))%>%
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
  theme_bw()+
  ggtitle(paste0('clusters - Differential TE only Genes\n Stepwise time factors\nClustering: HDBscan w/ min_cluster_size=10'))

#And on all coefficients


#on only differential TE 
modelpreddf%>%left_join(data_frame(gene_name= allcoefftbl$gene_name,cluster=coefhdb$cluster))%>%
  select(gene_name,time,assay,predicted_signal_full_step,cluster)%>%
  filter(!is.na(predicted_signal_full_step))%>%
  filter(!is.na(cluster))%>%
  group_by(cluster)%>%mutate(clustern=paste0('Cluster_',LETTERS[cluster+1],' n = ',n_distinct(gene_name)))%>%
  group_by(gene_name,assay)%>%
  mutate(predicted_signal_full_step = predicted_signal_full_step-mean(predicted_signal_full_step,na.rm=T))%>%
  ggplot(aes(x=as_factor(time),y=predicted_signal_full_step,group=gene_name,color=as.factor(clustern)))+
  #    geom_line(alpha=I(0.1))+
  geom_line(alpha=I(0.1))+
  facet_grid(assay~clustern)+
  coord_cartesian(ylim=c(-2,2))+
  scale_y_continuous('Centered log2-signal')+
  theme_bw()+
  ggtitle(paste0('clusters - All Genes\n Stepwise time factors\nClustering: HDBscan w/ min_cluster_size=10'))


data_frame(gene_name= difftecoeffs%>%rownames,cluster=difftecoefhdb$cluster)%>%filter(cluster=='1')%>%as.data.frame

'/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/xtail/xtail*.txt'%>%.[4]%>%Sys.glob()%>%map(fread)%>%bind_rows%>%filter(adj_p_value<0.05)%>%nrow


#TODO heatmaps with expr data not fold changes

#let's first just make some trajectory plots

#trajectory plot with raw data from heatamaps


dist(allcoefftbl[,-1],'cosine')
#coeff_deconstyle dbscan::hdbscan()



####Now, decon style, nearest neghbor clustering
#on only differential TE 
modelpreddf%>%left_join(data_frame(gene_name=rownames(difftecoeffs),cluster=difftecoefhdb$cluster))%>%
  select(gene_name,time,assay,predicted_signal_full_step,cluster)%>%
  filter(!is.na(predicted_signal_full_step))%>%
  filter(!is.na(cluster))%>%
  group_by(cluster)%>%mutate(clustern=paste0('Cluster_',LETTERS[cluster+1],' n = ',n_distinct(gene_name)))%>%
  group_by(gene_name,assay)%>%
  mutate(predicted_signal_full_step = predicted_signal_full_step-mean(predicted_signal_full_step,na.rm=T))%>%
  ggplot(aes(x=as_factor(time),y=predicted_signal_full_step,group=gene_name,color=as.factor(clustern)))+
  #    geom_line(alpha=I(0.1))+
  geom_line(alpha=I(0.1))+
  facet_grid(assay~clustern)+
  coord_cartesian(ylim=c(-2,2))+
  scale_y_continuous('Centered log2-signal')+
  theme_bw()+
  ggtitle(paste0('clusters - Differential TE only Genes\n Stepwise time factors\nClustering: HDBscan w/ min_cluster_size=10'))

X <- allcoefftbl[,-1]
cos.sim <- function(ix) 
{
  A = X[ix[1],]
  B = X[ix[2],]
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}   
n <- nrow(X) 
cmb <- expand.grid(i=1:n, j=1:n) 
C <- matrix(apply(cmb,1,cos.sim),n,n)

modelpreddf%>%head


gnames <- modelpreddf$gene_name%>%unique

clutomat<-modelpreddf%>%select(gene_name,dataset,time,assay,predicted_signal_full_step)%>%
    # filter(assay=='ribo')%>%
    filter(!is.na(predicted_signal_full_step))%>%
    group_by(gene_name,assay)%>%
    mutate(signal = predicted_signal_full_step - mean(na.rm=T,predicted_signal_full_step))%>%
    ungroup%>%select(dataset,signal,gene_name)%>%
    mutate(dataset=str_replace(dataset,'_\\d$',''))%>%
    spread(dataset,signal)%>%
    {set_rownames(as.matrix(.[,-1]),.$gene_name)}%>%
    as.data.frame




get_cluto_clusts<-function(clutomat,clustnum,
                           clutoinput= 'pipeline/clustering/cluto_clustering/tmp_clusterdata.txt'
){
  #create folder for cluto
  clutoinput%>%dirname%>%dir.create(rec=T,showW=F)
  clutoinput<-clutoinput%>%{file.path(normalizePath(dirname(.),mustWork = F),basename(.))}
  #put our table with a 'dims' line before it
  cat(file=clutoinput,x=str_interp('${nrow(clutomat)} ${ncol(clutomat)}\n'))
  as.data.frame(clutomat)%>%write_tsv(append = T, clutoinput,col_names = F)
  
  system('cd pipeline/clustering/cluto_clustering ;~/Downloads/cluto-2.1.2/Darwin-i386/vcluster ${clutoinput} ${clustnum} -agglofrom=50 -sim=cos -clmethod=graph -showtree -fulltree -showfeatures -showsummaries=cliques -labeltree'%>%str_interp%T>%message)
  Sys.sleep(2)
  clusteringfile <- clutoinput%>%{paste0(.,'.clustering.',clustnum)}
  stopifnot(file.exists(clusteringfile))
  clutoclusts <- clusteringfile%>%scan%>%setNames(rownames(clutomat))
  list(cluster=clutoclusts)
  
}

'pipeline/clustering/cluto_clustering/*'%>%Sys.glob%>%file.remove
clutoclustlist <- lapply(1:20%>%setNames(.,.),get_cluto_clusts,clutomat=clutomat,clutoinput='pipeline/clustering/testcluster_exp.txt')
clutoclusts <- clutoclustlist[['20']]$cluster
clutoclustlist%>%saveRDS('pipeline/clustering/cluto_clusters.rds')

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
source('src/R/go_term_funcs.R')

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




distgodistmatrix[1:10,1:10]

lapply(goenrichmentsbp,function(x) allterms %in% x)





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


save.image('data/cluto_clusts.Rdata')
