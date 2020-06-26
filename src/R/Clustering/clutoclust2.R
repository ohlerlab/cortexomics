library(here)
library(tidyverse)
library(ggplot2)
library(magrittr)
if(!exists('mscountvoom'))load('data/integrate_exprdata2.Rdata')
isdarwin = grepl(x=system(intern=T,'uname -a'),pattern='Darwin')
os = if(isdarwin) 'Darwin-i386' else 'Linux-x86_64' 
CLUTOBIN <- here(str_interp('Applications/cluto-2.1.2/${os}/vcluster'))

get_cluto_clusts<-function(clutomat,clustnum,
                           clutoinputfold= 'pipeline/clustering/cluto_clustering/',
                           clutodir = 'pipeline/clustering/cluto_clustering/',
                           clustfnpar = '',
                           agglofrom=50,
                           bisectstr='best',
                           clutobin=CLUTOBIN
){
  #create folder for cluto
  qs(clutomat,'M+')
  stopifnot	(bisectstr %in% c('best','large','largess'))
  clutoinputsig = paste0(attr(clutomat,'obname'),substr(digest::digest(mtcars),1,10))
  clutoinputfile = file.path(clutoinputfold,paste0('clutodata_',clutoinputsig,'.tsv'))
  clutoinputfile<-clutoinputfile%>%{file.path(normalizePath(dirname(.),mustWork = F),basename(.))}
  #put our table with a 'dims' line before it
  clutoinputfile%>%dirname%>%dir.create(rec=T,showW=F)
  cat(file=clutoinputfile,x=str_interp('${nrow(clutomat)} ${ncol(clutomat)}\n'))
  as.data.frame(clutomat)%>%write_tsv(append = T, clutoinputfile,col_names = F)
  #  
  system('cd ${clutodir} ;${clutobin} ${clutoinputfile} ${clustnum} -agglofrom=${agglofrom} -sim=cos ${clustfnpar} -cstype=${bisectstr} -clmethod=graph -showtree -fulltree -showfeatures -showsummaries=cliques -labeltree'%>%str_interp%T>%message)
  Sys.sleep(2)
  clusteringfile <- clutoinputfile%>%{paste0(.,'.clustering.',clustnum)}
  stopifnot(file.exists(clusteringfile))
  clutoclusts <- clusteringfile%>%scan%>%setNames(rownames(clutomat))
  list(cluster=clutoclusts,
  	data = clutomat,
  	name = str_interp('Cluto_K${clustnum}_Bis_${bisectstr}_${attr(clutomat,"obname")}')
  )
}

datamats=list()
datamats['norm_raw'] = list(mscountvoom$E)

# #datamats['proDApred'] =

# exprmat <- datamats[['norm_raw']]
# sample_info_df <- data.frame(name = colnames(exprmat),
#                              stringsAsFactors = FALSE)
# sample_info_df$time <- sample_info_df$name%>%str_extract('[^_]+')
# sample_info_df$assay <- sample_info_df$name%>%str_extract('(?<=_)[^_]+')

# proDAfit <- proDA(exprmat[1:100,], design = ~ condition, 
#              col_data = sample_info_df, reference_level = "E13")


# attr(proDAmat,'obname') <- proDAdata



###create CLUTO clusters
for(i in seq_along(datamats)){
	datamat=datamats[[i]]
	attr(datamats,'obname') = names(datamats)[[i]]
	kvals = 5:14
	clutoclustlist <- lapply(kvals%>%setNames(.,.),get_cluto_clusts,clutomat=datamat)
	# clutoclustlist_filt <- clutoclustlist%>%map(function(x){
	# 	x$cluster=x$cluster+1;x$cluster%<>%keep(`>`,0) ;x
	# })
	clutoclustlist_filt<-clutoclustlist
	names(clutoclustlist) = map_chr(clutoclustlist,'name')
	#enter into list of clusters to plot
	clustlist[names(clutoclustlist_filt)] = clutoclustlist_filt
}


make_cluster_trajplots<-function(clusts){
	indata = clusts$data %||% stop()
	plotname = clusts$name %||% stop()

  clutplot = indata%>%left_join(data_frame(gene_name=names(clusts),cluster=clusts))%>%
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
    ggtitle(paste0(plotname))

    here('plots/clusters/')%>%dir.create(showWarn=F)
    plotfile = here(paste0('plots/clusters/',plotname,'.pdf'))

	cairo_pdf(w=16,h=12,plotfile%T>%{normalizePath(.)%>%message})
	print(clutplot)
	dev.off()

}

lapply(clustlist,make_cluster_trajplots)
lapply(clustlist,make_cluster_goplots)

clustlist%>%saveRDS('data/clustlist.rds')

#mark satb2