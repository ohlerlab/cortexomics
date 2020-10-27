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
                           bisectstr='large',
                           clutobin=CLUTOBIN
){
  #create folder for cluto
  qs(clutomat,'M+')
  stopifnot	(bisectstr %in% c('best','large','largess'))
  clutoinputsig = paste0(attr(clutomat,'obname'),substr(digest::digest(mtcars),1,10))
  clutoinputfile = file.path(clutoinputfold,paste0('clutodata_',clutoinputsig,'K_',clustnum,'_Bisect',bisectstr,'.tsv'))
  clutoinputfile<-clutoinputfile%>%{file.path(normalizePath(dirname(.),mustWork = F),basename(.))}
  Sys.glob(paste0(clutoinputfile,'*'))%>%file.remove
  #put our table with a 'dims' line before it
  clutoinputfile%>%dirname%>%dir.create(rec=T,showW=F)
  cat(file=clutoinputfile,x=str_interp('${nrow(clutomat)} ${ncol(clutomat)}\n'))
  as.data.frame(clutomat)%>%write_tsv(append = T, clutoinputfile,col_names = F)
  #  
  clusteringfile = paste0(clutoinputfile,'.clustering.',clustnum)
  if(clutoinputsig%>%str_detect('effect')){
  clutocmd = 'cd ${clutodir} &&  
  ${clutobin} -agglofrom=${agglofrom} -sim=dist ${clustfnpar} -cstype=${bisectstr} -clmethod=graph -showtree -fulltree -nfeatures=5  -labeltree -clustfile ${clusteringfile}  ${clutoinputfile} ${clustnum}'%>%
    str_interp%T>%message

  }else{
   clutocmd = 'cd ${clutodir} &&  
    ${clutobin} -plotformat svg -plottree=${clutoinputfile}.treeplot.svg -plotmatrix=${clutoinputfile}.matplot.svg -plotsmatrix=${clutoinputfile}.matplots.svg -plotcluster=${clutoinputfile}.clplot.svg -plotsclusters=${clutoinputfile}.clplots.svg -agglofrom=${agglofrom} -sim=cos ${clustfnpar} -cstype=${bisectstr} -clmethod=graph -showtree -fulltree ${featshow} -nfeatures=5 ${showsumm} -labeltree -clustfile ${clusteringfile}  ${clutoinputfile} ${clustnum}'%>%str_interp%T>%message

  }
  clutostdout = system(clutocmd,intern=TRUE)
  cat(clutostdout,file=str_interp('${clutoinputfile}.stdout'))
  stopifnot(file.exists(clusteringfile))
  clusteringfile <- clutoinputfile%>%{Sys.glob(paste0(.,'.clustering.*'))}
  clutoclusts <- clusteringfile%>%scan%>%setNames(rownames(clutomat))
  dataname = attr(clutomat,"obname")%>%qs('s')
  list(cluster=clutoclusts,
  	data = clutomat,
  	name = str_interp('Cluto_K${clustnum}_Bis_${bisectstr}_${dataname}')
  )

}

get_hdbscan_clusts <- function(clutomat,minppoints,
  ){
  #create folder for cluto
  qs(clutomat,'M+')
  if(clutoinputsig%>%str_detect('effect')){

  }else{
   clutocmd = 'cd ${clutodir} &&  
    ${clutobin} -plotformat svg -plottree=${clutoinputfile}.treeplot.svg -plotmatrix=${clutoinputfile}.matplot.svg -plotsmatrix=${clutoinputfile}.matplots.svg -plotcluster=${clutoinputfile}.clplot.svg -plotsclusters=${clutoinputfile}.clplots.svg -agglofrom=${agglofrom} -sim=cos ${clustfnpar} -cstype=${bisectstr} -clmethod=graph -showtree -fulltree ${featshow} -nfeatures=5 ${showsumm} -labeltree -clustfile ${clusteringfile}  ${clutoinputfile} ${clustnum}'%>%str_interp%T>%message

  }
  clutostdout = system(clutocmd,intern=TRUE)
  cat(clutostdout,file=str_interp('${clutoinputfile}.stdout'))
  stopifnot(file.exists(clusteringfile))
  clusteringfile <- clutoinputfile%>%{Sys.glob(paste0(.,'.clustering.*'))}
  clutoclusts <- clusteringfile%>%scan%>%setNames(rownames(clutomat))
  dataname = attr(clutomat,"obname")%>%qs('s')
  list(cluster=clutoclusts,
    data = clutomat,
    name = str_interp('Cluto_K${clustnum}_Bis_${bisectstr}_${dataname}')
  )

}



nomissingids <- metainfo%>%filter(n_stagemissing==0)%>%.$uprotein_id
upid2gnm = metainfo%>%filter(isbest)%>%distinct(uprotein_id,gene_name)%>%{safe_hashmap(.[[2]],.[[1]])}

voomeffdf<-timeff_ciddf%>%mutate(effect=paste0(time,'_',assay))%>%filter(time!='E13')%>%
  filter(uprotein_id%in%upid2gnm$keys())%>%
  select(uprotein_id,effect,logFC)%>%
  spread(effect,logFC)%>%
  {set_rownames(as.matrix(.[,-1]),.[[1]])}

voomnomissing = voomeffdf%>%{voomeffdf[intersect(rownames(.),nomissingids),]}

datamats=list()
datamats['norm_proDD'] = list(bestmscountvoom$E)
datamats['voom_time_effects'] = list(voomeffdf)
datamats['voom_time_effects_nomiss'] = list(voomnomissing)
for(i in seq_along(datamats)){rownames(datamats[[i]]) %<>% {upid2gnm[[.]]}}
clustlist = list()
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
	attr(datamat,'obname') = names(datamats)[[i]]
	kvals = c(5,10,15,20)
	clutoclustlist <- lapply(kvals%>%setNames(.,.),projmemoise(get_cluto_clusts),clutomat=datamat)
	# clutoclustlist_filt <- clutoclustlist%>%map(function(x){
	# 	x$cluster=x$cluster+1;x$cluster%<>%keep(`>`,0) ;x
	# })
  names(clutoclustlist) = map_chr(clutoclustlist,'name')
	clutoclustlist_filt<-clutoclustlist
	#enter into list of clusters to plot
  #mark the unclustered ones, and then have them start at 1
  for(i in seq_along(clutoclustlist_filt)) clutoclustlist_filt[[i]]$cluster%<>%{.[.== -1] = NA;.=.+1;.}
	clustlist[names(clutoclustlist_filt)] = clutoclustlist_filt

}
for(i in seq_along(clustlist)) stopifnot(clustlist[[i]][[1]]%>%unique%>%Negate(is_in)(.,0)%>%all)
clustlist%>%map('cluster')%>%map(table)
  stopifnot(clutoclustlist_filt[[1]][[1]]%>%unique%>%Negate(is_in)(.,0)%>%all)


