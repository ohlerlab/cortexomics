library(here)
################################################################################
################################################################################
base::source(here::here('src/Rprofile.R'))

trid2gid<-load_hashmap('trid2gid.hmp')
gid2gnm<-load_hashmap('gid2gnm.hmp')
gnm2gid<-load_hashmap('gnm2gid.hmp')
trid2gnm<-load_hashmap('trid2gnm.hmp')
gnm2gid<-load_hashmap('gnm2gid.hmp')
gnm2trid<-load_hashmap('gnm2trid.hmp')
library(tidyverse)
library(ggplot2)
library(magrittr)
# if(!exists('mscountvoom'))load('data/integrate_exprdata2.Rdata')
isdarwin = grepl(x=system(intern=T,'uname -a'),pattern='Darwin')
os = if(isdarwin) 'Darwin-i386' else 'Linux-x86_64' 
CLUTOBIN <- here(str_interp('Applications/cluto-2.1.2/${os}/vcluster'))

featshow=''
showsumm=''

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

get_hdbscan_clusts <- function(clutomat,minppoints){
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



# nomissingids <- metainfo%>%filter(n_stagemissing==0)%>%.$uprotein_id
# upid2gnm = metainfo%>%filter(isbest)%>%distinct(uprotein_id,gene_name)%>%{safe_hashmap(.[[2]],.[[1]])}

# voomeffdf<-timeff_ciddf%>%mutate(effect=paste0(time,'_',assay))%>%filter(time!='E13')%>%
#   filter(uprotein_id%in%upid2gnm$keys())%>%
#   select(uprotein_id,effect,logFC)%>%
#   spread(effect,logFC)%>%
#   {set_rownames(as.matrix(.[,-1]),.[[1]])}

# voomnomissing = voomeffdf%>%{voomeffdf[intersect(rownames(.),nomissingids),]}




sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
sel_ms_mat<-readRDS('data/sel_ms_mat.rds')
countpred_df<-readRDS('data/countpred_df.rds')
tx_countdata<-readRDS('data/tx_countdata.rds')
allvoom <- readRDS(here('data/allvoom.rds'))
prediction_df = bind_rows(
  sel_prodpreds%>%
    mutate(assay='MS')%>%
    select(gene_id,time,assay,estimate,CI.L,CI.R),
  countpred_df%>%
    separate(contrast,c('time','assay'))%>%
    select(gene_id,time,assay,estimate=logFC,CI.L,CI.R)
)

if(!file.exists(here('data/sel_ms_impute.rds'))){
  sel_ms_impute <- sel_ms_mat[,]%>%
    as.data.frame%>%
    rownames_to_column('gene_id')%>%
    filter(!is.na(gene_id))%>%
    gather(dataset,signal,-gene_id)%>%
    separate(dataset,into=c('time','assay','replicate'))%>%
    safe_left_join(prediction_df%>%select(gene_id,time,assay,estimate,CI.L,CI.R))%>%
    mutate(signal = ifelse(is.na(signal),rnorm(length(signal),estimate,(CI.R-CI.L/3.92)),signal))%>%
    select(-estimate,-CI.R,-CI.L)%>%
    unite(dataset,time,assay,replicate)%>%
    spread(dataset,signal)
  saveRDS(sel_ms_impute,here('data/sel_ms_impute.rds'))
}else{
  sel_ms_impute<-readRDS(here('data/sel_ms_impute.rds'))
}

countcontr_df <- readRDS(here('data/countcontr_df.rds'))
stepcountcontrdf <- readRDS(here('data/stepcountcontrdf.rds'))
contrdf<-readRDS('data/contrdf.rds')
stepstepprotcontrdf<-readRDS('data/stepcontrdf.rds')

#Get expression levels, with imputations
exprdf = inner_join( 
  tx_countdata$abundance%>%
    as.data.frame%>%
    rownames_to_column('gene_id'),
  sel_ms_impute[,]
)
exprmat = as.matrix(exprdf[,-1])%>%set_rownames(exprdf[[1]])
rownames(exprmat) = gid2gnm[[rownames(exprmat)]]
{
#get t0 contrasts
t0_contrastdf = bind_rows(countcontr_df%>%filter(!is.na(time))%>%select(assay,gene_id,time,logFC),
  contrdf%>%filter(!is.na(time))%>%mutate(assay='MS')%>%select(assay,gene_id,time,logFC=diff))%>%
  filter(!str_detect(assay,'ribo'))%>%
  group_by(gene_id)%>%
  filter(!any(is.na(logFC)))%>%
  filter(n()==12)%>%
  unite(contrast,time,assay)%>%
  spread(contrast,logFC)
t0_contrastdf = as.matrix(t0_contrastdf[,-1])%>%set_rownames(t0_contrastdf[[1]])
rownames(t0_contrastdf) = gid2gnm[[rownames(t0_contrastdf)]]

#and stepwise
stepcontrastdf = bind_rows(stepcountcontrdf%>%filter(!is.na(time))%>%select(assay,gene_id,time,logFC),
  stepstepprotcontrdf%>%filter(!is.na(time))%>%mutate(assay='MS')%>%select(assay,gene_id,time,logFC=diff))%>%
  filter(!str_detect(assay,'ribo'))%>%
  group_by(gene_id)%>%filter(!any(is.na(logFC)))%>%
  filter(n()==12)%>%
  unite(contrast,time,assay)%>%
  spread(contrast,logFC)
stepcontrastdf = as.matrix(stepcontrastdf[,-1])%>%set_rownames(stepcontrastdf[[1]])
rownames(stepcontrastdf) = gid2gnm[[rownames(stepcontrastdf)]]

#get t0 contrasts
t0_minmax = prediction_df%>%
  filter(!str_detect(assay,'ribo'))%>%
  group_by(gene_id)%>%
  filter(n_distinct(assay)==3)%>%
  group_by(gene_id,assay)%>%
  mutate(estimate = (2^estimate)%>%{.-min(.)}%>%divide_by(max(.)))%>%
  group_by(gene_id)%>%filter(!any(is.na(estimate)))%>%
  unite(contrast,time,assay)%>%
  distinct(contrast,gene_id,estimate)%>%
  spread(contrast,estimate)
t0_minmax = as.matrix(t0_minmax[,-1])%>%set_rownames(t0_minmax[[1]])
rownames(t0_minmax) = gid2gnm[[rownames(t0_minmax)]]


# predictionmat = prediction_df%>%unite(contrast,time,assay)%>%
#   select(contrast,estimate,gene_name)
#   spread(contrast,estimate)
# }
#We need to perform cluto clustering with filled in data.

# exprdf$gene_name = gid2gnm[[exprdf$gene_id]]
# prediction_df$gene_name = gid2gnm[[prediction_df$gene_id]]

stop()

{
datamats=list()
datamats['exprlevels'] = list(exprmat)
datamats['t0_effect'] = list(t0_contrastdf)
datamats['stepwise_effect'] = list(stepcontrastdf)


# for(i in seq_along(datamats)){rownames(datamats[[i]]) %<>% {upid2gnm[[.]]}}

clustlist = list()
}


{
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
}

