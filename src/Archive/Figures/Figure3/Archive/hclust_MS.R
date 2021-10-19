################################################################################
################################################################################
base::source(here::here('src/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("src/Figures/load_annotation.R")
}

sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
proDAfitms<-readRDS('data/proDAfitms.rds')
sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
sel_ms_mat<-readRDS('data/sel_ms_mat.rds')

countpred_df<-readRDS('data/countpred_df.rds')
countcontr_df<-readRDS('data/countcontr_df.rds')
tx_countdata<-readRDS('data/tx_countdata.rds')


X=data_scaled
cordist <- function(X){
  
  require(matrixStats)
  
  X = X - (rowMeans(X))
  
  D = X %*% t(X)
  
  mags = sqrt(diag(D))
  
  Dc = t(t ( D / mags  ) / mags )
  
  as.dist(1 - Dc)
  
}

# countpred_df$gnm = gid2gnm[[countpred_df$gene_id]]

# clustdata = countpred_df%>%distinct(contrast,logFC,gnm)%>%spread(contrast,logFC)
# clustdata = tx_countdata$abundance
# data_scaled = clustdata[,colnames(clustdata)%>%str_detect('total')]
# colnames(data_scaled)
# rownames(data_scaled) = clustdata$gnm
# if(str_detect(rownames(data_scaled),'ENSMUSG')) gid2gnm[[rownames(data_scaled)]]
# data_scaled %<>% as.matrix
# col.dd <- hclust(cordist(data_scaled), method="complete")

################################################################################
########
################################################################################
  

require(ComplexHeatmap)
require(fastcluster)
  
dim(datamats$t0_effect)

tmp <- dbscan::hdbscan(datamats$exprlevels,20)
tmp <- dbscan::hdbscan(datamats$t0_effect,20)
tmp <- dbscan::hdbscan(datamats$stepwise_effect[changinggenes,],10)
tmp <- dbscan::hdbscan(t0_minmax,10)

tmp <- dbscan::hdbscan(cordist(t0_minmax),20)
nrow(t0_minmax)
row.hc  <- hclust(cordist(t0_minmax), method="complete")

cutree(row.hc,3)%>%table

clustdata<-t0_minmax
kmax=20
row.hc=row.dd
rowordering = names(cutree(row.hc,1))%>%match(rownames(clustdata))
{
cutreelistorig = (2:kmax)%>%setNames(.,paste0('k_',.))%>%lapply(function(k)letters[cutree(row.hc,k=k)])
cutreelist=cutreelistorig
#RECODE the bottom level to alphabetic ordering
unsort2sort_clustname <- cutreelistorig[[kmax-1]][rowordering]%>%unique%>%setNames(sort(.),.)
cutreelist[[kmax-1]] %<>% unsort2sort_clustname[.]%>%setNames(NULL)
#First of all, recode the hierarchical clusterings to give consistent colorings between layers of the plot
#so  a,b,c,d -> a,b,d -> a,d etc.
for(ind in kmax:3){
  nextvals <- split(
    cutreelist[[ind-1]],
    cutreelist[[ind-2]]
  )%>%map_chr(1)
  cutreelist[[ind-2]] %<>% recode_factor(!!!nextvals)%>%as.character
}
#createa  list of colors for the annotation bars
#collist = list(setNames(rainbow(kmax),letters[1:(kmax)]) )%>%rep(kmax-1)%>%setNames(names(cutreelist))
collist = list(setNames(rep('white',kmax),letters[1:(kmax)]) )%>%rep(kmax-1)%>%setNames(names(cutreelist))
#
#Now color in the new cluster in each of the sequence with balck
newval=cutreelist[[1]]%>%unique%>%head(1)
for(i in seq(1,length(collist),1)){
  if(i!=1) newval = cutreelist[[i]]%>%setdiff(cutreelist[[i-1]])%>%unique
  collist[[i]][]='white'
  collist[[i]][newval]='black'
  # names(cutreelist)[i]=newval
  # names(collist)[i]=newval

  #also color the borders of each cluster, so we can see where they are
  borderinds <- data.frame(cat = cutreelist[[i]],ind =rowordering)%>%
  group_by(cat)%>%slice(c(1,n()))%>%.$ind

  grpvect = cutreelist[[i]][rowordering]
  cutreelist[[i]][rowordering[cumsum(head(Rle(grpvect)@lengths,-1))]]=NA
}
}

# row_annotation = rowAnnotation( width = 0.2,df=as.data.frame(cutreelist),show_legend=F,col = collist,border=TRUE,na_col='black')
row_annotation = rowAnnotation( width = 0.2,df=as.data.frame(cutreelistorig)[1:10],show_legend=F,border=TRUE,na_col='black')

colorder = t0_minmax%>%colnames%>%.[order(str_detect(.,'TE'))]%>%.[order(str_detect(.,'MS'))]
#
colorder%>%split(.,str_extract(['']))%>%datahm <- Heatmap(
  #name='Min/Max Normed Gene Expression',
      t0_minmax[,colorder],
      row_order=NULL,
      cluster_rows = row.dd,
      cluster_columns = FALSE,
      show_row_names = FALSE,
      left_annotation = row_annotation,
      # row_dend_width = unit(4,'cm'),
      # na_col = 'black'
)
datahm

print(datahm)

groupfreqsort <- function(.)arrange(.,
		desc(txngroup%>%str_detect('Change')),desc(txngroup%>%str_detect('Both')),desc(txngroup%>%str_detect('Down')),
		desc(msgroup%>%str_detect('Change')),desc(msgroup%>%str_detect('Both')),desc(msgroup%>%str_detect('Down')),
		desc(tegroup%>%str_detect('Change')),desc(tegroup%>%str_detect('Both')),desc(tegroup%>%str_detect('Down')))%>%
	filter(!str_detect(txngroup,'Both'),!str_detect(tegroup,'Both'),!str_detect(msgroup,'Both'))%>%
	as.data.frame


clustdf = cutree(col.dd,3)%>%enframe('gene_name','clust')


txnte_trajdf %>%left_join(clustdf)%>%group_by(txngroup,tegroup,msgroup)%>%summarise(table(clust)%>%data.frame)%>%
	spread(clust,Freq)%>%groupfreqsort

mschangedf%>%left_join(techangedf,by='gene_id',suffix=c('_MS','_TE'))%>%
	group_by(up_MS,down_MS,up_TE,down_TE)%>%
	left_join(clustdf)%>%
	filter(! (up_MS&down_MS),!(up_TE & down_TE))%>%
	group_by(up_MS,down_MS,up_TE,down_TE)%>%
	summarise(n=n(),table(clust)%>%as.data.frame)%>%
	mutate(Freq = Freq / sum(Freq))%>%
	spread(clust,Freq)


stopifnot('clust' %in% colnames(hclustcatdf))

clustexprdf = clustdata%>%inner_join(hclustcatdf)

if(!'TE'%in%clustexprdf$assay & ('ribo'%in%clustexprdf$assay)){
	clustexprdf%<>%
	clustexprdf%<>%bind_rows(clustexprdf%>%filter(assay=='ribo')%>%mutate(assay='TE',signal = signal - (clustexprdf%>%filter(assay=='total')%>%.$signal)))
}

ggdf = clustexprdf%>%
  arrange(time,assay=='MS',assay=='TE',assay=='ribo')%>%
  mutate(assay=as_factor(assay))%>%
  group_by(gene_name,assay,time)%>%
  summarise(signal=mean(signal,na.rm=TRUE))%>%
  group_by(gene_name,assay,cat)%>%
  mutate(signal = signal - median(signal))

medggdf = ggdf%>%group_by(assay,cat,time)%>%summarise(signal = median(signal))

#now plot
plotfile<- here('plots/figures/figure1/hclust_traj.pdf')
pdf(plotfile,w=14,h=7)
ggdf%>%
    # filter(signal < -3)
  ggplot(.,aes(y=signal,x=as.numeric(as.factor(time)),group=gene_name))+
  geom_line(alpha=I(.8),color=I("grey"))+
  # scale_color_discrete(name='colorname',colorvals)+
  scale_x_continuous(paste0('Time'),labels = hclustexprdf$time%>%unique)+
  scale_y_continuous(paste0('log2(CPM/iBAQ) relative to Mean'),limits=c(-1,1))+
  ggtitle(paste0('Ribosomal Proteins - Expression Trajectory'))+
  facet_grid(clust~assay,scales='free')+
  geom_line(data=medggdf,aes(group=cat),color=I('black'))+
  theme_bw()
dev.off()
normalizePath(plotfile)