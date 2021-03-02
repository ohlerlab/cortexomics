################################################################################
################################################################################

{
base::source(here::here('src/R/Rprofile.R'))

source(('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure4/0_plotclustfuns.R'))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gid2gnm<-load_hashmap('data/gid2gnm.hmp')
gnm2gid<-load_hashmap('data/gnm2gid.hmp')

sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
proDAfitms<-readRDS('data/proDAfitms.rds')
sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
sel_ms_mat<-readRDS('data/sel_ms_mat.rds')

countpred_df<-readRDS('data/countpred_df.rds')
countcontr_df<-readRDS('data/countcontr_df.rds')
tx_countdata<-readRDS('data/tx_countdata.rds')

allxtail = Sys.glob('pipeline/xtail/*')%>%map_df(.id='time',fread)%>%group_by(gene_name)
allxtail$gene_id = gnm2gid[[ allxtail$gene_name]]
techangedf <- allxtail%>%group_by(gene_id,gene_name)%>%
  mutate(sig = (adj_p_value < 0.05)& (abs(log2fc)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (log2fc > 0))),
    down = as.numeric(any(sig & (log2fc < 0)))
  )
 techangedf%>%write_tsv('tables/xtailTEchange.tsv') 
techangegenes = techangedf%>%filter(up==1|down==1)%>%.$gene_name
teupgenes = techangedf%>%filter(up==1)%>%.$gene_name
tedowngenes = techangedf%>%filter(down==1)%>%.$gene_name

cordist <- function(X){
  
  require(matrixStats)
  
  X = X - (rowMeans(X))
  
  D = X %*% t(X)
  
  mags = sqrt(diag(D))
  
  Dc = t(t ( D / mags  ) / mags )
  
  as.dist(1 - Dc)
  
}


}
techangedf%>%filter(down==1)

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
########contrast data matrices
################################################################################
{ 
countcontr_df <- readRDS(here('data/countcontr_df.rds'))
stepcountcontrdf <- readRDS(here('data/stepcountcontrdf.rds'))
contrdf<-readRDS('data/contrdf.rds')
stepstepprotcontrdf<-readRDS('data/stepcontrdf.rds')

stepcountcontrdf%>%group_by(gene_id)%>%group_slice(1); countcontr_df%>%group_by(gene_id)%>%group_slice(1)

prediction_df = bind_rows(
  sel_prodpreds%>%
    mutate(assay='MS')%>%
    select(gene_id,time,assay,estimate,CI.L,CI.R),
  countpred_df%>%
    separate(contrast,c('time','assay'))%>%
    select(gene_id,time,assay,estimate=logFC,CI.L,CI.R)
)


#get t0 contrasts
t0_contrastdf = bind_rows(countcontr_df%>%filter(!is.na(time))%>%select(assay,gene_id,time,logFC),
  contrdf%>%filter(!is.na(time))%>%mutate(assay='MS')%>%select(assay,gene_id,time,logFC=diff))%>%
  filter(!str_detect(assay,'ribo'))%>%
  group_by(gene_id)%>%
  filter(!any(is.na(logFC)))%>%
  # filter(n()==16)%>%
  filter(n()==12)%>%
  unite(contrast,time,assay)%>%
  spread(contrast,logFC)
t0_contrastdf = as.matrix(t0_contrastdf[,-1])%>%set_rownames(t0_contrastdf[[1]])
rownames(t0_contrastdf) = gid2gnm[[rownames(t0_contrastdf)]]

t0_contrastdf

#get t0 contrasts
ribo_t0_contrastdf = bind_rows(countcontr_df%>%filter(!is.na(time))%>%select(assay,gene_id,time,logFC),
  contrdf%>%filter(!is.na(time))%>%mutate(assay='MS')%>%select(assay,gene_id,time,logFC=diff))%>%
  # filter(!str_detect(assay,'TE'))%>%
  group_by(gene_id)%>%
  filter(!any(is.na(logFC)))%>%
  filter(n()==16)%>%
  unite(contrast,time,assay)%>%
  spread(contrast,logFC)
ribo_t0_contrastdf = as.matrix(ribo_t0_contrastdf[,-1])%>%set_rownames(ribo_t0_contrastdf[[1]])
rownames(ribo_t0_contrastdf) = gid2gnm[[rownames(ribo_t0_contrastdf)]]

#and stepwise
stepcontrastdf = bind_rows(stepcountcontrdf%>%filter(!is.na(time))%>%select(assay,gene_id,time,logFC),
  stepstepprotcontrdf%>%filter(!is.na(time))%>%mutate(assay='MS')%>%select(assay,gene_id,time,logFC=diff))%>%
  # filter(!str_detect(assay,'ribo'))%>%
  group_by(gene_id)%>%filter(!any(is.na(logFC)))%>%
  filter(n()==12)%>%
  unite(contrast,time,assay)%>%
  spread(contrast,logFC)
stepcontrastdf = as.matrix(stepcontrastdf[,-1])%>%set_rownames(stepcontrastdf[[1]])
rownames(stepcontrastdf) = gid2gnm[[rownames(stepcontrastdf)]]

# #get t0 contrasts
# t0_minmax = prediction_df%>%
#   filter(!str_detect(assay,'ribo'))%>%
#   group_by(gene_id)%>%
#   filter(n_distinct(assay)==3)%>%
#   group_by(gene_id,assay)%>%
#   mutate(estimate = (2^estimate)%>%{.-min(.)}%>%divide_by(max(.)))%>%
#   group_by(gene_id)%>%filter(!any(is.na(estimate)))%>%
#   unite(contrast,time,assay)%>%
#   distinct(contrast,gene_id,estimate)%>%
#   spread(contrast,estimate)
# t0_minmax = as.matrix(t0_minmax[,-1])%>%set_rownames(t0_minmax[[1]])
# rownames(t0_minmax) = gid2gnm[[rownames(t0_minmax)]]


# predictionmat = prediction_df%>%unite(contrast,time,assay)%>%
#   select(contrast,estimate,gene_name)
#   spread(contrast,estimate)
#We need to perform cluto clustering with filled in data.
genesets <- readRDS(here('data/genesets.rds'))
changegenes = genesets[-1]%>%unlist%>%unique%>%gid2gnm[[.]]%>%intersect(rownames(t0_contrastdf))

}

contrdf%>%mutate(gene_name=gid2gnm[[gene_id]])%>%
  # rename('estimate':=diff)%>%
  ungroup%>%
  # filter(assay=='MS')%>%
  mutate(width=CI.R-CI.L)%>%
  group_by(gene_name)%>%
  summarise(meanw = mean(width))%>%
  mutate(rwidth = rank(meanw)/length(meanw))%>%
  # mutate(estimate = estimate - median(estimate))%>%
  filter(gene_name%>%str_detect(regex(ignore_case=T,'bcl11b')))


# exprdf%>%
#   # mutate(gene_name=gid2gnm[[gene_id]])%>%
#   # rename('estimate':=diff)%>%
#   ungroup%>%
#   # filter(assay=='MS')%>%
#   mutate(sd=sd(signal,na.rm))%>%
#   group_by(gene_name)%>%
#   mutate(estimate = estimate - median(estimate))%>%filter(gene_name%>%str_detect(regex(ignore_case=T,'bcl11b')))


# clustdata<-stepcontrastdf[techangegenes%>%intersect(rownames(stepcontrastdf)),]
# clustdata<-stepcontrastdf[changegenes,]
# clustdata<-stepcontrastdf[,]
# clustdata<-t0_contrastdf[,]
# clustdata<-t0_contrastdf[changegenes,]
# clustdata <- ribo_t0_contrastdf[changegenes,]

{
require(ComplexHeatmap)
require(fastcluster)
 
library(factoextra)
library(preprocessCore)


qnormwithdims = function(x){
  xn = preprocessCore::normalize.quantiles(x)
  rownames(xn)=rownames(x)
  colnames(xn)=colnames(x)
  xn
}

conttype='stepwise'
conttype='pca_t0'
conttype='ribot0'
conttype='pca_t0_noribo'
conttype='allPca_t0_noribo'
# conttype='pca_t0_jribo_changegenes'
# conttype='pca_step'

stop()

if(conttype%>%str_detect('step')){
  clustdata = t0_contrastdf
  clustdata_clust <- stepcontrastdf
}else if(conttype%>%str_detect('t0')){
  clustdata=t0_contrastdf
  clustdata_clust <- clustdata
}

if(conttype %>%str_detect('noribo')){
  clustdata%>%colnames
  clustdata_clust%>%colnames
  clustdata = clustdata%>%{.[,colnames(.)%>%str_detect(neg=T,'ribo|TE')]}
  clustdata_clust = clustdata_clust%>%{.[,colnames(.)%>%str_detect(neg=T,'ribo|TE')]}
}
if(conttype %>%str_detect('jribo')){
  stopifnot(clustdata%>%colnames%>%str_detect('ribo')%>%any)  
  clustdata = clustdata%>%{.[,colnames(.)%>%str_detect(neg=T,'TE')]}
  clustdata_clust = clustdata_clust%>%{.[,colnames(.)%>%str_detect(neg=T,'TE')]}
}

if(conttype%>%str_detect('pca')){
#as opposed to for the display
  colassays = colnames(clustdata_clust)%>%str_extract('(?<=_).*')
  pcas = lapply(unique(colassays)%>%setNames(.,.),function(cassay){
    # cassay<<-cassay
    pca <- clustdata_clust[,colassays==cassay]%>%princomp
  })
  clustdata_clust <- pcas%>%map(~.$scores[,1:3])%>%imap(.,function(x,i)set_colnames(x,paste0('PC',1:ncol(x),'_',i))%>%as.matrix)%>%do.call(what=cbind)
}

if(conttype%>%str_detect('allPca')){
#as opposed to for the display
  colassays = colnames(clustdata_clust)%>%str_extract('(?<=_).*')
  pca = clustdata_clust%>%princomp

  clustdata_clust <- pca$scores[,1:8]%>%{set_colnames(.,paste0('PCall',1:ncol(.)))%>%as.matrix}

}

#now plot
# plotfile<- here(paste0('plots/','pcaplot','.pdf'))
# pdf(plotfile)
#perform pca
pcafit <- princomp(ribo_t0_contrastdf%>%as.data.frame%>%select(-matches('ribo'))%>%as.matrix)
#make plots of the pcas
pcaplot<-here('plots/pcafit_limmafcs_stepwise.pdf')
pdf(pcaplot)
# plot(pcafit,main='Fold Change Over Time - PCA')
plot(pcafit$scores[,1:2],ylim=pcafit$scores[,1:2]%>%range,xlim=pcafit$scores[,1:2]%>%range)
plot(pcafit$scores[,2:3],ylim=pcafit$scores[,2:3]%>%range,xlim=pcafit$scores[,2:3]%>%range)
plot(pcafit,main='Fold Change Over Time Seq Only - PCA')
plot(pcafit$scores[,1:2],ylim=pcafit$scores[,1:2]%>%range,xlim=pcafit$scores[,1:2]%>%range)
plot(pcafit$scores[,2:3],ylim=pcafit$scores[,2:3]%>%range,xlim=pcafit$scores[,2:3]%>%range)
dev.off()
pcaplot%>%normalizePath%>%message
plotpcafit<-function(i=1,pcafit){
  pcafit$loading[,i]%>%
  enframe('dimension','loading')%>%
  arrange(str_detect(dimension,'MS'),str_detect(dimension,'ribo|TE'),str_detect(dimension,'P0'),str_detect(dimension,'E17'),str_detect(dimension,'E16'),str_detect(dimension,'E14'))%>%
  mutate(color = case_when(
     str_detect(dimension,'total|rna|all') ~ 'blue',
     str_detect(dimension,'ribo|TE') ~ 'green',
     str_detect(dimension,'MS') ~ 'orange'
    ))%>%
  mutate(dimension=factor(dimension,unique(dimension)))%>%
  ggplot(aes(x=dimension,y=loading,fill=I(color)))+stat_identity(geom='bar')+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,vjust=0.5))+
  ggtitle(str_interp('Developmental Fold Changes \n PCA${i}'))
}
#
pcaloadingsplot<-here('plots/hclust_pca_loadings.pdf')
pdf(h=5,w=8,pcaloadingsplot%T>%{normalizePath(.)%>%message})
ggpubr::ggarrange(ncol=2,plotlist=lapply(1:2,FUN=plotpcafit,pcafit))
dev.off()
pcaloadingsplot%>%normalizePath%>%message
#
pcaloadingsplot<-here('plots/hclust_pca_loadings_2.pdf')
pdf(h=5,w=8,pcaloadingsplot%T>%{normalizePath(.)%>%message})
ggpubr::ggarrange(ncol=2,plotlist=lapply(1:4,FUN=plotpcafit,pcafit))
dev.off()
pcaloadingsplot%>%normalizePath%>%message
pcafit%>%summary

# %>%
#   ggplot(.,aes())+
#   scale_color_discrete(name='colorname',colorvals)+
#   scale_x_continuous(paste0('xname'))+
#   scale_y_continuous(paste0('yname'))+
#   ggtitle(paste0('title'))+
#   theme_bw()
# dev.off()
# message(normalizePath(plotfile))



if(conttype%>%str_detect('changegenes'))clustdata = clustdata[changegenes,]
if(conttype%>%str_detect('changegenes'))clustdata_clust = clustdata_clust[changegenes,]

coltitle=paste0('Contrasts ',conttype)
#quantile normalisation
rnames = rownames(clustdata)
# stopifnot(colnames(clustdata_clust)==colnames(clustdata))
cnames = colnames(clustdata)

clustdata_clust%>% as.data.frame%>%rownames_to_column('gene.name')%>%write_tsv(here(paste0('data/',conttype,'.tsv')))


}



colorder = t0_contrastdf%>%colnames%>%.[order(str_detect(.,'ribo'))]%>%.[order(str_detect(.,'TE'))]%>%.[order(str_detect(.,'MS'))]

clustdata_clust %>% as.data.frame%>%rownames_to_column('gene.name')%>%write_tsv(here('data/clustdata_clust.tsv'))
t0_contrastdf[,colorder] %>% as.data.frame%>%rownames_to_column('gene.name')%>%write_tsv(here('data/t0_clustdat.tsv'))

library(cluster)

{

distfun=cordist
distfun=dist
clustfun=hclust
# clustfun=agnes


kmax=13
dists = distfun(clustdata_clust)
row.hc <- clustfun(dists, method="ward")
rowordering = names(cutree(row.hc,1))%>%match(rownames(clustdata_clust))

cutreelistorig = (2:kmax)%>%setNames(.,paste0('k_',.))%>%lapply(function(k)cutree(row.hc,k=k)%>%{setNames(letters[.],names(.))})
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

KMAX=13
MAXCLUSTNUM=12
{

hmclust=row.hc

  clusterorder = cutreelistorig[[MAXCLUSTNUM-1]][labels(hmclust)]%>%unique

  # grpcolors['a']='black'
genesofinterest=c('Nes','Tle4','Flna','Satb2','Bcl11b')
library(circlize)
library(dendextend)
# row.dend=as.dendrogram(row.hc)
# plot(color_branches(cutree(row.dend,17),17))
# cutheights = dendextendRcpp::dendextendRcpp_heights_per_k.dendrogram(row.dend)
# clustercutdend = cutheights[as.character(MAXCLUSTNUM)]
# cdend = cut(row.dend,h=clustercutdend)
#
#
hmcols = clustdata%>%colnames%>%.[order(str_detect(.,'ribo'))]%>%.[order(str_detect(.,'TE'))]%>%.[order(str_detect(.,'MS'))]
gs2plot = genesofinterest
#
selclusts = cutreelistorig%>%last
selclustsn = selclusts%>%unique%>%as.factor%>%as.numeric%>%setNames(selclusts%>%unique)
grpcolors = gg_color_hue(kmax)%>%setNames(names(selclustsn))
# row_annotation = rowAnnotation( width = 0.2,df=as.data.frame(cutreelist),show_legend=F,col = collist,border=TRUE,na_col='black')
row_annotation = rowAnnotation(
    width = unit(4,'cm'),
    df=as.data.frame(selclusts)%>%set_colnames('grps'),
    col=list(grps=grpcolors),
    show_legend=F,border=TRUE,na_col='black')

colorder = clustdata%>%colnames%>%.[order(str_detect(.,'ribo'))]%>%.[order(str_detect(.,'TE'))]%>%.[order(str_detect(.,'MS'))]
#
teup2plot = teupgenes%>%intersect(rownames(clustdata))
disprownames = rownames(clustdata)%>%setNames(.,.)
disprownames[!disprownames%in%genesofinterest]=''
datahm <- colorder%>%split(.,str_extract(.,'(?<=_).*$')) %>%rev%>%imap(function(hmcols,i){
   Heatmap(
    column_title=,
    # name='Min/Max Normed Gene Expression',
      clustdata%>%.[,hmcols]%>%set_rownames(disprownames),
      row_order=NULL,
      cluster_rows = hmclust,
      col  = colorRamp2(c(-4, 0, 4), c('#8904B1','white','#FF8000')),
      # col  =  c('#8904B1','white','#FF8000'),
      # row_labels = gt_render(rownames(clustdata), col = rownames(clustdata)%>%as.factor%>%rainbow(5)[.]),
      show_row_dend=T,
      cluster_columns = FALSE,
      show_row_names = hmcols%>%str_detect('MS'),
      row_dend_width = unit(8,'cm'),
      # na_col = 'black'
    )
})
teribostring = names(datahm)%>%str_extract('TE|ribo')%>%na.omit%>%unique
#now plot
plotfile<- here(paste0('plots/','contrasts_heatmap',str_replace_all(coltitle,' ','_'),'.pdf'))
pdf(plotfile,w=10)
draw(datahm[['all']]+datahm[['ribo']]+datahm[['TE']]+datahm[['MS']]+row_annotation,column_title=coltitle,auto_adjust=F)
# draw(row_annotation + datahm[['all']]+purrr::reduce(.f='+',.x=datahm[teribostring])+datahm[['MS']],column_title=coltitle)
dev.off()
message(normalizePath(plotfile))


}

{
teribostring
capture.output(cutreelistorig[[KMAX-1]][genesofinterest])%>%paste0(collapse='\n')%>%message
cutreelistorig[[KMAX-1]][genesofinterest]%>%n_distinct
message(map_dbl(1:length(cutreelistorig),function(i)cutreelistorig[[i]][genesofinterest]%>%n_distinct)%>%`==`(4)%>%which%>%first)
}


hclusterlist = c(list(rep('a',length(cutreelistorig[[1]]))),cutreelistorig)
MAXCLUSTNUM = 13
hclustob = list(data=clustdata,names=paste0(str_replace_all(coltitle,' ','_'),'_effect'),cluster = hclusterlist[[MAXCLUSTNUM]])
hclustob$cluster%<>%.[order(match(.,clusterorder))]
hclustob$clustercols = grpcolors[clusterorder]







# {a
# #now plot
# # cdistmat = cordist(clustdata)%>%as.matrix
# nbclustplot2_dist<-fviz_nbclust(as.matrix(dists), function(mat,n) list(cluster=hclusterlist[[n]]), method = "wss",k.max = 26) +
#   labs(subtitle = "Elbow method - cosine distances")
# plotfile<- here(paste0('plots/','nb_clust_cors','.pdf'))
# pdf(plotfile)
# print(nbclustplot2_dist)
# dev.off()
# normalizePath(plotfile)
# }
{
clusterdteglm = hclustob$cluster%>%
  split(.,.)%>%
  map(~setNames(names(.)%in%teupgenes,names(.))%>%enframe('gene_name','istechange'))%>%
  bind_rows(.id='cluster')%>%
  mutate(tmp=TRUE)%>%
  spread(cluster,tmp)%>%
  mutate_all(replace_na,FALSE)%>%
  glm(data=.,as.formula(paste0('istechange ~ 0+ ',paste0(collapse='+',letters[1:MAXCLUSTNUM]))))
clusterdteglmdown = hclustob$cluster%>%
  split(.,.)%>%
  map(~setNames(names(.)%in%tedowngenes,names(.))%>%enframe('gene_name','istechange'))%>%
  bind_rows(.id='cluster')%>%
  mutate(tmp=TRUE)%>%
  spread(cluster,tmp)%>%
  mutate_all(replace_na,FALSE)%>%
  glm(data=.,as.formula(paste0('istechange ~ 0+ ',paste0(collapse='+',letters[1:MAXCLUSTNUM]))))
 # 
clusterdteglm%>%summary
clusterdteglmdown%>%summary

allgenes = names(hclustob$cluster)

teupenrichdf = map_df(.id='cluster',letters[1:MAXCLUSTNUM]%>%setNames(.,.),function(l)table(
  allgenes %in% teupgenes,allgenes %in%(hclustob$cluster%>%split(.,.)%>%.[[l]]%>%names)
)%>%fisher.test%>%tidy)

tedownenrichdf = map_df(.id='cluster',letters[1:MAXCLUSTNUM]%>%setNames(.,.),function(l)table(
  allgenes %in% tedowngenes,allgenes %in%(hclustob$cluster%>%split(.,.)%>%.[[l]]%>%names)
)%>%fisher.test%>%tidy)

enrichdf = teupenrichdf%>%left_join(tedownenrichdf,suffix=c('_up','_down'),by='cluster')

clusters = letters[1:MAXCLUSTNUM]%>%setNames(.,.)
dteenrichdf = enrichdf %>%left_join(map_df(.id='cluster',clusters,function(cl){
  clgenes = names(hclustob$cluster)[hclustob$cluster%in%cl]
  out=c(up=sum(clgenes%in%teupgenes),down=sum(clgenes%in%tedowngenes))
  out['nodte']=length(clgenes)-sum(out)
  out
}))
}


hclustob%>%make_cluster_trajplots(dteenrichdf)

{
library(dendextend)
row.dend=as.dendrogram(row.hc)
# plot(color_branches(cutree(row.dend,17),17))
cutheights = dendextendRcpp::dendextendRcpp_heights_per_k.dendrogram(row.dend)
clustercutdend = cutheights[as.character(MAXCLUSTNUM)]
cdend = cut(row.dend,h=clustercutdend)

cdend$lower%<>%setNames(cdend$lower%>%map(labels)%>%map(~ hclusterlist[[MAXCLUSTNUM]][.])%>%map_chr(unique))
labels(cdend$upper) <- cdend$lower%>%map(labels)%>%map(~ hclusterlist[[MAXCLUSTNUM]][.])%>%map_chr(unique)
labelcols = labels(cdend$upper)%>%as.factor%>%as.numeric%>%{gg_color_hue(max(.))[.]}

#now plot
plotfile<- here(paste0('plots/','clusterdend','.pdf'))
pdf(plotfile)
print({
color_labels(cdend$upper,labels=labels(cdend$upper),col=labelcols)%>%set('labels_cex',2)%>%plot()
})
dev.off()
normalizePath(plotfile)

}


# hclustob$cluster%>%enframe('gene_name','cluster')%>%inner_join(mcshanethalfs)%>%mutate(is_ned=McShane_deg_cat=='NED')%>%
#   group_by(cluster)%>%
#   summarise(mean(is_ned))%>%
#   arrange(match(cluster,clusterorder))

#now plot
source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure4/2_clustergos.R")
clustergos<-get_cluster_gos(hclustob$cluster)

plotfile<- here(paste0('plots/','cluster_go_bp',hclustob$name,'.pdf'));cairo_pdf(h=21,w=21,plotfile)
go_comparison_plot(clustergos%>%filter(ontology=='BP')%>%{split(.,as_factor(.$cluster))})+ggtitle('Biological Process')
dev.off()
normalizePath(plotfile)%>%message
clustergos%>%write_tsv(here('tables/cluster_go.tsv'))

################################################################################
########Now make split heatmaps with TE classes
################################################################################
{
teclasses = list(te_up = teupgenes,te_down = tedowngenes, te_nochange=rownames(clustdata)%>%setdiff(c(teupgenes,tedowngenes)))
teclasses%<>%map(intersect,rownames(clustdata))
teclass = teclasses[[1]]
teclassname = names(teclasses)[1]

hmaps = lapply(teclasses,function(teclass){
    teclass = teclass[order(hclustob$cluster[teclass])]
    selclusts = hclustob$cluster[teclass]
    selclustsn = selclusts%>%unique%>%as.factor%>%as.numeric%>%setNames(selclusts%>%unique)
    #
    row_annotation = rowAnnotation(
    width = unit(4,'cm'),
    df=as.data.frame(selclusts)%>%set_colnames('grps'),
    col=list(grps=gg_color_hue(kmax)[selclustsn]%>%setNames(names(selclustsn))),
    show_legend=F,border=TRUE,na_col='black')
    genesofinterest=c('Nes','Tle4','Flna','Satb2')
    library(circlize)
    library(dendextend)
    hmcols = clustdata%>%colnames%>%.[order(str_detect(.,'ribo'))]%>%.[order(str_detect(.,'TE'))]%>%.[order(str_detect(.,'MS'))]
    gs2plot = genesofinterest
    #
    colorder = clustdata%>%colnames%>%.[order(str_detect(.,'ribo'))]%>%.[order(str_detect(.,'TE'))]%>%.[order(str_detect(.,'MS'))]
    #
    teup2plot = teupgenes%>%intersect(rownames(clustdata))
    disprownames = rownames(clustdata)%>%setNames(.,.)
    disprownames[!disprownames%in%genesofinterest]=''
    datahm <- colorder%>%split(.,str_extract(.,'(?<=_).*$')) %>%rev%>%imap(function(hmcols,i){

     datahm=  Heatmap(
        column_title=,
        # name='Min/Max Normed Gene Expression',
          clustdata%>%.[teclass,hmcols]%>%set_rownames(disprownames[teclass]),
          row_order=NULL,
          # cluster_rows = hmclust,
          col  = colorRamp2(c(-4, 0, 4), c('#8904B1','white','#FF8000')),
          # col  =  c('#8904B1','white','#FF8000'),
          # row_labels = gt_render(rownames(clustdata), col = rownames(clustdata)%>%as.factor%>%rainbow(5)[.]),
          show_row_dend=T,
          cluster_columns = FALSE,
          show_row_names = TRUE,
          row_dend_width = unit(8,'cm'),
          # na_col = 'black'
        )
      if(any(str_detect(hmcols,'MS')))datahm = datahm + row_annotation
      datahm
    })
})
}

{
teribostring = names(datahm)%>%str_extract('TE|ribo')%>%na.omit%>%unique
#now plot
plotfile<- here(paste0('plots/','contrasts_heatmap',str_replace_all(coltitle,' ','_'),'split..pdf'))
pdf(plotfile,w=10,h=10)
pushViewport(viewport(layout = grid.layout(nrow=3,ncol=1)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
teclassname='te_up'
draw(hmaps[[teclassname]][['all']]+hmaps[[1]][['MS']],column_title=teclassname,auto_adjust=F,newpage=F)
popViewport()
teclassname='te_down'
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(hmaps[[teclassname]][['all']]+hmaps[[2]][['MS']],column_title=teclassname,auto_adjust=F,newpage=F)
popViewport()
pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1))
p = hclustob$cluster%>%enframe%>%filter(name%in%unlist(teclasses[1:2]))%>%
  mutate(te_class=ifelse(name%in%teclasses[[1]],'up','down'))%>%group_by(te_class,value)%>%
  summarise(n=n())%>%
  mutate(value = factor(value,hclustob$cluster%>%unique))%>%
  ggplot(data=.,aes(y=n,fill=value,x=value,color=I('black'),alpha=te_class))+
    geom_bar(stat='identity',position='dodge',width=0.5)+
    xlab('cluster')+ylab('n_genes')+theme_minimal()+
    scale_fill_manual(values = grpcolors )
print(p,  vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
dev.off()
message(normalizePath(plotfile))
}

hclustob %>% saveRDS(here('data/hclustob.rds'))


{
#now plot
plotfile<- here(paste0('plots/','nb_clust_data','.pdf'))
pdf(plotfile)
nbclustplot2<-factoextra::fviz_nbclust(clustdata_clust, function(mat,n) list(cluster=hclusterlist[[n]]), method = "wss",k.max = 13) +
  labs(subtitle = "Elbow method - processed data")
print(nbclustplot2)
dev.off()
normalizePath(plotfile)%>%message
}


################################################################################
########Legacy code
################################################################################
  
if(F){

stop()

nbclustplot_gap<-fviz_nbclust(as.data.frame(clustdata), function(mat,n) clutoclustlist_filt[[n]],  method = "gap_stat", nboot = 200, k.max=20) +
  labs(subtitle = "Gap_stat Method")

#now plot
plotfile<- here(paste0('plots/','gap_stat_hclust_cos_t0change','.pdf'))
pdf(plotfile)
nbclustplot_gap_dist<-fviz_nbclust(cdistmat, function(mat,n) hclusterlist[[n]],  method = "gap_stat", nboot = 10, k.max=5) +
  labs(subtitle = "Gap_stat Method - cosine distances")
dev.off()
normalizePath(plotfile)

nbclustplot_gap_dist


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


 
countcontr_df%>%group_by(gene_id)%>%group_slice(1)
stepcountcontrdf%>%group_by(gene_id)%>%group_slice(1)
techangers = techangegenes%>%intersect(rownames(stepcontrastdf))
head(t0_contrastdf==stepcontrastdf)
tmp <- dbscan::hdbscan(cordist(t0_contrastdf[changegenes,]),20)

################################################################################
########DBscan experiments
################################################################################
  
dbscanclustob = list(data=t0_contrastdf[changegenes,],name='HDBSCAN mp20 t0 changegenes',cluster=tmp$cluster)

# tmp <- dbscan::hdbscan(datamats$stepwise_effect[changinggenes,],10)
# tmp <- dbscan::hdbscan(t0_minmax,10)
tecols <- colnames(t0_contrastdf)%>%str_detect('TE')
# tmp <- dbscan::hdbscan(cordist(t0_contrastdf),50)
t0tedists = cordist(t0_contrastdf[techangers,tecols])
steptedists = cordist(stepcontrastdf[techangers,tecols])
t0tedists%>%head
steptedists%>%head
t0db<- dbscan::hdbscan(t0tedists,20)
stepdb<- dbscan::hdbscan(steptedists,20)
t0db
stepdb

# tmp <- dbscan::hdbscan(cordist(t0_minmax),20)
nrow(t0_minmax)

}
