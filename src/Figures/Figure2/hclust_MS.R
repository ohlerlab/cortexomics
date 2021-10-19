################################################################################
################################################################################
library(ComplexHeatmap)

{
base::source(here::here('src/Rprofile.R'))

source(('src/Figures/Figure4/0_plotclustfuns.R'))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
# proDAfitms<-readRDS('data/proDAfitms.rds')
# sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
# sel_ms_mat<-readRDS('data/sel_ms_mat.rds')

# # countpred_df<-readRDS('data/countpred_df.rds')
# countcontr_df<-readRDS('data/countcontr_df.rds')

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

allxtail%>%group_by(gene_id,gene_name)%>%
  mutate(sig = (adj_p_value < 0.05)& (abs(log2fc)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (log2fc > 0))),
    down = as.numeric(any(sig & (log2fc < 0)))
  )%>%{list(sum(na.omit(.$up)),sum(na.omit(.$down)),sum(na.omit(.$up+.$down)==2))}

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
# stepcountcontrdf <- readRDS(here('data/stepcountcontrdf.rds'))
prot_contrdf<-readRDS('data/contrdf.rds')
# stepstepprotcontrdf<-readRDS('data/stepcontrdf.rds')

# prediction_df = bind_rows(
#   sel_prodpreds%>%
#     mutate(assay='MS')%>%
#     select(gene_id,time,assay,estimate,CI.L,CI.R),
#   countpred_df%>%
#     separate(contrast,c('time','assay'))%>%
#     select(gene_id,time,assay,estimate=logFC,CI.L,CI.R)
# )


#get t0 contrasts
t0_contrastdf = bind_rows(countcontr_df%>%filter(!is.na(time))%>%select(assay,gene_id,time,logFC),
  prot_contrdf%>%filter(!is.na(time))%>%mutate(assay='MS')%>%select(assay,gene_id,time,logFC=diff))%>%
  filter(!str_detect(assay,'ribo'))%>%
  group_by(gene_id)%>%
  filter(!any(is.na(logFC)))%>%
  # filter(n()==16)%>%
  filter(n()==12)%>%
  unite(contrast,time,assay)%>%
  spread(contrast,logFC)
t0_contrastdf = as.matrix(t0_contrastdf[,-1])%>%set_rownames(t0_contrastdf[[1]])
rownames(t0_contrastdf) = gid2gnm[[rownames(t0_contrastdf)]]

#get t0 contrasts
ribo_t0_contrastdf = bind_rows(countcontr_df%>%filter(!is.na(time))%>%select(assay,gene_id,time,logFC),
  prot_contrdf%>%filter(!is.na(time))%>%mutate(assay='MS')%>%select(assay,gene_id,time,logFC=diff))%>%
  # filter(!str_detect(assay,'TE'))%>%
  group_by(gene_id)%>%
  filter(!any(is.na(logFC)))%>%
  filter(n()==16)%>%
  unite(contrast,time,assay)%>%
  spread(contrast,logFC)
ribo_t0_contrastdf = as.matrix(ribo_t0_contrastdf[,-1])%>%set_rownames(ribo_t0_contrastdf[[1]])
rownames(ribo_t0_contrastdf) = gid2gnm[[rownames(ribo_t0_contrastdf)]]

# #and stepwise
# stepcontrastdf = bind_rows(stepcountcontrdf%>%filter(!is.na(time))%>%select(assay,gene_id,time,logFC),
#   stepstepprotcontrdf%>%filter(!is.na(time))%>%mutate(assay='MS')%>%select(assay,gene_id,time,logFC=diff))%>%
#   # filter(!str_detect(assay,'ribo'))%>%
#   group_by(gene_id)%>%filter(!any(is.na(logFC)))%>%
#   filter(n()==12)%>%
#   unite(contrast,time,assay)%>%
#   spread(contrast,logFC)
# stepcontrastdf = as.matrix(stepcontrastdf[,-1])%>%set_rownames(stepcontrastdf[[1]])
# rownames(stepcontrastdf) = gid2gnm[[rownames(stepcontrastdf)]]

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
# library(preprocessCore)


# qnormwithdims = function(x){
#   xn = preprocessCore::normalize.quantiles(x)
#   rownames(xn)=rownames(x)
#   colnames(xn)=colnames(x)
#   xn
# }

conttype='stepwise'
conttype='pca_t0'
conttype='ribot0'
conttype='pca_t0_noribo'
conttype='allPca_t0_noribo'
# conttype='pca_t0_jribo_changegenes'
# conttype='pca_step'


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
# library(circlize)
# library(dendextend)
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
dir.create('plots/Figures/Figure4/',rec=T,showWarn=F)
plotfile<- here(paste0('plots/Figures/Figure4/','contrasts_heatmap',str_replace_all(coltitle,' ','_'),'.pdf'))
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


# hclustob%>%make_cluster_trajplots(dteenrichdf)

# {
# library(dendextend)
# row.dend=as.dendrogram(row.hc)
# # plot(color_branches(cutree(row.dend,17),17))
# cutheights = dendextendRcpp::dendextendRcpp_heights_per_k.dendrogram(row.dend)
# clustercutdend = cutheights[as.character(MAXCLUSTNUM)]
# cdend = cut(row.dend,h=clustercutdend)

# cdend$lower%<>%setNames(cdend$lower%>%map(labels)%>%map(~ hclusterlist[[MAXCLUSTNUM]][.])%>%map_chr(unique))
# labels(cdend$upper) <- cdend$lower%>%map(labels)%>%map(~ hclusterlist[[MAXCLUSTNUM]][.])%>%map_chr(unique)
# labelcols = labels(cdend$upper)%>%as.factor%>%as.numeric%>%{gg_color_hue(max(.))[.]}

# #now plot
# plotfile<- here(paste0('plots/','clusterdend','.pdf'))
# pdf(plotfile)
# print({
# color_labels(cdend$upper,labels=labels(cdend$upper),col=labelcols)%>%set('labels_cex',2)%>%plot()
# })
# dev.off()
# normalizePath(plotfile)

# }

#now plot
source("src/Figures/Figure4/1_go_term_funcs.R")
library(topGO)
clustergos<-get_cluster_gos(hclustob$cluster)

#
plotfile<- here(paste0('plots/Figures/Figure4/','cluster_go_bp',hclustob$name,'.pdf'));cairo_pdf(h=21,w=21,plotfile)
go_comparison_plot(clustergos%>%filter(ontology=='BP')%>%{split(.,as_factor(.$cluster))})+ggtitle('Biological Process')
dev.off()
normalizePath(plotfile)%>%message
#
clustergos%>%write_tsv(here('tables/cluster_go.tsv'))
hclustob$cluster%>%enframe('gene_name','cluster')%>%mutate(gene_id=gnm2gid[[gene_name]])%>%write_tsv('tables/gene_clusters.tsv')

if(FALSE){

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
      # library(dendextend)
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

}

# {
# #now plot
# plotfile<- here(paste0('plots/','nb_clust_data','.pdf'))
# pdf(plotfile)
# nbclustplot2<-factoextra::fviz_nbclust(clustdata_clust, function(mat,n) list(cluster=hclusterlist[[n]]), method = "wss",k.max = 13) +
#   labs(subtitle = "Elbow method - processed data")
# print(nbclustplot2)
# dev.off()
# normalizePath(plotfile)%>%message
# }

