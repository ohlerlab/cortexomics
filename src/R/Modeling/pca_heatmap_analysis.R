#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
# suppressMessages(library(DESeq2))
message('...done')


argv <- c(
  foldchangesfile='pipeline/exprdata/limma_fold_changes.txt',
  transformdexprfile=file.path('exprdata/transformed_data.txt')
)

if(length(commandArgs,trailing=T)) argv[]<- commandArgs(trailingOnly=TRUE)

for(i in names(argv)) assign(i, argv[i])

#read in expr data
exprtbl <- read_tsv(transformdexprfile) 
exprtbl %<>% select(gene_name, everything())

#read in the data
allcoefdf <- ebayes_stepwise%>%topTable(number=Inf)%>%rownames_to_column('uprotein_id')%>%select(-AveExpr,-F,-P.Value,-adj.P.Val)
allcoefmat <- allcoefdf[,-1]%>%set_rownames(allcoefdf[,1])%>%as.matrix

#perform pca
pcafit <- princomp(allcoefmat)
pcafitSeq <- princomp(allcoefmat%>%select(-dplyr::matches('MS')))

#make plots of the pcas
pcaplot<-here('plots/pcafit_limmafcs_stepwise.pdf')
pdf(pcaplot)
plot(pcafit,main='Fold Change Over Time - PCA')
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


pcaloadingsplot<-here('plots/limmafc_pca_loadings_spline.pdf')
pdf(h=5,w=8,pcaloadingsplot%T>%{normalizePath(.)%>%message})
ggpubr::ggarrange(ncol=2,plotlist=lapply(1:4,FUN=plotpcafit,pcafit))
dev.off()
pcaloadingsplot%>%normalizePath%>%message

pcaloadingsplot<-here('plots/limmafc_pca_loadings_spline_seq.pdf')
pdf(h=5,w=8,pcaloadingsplot%T>%{normalizePath(.)%>%message})
ggpubr::ggarrange(ncol=2,plotlist=lapply(1:4,FUN=plotpcafit,pcafitSeq))
dev.off()
pcaloadingsplot%>%normalizePath%>%message


hmapfile <- file.path(paste0('heatmaps/',basename(transformdexprfile),'_limma_cs.pdf'))
hmapfile%>%dirname%>%dir.create



#
hmapdata <- allcoefmat%>%pmax(.,-3)%>%pmin(.,3)
#
cairo_pdf(normalizePath(hmapfile) %T>% message,w = 10, h = 10)
#
cbreaks <-
  cut_number(hmapdata,100,ordered=TRUE)%>%levels%>%str_split(',')%>%map(str_replace_all,'[^0-9\\-\\.]','')%>%{c(.[[1]][[1]],map_chr(.,2))}%>%as.numeric
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
hmapfile%>%normalizePath%>%message


#TODO heatmaps with expr data not fold changes

#let's first just make some trajectory plots

#trajectory plot with raw data from heatamaps

geneofinterest='Shank2'

anticorgenes <- allcoefftbl%>%filter(abs(`timeP0`)>2,abs(`timeP0`+`timeP0:assayMS`)<0.25)%>%.$gene_name
selgenes <- sample(anticorgenes,4)

pdf('tot_ms_fc_anticorgenes_trajectory.pdf')
#
for(geneofinterest in selgenes[1]){
  par(mfrow=c(1,3))
  #get transformed expression data
  exprplotdata <-
    exprtbl%>%
    filter(gene_name==geneofinterest)%>%
    gather(key = etype, value = signal, -gene_name)%>%
    slice(-1) %>%
    separate(etype,into=c('time','assay','rep'))%>%
    # mutate(assay = ifelse(is.na(assay), 'common',assay))%>%
    mutate(signal = as.numeric(signal))
  #get limma fold changes
  foldchanges <-
    allcoefftbl%>%
    filter(gene_name==geneofinterest)%>%
    gather(key = vartype, value = signal, -gene_name)%>%
    # slice(-1) %>%
    separate(vartype,into=c('time','assay'))%>%
    mutate(assay = ifelse(is.na(assay), 'common',assay))%>%
    mutate(signal = as.numeric(signal))
  foldchanges$time%<>%str_replace('time','')
  foldchanges$assay%<>%str_replace('assay','')
    #now plot it all
    bind_rows(
      exprplotdata%>%mutate(level='transformed_expression'),
      foldchanges%>%mutate(level='foldchanges')
    )%>%
     { 
        ggplot(data=.,aes(y=signal,x = time))+
          geom_point(position='identity')+
          facet_grid(scale = 'free',level~assay )+
          # expand_limits(y=c(floor(min(.$signal)),ceiling(max(.$signal))))+
          ggplot2::labs( ylab =  'Limma Fold Changes')+
          ggtitle(str_interp("Limma Fold Change trajectory for ${geneofinterest}"))+
          theme_bw()+
          theme(text = element_text(size = 16),axis.text.x=element_text(angle=45,size=12,vjust= 0.5))
      }%>%print
  #
  message(geneofinterest)
}
getwd()%>%message
dev.off()