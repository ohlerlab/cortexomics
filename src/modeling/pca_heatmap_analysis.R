#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
# suppressMessages(library(DESeq2))
message('...done')


defaultargs <- c(
  foldchangesfile='exprdata/limma_fold_changes.txt',
  transformdexprfile=file.path('exprdata/transformed_data.txt')
)

if(!exists('cargs')) cargs<- commandArgs(trailingOnly=TRUE)

cargs <- coalesce(
  cargs[1:length(defaultargs)]%>%setNames(names(defaultargs)),
  defaultargs
)
for(i in names(cargs)) assign(i, cargs[i])

#read in expr data
exprtbl <- read_tsv(transformdexprfile) 
exprtbl %<>% select(gene_name, everything())

#read in the data
allcoefftbl<-read_tsv(foldchangesfile)
#perform pca
pcafit <- princomp(allcoefftbl[,-1])

#make plots of the pcas
pdf('tmp.pdf')
plot(pcafit)
plot(pcafit$scores[,1:2])
plot(pcafit$scores[,c(1:3)])
plot(pcafit$scores[,c(2,3)])
dev.off()
getwd()%>%message



# svdfit <- svd(t(allcoefftbl[,-1]))
# svdfit%>%str


pdf('tmp.pdf')
plot(pcafit)
dev.off()
getwd()%>%message


hmapfile <- file.path(paste0('heatmaps/',basename(transformdexprfile),'_limma_cs.pdf'))
hmapfile%>%dirname%>%dir.create

#
hmapdata <- as.matrix(allcoefftbl[,-1])%>%pmax(.,-3)%>%pmin(.,3)
#
cairo_pdf(normalizePath(hmapfile) %T>% message,w = 10, h = 10)
#
cbreaks <-
  cut_number(as.matrix(allcoefftbl[,-1]),100,ordered=TRUE)%>%levels%>%str_split(',')%>%map(str_replace_all,'[^0-9\\-\\.]','')%>%{c(.[[1]][[1]],map_chr(.,2))}%>%as.numeric
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

#TODO heatmaps with expr data not fold changes

#let's first just make some trajectory plots

#trajectory plot with raw data from heatamaps

geneofinterest='Shank2'

anticorgenes <- allcoefftbl%>%filter(abs(`timeP0`)>2,abs(`timeP0`+`timeP0:assayMS`)<0.25)%>%.$gene_name
selgenes <- sample(anticorgenes,20)

pdf('tot_ms_fc_anticorgenes_trajectory.pdf')
#
for(geneofinterest in selgenes){
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