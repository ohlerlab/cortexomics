stopifnot(exists('ribo_t0_contrastdf'))
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