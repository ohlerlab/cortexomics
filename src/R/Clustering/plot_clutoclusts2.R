
make_cluster_trajplots<-function(clusts){
	indata = clusts$data %||% stop()
	plotname = clusts$name %||% stop()
    #process data for plot
 
  ggdf =  indata%>%as.data.frame%>%rownames_to_column('gene_name')%>%
     gather(dataset,value,-gene_name)%>%
    separate(dataset,c('time','assay','rep'),extra='warn',fill='right')%>%
    left_join(tibble(gene_name=names(clusts$cluster),cluster=clusts$cluster))%>%
    select(gene_name,time,assay,value,cluster)%>%
    filter(!is.na(value))%>%
    filter(!is.na(cluster))%>%
    group_by(cluster)%>%
    mutate(clustern=paste0('Cluster_',LETTERS[cluster],' n = ',n_distinct(gene_name)))%>%
    group_by(clustern,gene_name,assay,time)%>%
    summarise(value = mean(value))%>%
    group_by(gene_name,assay)%>%
    mutate(value = value-mean(value,na.rm=T))
    
  clutplot =ggplot(ggdf,aes(x=as_factor(time),y=value,group=gene_name,color=as.factor(clustern)))+
    #    geom_line(alpha=I(0.1))+
    geom_line(alpha=I(0.1))+
    facet_grid(assay~clustern)+
    coord_cartesian(ylim=c(-2,2))+
    scale_y_continuous('Centered log2-signal')+
    stat_summary(aes(x=as_factor(time),group=clustern,color=as.factor(clustern)),alpha=I(1),color=I('black'),fun=median, fun.min=median, fun.max = median, geom='line',linetype=2)+
    theme_bw()+
    ggtitle(paste0(plotname))
    #now save
    here('plots/clusters/')%>%dir.create(showWarn=F)
    plotfile = here(paste0('plots/clusters/',plotname,'.pdf'))

	cairo_pdf(w=16,h=12,plotfile%T>%{normalizePath(.)%>%message})
	print(clutplot)
	dev.off()

}

make_cluster_goplots <- function(oclustgores,clusteringname,nterms=10)  {
  #go_comparison_plot <-
   onts = c('BP','MF','CC') 
   ont=onts[1]
   lapply(onts,function(ont){
        clusteringname <- paste0(clusteringname,'_',ont)
    fname = here(str_interp('plots/${clusteringname}.goplot.pdf'))
    #
    gplot <- oclustgores%>%group_by(cluster)%>%
        filter(ontology==ont)%>%
        filter(elimFisher<0.05)%>%
        slice(1:nterms)%>%
        select(Term,elimFisher)%>%
        arrange(cluster,elimFisher)%>%
        ungroup%>%
        mutate(Term = as_factor(Term))%>%
        mutate(cluster = as.factor(LETTERS[as.numeric(cluster)]))%>%
    ggplot(.,aes(x=cluster,color=-log10(elimFisher),
                 # size=-log10(elimFisher),
                    y=Term),
         
           )+
    geom_point()+
    theme_bw()+
    facet_grid(cluster~.,scale='free_y')+
    ggtitle(str_interp('Shared GO Term Plot - ${clusteringname}\n ${ont}'))
    #
    cairo_pdf(w=16,h=12,fname%T>%{normalizePath(.)%>%message})
    print(gplot)
    dev.off()
    fname
})
}

lapply(clustlist,make_cluster_trajplots)%>%unlist
stop()
imap(clustergos,make_cluster_goplots)
clustlist[['Cluto_K10_Bis_large_voom_time_effects']]$cluster%>%{names(.)[.==1]}%>%na.omit

clustlist

