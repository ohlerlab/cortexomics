# clusts<-hclustob
# clusts$data%>%head
# dtedf = dteenrichdf
make_cluster_trajplots<-function(clusts,dteenrichdf){
	indata = clusts$data %||% stop()
	plotname = clusts$name %||% stop()
    ccols = clusts$clustercols %||% stop()
    names(clusts$cluster) %||% stop('no names on the cluster vect')

    if(is.character(clusts$cluster)) clusts$cluster%<>%as_factor%>%setNames(names(clusts$cluster))
    clusters = clusts$cluster
    clustersizes = clusts$cluster%>%table
    clusternames = paste0('Cluster_',names(clustersizes),' n = ',clustersizes)%>%setNames(names(clustersizes))
    clusters = clusternames[clusters]%>%setNames(names(clusters))

    ccols = ccols[names(clusternames)]%>%setNames(clusternames)

    #process data for plot
  ggdf =  indata%>%as.data.frame%>%rownames_to_column('gene_name')%>%
     gather(dataset,value,-gene_name)%>%
    separate(dataset,c('time','assay','rep'),extra='warn',fill='right')%>%
    left_join(tibble(gene_name=names(clusters),clustern=clusters))%>%
    select(gene_name,time,assay,value,clustern)%>%
    filter(!is.na(value))%>%
    filter(!is.na(clustern))%>%
    group_by(clustern,gene_name,assay,time)%>%
    summarise(value = mean(value))%>%
    group_by(gene_name,assay)
  if(!'E13'%in%ggdf$time)ggdf=bind_rows(ggdf%>%distinct(clustern,gene_name,assay)%>%mutate(time='E13',value=0),ggdf)

  ggdf%<>%mutate(value = value-value[time=='E13'],na.rm=T)

  assayorder =  distinct(ungroup(ggdf),assay)%>%arrange(assay=='MS',assay=='TE',assay=='ribo')%>%.$assay

  ggdf$assay%<>%factor(assayorder)
  limwidth = if(str_detect(plotname,'effect')) 2 else 4

  ggdf$clustern%<>%factor(clusternames)

  ggdfsum<-ggdf%>%group_by(time,clustern,assay)%>%
    summarise(med=median(value),lower=quantile(value,0.25),upper=quantile(value,0.75))

  clutplot =ggplot(ggdf,aes(x=as_factor(time),y=value,group=gene_name,color=clustern))+
    # geom_line(alpha=I(0.1))+
    facet_grid(clustern ~ assay,scale='free')+
    scale_y_continuous('Centered log2-signal')+
    # coord_cartesian(ylim=c(-limwidth,limwidth))+
    # stat_summary(aes(x=as_factor(time),group=clustern,fill=clustern),fun=median, 
        # fun.min = function(x) quantile(x, 0.25), 
        # fun.max  = function(x) quantile(x, 0.75), 
        # geom='ribbon',alpha=I(0.5))+
    geom_ribbon(data=ggdfsum,aes(x=as_factor(time),group=clustern,fill=clustern,y=med,ymin=lower,ymax=upper))+
     stat_summary(aes(x=as_factor(time),group=clustern,fill=clustern),alpha=I(1),color=I('black'),fun=median, 
        geom=c('line'),linetype='dashed',alpha=I(1))+
     scale_fill_manual(values = ccols)
    theme_bw()+
    ggtitle(paste0(plotname))
    fracffdf=dteenrichdf%>%select(cluster,up,down,nodte)%>%gather(class,n,-cluster)%>%group_by(cluster)%>%mutate(n=n/sum(n))
    enrichdf=

    enrichggdf = dteenrichdf%>%select(cluster,estimate_up,estimate_down)%>%gather(class,estimate,-cluster)%>%mutate(class=str_replace(class,'estimate_',''))%>%
        left_join(dteenrichdf%>%select(cluster,conf.low_up,conf.low_down)%>%gather(class,conf.low,-cluster)%>%mutate(class=str_replace(class,'conf.low_','')))%>%
        left_join(dteenrichdf%>%select(cluster,conf.high_up,conf.high_down)%>%gather(class,conf.high,-cluster)%>%mutate(class=str_replace(class,'conf.high_','')))

 siglabels = dteenrichdf%>%select(cluster,p.value_down,p.value_up)%>%gather(class,val,-cluster)%>%mutate(siglabel=ifelse(p.adjust(val)<0.05,'*',''))%>%mutate(class=class%>%str_replace('p.value_',''))

  dteplot = ggplot(data=fracffdf%>%
    #filter(class!='nodte')%>%
    left_join(siglabels)%>%mutate(cluster = factor(cluster,unique(hclustob$cluster))),aes(x='',fill=class,y=n))+
    scale_fill_manual(values=c('up'='red','down'='blue','nodte'='grey'))+
    # scale_y_continuous('Fraction of Cluster in Class')+
    # stat_identity(geom='bar')+
    geom_bar(stat='identity',width=1,geom='bar')+
    facet_grid(cluster~.)+
    # geom_text(aes(label=siglabel))+
    coord_polar('y',start=0)+
    ggtitle('dTE class')+
     theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        # panel.grid  = element_blank()
    )

    enrichplot = enrichggdf%>%left_join(siglabels)%>%mutate(cluster = factor(cluster,unique(hclustob$cluster)))%>%ggplot(data=.,aes(x=class,y=log2(estimate),ymin=log2(conf.low),ymax=log2(conf.high),fill=class))+
     scale_fill_manual(values=c('up'='red','down'='blue','nodte'='grey'))+
     geom_bar(stat='identity',width=1,geom='bar')+
     facet_grid(cluster~.)+
     geom_text(aes(label=siglabel))+
     ggtitle('dTE Enrichment')
    # enrichplot
    #now save
    # here('plots/clusters/')%>%dir.create(showWarn=F)
    # plotfile = here(paste0('plots/clusters/',plotname,'.pdf'))
    # cairo_pdf(w=9,h=1*n_distinct(ggdf$clustern),plotfile%T>%{normalizePath(.)%>%message})
    # print(ggarrange(ncol=3,plotlist=list(clutplot,enrichplot,dteplot),widths=c(3,2,2)))
    # dev.off()
    ggarrange(ncol=3,plotlist=list(clutplot,enrichplot,dteplot),widths=c(3,2,2))

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

