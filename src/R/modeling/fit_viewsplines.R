


# stopifnot(nrow(l2fctable)>10e3)


#sum expr for isoforms
expression_table_gsum<-
  expression_table%>%group_by(gene_name,time,assay,replicate)%>%
    summarise(signal=sum(na.omit(signal)))
expression_table_gsum%<>%ungroup

#get genes with complete info
expression_table_gsum%<>% group_by(gene_name)%>%filter(sum(!is.na(signal))>1)%>%ungroup
expression_table_gsum%<>% group_by(gene_name)%>%filter(sum(signal>0)>1)%>%ungroup

#relative expression for plotting
expression_table_gsum%<>% group_by(gene_name,assay)%>%mutate(relative_signal = log2(signal/mean(na.omit(signal[time=='E13']))))

expression_table_gsum$assay[expression_table_gsum$assay=='MS']='MassSpec'



stop()


confinttbl%>%
  filter(gene_name %>% {. %in% unique(.)[1:100]})%>%
  


# modfunc <-  function(data){
  
#   lm(data=data,log2fc ~ ns(ntime,3)*assay)
#   # lm(data=data,log_signal ~ ntime*assay)
# }


# expression_table_models<-
#   l2fctable%>%group_by(gene_name)%>%group_slice(1:10) %>%
#   nest %>%
#   mutate(model = map(data,safely(modfunc) ))

# expression_table_models$model[1]

# expression_table_models%>%
#   mutate(
#     modelworked = map(model,'error')%>%map_lgl((is.null)),
#     model = map(model,'result')) %>%
#   filter(modelworked)

# expression_table_models<-expression_table_models%>%mutate(predicted = map2(model,data,safely(predict))) %>%
#   filter(!map(predicted,'error')%>%map_lgl(is.error)) %>%
#   filter(map(predicted,'warning')%>%map_lgl(is.null)) %>%
#   mutate(predicted = map(predicted,'result')) %>%
#   filter(map_dbl(predicted,length)%>%`==`(35))

# expression_table_models$model[[1]]%>%.$residuals%>%`^`(2)%>%sum
# expression_table_models%<>%mutate(sqe = map(model,'residuals')%>%map(`^`,2)%>%map_dbl(sum))

# expression_table_models$sqe%>%log10%>%hist(breaks=20,xlim = c(-2,3))
# expression_table_models$sqe%>%sum

# gene_spline_coeffs<-expression_table_models%>%mutate(coeffs = map(model,.%>%coefficients%>%stack%>%set_colnames(c('coef_value','coef_name'))))%>%
#   select(-model,-predicted,-data)%>%
#   unnest%>%
#   spread(coef_name,coef_value)

# gene_spline_mat<-gene_spline_coeffs%>%
#   select(matches('ns'))%>%
#   as.matrix%>%
#   set_rownames(gene_spline_coeffs$gene_name)%>%na.omit



#let's also calculate genes with very different trajectories
# txncoeffs=clusteringdata%>%as.data.frame%>%select(matches('ns\\(ntime.*?\\)\\d$'))
# mscoeffs=clusteringdata%>%as.data.frame%>%select(matches('ns\\(ntime.*?\\)\\d:assayMassSpec$'))
# diffs = abs(txncoeffs- mscoeffs)%>%rowSums
# tdiff_df<-diffs%>%stack%>%set_colnames(c('tdiff','gene_name'))


# #run kmeans on the number of clusters
gene_spline_dists <- dist(clusteringdata_noms[,-1])

# gene_spline_kmeans <- kmeans(gene_spline_dists,centers = 6)
# gene_spline_clustedf <- gene_spline_kmeans$cluster%>%stack%>%set_colnames(c('cluster','gene_name'))

# expression_table_predictions <- expression_table_models%>%select(-model)%>%unnest
# expression_table_predictions%<>%inner_join(gene_spline_clustedf)
# expression_table_predictions%<>%left_join(tdiff_df)




clusteringdata_noms$gene_name%>%n_distinct

clusteringdata_noms<-
  l2fctable%>%
  filter(assay!='MassSpec')%>%
  unite(dset,time,assay)%>%
  select(-ntime)%>%
  spread(dset,log2fc)

gene_spline_dists<- dist(clusteringdata_noms%>%{set_rownames(as.matrix(.[,-1]),.[[1]])})

havenas = clusteringdata_noms[,-1]%>%apply(1,is.na)%>%t%>%apply(1,any)
clusteringdata_noms%<>%.[!havenas,]
k.max <- 15

kmeanfits <- mclapply(1:k.max, function(k){
  clusteringdata_noms[,-1]%>%
  select(matches('TE'))%>%
  as.matrix%>%
  kmeans( k, nstart=50,iter.max = 8 )
})
kmeanfits%<>%simplify2array
plot(1:k.max, kmeanfits['tot.withinss',],
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# for(clustnum in 1:8){
for(clustnum in 6){

  kclusttable <-  
    data_frame(
      gene_name = clusteringdata_noms$gene_name,
      cluster = kmeanfits[1,clustnum][[1]] 
    )%>%
    left_join(l2fctable)%>%
    group_by(cluster)

  clustnames<-kclusttable%>%summarise(n=n_distinct(gene_name))%>%transmute(cluster,clustnamed=paste0(letters[cluster],' - ',n))
  kclusttable %<>% left_join(clustnames)
 kplot<- kclusttable%>%
  # group_slice(1:100)%>%
    # filter(!is.na(assay))%>%
    ggplot(aes(x=as.numeric(ntime),y=log2fc,color=assay,fill=assay))+
    stat_summary(fun.y = mean, geom='line')+
    # stat_summary(fun.ymax = partial(quantile,p=0.95) , fun.ymin = partial(quantile,p=0.05), geom='ribbon', alpha = I(0.5))+
    stat_summary(fun.data=mean_cl_boot, geom='ribbon', alpha = I(0.5))+
    scale_y_continuous(limits = c(-2,2))+
    scale_x_continuous(labels = timepoints[-1])+
    facet_grid(clustnamed~assay)
  ggsave(plot=kplot,w=9,file=str_interp(file.path(root,'plots/te_clustonly_cluster_lf2c_smoothplot_k_${clustnum}.pdf')))
  
  kclusttable%>%filter(gene_name=='Satb2')
}

stop()

#look at size of clusters...
kmeanfits[1,]%>%map(table)%>%map(~100*round(sort(./sum(.)),3))


#let's plot expression for sine utberesting genes.
expression_table<-bind_rows(
  totalmslfq%>%
      mutate(assay='MS')%>%ungroup%>%
      select(signal,isoform=Protein_IDs,gene_name,assay,time,replicate),
  kallistocounts%>%
    select(signal,isoform=transcript_id,gene_name,assay,time,replicate=rep)
)


#genes with translational regulation but no change in RNA?
techangegenes <- te_fc_table %>%
  group_by(gene_name) %>%
  mutate(sig_te = (assay=='TE') & (abs(log2fc)>log2(1.5) )) %>%
  filter(sum(sig_te)>1) %>%#has more than one stage te
  filter(n_distinct(log2fc[sig_te] > 0)==1 ) %>% #consistent direction of TE
  .$gene_name

techangegenes_nornachange <- totaltable%>%
  group_by(gene_name)%>%
  filter(!any(adj_p_value<0.05))%>%
  filter(gene_name %in% techangegenes)

#for filtering enes
totexpr<-expression_table%>%group_by(assay,gene_name,time)%>%summarise(signal=mean(na.omit(signal)))
noexpr_genes <- totexpr%>%group_by(gene_name,assay)%>%filter(all(signal==0))%>%.$gene_name


techangegenes_nornachange%<>%
  mutate(l2var = sd(log2fc))%>%
  arrange(l2var)
techangegenes_nornachange_hassig<- techangegenes_nornachange %>%filter(!gene_name %in% noexpr_genes)

expression_table_gsum$signal[expression_table_gsum$signal==0]<-NA

has_ms_genes<-  expression_table_gsum%>%
  group_by(gene_name)%>%
  filter(sum((assay=='MS') & (!is.na(signal) ))>2)%>%
  .$gene_name%>%unique

techangegenes_nornachange_hassig_hasms<- 
  techangegenes_nornachange %>%
  filter(gene_name %in% has_ms_genes)%>%
  .$gene_name%>%unique


'Satb2' %in% techangegenes

gene_spline_dists%>%colnames





# 
# 
# mod_model = expression_table_models%>%filter(gene_name=='Magohb')%>%.$model%>%.[[1]]
# 
# mod_model$coefficient[2:5]=gene_spline_kmeans$centers[1,]
# 
# #plot predicted clustering group
# expression_table_predictions%>%
#   ggplot(data=.,aes(y=log_signal,x = ntime))+
#   geom_smooth+
#   facet_grid( scale = 'free',.~assay )+
#   expand_limits(y=c((min(na.omit(.$log_signal))),(max(na.omit(.$log_signal)))))+
#   ggplot2::labs( ylab =  'T/T0 (TPM or LFQshares)')+
#   ggtitle(str_interp("TPM/IMS for ${name_i}"))+
#   theme_bw()+
#   theme(text = element_text(size = 16))
# 



#make fake data to illustrate the clusters

make_expr_dotplot<-function(name_i,expression_table){
  #
  assert_that(all(expression_table %has_name% c('signal','gene_name','assay','time')))
  assert_that(name_i %in% expression_table$gene_name)
  scale_y_log2_linlabels<-scale_y_continuous(
    breaks = function(lims) {
      nb = 7
      step = ceiling((lims[2] - lims[1])/7)
      lims[1] = nb*(floor(lims[1]/nb))
      lims[2] = nb*(ceiling(lims[2]/nb))
      seq(from=lims[1],to=lims[2],by=step);
    },
    labels = function(x) format((2^x),digits = 2)
  )
  names_i = c(name_i,str_replace(name_i,'\\-\\d+$',''))
  message(names_i)
  #
  expression_table %>% 
    filter(gene_name%in%names_i)%>%
    mutate(signal = signal+0.001)%>%
    group_by(gene_name)%>%
    mutate(signal = log2(signal/mean(na.omit(signal[time=='E13']))))%>%
    { 
      ggplot(data=.,aes(y=signal,x = time,color=ispredicted))+
        geom_point(position='identity')+
        facet_grid( scale = 'free',~assay )+
        expand_limits(y=c(floor(min(.$signal)),ceiling(max(.$signal))))+
        scale_y_log2_linlabels+
        ggplot2::labs( ylab =  'T/T0 (TPM or LFQshares)')+
        ggtitle(str_interp("TPM/IMS for ${name_i}"))+
        theme_bw()+
        theme(text = element_text(size = 16))
      
    }
}

expression_table

expression_table_predictionsplot <- expression_table_predictions%>%
  mutate(predicted = 2**predicted)%>%
  gather(ispredicted,signal,signal,predicted)

clusters = unique(expression_table_predictionsplot$cluster)
dir.create(file.path(root,paste0('exploration/clustplots/')),showWarnings = F)

icluster=2

mclapply( clusters,function(icluster){
  pdf(file.path(root,paste0('exploration/clustplots/cluster_',icluster,'_all.pdf'))%T>%message)
  genes = expression_table_predictionsplot%>%filter(cluster==icluster)%>%.$gene_name%>%unique
  # genes = genes[1:min(10,length(genes))]
  for(gene in genes){
    print(make_expr_dotplot(gene,expression_table_predictionsplot))
  }
  dev.off()
})
tightcoupledgenes <- expression_table_predictionsplot%>%arrange((tdiff))%>%.$gene_name%>%unique
decoupledgenes <- expression_table_predictionsplot%>%group_by(gene_name)%>%filter(!any(is.na(log_signal)))%>%arrange(desc(tdiff))%>%.$gene_name%>%unique

expression_table_predictionsplot$gene_name%>%n_distinct

print(make_expr_dotplot(tightcoupledgenes[5],expression_table_predictionsplot))

print(make_expr_dotplot(decoupledgenes[5],expression_table_predictionsplot))


expression_table_predictionsplot%>%filter(gene_name=='Pa2g4')%>%.$cluster

expression_table_predictionsplot%>%distinct(gene_name,cluster)%>%write_tsv(file.path(root,'exploration/gene_cluster_df.txt'))





