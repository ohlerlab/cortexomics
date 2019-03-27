

make_expr_dotplot1<-function(expression_table){
  

  assert_that(all(expression_table %has_name% c('signal','gene_name','assay','time')))
  assert_that(all(expression_table$gene_name%>%n_distinct%>%`==`(1)))
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
  name_i=expression_table$gene_name[1]
  expression_table$signal[expression_table$signal==0]<-NA
  #
  # browser()
  expression_table %>%
    # mutate(signal = signal+0.001)%>%
    group_by(gene_name,assay)%>%

    
    { 
      ggplot(data=.,aes(y=relative_signal,x = time,color=time))+
        geom_point(position='identity')+
        facet_grid( scale = 'free',~assay )+
        expand_limits(y=c(floor(min(.$relative_signal)),ceiling(max(.$relative_signal))))+
        scale_y_log2_linlabels+
        ggplot2::labs( ylab =  'T/T0 (TPM or LFQshares)')+
        ggtitle(str_interp("TPM/IMS for ${name_i}"))+
        theme_bw()+
        theme(text = element_text(size = 16))
      
    }
}




#Genes with no total rna regulation and te are MUCH more likely to be RPS!!
intgenes_is_rp<-techangegenes_nornachange_hassig_hasms%>%str_detect('^Rp[sl]')
nonintgenes_is_rp<-expression_table_gsum%>%
  group_by(gene_name)%>%
  filter(sum((assay=='MS') & (!is.na(signal) ))>2)%>%
  .$gene_name%>%unique%>%str_detect('^Rp[sl]')

fisher.test(matrix(c(table(nonintgenes_is_rp),table(intgenes_is_rp)),nrow=2))

table(nonintgenes_is_rp)
table(nonintgenes_is_rp)%>% {./sum(.)}


te_rpl_table<-te_fc_table%>%group_by(gene_name)%>%summarise(hastechange=  any(
      (adj_p_value[assay=='TE']<0.05)&(abs(log2fc)>log2(1.5))
      # (adj_p_value[assay=='TE']<0.05)
    ))%>%
  mutate(isrpl = str_detect(gene_name,'Rp[sl]\\d+[^k]'))%>%
  {table(.$hastechange,.$isrpl)}
te_rpl_table
te_rpl_table%>%fisher.test
te_rpl_table%>%apply(2,function(x) x/sum(x))

#but not for individual time stages? - barely signiicant, weirdly
te_rpl_table<-te_fc_table%>%group_by(gene_name)%>%filter(time=='E16')%>%summarise(hastechange=  any(
      (adj_p_value[assay=='TE']<0.05)&(abs(log2fc)>log2(1.5))
      # (adj_p_value[assay=='TE']<0.05)
    ))%>%
  mutate(isrpl = str_detect(gene_name,'Rp[sl]\\d+[^k]'))%>%
  {table(.$hastechange,.$isrpl)}
te_rpl_table
te_rpl_table%>%fisher.test
te_rpl_table%>%apply(2,function(x) x/sum(x))





# #for satb2 - we have consistent translational amplificaiton of upreg.
# make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Satb2'))
# totexpr%>%filter(gene_name=='Satb2')


# #this one is odd - not ottal RNA signal... 
# make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Sumo2'))
# totexpr%>%filter(gene_name=='Sumo2')
# ]

# #Nice - translational regulation, but no MS
# make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Elmod3'))
# totexpr%>%filter(gene_name=='Elmod3')

techangegenes_nornachange_hassig_hasms

techangegenes_nornachange_hassig_hasms_rps = techangegenes_nornachange_hassig_hasms%>%str_subset('^[rR]p[ls]')
techangegenes_nornachange_hassig_hasms_norps = techangegenes_nornachange_hassig_hasms%>%setdiff(techangegenes_nornachange_hassig_hasms_rps)
# #Nice - translational reuglation it seems - but the total has something wrong, missing data....
make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Rps26'))
# totexpr%>%filter(gene_name=='Rps26')

#This one is good - translational regulation only
make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Mpp6'))
totexpr%>%filter(gene_name=='Mpp6')

#this one kind of looksl ike it should have total RNA upregulation
make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Armcx3'))
totexpr%>%filter(gene_name=='Armcx3')

#beautiful upreg of ribo MS without total
make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Fam160a2'))
totexpr%>%filter(gene_name=='Fam160a2')

#Now that's weird - look at the pick up total is constant enough, ribo eratic, ms consistent drop
make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Rpl13'))
totexpr%>%filter(gene_name=='Rpl13')

#Now that's weird - MS kiiiinda looks like it's going up, but the ribo is giong down???
make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Ssr4'))
totexpr%>%filter(gene_name=='Ssr4')

#More weirdnes....
make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Sdf2'))
totexpr%>%filter(gene_name=='Sdf2')

#rps 9  - nice consistent translational down regulation.
make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Rps9'))
totexpr%>%filter(gene_name=='Rps9')


make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Rps9'))
totexpr%>%filter(gene_name=='Rps9')

#okay lots plot all the RPs.
pdf('rp_te_dotplots.pdf')
for(rpname in techangegenes_nornachange_hassig_hasms_rps){
  print(make_expr_dotplot1(expression_table_gsum%>%filter(gene_name==rpname)))
}
dev.off()


allrp_genes <- expression_table_gsum$gene_name%>%str_subset('Rp[sl]')%>%unique%>%grep(value=TRUE,patt='\\-ps',inver=TRUE,x=.)
pdf('rp_dotplots.pdf')
for(rpname in allrp_genes){
  message(rpname)
  print(make_expr_dotplot1(expression_table_gsum%>%filter(gene_name==rpname)))
}
dev.off()


make_agg_dotplot<-function(expression_table){
  

  assert_that(all(expression_table %has_name% c('signal','gene_name','assay','time')))
  assert_that(all(expression_table$gene_name%>%n_distinct%>%`==`(1)))
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
  name_i=expression_table$gene_name[1]
  expression_table$signal[expression_table$signal==0]<-NA
  #
  # browser()
  expression_table %>%
    # mutate(signal = signal+0.001)%>%
    group_by(gene_name,assay)%>%

    
    { 
      ggplot(data=.,aes(y=relative_signal,x = ntime))+
        geom_point(position='identity')+
        facet_grid( scale = 'free',~assay )+
        expand_limits(y=c(floor(min(.$relative_signal)),ceiling(max(.$relative_signal))))+
        scale_y_log2_linlabels+
        ggplot2::labs( ylab =  'T/T0 (TPM or LFQshares)')+
        ggtitle(str_interp("TPM/IMS for ${name_i}"))+
        scale_x_continuous(labels=timepoints,name='Stage')+
        stat_summary(fun.data=mean_cl_boot, geom='ribbon', alpha = I(0.5),c)+
        theme_bw()+
        theme(text = element_text(size = 16))
      
    }
}
expression_table_gsum$assay%<>%factor
levels(confinttbl$assay) <- c('MassSpec','RiboSeq','Total RNAseq',)

expression_table_gsum
expression_table_gsum$time%<>%str_replace('p','')
expression_table_gsum%<>%mutate(ntime = match(time,timepoints))
pdf('plots/agg_rpsplot',w=9,h=4) ;print( make_agg_dotplot(expression_table_gsum%>%filter(str_detect(gene_name,'Rp[sl]\\d+[^k]'))%>%ungroup%>%mutate(gene_name='All_Rps'))+
  coord_cartesian(ylim=c(-2,2)) ); dev.off()


# #also outlier
# make_expr_dotplot1(expression_table_gsum%>%filter(gene_name=='Armcx3'))
# totexpr%>%filter(gene_name=='Jund')

# totexpr%>%filter(gene_nm)