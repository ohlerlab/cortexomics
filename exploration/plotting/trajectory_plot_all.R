source(file.path(root,'exploration','make_expression_table.R'))
source(file.path(root,'exploration','modeling/make_fc_table.R'))

tetbl<- lfc_tbl%>%
  group_by(gene_name)%>%filter(assay=='TE')%>%
  mutate(adj_p_value = ifelse(is.na(adj_p_value),1,adj_p_value))%>%
  mutate(sig_lf2c = ifelse(adj_p_value<0.05,log2fc,0))

gene_classes <- tetbl %>%arrange(time)%>%summarise(class = case_when(
    all(sig_lf2c >= 0 ) & all(na.omit(sig_lf2c>lag(sig_lf2c))) ~ 'monotonic increasing (darkgreen)',
    all(sig_lf2c >= 0 ) & all(na.omit(sig_lf2c<lag(sig_lf2c))) ~ 'monotonic decreasing (darkred)',
    (sum(sig_lf2c!=0) == 1) & (last(sig_lf2c)!=0) ~ 'single TP TE change but P0 (yellow)', 
    (sum(sig_lf2c!=0) == 1) ~ 'single TP TE change (brown)', 
    ((any(sig_lf2c>0))&all(sig_lf2c>=0)) ~ 'inconsistent increasing (darkgreen)',
    ((any(sig_lf2c<0))&all(sig_lf2c<=0)) ~ 'inconsistent decreasing (darkred)',
    any(sig_lf2c!=0) ~ 'inconsistent TE change (blue)' ,
    all(sig_lf2c==0) ~ 'No sig TE change (grey)' ,
    TRUE ~ 'other (white)'
  ))


gene_classes <- tetbl %>%arrange(time)%>%summarise(class = case_when(
    # all(sig_lf2c>0) ~ 'allincreasing (darkgreen)',
    # all(sig_lf2c<0) ~ 'alldecreasing (darkred)',    
    any(sig_lf2c>0)&any(sig_lf2c<0) ~ 'inconsistent (pink)' ,
    any(sig_lf2c>0) ~ 'increasing (green)',
    any(sig_lf2c<0) ~ 'decreasing (red)',
    all(sig_lf2c==0) ~ 'No sig TE change (grey)' ,
    TRUE ~ 'other (white)'
  ))


# gene_classes <- tetbl %>%arrange(time)%>%summarise(class = case_when(
#     any(sig_lf2c!=0) ~ 'TE change (red)' ,
#     TRUE ~ 'No sig TE change (grey)'
#   ))

class_color_dict = gene_classes%>%distinct(class)%>%
  mutate(class_color = class%>%str_extract(regex('[^\\s]+?$'))%>%str_replace_all('[\\(\\)]',''))%>%
  mutate(class = class%>%str_replace_all(regex(' [^\\s]+?$'),''))%>%
  {setNames(.$class_color,.$class)}

gene_classes%<>%mutate(class = class%>%str_replace_all(regex(' [^\\s]+?$'),''))
stopifnot(all(gene_classes$class %in% names(class_color_dict)))

gene_classes$class%>%table

# gene_classes%>%group_by(class)%>%tally%>%ungroup%>%mutate(percentage = 100*(n / sum(n)) )%>%arrange(n)
# gene_classes$class%>%table

# lfc_tbl%>%group_by(gene_name)%>%mutate(sig_lf2c = ifelse(adj_p_value<0.05,log2fc,0))%>%

# gene_classes%>%filter(class=='inconsistent increasing')%>%filter(gene_name%>%{.==unique(.)[1]})%>%inner_join(tetbl)%>%
#   transmute(gene_name,adj_p_value,log2fc,sig_lf2c)
# gene_classes%>%filter(class=='single TP TE change')%>%filter(gene_name%>%{.==unique(.)%>%sample(1)})%>%inner_join(tetbl)%>%
#   transmute(gene_name,adj_p_value,log2fc,sig_lf2c)



# lfc_tbl_plot<-lfc_tbl%>%filter(gene_name=='Satb2')%>%filter(assay=='TE')#%>%make_expr_trajplot()
# lfc_tbl%>%filter(gene_name%>%{.%in%unique(.)[1:10]})%>%filter(assay=='TE')#%>%make_expr_trajplot()
# lfc_tbl_plot<-tetbl%>%filter(gene_name=='Satb2')%>%left_join(gene_classes)#%>%make_lfc_trajplot

#testrun
# trajplot <- tetbl%>%ungroup%>%filter((gene_name=='Satb2')|(gene_name%in%sample(unique(gene_name),10)))%>%left_join(gene_classes,by='gene_name')%>%make_lfc_trajplot



# stop()

make_lfc_trajplot<-function(lfc_tbl_plot){
  
  assert_that(all(  lfc_tbl_plot %has_name% c('log2fc','lmax','lmin','gene_name','assay','time')))
  assert_that(nrow(lfc_tbl_plot)!=0)

  lfc_tbl_plot %<>% filter(time %>% {.==unique(.)[1]}) %>%mutate(log2fc=0,time=factor('E13',levels=levels(time)))%>%rbind(lfc_tbl_plot)

  opac <- 1

    ggplot(data=lfc_tbl_plot,aes(y=log2fc,x = as.numeric(time),group=gene_name,color=class))+
      geom_line(data=lfc_tbl_plot%>%filter(class%in%c('other','No sig TE change')),position='identity',alpha = I(1)) +
      geom_line(data=lfc_tbl_plot%>%filter(!class%in%c('No sig TE change')),position='identity',alpha=I(opac)) +
      # geom_line(position='identity',alpha='0.5') +
      # expand_limits(y=c(floor(min(lfc_tbl_plot$log2fc)),ceiling(max(lfc_tbl_plot$log2fc))))+
      # coord_cartesian(y=c(-3,3))+
      ggtitle(str_interp("TE trajectory - all genes"))+
      theme_bw()+
      scale_color_manual(values=class_color_dict)+
      guides(colour = guide_legend(override.aes = list(alpha = 1)))+
      scale_x_continuous(labels=levels(time),name='Stage')+
      theme(text = element_text(size = 16))

      # geom_ribbon(data=lfc_tbl_plot,aes(ymax=lmax,ymin=lmin,x=ntime),alpha=0.5,fill=I('darkgreen'))+
        # geom_line(data=lfc_tbl_plot,aes(y=log2fc,x=ntime))+
        
}
tetbl%>%left_join(gene_classes,by='gene_name')%>%filter(sig_lf2c<0,class=='increasing')

gene_classes$class%>%table
#now for all genes
trajplot <- tetbl%>%left_join(gene_classes,by='gene_name')%>%make_lfc_trajplot
ggsave(plot=trajplot,file= file.path(root,'plots/trajectory_all_opaque.pdf')%T>%message)

gene_classes%>%qplot(data=.,x=class,fill=class,geom='bar')+scale_fill_manual(values=class_color_dict)+theme_bw()