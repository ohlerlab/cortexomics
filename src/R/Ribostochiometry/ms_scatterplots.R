SCATTERSIGTYPE='norm_iBAQ'
ms_tall_sig <- ms_tall%>%filter(sigtype==SCATTERSIGTYPE)

translationprotids%<>%keep(~ ! . %in% rids$Protein_IDs)
protcatlist<-list(
  'RPL +'=rids %>% filter((`RPL_+`))%>%.$`Protein_IDs`,
  'RPS +'=rids %>% filter((`RPS_+`))%>%.$`Protein_IDs`,
  'Translation Associated' = translationprotids,
  'Ebp1'=ebp1pids)

scattercoldf<-file.path(root,'exploration/tables/scattercolors.tsv')%>%read_tsv
gethex <- function(color){
  c<-(col2rgb(color)[,1])
  sprintf("#%02X%02X%02X", c[1],c[2],c[3])
}
# scattercoldf%>%mutate(col = map2(fill,border,averagehexes))%>%unnest%>%write_tsv('~/cortexomics/scattercolors.tsv')
protcatcols = setNames(c(scattercoldf$col,gethex('grey'),gethex('black')),c('RPL +','RPS +','Ebp1','Other','Translation Associated'))


#
ebp1pids<-ms_data_all%>%filter(gene_name.x=='Pa2g4')%>%.$Protein_IDs
#did we add any

#let's replicate Matt's scatter plots
lfqe13<-ms_tall_sig%>%
    filter(time=='E13')%>%
    group_by(Protein_IDs,fraction)%>%
    summarize(Signal_E13=median(na.omit(signal)))
lfqlater<-ms_tall_sig%>%
    filter(!time=='E13')%>%
    group_by(Protein_IDs,fraction,time)%>%
    summarize(Signal=median(na.omit(signal)))
scatterdf<-left_join(by=c('fraction','Protein_IDs'),lfqe13,lfqlater)


message('Labeling protein IDs with annotations')
pids = ms_tall%>%ungroup%>%distinct(Protein_IDs)
pids%<>%mutate(pcat = case_when(
  (Protein_IDs==ebp1pid) ~ "Ebp1",
  sep_element_in(Protein_IDs,ridssplit) ~ "Ribosomal",
  sep_element_in(Protein_IDs,translationprotids) ~ "Translation Associated",
  TRUE ~ "other"
))
stopifnot("Ebp1" %in% pids$pcat)
# onegeneprotgroups = protiddt%>%



ggdf<-scatterdf%>%
  # filter(time=='P0')%>%
  filter(fraction %in% c('poly','cyto','80S'))%>%
  mutate(Protein = Protein_IDs%>%catagorize(protcatlist,default='Other'))%>%
  mutate(label=ifelse(Protein=='Ebp1','Ebp1',''))%>%
  filter(Protein_IDs%>%str_detect(';'))

ggdf$Protein%>%unique
(ggdf$Protein=='Other')%>%table
table(ggdf$Protein %in% names(protcatcols))
    # facet_grid(fraction~time)+
scatterplots<-ggdf%>%
  group_by(fraction,time)%>%
  # group_slice(1)%>%
  do({
    # browser();
    list(
    ggplot(.,aes(label=label,color = Protein,y=Signal_E13,x=Signal))+

    scale_color_manual(values=protcatcols)+
    geom_point(alpha=I(0.2),data=filter(.,Protein=='Other'),size=I(3))+
    geom_point(alpha=I(1),data=filter(.,Protein!='Other'),size=I(3))+
    geom_point(alpha=I(1),data=filter(.,Protein=='Ebp1'),size=I(3))+
    geom_text(show.legend=FALSE,nudge_x = 0.5,nudge_y=-0.5)+
    scale_x_log10(name = paste0('LFQ ',.$time[1]),lim=c(1e5,1e11))+
    scale_y_log10(name = paste0('LFQ ','E13'    ),lim=c(1e5,1e11))+
    geom_abline(slope=1,linetype='dashed')+
    ggtitle(paste0('Fraction: ',.$fraction[1]))+
    guides(text=FALSE)+
    theme_bw()+theme(aspect.ratio=1)+ 
    theme(plot.title = element_text(hjust = 0.5))
  )%>%data_frame(plot=.)})

scatterplotfile<-file.path(root,str_interp('plots/ms_fraction_scatterplots_${SCATTERSIGTYPE}.pdf'))
pdf(scatterplotfile,useDingbats=FALSE,w=24,h=18)
do.call(gridExtra::grid.arrange,c(scatterplots$plot,ncol=4))
dev.off()
scatterplotfile%>%message
stop()



message('Labeling protein IDs with annotations')
pids = ms_tall%>%ungroup%>%distinct(Protein_IDs)
pids%<>%mutate(pcat = case_when(
  (Protein_IDs==ebp1pid) ~ "Ebp1",
  sep_element_in(Protein_IDs,ridssplit) ~ "Ribosomal",
  sep_element_in(Protein_IDs,translationprotids) ~ "Translation Associated",
  TRUE ~ "other"
))
stopifnot("Ebp1" %in% pids$pcat)