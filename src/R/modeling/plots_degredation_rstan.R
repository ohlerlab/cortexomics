
#` #### Plots showing distribution of overall fits
#Histogram of the distriubtion of sd for lrTE
pdf('../plots/modelling/sd_dist_rTE.pdf');
plot<-
  allsummary%>%filter(parameter=='lrTE')%>%
  qplot(data=.,x=sd,geom='histogram')+
    theme_bw()
print(plot)
dev.off()
message(normalizePath('../plots/modelling/sd_dist_rTE.pdf'))

#mark these genes

allsummary%>%colnames


pdf('../plots/modelling/indiv_rTEs_linerange.pdf');
plot<-
allsummary%>%filter(parameter=='lrTE')%>%
  mutate(nondet_rTE=sd>2.5)%>%
  arrange(nondet_rTE)%>%
  qplot(data=.,x=mean,xmin=`2.5%`,xmax=`97.5%`,y=as.numeric(as_factor(gene_name)),geom='blank')+
    ggstance::geom_linerangeh()+
    scale_x_continuous(name='Log rTE')+
    scale_y_continuous(name='Gene')+
    theme_bw()
print(plot)
dev.off()
message(normalizePath('../plots/modelling/indiv_rTEs_linerange.pdf'))



hmu_lrTE <- hallsummary%>%filter(parameter%in%c('hmu_lrTE'))%>%.$mean
hsig_lrTE <- hallsummary%>%filter(parameter%in%c('hsig_lrTE'))%>%.$mean%>%exp
groupparamdata <- tibble(mean=hmu_lrTE,`2.5%`=hmu_lrTE-(2*hsig_lrTE),`97.5%`=hmu_lrTE+(2*hsig_lrTE),gene_name='Group')
#Now on the joint model
pdf('../plots/modelling/hierarch_rTEs_linerange.pdf');
plot<-
hallsummary%>%
  filter(parameter%in%c('lrTE'))%>%
  # mutate(nondet_rTE=sd>2.5)%>%
  # arrange(nondet_rTE)%>%
  qplot(data=.,x=mean,xmin=`2.5%`,xmax=`97.5%`,y=as.numeric(as_factor(gene_name)),geom='blank')+
    ggstance::geom_linerangeh()+
    ggstance::geom_linerangeh(data=groupparamdata)+
    scale_x_continuous(name='Log rTE')+
    scale_y_continuous(name='Gene')+
    facet_grid(scale='free',ifelse(gene_name=='Group','Group','Indiv') ~ . )+
    theme_bw()
print(plot)
dev.off()
message(normalizePath('../plots/modelling/hierarch_rTEs_linerange.pdf'))

library(ggExtra)

pdf('../plots/modelling/indiv_rTEs_linerange_detonly.pdf',h=12);
plot<-
allsummary%>%filter(parameter=='lrTE')%>%
  mutate(nondet_rTE=sd>2.5)%>%
  filter(!nondet_rTE)%>%
  qplot(data=.,x=log10(exp(mean)),xmin=log10(exp(`2.5%`)),xmax=log10(exp(`97.5%`)),y=as.numeric(as_factor(gene_name)),fill=I('black'),geom='point')+
    ggstance::geom_linerangeh()+
    scale_x_continuous(name='Log10 rTE')+
    scale_y_continuous(name='Gene')+
    geom_point()+
    theme_bw()
print(ggMarginal(plot,type='histogram',margins='x'))
dev.off()
message(normalizePath('../plots/modelling/indiv_rTEs_linerange_detonly.pdf'))





pdf('../plots/modelling/indiv_ldeg_linerange.pdf');
plot<-
allsummary%>%filter(parameter=='ldeg')%>%
  mutate(nondet_deg=sd>2.5)%>%
  arrange(nondet_deg)%>%
  qplot(data=.,x=mean,xmin=`2.5%%`,xmax=`97.5%`,y=as.numeric(as_factor(gene_name)),geom='blank')+
    ggstance::geom_linerangeh()+
    scale_x_continuous(name='Log10 deg')+
    scale_y_continuous(name='Gene')+
    theme_bw()
print(plot)
dev.off()
message(normalizePath('../plots/modelling/indiv_ldeg_linerange.pdf'))

##Now the hierarchical model of ldeg...
hmu_ldeg <- hallsummary%>%filter(parameter%in%c('hmu_ldeg'))%>%.$mean
hsig_ldeg <- hallsummary%>%filter(parameter%in%c('hsig_ldeg'))%>%.$mean%>%exp
groupparamdata <- tibble(mean=hmu_ldeg,`2.5%`=hmu_ldeg-(2*hsig_ldeg),`97.5%`=hmu_ldeg+(2*hsig_ldeg),gene_name='Group')
#Now on the joint model
pdf('../plots/modelling/hierarch_degs_linerange.pdf');
plot<-
hallsummary%>%
  filter(parameter%in%c('ldeg'))%>%
  mutate(nondet_deg=sd>2.5)%>%
  arrange(nondet_deg)%>%
  qplot(data=.,x=mean,xmin=`2.5%`,xmax=`97.5%`,y=as.numeric(as_factor(gene_name)),geom='blank')+
    ggstance::geom_linerangeh()+
    ggstance::geom_linerangeh(data=groupparamdata)+
    scale_x_continuous(name='Log deg')+
    scale_y_continuous(name='Gene')+
    facet_grid(scale='free',ifelse(gene_name=='Group','Group','Indiv') ~ . )+
    theme_bw()
print(plot)
dev.off()
message(normalizePath('../plots/modelling/hierarch_degs_linerange.pdf'))

library(ggExtra)







#Generally Genes have either a determined rTE or a determined degredation rate - not both.
allsummary%>%
  mutate(nondet=sd>4.5)%>%
  group_by(gene_name)%>%
  filter(parameter%in%c('lrTE','ldeg'))%>%
  summarise(determined = paste0(collapse=' ',parameter[!nondet]))%>%group_by(determined)%>%tally
  


pdf('../plots/modelling/indiv_rTEs_linerange.pdf');
plot<-
allsummary%>%filter(parameter=='lrTE')%>%
  mutate(nondet_rTE=sd>2.5)%>%
  filter(!nondet_rTE)%>%
  qplot(data=.,x=log10(exp(mean)),xmin=log10(exp(`2.5%`)),xmax=log10(exp(`97.5%`)),y=as.numeric(as_factor(gene_name)),geom='blank')+
    ggstance::geom_linerangeh()+
    scale_x_continuous(name='Log10 rTE')+
    scale_y_continuous(name='Gene')+
    theme_bw()
print(plot)
dev.off()
message(normalizePath('../plots/modelling/indiv_rTEs_linerange.pdf'))



#Comparing similtaneous and individual fits
pdf('tmp.pdf');
print(allsummary%>%filter(parameter=='ldeg')%>%select(gene_name,mean)%>%
  left_join(allsummarysim%>%filter(parameter=='ldeg')%>%select(gene_name,mean),by='gene_name')%>%
  qplot(data=.,x=mean.x,y=mean.y))
dev.off()
message(normalizePath('tmp.pdf'))



