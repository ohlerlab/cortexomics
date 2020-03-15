
library(purrr)
library(ggExtra)
library(ggpubr)
library(tidyverse)
library(Rcpp)
library(doMC)
library(rstan)
#BiocManager::install('rstan')
library(tidyverse)
library(magrittr)
library(data.table)
library(stringr)
library(magrittr)
library(splines)
library(parallel)

source('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/R/cortexomics_myfunctions.R')

################################################################################
########Create trajectory plots from stan fits
################################################################################
  
testgeneset <- c("Acadvl", "Ak1", "Asna1", "Cald1", "Cst3", "Ctnnd2", "Cul2",
"Dclk1", "Dpysl3", "Epb41l1", "Flna", "Hspa12a", "Igsf21", "Nos1",
"Orc3", "Phactr4", "Rap1b", "Rasa3", "Rps6ka5", "Satb2", "Srcin1",
"Stx1b", "Tbc1d24", "Trmt1l", "Zc2hc1a", "Zmym3")

#plotting individual fits
gene_ind=1 # legacy, indiv objects hould only have 1 index now
# gnamei <-'Orc3'
# gnamei <- exprdata$gene_name%>%sample(1)
# gnamei <- increasinggenes%>%sample(1)
# gnamei <- decreasinggenes%>%sample(1)
gnamei <- testgeneset%>%sample(1)
fitfiles <- Sys.glob(paste0('stanfits/*nonh_vs_lin.stanfit.rds'))
# for(gnamei in testgeneset){



################################################################################
########reading stan csvs
################################################################################
  
testhierarchfit<-readRDS('../data/testhierarchfit.rds')



################################################################################
########
################################################################################
  

# for(gnamei in 'Satb2'){
gnamei = 'Satb2'

  fitfile<-fitfiles%>%str_subset(paste0(gnamei,'non'))

  stopifnot(length(fitfile)==1)
  stopifnot(fitfile%>%file.exists)

  message('read in model fit')
  fitstem <- here('pipeline/stansamples/*')
  fitcsvs<-Sys.glob(paste0(  fitstem))%T>%{stopifnot(length(.)>0)}
  fit <- %>%lapply()

  ###MOdel fits from individual RDS files
  fit<-fitfile%>%readRDS
  realstanfit <- fit%>%.$result%>%.$kinetic
  realstanfit_lin <- fit%>%.$result%>%.$linear
  fitdata<-fit%>%.$result%>%.$data

  realstanfit%>%summary



  message('Process')
  #get maximim likelihood (kinda) and mean values from our two models
  parsed_ml<-realstanfit%>%get_ml_stanfit
  parsed_ml_lin<-realstanfit_lin%>%get_ml_stanfit
  parsed_ml<-realstanfit%>%get_parsed_summary%>%mutate(val=`50%`)
  parsed_ml_lin<-realstanfit_lin%>%get_parsed_summary%>%mutate(val=`50%`)


  fread(cmd=paste("grep -ve '#' stansamples_RtmpbokGcP_degmodel_hierarch.stan_modsamples_4.csv"),sep=',', nrow=10)%>%apply(2,mean)
  
  #function to pull out the mcmc samples for plotting 
  fit <- realstanfit

  rdata2plot<-fitdata$MS[,,gene_ind]%>%as.data.frame%>%set_colnames(tps)%>%mutate(rep=1:3)%>%gather(time,signal,-rep)

  ml2plot<-parsed_ml%>%filter(parameter=='prot')%>%transmute(gene,time=as_factor(tps[time]),signal=(val))%>%filter(gene==gene_ind)
  ml2plot_lin<-parsed_ml_lin%>%filter(parameter=='prot')%>%transmute(gene,time=as_factor(tps[time]),signal=(val))%>%filter(gene==gene_ind)

  protsampleslin<- get_prot_samples(realstanfit_lin)%>%mutate(sample=factor(sample),signal=value,time=as_factor(tps[time]),model=as_factor('Linear'))%>%filter(gene==gene_ind)
  protsamples<- get_prot_samples(realstanfit)%>%mutate(sample=factor(sample),signal=value,time=as_factor(tps[time]),model=as_factor('Kinetic'))%>%filter(gene==gene_ind)

  # 
  # protsamples_lowrTEmode<- get_prot_samples(realstanfit%>%as.data.frame%>%filter(`rTE[1]`<0.01))%>%mutate(sample=factor(sample),signal=value,time=as_factor(tps[time]),model=as_factor('Kinetic'))%>%filter(gene==gene_ind)
  # protsamples_highrTEmode<- get_prot_samples(realstanfit%>%as.data.frame%>%filter(`rTE[1]`>0.1))%>%mutate(sample=factor(sample),signal=value,time=as_factor(tps[time]),model=as_factor('Kinetic'))%>%filter(gene==gene_ind)
  # 
  assay2model <- c('ribo'='Riboseq Data','total'='RNAseq Data')
  model2color <- c('Kinetic'='red','Linear'='blue','MS Data'='black','rTE: 0'='purple','Riboseq Data'='dark green','RNAseq Data'='purple')
  seqdata <- preddf%>%filter(gene_name==gnamei,assay%in%c('total','ribo'),!is.na(predicted_signal_spline_3))%>%mutate(signal=2^predicted_signal_spline_3,model=assay2model[assay])

  stopifnot(all(assay2model %in% seqdata$model))

  stopifnot(is.charcter(gnamei) )
  stopifnot(c('time','signal','model') %in% colnames(rdata2plot) )
  stopifnot(c('time','signal','model') %in% colnames(rdata2plot) )


  message('Plot')

  #plot showing trajectories and fits with MCMC samples 
  trajectoryplot<-ggplot(rdata2plot%>%mutate(model='MS Data'),aes(color=model,y=log2(signal),x=as.numeric(as_factor(time))))+
    stat_summary(geom='line',fun.y=mean)+
    geom_point()+
    geom_line(size=I(1),data=ml2plot%>%mutate(model='Kinetic'),linetype=1)+
    geom_line(size=I(1),data=ml2plot_lin%>%mutate(model='Linear'),linetype=1)+
    geom_line(alpha=I(0.005),data=protsampleslin,aes(group=sample))+
    geom_line(alpha=I(0.005),data=protsamples,aes(group=sample))+
    scale_x_continuous(name='Stage',labels=tps)+
    scale_y_continuous(name='Log2 LFQ / Log2 Normalized Counts')+
    facet_grid(scale='free',ifelse(model%in%names(assay2model),'Seq Data','MS') ~ . )+
    scale_color_manual(values = model2color)+
    theme_bw()+
    geom_line(data=seqdata)+
    ggtitle(label = str_interp('Linear vs. Kinetic Model - ${gnamei}'),
      sub="Faded Lines represent samples from Posterior\nSolid Line is Median Value")

    
  # trajectoryplot
  gglayout = c(1,1,1,2,3,4)%>%matrix(ncol=2)

  stopifnot('ldeg[1]' %in% (realstanfit%>%as.data.frame%>%colnames))

  arrangedplot<-gridExtra::grid.arrange(
    trajectoryplot,
    realstanfit%>%as.data.frame%>%.$`lrTE[1]`%>%qplot(bins=50,main=str_interp('Posterior Distribution log rTE\n${gnamei}'))+theme_bw(),
    realstanfit%>%as.data.frame%>%.$`ldeg[1]`%>%qplot(bins=50,main=str_interp('Posterior Distribution - log(degredation / 1.5 days) ${gnamei}'))+theme_bw(),
    realstanfit%>%as.data.frame%>%.$`ldeg[1]`%>%exp%>%qplot(bins=50,main=str_interp('Posterior Distribution - degredation / 1.5 days ${gnamei}'))+theme_bw(),
  layout_matrix=gglayout)

  # gnamei='Satb2'
  trjplotfile<-str_interp('../plots/modelling/stanmodelcomp_${gnamei}.pdf')
  ggsave(file=trjplotfile,arrangedplot,w=12,h=10)
  trjplotfile%>%normalizePath(.)%>%message
  
# }
  here('plots/modelling/stanmodelcomp_Satb2.pdf')



#' sab2 increasing, linear fit, kind of bimodal posterior now.
#' Orc3 weird zig zag in the riboseq that the model fits by just completely disregarding it.... This may be an instance of the riboseq signal itself being weird, rnaseq is slightly better...
#' being fucked up
#' Ewsr1 - another fucking weird, zig-zaggy riboseq signal, that the model just chucks away.
#' 'Ap1' nice clear signal for rTE, increasing, didn't go this time
#' 'Rasa3' - looks like shit. The linear fit is one step off - as if the MS just takes a single stage to pick up...
#' Rap1b looks fucked' - looks like shit. In this case the degredation model works great but only by fitting zero degredation to what is clearly just unrelated Ribo and RNA - ribo goes way down at P0, protein doesn't.
#' Zc2hc1a  looks kind of like an okay fit, although rTE is a little bimodal, still
#' #' Rap1b looks fucked' - looks like shit. In this case the degredation model works great but only by fitting zero degredation to what is clearly just unrelated Ribo and RNA - ribo goes way down at P0, protein doesn't.
#' Acadvl - rTE sucks, cos again, the Riboseq doesn't really bear much relationship to the mS...
#' Cul2 - once again, bananas!! wtf is happening here. RNAseq maaaayyyyyybe looks like it mtches better but... not that much
#' Zmym3 - again something weird...
#' Trmt1l - wtf???? So little correlation between the MS and the other two... I don't know what the fuck is happening...
#' 'Hspa12a' - is increasing, and looks correct.
#' Rps6ka5 - also increasing, also looks correct
#' Stx1b increasing - a case in which it still kind of looks like the RNAseq fits better, just like, compressed a bit - making me wonder if tp specific factors + linear would work best
#' Srcin1 Same phenomenon - I wonder if these guys tend to have like, 
#' Ak1 Same phenomenon - again
#' Tbc1d24 Same phenomenon - again
#' Cald1 decreases - again it's like, linear fit pretty good, but just underestimates the effect that the Riboseq has on the MS
#' Phactr4 too
