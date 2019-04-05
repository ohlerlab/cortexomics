#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(limma))
suppressMessages(library(statmod))
message('...done')


argv <- c(
  transformdexprfile=file.path('exprdata/transformed_data.txt'),
  designmatrixfile=file.path('exprdata/designmatrix.txt'),
  foldchangesfile='exprdata/limma_fold_changes.txt'
)

argv <- if(!interactive()) commandArgs(trailingOnly=TRUE) else argv

for(i in names(argv)) assign(i, argv[i])


# save.image();stop('imagesaved')

#and export
dir.create('exprdata',showWarnings = FALSE)
exprtbl <- read_tsv(transformdexprfile) 
exprtbl %<>% select(gene_name, everything())
assert_that(map_chr(exprtbl,class)[1] == 'character')
assert_that(all(map_chr(exprtbl,class)[-1] == 'numeric'))

exprmatrix <- exprtbl  %>% { set_rownames(as.matrix(.[,-1]),.[[1]]) }


designmatrix <- read_tsv(designmatrixfile)

levels(designmatrix$assay) <- c('total','ribo','MS')

design = model.matrix( ~ time + assay + time:assay , designmatrix, xlev = list(assay = c('total','ribo','MS')) )

design%>%colnames

limmafit = limma::lmFit(exprmatrix,design=design)

bayesfit = limma::eBayes(limmafit,trend=TRUE, robust=TRUE)

coefs<-bayesfit$p.value%>%colnames

allcoefs<-coefs%>%setNames(.,.)%>%
  map(~topTable(bayesfit,number=nrow(exprmatrix),coef=.,confint=0.95)%>%
  {cbind(gene=rownames(.),.)})%>%
  bind_rows(.id='coefficient')%>%
  as_data_frame

allcoefs %>%  write_tsv(paste0(foldchangesfile,'full.txt'))

coeffstoexport <- allcoefs %>%
  # filter(! coefficient=='(Intercept)')%>%
  filter(coefficient%>%str_detect('time'))
  identity

stopifnot(n_distinct(coeffstoexport$coefficient)>11)

coeffstoexport%>%
  select(gene_name=gene,logFC,coefficient)%>%
  spread(coefficient,logFC)%>%
  write_tsv(paste0(foldchangesfile)%T>%message)
  
  coefs<-bayesfit$p.value%>%colnames

allcoefs<-coefs%>%setNames(.,.)%>%
  map(~topTable(bayesfit,number=nrow(exprmatrix),coef=.,confint=0.95)%>%
  {cbind(gene=rownames(.),.)})%>%
  bind_rows(.id='coefficient')%>%
  as_data_frame


###Look at fits, plot with the data
gnames <- exprmatrix%>%rownames




####Now get data and predictions:
exprdata <- exprmatrix%>%
  set_colnames(designmatrix$dataset)%>%
  as.data.frame%>%
  rownames_to_column('gene_name')%>%
  gather(dataset,signal,-gene_name)%>%
  left_join(designmatrix)

exprdata%>%group_by(gene_name)%>%nest%>%slice(1)%>%unnest

get_limmafit_predvals <- function(limmafit,designmatrix){
  (limmafit$coefficients %*% t(limmafit$design))%>%
  set_colnames(designmatrix$dataset)%>%
  as.data.frame%>%
  rownames_to_column('gene_name')%>%
  gather(dataset,predicted_signal,-gene_name)%>%
  left_join(designmatrix)%>%
  distinct(gene_name,time,assay,.keep_all = TRUE)
}


get_limmafit_stdevs <- function(limmafit,designmatrix){
  (((limmafit$stdev.unscaled)^2) %*% t(limmafit$design))%>%
  set_colnames(designmatrix$dataset)%>%
  as.data.frame%>%
  rownames_to_column('gene_name')%>%
  gather(dataset,var_signal,-gene_name)%>%
  mutate(sd_signal = sqrt(var_signal))%>%
  left_join(designmatrix)%>%
  distinct(gene_name,time,assay,.keep_all = TRUE)
}

predictedvals <- get_limmafit_predvals(limmafit,designmatrix)
predictedstdevs <- get_limmafit_stdevs(limmafit,designmatrix)
predictedvals%<>%left_join(predictedstdevs)
exprdata %<>% left_join(predictedvals%>%filter(rep==1)%>%select(-rep))

#and view this
exprdata%>%filter(gene_name==gene_name[2])%>%
{ggplot(data=.,aes(x=as.numeric(as.factor(time)),y=signal,color=assay))+
    geom_point(size=3)+
    geom_line(data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal,color=assay,fill=assay,ymax=signal+sd_signal,ymin=signal-sd_signal),size=3)
  geom_ribbon(data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=signal,color=assay,fill=assay,ymax=signal+sd_signal,ymin=signal-sd_signal),size=3)
}



########Now, let's compare models of different complexity levels
designmatrix$ribo <- designmatrix$assay %in% c('ribo','MS')
designmatrix$MS <- designmatrix$assay %in% c('MS')

limmafit_nTE_change = limma::lmFit(exprmatrix,
                        design=model.matrix( ~ time + ribo + MS, designmatrix)
)

limmafit_TE_change = limma::lmFit(exprmatrix,
                                   design=model.matrix( ~ time + ribo*time + MS, designmatrix)
)

limmafit_MS_change = limma::lmFit(exprmatrix,
                                  design=model.matrix( ~ time + ribo + MS*time, designmatrix)
)

limmafit_full = limma::lmFit(exprmatrix,
                                  design=model.matrix( ~ time + ribo*time + MS*time, designmatrix)
)

rlimmafit_nMS_change = limma::lmFit(exprmatrix[,designmatrix%>%filter(!assay=='total')%>%pluck('dataset')],
                                   design=model.matrix( ~ time + MS, designmatrix%>%filter(!assay=='total'))
)

rlimmafit_MS_change = limma::lmFit(exprmatrix[,designmatrix%>%filter(!assay=='total')%>%pluck('dataset')],
                                   design=model.matrix( ~ time + MS + MS:time, designmatrix%>%filter(!assay=='total'))
)

rnalimmafit_nMS_change = limma::lmFit(exprmatrix[,designmatrix%>%filter(!assay=='ribo')%>%pluck('dataset')],
                                    design=model.matrix( ~ time + MS, designmatrix%>%filter(!assay=='ribo'))
)

rnalimmafit_MS_change = limma::lmFit(exprmatrix[,designmatrix%>%filter(!assay=='ribo')%>%pluck('dataset')],
                                   design=model.matrix( ~ time + MS + MS:time, designmatrix%>%filter(!assay=='ribo'))
)


exprdata$ribo <- exprdata$assay %in% c('ribo','MS')
exprdata$MS <- exprdata$assay %in% c('MS')

#variable for playing with pcas
allcoefftbl <- limmafit$coefficients%>%as.data.frame%>%rownames_to_column('gene_name')%>%select(-`(Intercept)`)
allcoefftbl <- limmafit_nTE_change$coefficients%>%as.data.frame%>%rownames_to_column('gene_name')%>%select(-`(Intercept)`)
allcoefftbl <- limmafit_TE_change$coefficients%>%as.data.frame%>%rownames_to_column('gene_name')%>%  select(-`(Intercept)`)
allcoefftbl <- limmafit_full$coefficients%>%as.data.frame%>%rownames_to_column('gene_name')%>%select(-`(Intercept)`)
allcoefftbl <- rlimmafit_MS_change$coefficients%>%as.data.frame%>%rownames_to_column('gene_name')%>%select(-`(Intercept)`)





getlimmapredictions <- function(modfit,modname,designmatrix){
  predictedvals <- get_limmafit_predvals(modfit,designmatrix)
  predictedstdevs <- get_limmafit_stdevs(modfit,designmatrix)
  predictedvals%<>%left_join(predictedstdevs)
  predictedvals%<>%select(-matches('var_'))
  predictedvals%<>%filter(rep==1)
  predictedvals%<>%select(-rep)
  colnames(predictedvals) %<>%str_replace('predicted_signal',paste0('predicted_signal_',modname))
  colnames(predictedvals) %<>%str_replace('sd_signal',paste0('sd_signal_',modname))
  predictedvals
}
  
exprdata_withpred <- exprdata %>%select(-matches('predicted|_signal'))%>%
  left_join(getlimmapredictions(limmafit_nTE_change,'TE_change',designmatrix))%>%
  left_join(getlimmapredictions(limmafit_full,'full',designmatrix))

exprdata_withpred%<>% mutate(residual_TE_change = signal - predicted_signal_TE_change)

exprdata_withpred%<>% mutate(residual_full = signal - predicted_signal_full)

rssdf <- exprdata_withpred%>%filter(assay=='MS')%>%group_by(gene_name)%>%summarise(
  rss_te = sum(residual_TE_change^2,na.rm=T),
  rss_full = sum(residual_full^2,na.rm=T))%>%
  mutate(rssdiff = log(rss_full) - log(rss_te))%>%
  arrange((rssdiff))

rssdf%>%head

#and view this
i=2
exprdata_withpred%>%filter(gene_name==rssdf$gene_name[i])%>%
{ggplot(data=.,aes(x=as.numeric(as.factor(time)),y=signal,color=assay))+
    geom_point(size=I(8))+
    geom_line(data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_TE_change,color=assay,fill=assay),size=3)+
    geom_line(linetype=2,data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_full,color=assay,fill=assay),size=3)+
    ggtitle(paste0('Comparing steady state model (solid) to optimal fit (dashed):\n',rssdf$gene_name[i]))
  # geom_ribbon(data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=signal,color=assay,fill=assay,ymax=signal+sd_signal,ymin=signal-sd_signal),size=3)
}

#####Now the ribo only model

exprdata_withpred_ronly <- exprdata %>%select(-matches('predicted|_signal'))%>%
  left_join(getlimmapredictions(rlimmafit_nMS_change,'nMS_change',designmatrix%>%filter(!assay=='total'))) %>%
  left_join(getlimmapredictions(rlimmafit_MS_change,'MS_change',designmatrix%>%filter(!assay=='total')))

exprdata_withpred_ronly%<>% mutate(residual_nMS_change = signal - predicted_signal_nMS_change)
exprdata_withpred_ronly%<>% mutate(residual_MS_change = signal - predicted_signal_MS_change)

#and the rna only for comparison
exprdata_withpred_rnaonly <- exprdata %>%select(-matches('predicted|_signal'))%>%
  left_join(getlimmapredictions(rnalimmafit_nMS_change,'nMS_change',designmatrix%>%filter(!assay=='ribo'))) %>%
  left_join(getlimmapredictions(rnalimmafit_MS_change,'MS_change',designmatrix%>%filter(!assay=='ribo')))

exprdata_withpred_rnaonly%<>% mutate(residual_nMS_change = signal - predicted_signal_nMS_change)
exprdata_withpred_rnaonly%<>% mutate(residual_MS_change = signal - predicted_signal_MS_change)


exprdata_withpred_ronly%>%colnames
rssdf <- exprdata_withpred_ronly%>%filter(assay=='MS')%>%group_by(gene_name)%>%summarise(
  rss_te = sum(residual_nMS_change^2,na.rm=T),
  rssMS_change = sum(residual_MS_change^2,na.rm=T))%>%
  mutate(rssdiff = log(rssMS_change) - log(rss_te))%>%
  arrange((rssdiff))

rssdf%>%head

#and view this
i=2

exprdata_withpred_ronly%>%filter(gene_name==rssdf$gene_name[i])%>%
{ggplot(data=.,aes(x=as.numeric(as.factor(time)),y=signal,color=assay))+
    geom_point(size=I(8))+
    geom_line(data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_nMS_change,color=assay,fill=assay),size=3)+
    geom_line(linetype=2,data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_MS_change,color=assay,fill=assay),size=3)+
    ggtitle(paste0('Comparing steady state model (solid) to optimal fit (dashed):\n',rssdf$gene_name[i]))
  # geom_ribbon(data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=signal,color=assay,fill=assay,ymax=signal+sd_signal,ymin=signal-sd_signal),size=3)
}
#coefs same after eBayes
ebrlimmafit_MS_change <- eBayes(rlimmafit_MS_change)


#We specify a 1 column matrix, whose elements select the time:MS interactions
coefcomb = str_detect(colnames(rlimmafit_MS_change$coefficients),'time.*MS')
?contrasts.fit

mstimetoptable<-rlimmafit_MS_change%>%
  contrasts.fit(coefficients = coefcomb)%>%
  eBayes%>%
  topTable(number=1e12)
  
mstimetoptable$adj.P.Val%>%hist(40,main='Pvalue Histogram: Time,MS specific effects')

markergenes <- '/fast/groups/ag_ohler/work/dharnet_m/cortexomics/ext_data/marker_genes.tsv'%>%fread%>%.$gene_name

mstimetoptable%>%as.data.frame%>%rownames_to_column('gene_name')%>%mutate(rank = paste0(rank(1:nrow(.)),'/',nrow(.)))%>%
    filter(tolower(gene_name)%in%tolower(markergenes))%>%
    select(gene_name,adj.P.Val,rank)

mstimetoptable%>%head(30)
i=26
# exprdata_withpred_ronly%>%filter(gene_name==mstimetoptable%>%rownames%>%tail(26)%>%head(1))%>%
satb2lineardevplot<-exprdata_withpred_ronly%>%filter(gene_name=='Satb2')%>%
{ggplot(data=.,aes(x=as.numeric(as.factor(time)),y=signal,color=assay))+
    geom_point(size=I(8))+
    geom_line(data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_nMS_change,color=assay,fill=assay),size=3)+
    geom_line(linetype=2,data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_MS_change,color=assay,fill=assay),size=3)+
    ggtitle(paste0('Comparing steady state model (solid) to optimal fit (dashed):\n',.$gene_name[1]))+
    theme_bw()
  # geom_ribbon(data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=signal,color=assay,fill=assay,ymax=signal+sd_signal,ymin=signal-sd_signal),size=3)
}

cairo_pdf('plots/modelling/satb2lineardevplot.pdf'%T>%{normalizePath(.)%>%message})
print(satb2lineardevplot,);dev.off()

rlimmafit_MS_change$coefficients%>%as.data.frame%>%
    rownames_to_column('gene_name')%>%select(gene_name,matches('time'))%>%gather(coef,val,-gene_name)%>%
    mutate(ms=str_detect(coef,'MSTRUE'))%>%
    mutate(time=str_replace(coef,'time','')%>%str_replace(':MSTRUE',''))%>%
    filter(gene_name==unique(gene_name)[2])%>%
    qplot(data=.,x=time,y=val,color=ms)+geom_point(size=10)

####Calculate how much mass spec variance is explained by the 
vardf<-exprdata_withpred_ronly%>%
  filter(gene_name%in%gene_name)%>%group_by(gene_name)%>%
  filter(assay=='MS')%>%  
  filter(rep==1)%>%arrange(gene_name)
  
vardfrna<-exprdata_withpred_rnaonly%>%
  filter(gene_name%in%gene_name)%>%group_by(gene_name)%>%
  filter(assay=='MS')%>%  
  filter(rep==1)%>%arrange(gene_name)

#for one gene
vardf%>%filter(gene_name=='Satb2')%>%colnames
anova(lm(data=vardf%>%filter(gene_name=='Satb2'),predicted_signal_MS_change ~ predicted_signal_nMS_change))$"Sum Sq"%>%{.[1]/sum(.)}

#now do this for all genes
varexpldf<-vardf%>%group_by(gene_name)%>%
  nest%>%
  mutate(varexplained = map(.$data,~lm(data = ., predicted_signal_MS_change ~ predicted_signal_nMS_change)%>%anova%>%.$"Sum Sq"%>%{.[1]/sum(.)} ))%>%
  identity%>%
  select(-data)%>%
  unnest

anova(lm(data=vardfrna%>%filter(gene_name=='Satb2'),predicted_signal_MS_change ~ predicted_signal_nMS_change))$"Sum Sq"%>%{.[1]/sum(.)}



cor(
  vardf%>%filter(gene_name=='Satb2')%>%.$predicted_signal_MS_change,
  vardf%>%filter(gene_name=='Satb2')%>%.$predicted_signal_nMS_change  
)

cor(
  vardfrna%>%filter(gene_name=='Satb2')%>%.$predicted_signal_MS_change,
  vardfrna%>%filter(gene_name=='Satb2')%>%.$predicted_signal_nMS_change  
  )


#and for rna
varexpldfrna<-vardfrna%>%group_by(gene_name)%>%
  nest%>%
  mutate(varexplained = map(.$data,~lm(data = ., predicted_signal_MS_change ~ predicted_signal_nMS_change)%>%anova%>%.$"Sum Sq"%>%{.[1]/sum(.)} ))%>%
  identity%>%
  select(-data)%>%
  unnest


#These pretty much look the same
par(mfrow=c(2,1))
varexpldf$varexplained%>%hist(40,main='Proportion of Variance in Mass Spec Fold Changes\nExplained by RiboSeq')
varexpldfrna$varexplained%>%hist(40,main='Proportion of Variance in Mass Spec Fold Changes\nExplained by RNA')

qplot(varexpldf$varexplained,varexpldfrna$varexplained,color=varexpldf$gene_name %in% diffTEnames,alpha=I(0.1))

diffTEnames<-'/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/xtail/xtail_P0.txt'%>%read_tsv%>%filter(adj_p_value < 0.05)%>%
  select(gene_id=feature_id)%>%
  left_join(by='gene_id',read.table(header=T,'/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/ids.txt'))%>%
  .$gene_name

#plot the log time spec coefficients of the prot with the mass spec
#often roughly linear!
i=9

rlimmafit_MS_change$coefficients%>%as.data.frame%>%rownames_to_column('gene_name')%>%select(gene_name,matches('time'))%>%gather(coef,val,-gene_name)%>%
  mutate(ms=str_detect(coef,'MSTRUE'))%>%
  mutate(time=str_replace(coef,'time','')%>%str_replace(':MSTRUE',''))%>%
  filter(gene_name==unique(gene_name)[i])%>%
  mutate(ms=ifelse(ms,'mass_spec_coef','ribo_coef'))%>%
  select(ms,val,time)%>%
  spread(ms,val)%>%
  qplot(data=.,x=ribo_coef,y=mass_spec_coef,size=time)



i=i+1



#function simulating ms given a degredation constant, a constant value for rTE, starting ms, the riboseq,
simulate_ms <- function(ldeg,rTE,ribo,ms0){
  deg = exp(ldeg)
  ms = rep(NA,length(ribo))
  ms[1] = ms0
  for (i in 2:length(ms)){
    ms[i] = ms[i-1]+ (rTE*ribo[i])-(ms[i-1]*deg) 
  }
  ms
}


toyfit <- data_frame(g=1,time=as.factor(rep(1:5)),ribo=(c(1,10,29,50,100)))%>%
    mutate(ms=simulate_ms(ldeg=log(0.8),rTE=10,ms0=10,ribo=ribo))%>%
    mutate(lribo=log(ribo),lms=log(ms))%>%
    select(-lribo,-lms)%>%
    gather(assay,signal,-time,-g)%>%
    mutate(assay=factor(assay,unique(assay)))%>%
    group_by(g)%>%nest%>%mutate(fit=map(data,~lm(data=.,signal ~ time + time:assay + assay)))%>%
    mutate(slopes = map(fit, ~ .$coefficients %>%.[c(-1,-6)]%>%{.[1:4]/.[5:8]} ))


rstan::stan(model_code = "

            
            
",data=stan_data)






pcafit <- princomp(allcoefftbl%>%select(matches('time')))

mostloadedgenespca1<-allcoefftbl$gene_name[pcafit$scores%>%order%>%tail(10)]





####OKay so I should figure out whether i am doing confidence intervals on the linear fit correctly - no. I need to do msm
####I also need to implement splines for riboseq read removal
####And probably also simply implement linear modelling with ribo and MS alone, to see how well they fit
#####Also, why are is rTE consistently lower???
####Do I model the MS and ribo as responses both, or with MS as response and ribo as prector?
####The problem is, genes differ in how exponential their riboseq looks
####What if I do something like, subtract the the exponential component from the Riboseq
####We could limit the analysis only to increasing genes but that then suffers from an issue - we don't measure increase, only relative increase....
###I need to think about how I'm transforming thins... rlog might be best.


##Check discrepency in fold changes
library(data.table)

fread('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/ids.txt')




####Try to illustrate fit for MS being better in Riboseq than for RNAseq
# exprdata_withpred_ronly%>%filter(gene_name==mstimetoptable%>%rownames%>%tail(26)%>%head(1))%>%
satb2lineardevplotrna<-exprdata_withpred_rnaonly%>%filter(gene_name=='Satb2')%>%
  filter(assay=='MS')%>%
  {ggplot(data=.,aes(x=as.numeric(as.factor(time)),y=signal,color=assay))+
   # geom_point(size=I(2))+
    geom_line(linetype=2,data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_nMS_change,color=assay,fill=assay),size=1.5)+
    geom_line(linetype=1,data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_MS_change,color=assay,fill=assay),size=1.5)+
    ggtitle(paste0('Mass Spec Fit - RNA based Model:\n',.$gene_name[1]))+
    theme_bw()
  # geom_ribbon(data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=signal,color=assay,fill=assay,ymax=signal+sd_signal,ymin=signal-sd_signal),size=1.5)
}

satb2lineardevplotrna_nuc<-exprdata_withpred_rnaonly%>%filter(gene_name=='Satb2')%>%
  filter(assay=='total')%>%
  {ggplot(data=.,aes(x=as.numeric(as.factor(time)),y=signal,color=assay))+
      # geom_point(size=I(2))+
      geom_line(linetype=2,data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_nMS_change,color=assay,fill=assay),size=1.5)+
      geom_line(linetype=1,linetype=2,data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_MS_change,color=assay,fill=assay),size=1.5)+
      ggtitle(paste0('RNA Fit - RNA based Model:\n',.$gene_name[1]))+
      theme_bw()
    # geom_ribbon(data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=signal,color=assay,fill=assay,ymax=signal+sd_signal,ymin=signal-sd_signal),size=1.5)
  }

satb2lineardevplot<-exprdata_withpred_ronly%>%filter(gene_name=='Satb2')%>%
  filter(assay=='MS')%>%
  {ggplot(data=.,aes(x=as.numeric(as.factor(time)),y=signal,color=assay))+
         # geom_point(size=I(2))+
      geom_line(linetype=2,data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_nMS_change,color=assay,fill=assay),size=1.5)+
      geom_line(linetype=1,data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_MS_change,color=assay,fill=assay),size=1.5)+
      ggtitle(paste0('Mass Spec Fit - Ribo based Model:\n',.$gene_name[1]))+
      theme_bw()
    # geom_ribbon(data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=signal,color=assay,fill=assay,ymax=signal+sd_signal,ymin=signal-sd_signal),size=1.5)
  }
satb2lineardevplotnuc<-exprdata_withpred_ronly%>%filter(gene_name=='Satb2')%>%
  filter(assay=='ribo')%>%
  {ggplot(data=.,aes(x=as.numeric(as.factor(time)),y=signal,color=assay))+
      # geom_point(size=I(2))+
      geom_line(linetype=2,data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_nMS_change,color=assay,fill=assay),size=1.5)+
      geom_line(linetype=1,data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=predicted_signal_MS_change,color=assay,fill=assay),size=1.5)+
      ggtitle(paste0('Ribo Fit - Ribo based Model:\n',.$gene_name[1]))+
      theme_bw()
    # geom_ribbon(data=filter(.,rep==1),aes(x=as.numeric(as.factor(time)),alpha=I(.5),y=signal,color=assay,fill=assay,ymax=signal+sd_signal,ymin=signal-sd_signal),size=1.5)
  }

library(ggpubr)
compfitplot <- '~/projects/cortexomics/plots/modelling/rna_vs_ribofit_satb2.pdf'
pdf(compfitplot,w=12,h=12)
ggpubr::ggarrange(ncol=2,nrow=2,plotlist = list(satb2lineardevplotrna,satb2lineardevplot,satb2lineardevplotrna_nuc,satb2lineardevplotnuc))
dev.off()
compfitplot%>%normalizePath%>%message


