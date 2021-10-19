if(!exists('get_dp_standata_withpriors')) source('src/R/Modeling/allgene_sampling_clust.R')

base::source('../cortexomics/src/Archive/R/Functions/rstan_functions.R')
base::source('src/Rprofile.R')

library(tidyverse)

dp_model_ribnorm = rstan::stan_model('src/Stan/mod_proDD_ribnorm.stan')

library(rstan)

dp_model_ribprotnormfix = fix_param(dp_model_ribnorm,vars2fix = c('ribnorm','protnorm'))%>%{f='src/Stan/mod_proDD_ribnorm_ribnormfix.stan';cat(.,file=f);f}%>%stan_model
dp_model_pihalffix = fix_param(dp_model_ribnorm,vars2fix = c('l_pihalf','protnorm'))%>% {f='src/Stan/mod_proDD_ribnorm_lpihalffix.stan';cat(.,file=f);f}%>%stan_model(.)
dp_model_bothfix = fix_param(dp_model_ribnorm,vars2fix = c('l_pihalf','ribnorm','protnorm'))%>%{f='src/Stan/mod_proDD_ribnorm_bothfix.stan';cat(.,file=f);f}%>%stan_model

metainfo<-suppressMessages({read_tsv(here('data/metainfo.tsv'))})

gname2pid = metainfo%>%
  distinct(gene_name,protein_id)%>%
  {safe_hashmap(.[['gene_name']],.[['protein_id']])}

gnm2uid = metainfo%>%
  filter(!is.na(uprotein_id))%>%filter(!uprotein_id=='NA')%>%
  filter(isbest)%>%
  distinct(gene_name,uprotein_id)%>%
  {safe_hashmap(.[['gene_name']],.[['uprotein_id']])}

stopifnot(!gnm2uid[["Satb2"]]=='NA')

################################################################################
########Attempt sampling with multiple genes
################################################################################

edist <- dist(matchedms_mat_rscl)
{
satb2like_uids <- as.matrix(edist)[gnm2uid[['Satb2']],]%>%sort%>%head(100)%>%names%>%sample(10)
flnalike_uids <- as.matrix(edist)[gnm2uid[['Flna']],]%>%sort%>%head(100)%>%names%>%sample(10)
# mgids%T>%{assert_that(length(.)>0,msg='satb2 not in the uids')}

mgids <- c(satb2like_uids,flnalike_uids)
mginits <- bestfitinits%>%map(~ .[mgids])%>%lapply(get_comb_initvals)%>%.$ribo
mgdata = get_dp_standata_withpriors(c(mgids),ribomatrscl=countmats$ribo,ribo_sigma=sigmas$ribo)
mgdata$ribnorm = rep(0,5) %>%array(.,dim=5)
mgdata$protnorm = rep(0,5)%>%array(.,dim=5) 
mgdata$l_pihalf = rep(-4,length(mgids))

mgopt<-rstan::optimizing(dp_model_bothfix,data=mgdata,
                              init=mginits,iter=20e3,verbose=F,as_vector=F)
mginits$ribnorm = rep(0,5) %>%array(.,dim=5)
mginits$protnorm = rep(0,5)%>%array(.,dim=5) 
mginits$l_pihalf = rep(-4,length(mgids))
mgopt_normvar<-rstan::optimizing(dp_model_pihalffix,data=mgdata,
                         init=mginits,iter=20e4,verbose=F,as_vector=F)
mgopt_allvar <- rstan::optimizing(dp_model_ribnorm,data=mgdata,
                  init=mginits,iter=20e4,verbose=F,as_vector=F)
#now we want to plot the different fit types on top of one another and the actual data
# flnanum = which(mgids==gnm2uid[['Satb2']])%T>%{assert_that(length(.)>0,msg='satb2 not in the uids')}
visgene = 'Flna'
flnanum = which(mgids==gnm2uid[['Flna']])
flnanum= sample((1:length(mgids)),1)

ggdf <- mgdata$lMS[flnanum,]%>%enframe%>%separate(name,c('time','model','rep'))%>%
    mutate(time=as.numeric(as.factor(time)))%>%
    mutate(model='data')%>%
    rbind(data.frame(time=1:5,model='linear',rep=NA,value=mgopt$par$prot[flnanum,]))%>%
    rbind(data.frame(time=1:5,model='linear_normvar',rep=NA,value=mgopt_normvar$par$prot[flnanum,]))%>%
    rbind(data.frame(time=1:5,model='nonlinear_normvar',rep=NA,value=mgopt_allvar$par$prot[flnanum,]))

ggpubr::ggarrange(plotlist=list(
ggdf%>%ggplot(data=.,aes(x=time,y=value,color=model,group=model))+
    geom_line(data=ggdf%>%group_by(time,model)%>%summarise(value=mean(value)))+
  geom_point()+ggtitle(str_interp('model comparison ${mgids[flnanum]}')),

data.frame(x=mgopt_normvar$par$ribnorm)%>%qplot(data=.,y=x,x=seq_along(x),main='ribonorm')
))

}

#lolwhat mgids[flnanum]

###Goal get the best fit norm values for FLna like genes, Then fix these values, and optimize on everything, see if our fit is better
#first let's try the blank linear fit on these.




#use the individual estimates to get joint estimates
mgids <- c(satb2like_uids,flnalike_uids)
mginits <- bestfitinits%>%map(~ .[mgids])%>%lapply(get_comb_initvals)


mgsampling<-rstan::optimizing(dp_model_ribnorm,data=get_dp_standata_withpriors(c(mgids),ribomatrscl=countmats$ribo,ribo_sigma=sigmas$ribo),
                              init=mginits$ribo,iter=10)

mgsampling<-rstan::sampling(dp_model_ribnorm,data=get_dp_standata_withpriors(c(mgids),ribomatrscl=countmats$ribo,ribo_sigma=sigmas$ribo),chains=4,
                            init=function(){mginits$ribo},control=list(adapt_delta=.98,max_treedepth=15),iter=2e3)
