library(purrr)
library(ggExtra)
library(ggpubr)
library(tidyverse)
library(Rcpp)
library(rstan)



#function simulating ms given a degredation constant, a constant value for rTE, starting ms, the riboseq,
simulate_data <- function(ldeg,rTE,ribo,prot0,ms_sd=1,n_reps=K){
  deg = exp(ldeg)
  prot = rep(NA,tps)
  prot[1] = prot0
  for (i in 2:length(prot)){

    prot[i] = rTE*ribo[i] + (prot[i-1]*(1-deg)) 
  }
  prot

  # MS=replicate(n_reps,{exp(rnorm(n=seq_along(prot),mean=logprot),sd=ms_sd))})
  MS=replicate(n_reps,{(rnorm(n=seq_along(prot),mean=prot,sd=ms_sd))})

  data_frame(ldeg=ldeg,rTE=rTE,ribo=list(ribo),prot0=prot0,ms_sd=ms_sd,MS=list(MS))
}

#some basic patterns for the synthesis change pattern
ribopatterns = list(
  increasing=c(100,100,1000,1000,8000),
  stable=c(200,200,200,200,200),
  decreasing=c(3000,3000,100,100,50)
)

#constants describing data shape
tps = 5
n_genes = 1
K = 3

#mostly increasing, a few stable/decreasing genes
ribopatinds <- sample(1:3,n_genes,rep=T,prob=c(.9,.05,.05))
ribopatinds <- c(1)

ribopatinds[1]<-1
message(capture.output(table(ribopatinds))%>%paste0(collapse='\n'))
#choose a pattern
ribo = ribopatterns[ribopatinds]%>%simplify2array()%>%t
ribo = ribo + rpois(length(ribo),10)
#simulaterTE values
rTEs = rnorm(n_genes,1000,10)
#simulate degredations based on a half life distribution with median 48, sd~30 in the log10 scale 
k= - log(0.5) / (10^rnorm(n_genes,log10(48),log10(30)))
degs <- 1-exp(-k*36)
#simualte the ms0 starting values
ms0=ribo[,1]*exp(rnorm(n_genes,0,3))*rTEmu
#
simdata <- map_df(1:n_genes,~simulate_data(
  log(degs[.]),
  rTEs[.],
  ribo[.,T],
  ms_sd=10,
  prot0=(ribo[,1]*rTE)/2)
)

ms_array <- simdata$MS%>%simplify2array%>%aperm(c(1,2,3))
ribo_mat <- simdata$ribo%>%simplify2array()
simdata$prot0

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#get the names of some increasing rna dudes
exprdata%>%
  group_by(gene_name)%>%
  filter(rep==1,assay=='total')%>%
  filter(signal[5] > 3+signal[1])%>%.$gene_name%>%
  sample(10)

#we need a 
library(zoo)
# 
tps <- c('E13','E145','E16','E175','P0')
gnamei <- 'Satb2'
realdata <- list(
   G=1,T=length(tps),K=3,
   # gene_name=gnamei,
   MS=exprdata%>%filter(gene_name==gnamei)%>%filter(assay=='MS')%>%.$signal%>%array(dim=c(3,5,1))%>%aperm(c(1,2,3))%>%{2^.},
   ribo =matrix(exprdata%>%filter(gene_name==gnamei)%>%filter(assay=='ribo',rep==1)%>%.$predicted_signal%>%{2^.})
)
matrix(
  exprdata%>%filter(gene_name==gnamei)%>%filter(assay=='ribo',rep==1)%>%.$predicted_signal%>%{2^.}%>%rollmean(k=2)
       )%>%
  
  
# 
# stanfit <- rstan::stan(file='src/Stan/degmodel_simple.stan',data=list(G=n_genes,T=tps,K=K,
#                                                           MS=ms_array,ribo=ribo_mat),
#                        control=list(adapt_delta=0.95,max_treedepth=20),
#                        chains=4,iter=2e3,
#                        # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
#                        verbose=TRUE)

realstanfit <- rstan::stan(file='src/Stan/degmodel_simple.stan',data=realdata,
            control=list(adapt_delta=0.95,max_treedepth=20),
            chains=4,iter=2e3,
            # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
            verbose=TRUE)

#now fit linear model
realstanfit_lin <- rstan::stan(file='src/Stan/degmodel_simple_linear.stan',data=realdata,
                           control=list(adapt_delta=0.95,max_treedepth=20),
                           chains=4,iter=2e3,
                           # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
                           verbose=TRUE)
##Now with fixed values
realstanfit_test <- rstan::stan(file='src/Stan/degmodel_simple_linear.stan',data=realdata,
                               control=list(adapt_delta=0.95,max_treedepth=20),
                               init = list(list(lrTE=array(0,dim=c(1)))),
                               chains=1,iter=1,warmup = 0,
                               # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
                               verbose=TRUE)

stanpars <- colnames(as.data.frame(realstanfit))

realdata$MS

parse_stan_pars<-function(stanpars,indnames=c()){
  parsedpars<-stanpars%>%str_match('([^\\[]+)\\[?(\\d*),?(\\d*)\\]?')%>%as.data.frame%>%
    .[,-1]
  n_inds <- length(colnames(parsedpars))-1

    parsedpars[,-1]%>%.[,]%>%apply(1,function(x){k = keep(x,~ . !='');c(rep(NA,n_inds-length(k)),k)})%>%t%>%
            set_colnames(c('time','gene'))%>%
            as.data.frame%>%
            map_df(.,as.integer)%>%
          mutate(parameter=parsedpars[,1])%>%
      select(parameter,time,gene)%>%
      split(.,seq_len(nrow(.)))

}
#get maximim likelihood (kinda) and mean values from our two models, with these functions
get_ml_stanfit <- function(fit){fit%>%as.data.frame%>%slice(which.max(lp__))%>%t%>%as.data.frame%>%rownames_to_column('par')%>%mutate(ppars=parse_stan_pars(par))%>%unnest%>%select(parameter,val=V1,time,gene)}
get_parsed_summary<-function(fit) fit %>%summary%>%.$summary%>%as.data.frame%>%rownames_to_column('par')%>%mutate(ppars=parse_stan_pars(par))%>%unnest
#get_samples_stanfit <- function(fit){fit%>%as.data.frame%>%slice(which.max(lp__))%>%t%>%as.data.frame%>%rownames_to_column('par')%>%mutate(ppars=parse_stan_pars(par))%>%unnest%>%select(parameter,val=V1,time,gene)}

#get maximim likelihood (kinda) and mean values from our two models
parsed_ml<-realstanfit%>%get_ml_stanfit
parsed_ml_lin<-realstanfit_lin%>%get_ml_stanfit
parsed_ml<-realstanfit%>%get_parsed_summary%>%mutate(val=mean)
parsed_ml_lin<-realstanfit_lin%>%get_parsed_summary%>%mutate(val=mean)
#function to pull out the mcmc samples for plotting 
get_prot_samples<-function(fit) fit %>%as.data.frame%>%select(matches('prot'))%>%mutate(sample=1:nrow(.))%>%gather(par,value,-sample)%>%mutate(ppars=parse_stan_pars(par))%>%unnest%>%
  filter(parameter=='prot')%>%select(time,value,sample)
#e.g.
#realstanfit %>% get_prot_samples


rdata2plot<-realdata$MS%>%as.data.frame%>%set_colnames(tps)%>%mutate(rep=1:3)%>%gather(time,signal,-rep)

ml2plot<-parsed_ml%>%filter(parameter=='prot')%>%transmute(time=as_factor(tps[time]),signal=(val))
ml2plot_lin<-parsed_ml_lin%>%filter(parameter=='prot')%>%transmute(time=as_factor(tps[time]),signal=(val))
parsed_summ%>%filter(gene%in%c(1,NA))%>%filter(parameter=='prot',gene==1)

protsampleslin<- get_prot_samples(realstanfit_lin)%>%mutate(sample=factor(sample),signal=value,time=as_factor(tps[time]),model=as_factor('Linear'))
protsamples<- get_prot_samples(realstanfit)%>%mutate(sample=factor(sample),signal=value,time=as_factor(tps[time]),model=as_factor('Kinetic'))

#plot showing trajectories and fits with MCMC samples 
trajectoryplot<-ggplot(rdata2plot%>%mutate(model='Data'),aes(color=model,y=log2(signal),x=as.numeric(as_factor(time))))+geom_point()+
  geom_line(size=I(1),data=ml2plot%>%mutate(model='Kinetic'),linetype=1)+
  geom_line(size=I(1),data=ml2plot_lin%>%mutate(model='Linear'),linetype=1)+
  geom_line(alpha=I(0.005),data=protsampleslin,aes(group=sample))+
  geom_line(alpha=I(0.005),data=protsamples,aes(group=sample))+
  
  # geom_line(size=I(2),data=get_prot_samples(ml2plot_lin)%>%mutate(model='Linear'),linetype=1)+
  scale_x_continuous(name='Stage',labels=tps)+
  scale_y_continuous(name='Log2 LFQ')+
  scale_color_manual(values = c('Kinetic'='red','Linear'='blue','Data'='black'))+
  theme_bw()+
  ggtitle(label = str_interp('Linear vs. Kinetic Model - ${gnamei}'),sub="Faded Lines represent samples from Posterior")
  

trajectoryplot

gglayout = c(1,1,1,2,3,4)%>%matrix(ncol=2)
arrangedplot<-gridExtra::grid.arrange(
  trajectoryplot,
    realstanfit%>%as.data.frame%>%.$`rTE[1]`%>%qplot(bins=50,main='Posterior Distribution - rTE Satb2')+theme_bw(),
  realstanfit%>%as.data.frame%>%.$`ldeg[1]`%>%qplot(bins=50,main='Posterior Distribution - log(degredation / 1.5 days) Satb2')+theme_bw(),
  realstanfit%>%as.data.frame%>%.$`deg[1]`%>%qplot(bins=50,main='Posterior Distribution - degredation / 1.5 days Satb2')+theme_bw()
,layout_matrix=gglayout)

ggsave(str_interp('plots/modelling/stanmodelcomp_${gnamei}.pdf')%T>%{normalizePath(.)%>%message},arrangedplot,w=12,h=10,)



# 
# stanfit <- rstan::stan(file='src/degmodel.stan',data=list(G=n_genes,T=tps,K=K,
#                         MS=ms_array,ribo=ribo_mat),
#                        control=list(adapt_delta=0.90,max_treedepth=10),
#                        chains=4,iter=2e3,
#                        # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
#                        verbose=TRUE)


# 

samples=stanfit%>%as.data.frame

#Look at how rTE is fitting across genes
simdata$rTE

comp_paramval<-function(param,realval){
summary(stanfit)%>%as.data.frame%>%rownames_to_column('parameter')%>%filter(parameter%>%str_detect(paste0('^',param,'\\['))) %>%
  select(summary.mean,summary.2.5.,summary.97.5.) %>%
  mutate(actual = realval) %>%
  rowwise%>%
  mutate(accurate=between(actual,summary.2.5.,summary.97.5.))%>%
  ungroup%>%
  mutate(type=names(ribopatterns)[ribopatinds])%>%
  mutate(dit = summary.mean-actual)%>%
  mutate(parameter = param)%>%select(parameter,everything())
}
comp_paramval(param='rTE',  realval=simdata$rTE)
comp_paramval('MS0',simdata$prot0)
comp_paramval('deg',exp(simdata$ldeg))


summary(stanfit)%>%as.data.frame%>%rownames_to_column('parameter')%>%select(parameter,summary.Rhat)




#
g2plot=1
allpars <- stanfit%>%as.data.frame%>%colnames
parstoplot <- allpars%>%
  # grep(inv=F,val=T,patt=str_interp('\\[${g2plot}\\]|\\,${g2plot}\\]'))%>%
  grep(inv=T,val=T,patt='prot')
parstoplot%<>%append(allpars%>%{.[!str_detect(.,'\\[')]})
parstoplot%<>%setdiff(c('deg','ms0logratio[1]'))
print(pairs(stanfit,pars=parstoplot))
plot(stanfit,pars=parstoplot%>%str_subset('prot'))


pairs(realstanfit,pars=allpars%>%str_subset('prot|_|tau'))

pdf('tmp.pdf');pairs(stanfit,pars=parstoplot);dev.off()
pdf('tmp.pdf');pairs(stanfit,pars=parstoplot%>%str_subset('prot|_'));dev.off()

samples[['rTE[1]']]%>%hist(40)


#now display relationshiop of prot latent vars to the parameters
allpars <- stanfit%>%as.data.frame%>%colnames
parstoplot <- allpars%>%
  grep(inv=F,val=T,patt=str_interp('\\[${g2plot}\\]|\\,${g2plot}\\]'))
# parstoplot%<>%append(allpars%>%{.[!str_detect(.,'\\[')]})
# parstoplot%<>%setdiff(c('deg','ms0logratio[1]'))
print(pairs(stanfit,pars=parstoplot))


allpars <- stanfit%>%as.data.frame%>%colnames
parstoplot <- allpars%>%
  grep(inv=F,val=T,patt=str_interp('\\[${g2plot}\\]|\\,${g2plot+1}\\]'))
# parstoplot%<>%append(allpars%>%{.[!str_detect(.,'\\[')]})
# parstoplot%<>%setdiff(c('deg','ms0logratio[1]'))
print(pairs(stanfit,pars=parstoplot))


####So... parametrizing in terms of half life doesn't really seem to have helped. Makes everything pretty bimodal in fact.

##back to ldeg parametrization and things ar e

#What does this distribution of degs look like??

#What if i simulate values with the stan model?

map_df(1:10,~simulate_data(log(0.5),100,c(100,200,200,500,500),5e4,n_reps=3))%>%
  .$MS%>%simplify2array%>%aperm(c(1,2,3))%>%
  identity()%>%
  

k = - log(0.5) / (10^rnorm(10e3,log10(48),log10(30)))
degs <- 1-exp(-k*36)
degs = degs - min(degs/2)


MASS::fitdistr(degs,"beta",start=list(shape1=3,shape2=3))

1 - degs = exp(-k*36)
log(1-degs) = -k*36
 - log(1-degs) / 36 = k

rnorm(10e3,log10(48),log10(30))%>%hist(20)

10^(log10(12))

exp(log10(12))*exp(log10(10))

2^4
3^4

hist(degs,50)





library(rstan)
funnel <- stan_demo("funnel", seed = 12345)   # has 5 divergent transitions
pairs(funnel, pars = c("y", "x[1]", "lp__"), las = 1) # below the diagonal
funnel_reparam <- stan_demo("funnel_reparam") # has no divergent transitions
pairs(funnel_reparam, pars = c("y", "x[1]", "lp__"), las = 1) # below the diagonal


#ldeg_rTE parameter did NOT help, it just blows up to minus infinity
#' Looks like scaling of the standard deviation makes things work just fine
#' Now I"ll try putting things on a log scale