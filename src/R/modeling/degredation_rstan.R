library(purrr)
library(ggExtra)
library(ggpubr)
library(tidyverse)
library(Rcpp)
library(rstan)

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

#increasing genes
gnamesi <- c('Satb2')
#decreasing genes
'Isoc1'
gnamesi <- c('')
# gnamesi <- c('Satb2')
gnamesi <- c('Isoc1')

gnamesi <- exprdata%>%.$gene_name%>%sample(1)

gnamesi <- 'Orc3'

getstanfits <- function(gnamesi,...){

realdata <- list(
   G=length(gnamesi),T=length(tps),K=3,
   # gene_name=gnamei,
   MS=exprdata%>%filter(gene_name%in%gnamesi)%>%filter(assay=='MS')%>%arrange(match(gene_name,gnamesi),time,rep)%>%.$signal%>%array(dim=c(3,5,length(gnamesi)))%>%aperm(c(1,2,3))%>%{2^.},
   # ribo =exprdata%>%filter(gene_name==gnamei)%>%filter(assay=='ribo',rep==1)%>%.$predicted_signal%>%{2^.}%>%rollmean(k=2)%>%matrix
   ribo =exprdata%>%filter(gene_name%in%gnamesi)%>%filter(assay=='ribo',rep==1)%>%arrange(match(gene_name,gnamesi),time)%>%.$predicted_signal%>%{2^.}%>%matrix(ncol=length(gnamesi))
)
realdata
realdataold
# realdataold <- realdata

#Fit non-hierarchical model 
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

list(realstanfit,realstanfit_lin)
}
allstanfits <- lapply(exprdata$gene_name%>%unique,safely(getstanfits))
allstanfits %<>% setNames(exprdata$gene_name%>%unique)
##Now with fixed values
# realstanfit_test <- rstan::stan(file='src/Stan/degmodel_simple.stan',data=realdata,
#                                control=list(adapt_delta=0.95,max_treedepth=20),
#                                init = list(list(lrTE=array(0,dim=c(1)))),
#                                chains=1,iter=1,warmup = 0,
#                                # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
#                                verbose=TRUE)

stanpars <- colnames(as.data.frame(realstanfit))

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
parsed_ml<-realstanfit%>%get_parsed_summary%>%mutate(val=`50%`)
parsed_ml_lin<-realstanfit_lin%>%get_parsed_summary%>%mutate(val=`50%`)
#function to pull out the mcmc samples for plotting 
fit <- realstanfit
get_prot_samples<-function(fit) fit %>%as.data.frame%>%select(matches('prot'))%>%mutate(sample=1:nrow(.))%>%gather(par,value,-sample)%>%mutate(ppars=parse_stan_pars(par))%>%unnest%>%
  filter(parameter=='prot')%>%select(time,value,sample,gene)
#e.g.
#realstanfit %>% get_prot_samples

gene_ind <- 1
gnamei <- gnamesi[gene_ind]

rdata2plot<-realdata$MS[,,gene_ind]%>%as.data.frame%>%set_colnames(tps)%>%mutate(rep=1:3)%>%gather(time,signal,-rep)

ml2plot<-parsed_ml%>%filter(parameter=='prot')%>%transmute(gene,time=as_factor(tps[time]),signal=(val))%>%filter(gene==gene_ind)
ml2plot_lin<-parsed_ml_lin%>%filter(parameter=='prot')%>%transmute(gene,time=as_factor(tps[time]),signal=(val))%>%filter(gene==gene_ind)

protsampleslin<- get_prot_samples(realstanfit_lin)%>%mutate(sample=factor(sample),signal=value,time=as_factor(tps[time]),model=as_factor('Linear'))%>%filter(gene==gene_ind)
protsamples<- get_prot_samples(realstanfit)%>%mutate(sample=factor(sample),signal=value,time=as_factor(tps[time]),model=as_factor('Kinetic'))%>%filter(gene==gene_ind)

# 
# protsamples_lowrTEmode<- get_prot_samples(realstanfit%>%as.data.frame%>%filter(`rTE[1]`<0.01))%>%mutate(sample=factor(sample),signal=value,time=as_factor(tps[time]),model=as_factor('Kinetic'))%>%filter(gene==gene_ind)
# protsamples_highrTEmode<- get_prot_samples(realstanfit%>%as.data.frame%>%filter(`rTE[1]`>0.1))%>%mutate(sample=factor(sample),signal=value,time=as_factor(tps[time]),model=as_factor('Kinetic'))%>%filter(gene==gene_ind)
# 

#plot showing trajectories and fits with MCMC samples 
trajectoryplot<-ggplot(rdata2plot%>%mutate(model='MS Data'),aes(color=model,y=log2(signal),x=as.numeric(as_factor(time))))+geom_point()+
  geom_line(size=I(1),data=ml2plot%>%mutate(model='Kinetic'),linetype=1)+
  geom_line(size=I(1),data=ml2plot_lin%>%mutate(model='Linear'),linetype=1)+
  geom_line(alpha=I(0.005),data=protsampleslin,aes(group=sample))+
  geom_line(alpha=I(0.005),data=protsamples,aes(group=sample))+
  # geom_line(size=I(2),data=get_prot_samples(ml2plot_lin)%>%mutate(model='Linear'),linetype=1)+
  scale_x_continuous(name='Stage',labels=tps)+
  scale_y_continuous(name='Log2 LFQ / Log2 Normalized Counts')+
  scale_color_manual(values = c('Kinetic'='red','Linear'='blue','MS Data'='black','rTE: 0'='purple','Riboseq Data'='dark green'))+
  theme_bw()+
  geom_line(data=exprdata%>%filter(gene_name==gnamesi[gene_ind],assay=='ribo',!is.na(predicted_signal))%>%mutate(signal=2^predicted_signal,model='Riboseq Data'))+
  ggtitle(label = str_interp('Linear vs. Kinetic Model - ${gnamei}'),sub="Faded Lines represent samples from Posterior\nSolid Line is Median Value")
  
# trajectoryplot


gglayout = c(1,1,1,2,3,4)%>%matrix(ncol=2)

dev.off()

stop()

arrangedplot<-gridExtra::grid.arrange(
  trajectoryplot,
    realstanfit%>%as.data.frame%>%.$`lrTE[1]`%>%qplot(bins=50,main=str_interp('Posterior Distribution log rTE\n${gnamei}'))+theme_bw(),
  realstanfit%>%as.data.frame%>%.$`ldeg[1]`%>%qplot(bins=50,main=str_interp('Posterior Distribution - log(degredation / 1.5 days) ${gnamei}'))+theme_bw(),
  realstanfit%>%as.data.frame%>%.$`deg[1]`%>%qplot(bins=50,main=str_interp('Posterior Distribution - degredation / 1.5 days ${gnamei}'))+theme_bw(),
layout_matrix=gglayout)


#' sab2 increasing
#' Orc3 weird zig zag in the riboseq that the model fits by just completely disregarding it.... This may be an instance of the riboseq signal itself
#' being fucked up
#' Ewsr1 - another fucking weird, zig-zaggy riboseq signal, that the model just chucks away.
#' 'Ap1' nice clear signal for rTE, increasing

ggsave(str_interp('plots/modelling/stanmodelcomp_${gnamei}.pdf')%T>%{normalizePath(.)%>%message},arrangedplot,w=12,h=10)


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





stanfit<-realstanfit
#
g2plot=1
allpars <- stanfit%>%as.data.frame%>%colnames
parse
  
parstoplot%<>%append(allpars%>%{.[!str_detect(.,'\\[')]})
parstoplot%<>%setdiff(c('deg','ms0logratio[1]'))
print(pairs(stanfit,pars=parstoplot))
plot(stanfit,pars=parstoplot%>%str_subset('prot'))


pairs(realstanfit,pars=allpars%>%str_subset('prot|_|tau'))

pdf('tmp.pdf');pairs(stanfit,pars=parstoplot);dev.off()
pdf('tmp.pdf');pairs(stanfit,pars=parstoplot%>%str_subset('prot|_'));dev.off()

samples[['rTE[1]']]%>%hist(40)


#now display relationshiop of prot latent vars to the parameters
parstoplot <- stanfit%>%as.data.frame%>%colnames%>%data_frame(par=.)%>%mutate(ppar=parse_stan_pars(par))%>%unnest%>%
  filter(!parameter%in%c('dP','sP','ltau','ldeg','lrTE','ms0logratio','degfact'))%>%
  filter(parameter!='prot')%>%
  .$par

pairs(stanfit,pars=parstoplot)
dev.off()
plot(1)
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
#' Okay so applying a prior like that to transformed parameters is bad (warning about needing to add in the jacobian of the transform)
#' so I now have rTE and tau getting initilialized on a log scale.
#' Combining the two genes Isoc1 seems to make the second one kind of divergent... (low Rhat, neff)