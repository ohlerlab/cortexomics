################################################################################
################################################################################
base::source(here::here('src/R/Rprofile.R'))
library(rstan)

base::source('src/Archive/R/Functions/rstan_functions.R')

if(!exists("cdsgrl")) {
	# base::source("src/Figures/Figure0/0_load_annotation.R")
	# load("src/Figures/Figure0/0_load_annotation.R")
}
#Things we want ot sim - production slow half life, production high half life, 
bmodel <- rstan::stan_model('src/Archive/Stan/becker_proda.stan')
bmodel_stationary <- rstan::stan_model('src/Archive/Stan/becker_proda_stationary.stan')
bmodel_degonly <- rstan::stan_model('src/Archive/Stan/bmodel_degonly.stan')
bmodel_linear <- rstan::stan_model('src/Archive/Stan/becker_proda_linear.stan')
# bmodel_stationary = fix_param(bmodel,vars2fix = c('l_st','l_pihalf'))%>%{f='src/Stan/bmodel_stationary.stan';cat(.,file=f);f}%>%stan_model
#bmodel_degonly = fix_param(bmodel,vars2fix = c('l_st'))%>%{f='src/Stan/bmodel_degonly.stan';cat(.,file=f);f}%>%stan_model

sim_prot_withopt <- function(
	model,
	l_pihalf=log(.5),
	l_st=2,
	ribo = 	ribo <- c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 100, 100, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30)
	){
	require(abind)
	# ribo = ribo[10:18]
	P = 2*ribo#opt will fill this in
	nT=length(ribo)
	sampdata=list()
	sampdata$G=1
	sampdata$T=nT
	sampdata$l_st_priorsd=3
	sampdata$l_ribo_priorsd=3
	sampdata$l_pihalf_priormu=0.5
	sampdata$l_pihalf_priorsd=3
	sampdata$lMSmu = t(log(P)) #- median(log(P))
	sampdata$lSeqmu = t(log(ribo)) #- median(log(ribo))
	sampdata$lMSsigma = array(diag(rep(0.1,nT)),c(1,nT,nT))
	sampdata$lMSsigma = t(rep(0.1,nT))
	sampdata$lSeqsigma = array(diag(rep(0.1,nT)),c(1,nT,nT))
	sampdata$lSeqsigma = t(rep(0.1,nT))
	
	initvals=list()
	initvals$l_pihalf <- array(l_pihalf,1)
	initvals$lribo = array(log(ribo),c(1,length(ribo)))
	initvals$l_st = array(l_st,1)
	initvals$lprot0 = array(initvals$l_st+initvals$lribo[1],1)

	pars = get_stanpars(model)
	initvals = initvals[pars]

	# opt <- rstan::optimizing(model,
	# data=sampdata,
	# verbose=TRUE,
	# init=initvals,
	# as_vector=F,
	# iter=1,
	# # draws=2,
	# save_iterations=TRUE)

	s = sampling(model,data=sampdata,init=function()initvals,chains=1,iter=1)
	sampdata$lMSmu = s%>%{rstan::extract(.)}%>%.$prot%>%log%>%adrop(1)
	sampdata
}

models = list(
	production = bmodel,
	stationary = bmodel_stationary,
	degredation = bmodel_degonly,
	linear = bmodel_linear
)

sims = lapply(models,sim_prot_withopt)
sims$random = sims$production%>%{.$lMSmu[]=.$lMSmu%>%sample;.}
library(txtplot)
# sims$production$lMSmu%>%txtplot
# sims$linear$lMSmu%>%txtplot
# sims$degredation$lMSmu%>%txtplot
# sims$stationary$lMSmu%>%unlist%>%{txtplot(.,ylim=c(min(.)*.9,max(.)*1.1))}
# sims$random$lMSmu%>%unlist%>%{txtplot(.,ylim=c(min(.)*.9,max(.)*1.1))}
simdata = sims[[1]]
simopts = lapply(sims,function(simdata){
	lapply(models,function(model){
		initvals=list()
		initvals$l_pihalf <- array(log(0.5),1)
		initvals$lribo = array(simdata$lSeqmu,c(1,length(simdata$lSeqmu)))
		initvals$l_st = array(simdata$lMSmu[[1]]-simdata$lSeqmu,1)
		initvals$lprot0 = array(simdata$lMSmu[[1]],1)
		reopt <- rstan::optimizing(model,data=simdata,init=initvals,as_vector=F,save_iterations=TRUE)
		reopt
	})
})

opt_simdata <- function(model,simdata){
		initvals=list()
		initvals$l_pihalf <- array(log(0.5),1)
		initvals$lribo = array(simdata$lSeqmu,c(1,length(simdata$lSeqmu)))
		initvals$l_st = array(simdata$lMSmu[[1]]-simdata$lSeqmu,1)
		initvals$lprot0 = array(simdata$lMSmu[[1]],1)
		reopt <- rstan::optimizing(model,data=simdata,init=initvals,as_vector=F,save_iterations=TRUE)
		reopt
}

simopts%>%map_depth(2,~.$value)
simopts%>%map_depth(2,~.$value)%>%map(.%>%unlist%>%which.max)

#now let's do a chi squared test
datnames = names(sims)%>%setNames(.,.)
modnames = names(models)%>%setNames(.,.)
model_tests = lapply(datnames,function(datname){
	lapply(modnames,function(modname){
		perrors =(log(simopts[[datname]][[modname]]$par$prot) - sims[[datname]]$lMSmu)
		w_perrors = perrors/sims[[datname]]$lMSsigma
		rerrors = (simopts[[datname]][[modname]]$par$lribo - sims[[datname]]$lSeqmu)
		w_rperrors = rerrors/sims[[datname]]$lSeqsigma
		n_df = simopts[[datname]][[modname]]$par[get_stanpars(models[[modname]])]%>%unlist%>%length
		errorsum = sum(c(w_perrors,rerrors)^2)
		BIC = -2*simopts[[datname]][[modname]]$value+(log(length(c(w_perrors,rerrors))))*n_df
		pval = pchisq(errorsum,n_df)
		c(BIC=BIC,pval=pval)
		# errorsum
		# errorsum
	})
})

model_tests%>%map_df(.id='data',.%>%bind_rows(.id='model'))%>%mutate(passtest = pval < 0.05)%>%
	group_by(data)%>%mutate(best=BIC==min(BIC))







####Lets plot some simulated data

modelplot <- function(simopts,gdatas,fname){
	stopifnot(!is.null(names(simopts)))
	stopifnot(!is.null(names(gdatas)))
modelvals <- map_df(.id='gene',names(simopts)%>%setNames(.,.),function(exgn){
	map_df(.id='model',names(simopts[[exgn]])%>%setNames(.,.),function(model){
	optres = simopts[[exgn]][[model]]
	# parsds = 1.96*sqrt(diag(solve(-optres$hessian)))
	ntps = optres$par$prot%>%length
	modelests = list(
		prot=log(optres$par[['prot']])%>%as.vector,
		ribo=optres$par$lribo%>%as.vector)%>%
		enframe('assay','value')%>%
		unnest(value)%>%
		mutate(time=c(1:ntps,1:ntps))
	modelests
	})
})
modeldata <- map_df(.id='gene',names(gdatas)%>%setNames(.,.),function(exgn){
	gdata <- gdatas[[exgn]]
	datavals = list(prot = gdata$lMSmu%>%as.vector,ribo = gdata$lSeqmu%>%as.vector)%>%
	enframe('assay','value')%>%unnest(value)
	datasigs = list(prot = gdata$lMSsig%>%as.vector,ribo = gdata$lSeqsig%>%as.vector)%>%
	enframe('assay','value')%>%unnest(value)
	datavals$lower = datavals$value - datasigs$value*1.96
	datavals$upper = datavals$value + datasigs$value*1.96
	ntps = gdata$lMSmu%>%length
	datavals%>%mutate(time=c(1:ntps,1:ntps))
})
modeldata$model='data'
gplotdf = bind_rows(modeldata,modelvals)
# g2plotttile = paste0(exampgenes,'\n',names(exampgenes))%>%setNames(exampgenes)
# modelcols2 <- c(modelcols,data='black')
modelcols2<-c('data'='black','degredation'='green','production'='blue','stationary'='purple','linear'='lightblue',msdev='red')
#
gplotdf$model%>%is_in(modelcols2%>%names)
# browser()
#now plot
plotfile<- here(paste0('plots/',fname))
pdf(plotfile,w=7,h=4)
plotlist = gplotdf%>%
	split(.,.$gene)%>%
	lapply(function(x)x%>%
		filter(assay=='prot' | model=='data')%>%
		# mutate(gene=g2plotttile[gene])%>%
		# filter(gene%>%str_detect(model)|model=='data')%>%
		# group_slice(1)%>%
		ggplot(aes(y=value,x=time,color=model))+
		geom_line()+
		ylab('')+
		scale_color_manual(values=modelcols2)+
		facet_grid(scale='free',assay~.)+
		# geom_errorbar(width=(0.2),aes(ymin=lower,ymax=upper))+
		theme_bw()
	)
# plotlist[[5]] %<>% add(scale_color_manual(values=modelcols2))
print(ggarrange(ncol=length(plotlist),plotlist=plotlist))
dev.off()
message(normalizePath(plotfile))
}
simribo =  c( 10, 10, 10, 10, 10, 100, 100, 30, 30, 30,30,30)
simdata <- sim_prot_withopt(bmodel, ribo=simribo, l_pihalf = log(2))
simopts <- map(models[c('production','linear')],opt_simdata,simdata)
pdf<-grDevices::pdf
modelplot(list(test=simopts),list(test=simdata),'kinetics_redo/demoplot.pdf')

stop()
# make_simdata
# #data looks like this
# list(
# 	G;// number of proteins
#   int T;// info on the number of conditions
#   matrix[G,T] lMSmu;
#   matrix[G,T] lSeqmu;
#   matrix[T,T] lMSsigma[G];
#   matrix[T,T] lSeqsigma[G];
#   real l_st_priorsd;
#   real l_ribo_priorsd;
#   real l_pihalf_priormu;
#   real l_pihalf_priorsd;

#goal - a chi-squared test to distinguish the different types
#generate the different types

################################################################################
########This mess
################################################################################
if(F){

	library(tidyverse)
library(here)
library(magrittr)
library(splines2)
library(splines)


################################################################################
########
################################################################################
library(tidyverse)	
library(here)	
sampdata <- here('data/sampdata.rds')%>%read_rds
ribo <- c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 100, 100, 30, 
          30, 30, 30, 30, 30, 30, 30, 30, 30)
nT = length(ribo)
Pe = ribo*3
P = Pe
stopifnot(length(Pe) == length(ribo))


sampdata$T=nT
sampdata$lMSmu = t(log(P)) #- median(log(P))
sampdata$lSeqmu = t(log(ribo)) #- median(log(ribo))
sampdata$lMSsigma = array(diag(rep(0.1,nT)),c(1,nT,nT))
sampdata$lSeqsigma = array(diag(rep(0.1,nT)),c(1,nT,nT))

time = 1:length(ribo)
timeknots <- time[c(-1,-length(time))]
sampdata$mybs =  cbind(1,ibs(time, knots = timeknots,degree = 1, intercept = TRUE))
sampdata$mydbs =bs(time, knots = timeknots,degree = 1, intercept = TRUE)

initvals=list()
initvals$prot0 <- array(sampdata$lMSmu[1],1)

d_0 = 0.6
Kd = d_0
lKd = log(Kd)
lpihalf = log(log(2)) - lKd
pihalf = log(2) / Kd
lpihalf = log(log(2)) - lKd
log(log(2)) - lpihalf
lKd
#log(log(2)) - lKd




sampdata[names(sampdata)%>%str_subset('sigma$')] %<>% map(multiply_by,0.001)
# 
# sampdata%>%map(dim)
# sampdata
# proDAsigmastan2

d_0 = 0.6;
Ks = 40;
ribo= c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 100, 100, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30);
P = rep(NA,length(ribo))
P[1] = 3* (ribo[1] * Ks)

for(i in 2:length(ribo)){
  P[i] = stepsynthdeg(d_0=d_0,Ks=Ks,T=1, s_0 = ribo[i-1],s_1 =  ribo[i],P_0 = P[i-1])
}

lribo = log(ribo)
l_st = log(Ks)-log(d_0)
sampdata$prot0 = array(initvals$l_st+initvals$lribo[1],1)
initvals$l_pihalf <- array(log(1),1)
initvals$lribo = array(log(ribo),c(1,length(ribo)))
initvals$l_st = array(log(2),1)
initvals$lprot0 = 1+array(initvals$l_st+initvals$lribo[1],1)
# initvals$lprot0 = NA
initvals$prot0=NULL
Kd = log(log(2)) -  initvals$l_pihalf;
Ks = exp(initvals$l_st - lKd);
ribo = exp(lribo);
exp(initvals$lprot0);
beckermodel <- rstan::stan_model('src/Stan/becker_proda.stan')
opt <- rstan::optimizing(beckermodel,data=sampdata,verbose=TRUE,init=initvals,as_vector=F,iter=0,save_iterations=TRUE)
library(txtplot)
txtplot(log(opt$par$prot))
txtplot(opt$par$lribo)
opt$par$l_st
opt$par$lprot0
opt$par$l_st
log(opt$par$prot)
sampdata$lMSmu
initvals$lprot0
initvals$l_st

sampdata$lMSmu
sampdata2 <- sampdata
sampdata2$lMSmu <- log(opt$par$prot)
opt2 <- rstan::optimizing(beckermodel,data=sampdata2,verbose=TRUE,as_vector=F,save_iterations=TRUE)
txtplot(log(opt$par$prot))
opt2$par$l_pihalf
opt$par$l_pihalf
initvals$l_pihalf

opt2$par$l_st
opt$par$l_st

opt2$par$prot/opt$par$prot
txtplot(log(opt2$par$prot),log(opt$par$prot))


}