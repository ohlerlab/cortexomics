################################################################################
################################################################################
base::source(here::here('src/R/Rprofile.R'))
library(rstan)

base::source('src/R/Functions/rstan_functions.R')

if(!exists("cdsgrl")) {
	# base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
	# load("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
}
#Things we want ot sim - production slow half life, production high half life, 
bmodel <- rstan::stan_model('src/Stan/becker_proda.stan')
bmodel_stationary <- rstan::stan_model('src/Stan/becker_proda_stationary.stan')
bmodel_degonly <- rstan::stan_model('src/Stan/bmodel_degonly.stan')
# bmodel_stationary = fix_param(bmodel,vars2fix = c('l_st','l_pihalf'))%>%{f='src/Stan/bmodel_stationary.stan';cat(.,file=f);f}%>%stan_model
#bmodel_degonly = fix_param(bmodel,vars2fix = c('l_st'))%>%{f='src/Stan/bmodel_degonly.stan';cat(.,file=f);f}%>%stan_model

sim_prot_withopt <- function(model,l_pihalf=log(.5),l_st=2){
	require(abind)
	ribo <- c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 100, 100, 30, 
          30, 30, 30, 30, 30, 30, 30, 30, 30)
	ribo = ribo[10:18]
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
	degredation = bmodel_degonly
)

sims = lapply(models,sim_prot_withopt)
sims$random = sims$production%>%{.$lMSmu[]=.$lMSmu%>%sample;.}
library(txtplot)
sims$production$lMSmu%>%txtplot
sims$degredation$lMSmu%>%txtplot
sims$stationary$lMSmu%>%unlist%>%{txtplot(.,ylim=c(min(.)*.9,max(.)*1.1))}
sims$random$lMSmu%>%unlist%>%{txtplot(.,ylim=c(min(.)*.9,max(.)*1.1))}
simdata = sims[[1]]
opts = lapply(sims,function(simdata){
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

opts%>%map_depth(2,~.$value)
opts%>%map_depth(2,~.$value)%>%map(.%>%unlist%>%which.max)

#now let's do a chi squared test
datnames = names(sims)%>%setNames(.,.)
modnames = names(models)%>%setNames(.,.)
model_tests = lapply(datnames,function(datname){
	lapply(modnames,function(modname){
		perrors =(log(opts[[datname]][[modname]]$par$prot) - sims[[datname]]$lMSmu)
		w_perrors = perrors/sims[[datname]]$lMSsigma
		rerrors = (opts[[datname]][[modname]]$par$lribo - sims[[datname]]$lSeqmu)
		w_rperrors = rerrors/sims[[datname]]$lSeqsigma
		n_df = opts[[datname]][[modname]]$par[get_stanpars(models[[modname]])]%>%unlist%>%length
		errorsum = sum(c(w_perrors,rerrors)^2)
		BIC = -2*opts[[datname]][[modname]]$value+(log(n_df))*n_df
		pval = pchisq(errorsum,n_df)
		c(BIC=BIC,pval=pval)
		# errorsum
		# errorsum
	})
})

model_tests


