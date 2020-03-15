
#+ setup, include=FALSE, echo=FALSE, eval=T
library(rmarkdown)
library(knitr)
library(here)
knitr::opts_chunk(root.dir = here::here(),eval=FALSE,cache=FALSE,echo=FALSE,warning = FALSE,message = FALSE,include=FALSE)
isknitr<-isTRUE(getOption('knitr.in.progress'))
if(!isknitr) rmarkdown::render(knitr::spin(here('src/R/Modelling/run_degmodel_dropout.R'),knit=F),output_dir=here('Reports'),knit_root_dir=here())

try(silent = T,{library(colorout)})
message('loading libraries')
library(checkmate)
library(here)
library(assertthat)
library(stringr)
library(tidyverse)
library(abind)
library(ggpubr)
library(data.table)
library(zeallot)
library(splines)
library(parallel)
library(broom)
# #
library(rstan)
library(txtplot)

tryCatch(library(proDD),error=function(e){
  options(repos = getOption("repos")["CRAN"])
  BiocManager::install(c('MSnbase'))
  #This package uses a similiar stan model to mine
  #I'm using their EM function to get some global hyperparameters (e.g. global mean, variance, dropout frequency of MS in each library)
  devtools::install_github("const-ae/proDD")
})


#set up stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

source('src/R/Functions/rstan_functions.R')
library(bayesplot)

#set up stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#simplest model for single genes, but with dropout
dp_stanfile = here('src/Stan/mod_proDD.stan')%T>%{stopifnot(file.exists(.))}
dp_model = rstan::stan_model(dp_stanfile)


# dp_stanfile_protonly = here('src/Stan/dropout_protonly.stan')%T>%
# 	{stopifnot(file.exists(.))}
# dp_model_protonly = rstan::stan_model(dp_stanfile_protonly)

#Model with additional normalization parameters shared between genes, for each ribo and prot timepoint
dp_stanfile_ribnorm = here('src/Stan/mod_proDD_ribnorm.stan')%T>%{stopifnot(file.exists(.))}
dp_model_ribnorm = rstan::stan_model(dp_stanfile_ribnorm)#
#experimental model with deviances from the kinetic trajectory
dp_stanfile_dev = here('src/Stan/mod_proDD_dev.stan')%T>%{stopifnot(file.exists(.))}
dp_model_dev = rstan::stan_model(dp_stanfile_dev)#



#+ notinreportyet, include=F,eval=F

#Get a hold of the data
#load the metadata on all our genes and the many IDs involved
metainfo<-read_tsv(here('data/metainfo.tsv'))
#Pull out the gene names we want to analyze
uids4stan <- metainfo%>%
	filter(isbest)%>%#these are thee final pairs of gene/cds/mass spec ids that we use for modeling MS
	filter(sig_MS_change)%>%
	filter(n_stagemissing<=2)%>%#now filter them further - those with too many missing values will crash Rstan
	.$uprotein_id
##Should match the input objects used above
best_uprotein_ids <- uids4stan
best_ms_ids <- metainfo%>%{.$ms_id[match(best_uprotein_ids,.$uprotein_id)]}
best_protein_ids <- metainfo%>%{.$protein_id[match(best_uprotein_ids,.$uprotein_id)]}



# ribomat <- mscountvoom$E%>%.[str_subset(colnames(.),'ribo')]
# ribomat %>% saveRDS('data/ribomat.rds')
# matchedms_mat %>% saveRDS('data/matchedms_mat.rds')
# library(zeallot)
# c(matched_ms_mat,mscountvoom$E) %<-% with(new.env(),
	# load('data/integrate_exprdata2.Rdata')
	# list(matched_ms_mat,mscountvoom$E)	
# })

matchedms_mat <- memoise(readRDS)(here('data/matched_ms_mat.rds'))
ribomat <- memoise(readRDS)(here('data/ribomat.rds'))

# matchedms_mat %>%saveRDS('data/matched_ms_mat.rds')
# ribomat %>%saveRDS('data/ribomat.rds')

# matchedms_mat%>%colMedians(na.rm=T)%>%txtplot





################################################################################
########get our data, scale it, get prior parameters use proDD
################################################################################
	
if(!file.exists('data/proddparams.Rdata')){
	library(proDD)
	#get matrices of riboseq and ms data, median norm them, then subtract again so the values center on zero
	#(the former is necessary, the second is simply for numeric reasons - stan will freak out if values are high)
	matchedms_mat_rscl <- matchedms_mat[best_ms_ids,]
	matchedms_mat_rscl <- matchedms_mat_rscl%>%{proDD::median_normalization(.)}
	msmed <- matchedms_mat_rscl%>%apply(2,median,na.rm=T)%>%median(.,na.rm=T)
	matchedms_mat_rscl %<>% subtract(msmed)
	matchedms_mat_rscl%<>% set_rownames(best_uprotein_ids)
	length(best_protein_ids)
	length(best_uprotein_ids)

	best_protein_ids<-best_uprotein_ids%>%str_extract('[^_]+')
	ribomed <- ribomat[]%>%median(.,na.rm=T)
	ribomatrscl <- ribomat[best_uprotein_ids,]%>%
		set_rownames(best_uprotein_ids)
	ribomatrscl <- ribomatrscl%>%{proDD::median_normalization(.)}
	ribomed <- ribomatrscl%>%apply(2,median,na.rm=T)%>%median(.,na.rm=T)
	ribomatrscl %<>% subtract(ribomed)

	allvoom_sigma <- countvoom$weights%>%{1/.}%>%
		set_rownames(rownames(countvoom$E))%>%
		set_colnames(colnames(countvoom$E))%>%
		.[best_protein_ids,,drop=F]%>%
		set_rownames(best_uprotein_ids)
	tpvect<-matchedms_mat_rscl%>%colnames%>%str_extract('[^_]+')
	#
	proddparam <- proDD::fit_parameters(matchedms_mat_rscl,tpvect)
	save(matchedms_mat_rscl,ribomatrscl,ribomed,msmed,proddparam,
		tpvect,allvoom_sigma,file='data/proddparams.Rdata')
}else{
	library(zeallot)
	# c(matchedms_mat_rscl,ribomatrscl,ribomed,msmed,proddparam,tpvect) %<-% get()
	# file.remove(file='data/proddparams.Rdata')
	load(file='data/proddparams.Rdata')
}

stopifnot(all(uids4stan%in%rownames(matchedms_mat_rscl)))#the mass spec values, with missing values
stopifnot(all(uids4stan%in%rownames(ribomatrscl)))#the riboseq values, with 
stopifnot(all(uids4stan%in%rownames(allvoom_sigma)))

#function that pulls out the stan data we need in the right shape, and pairs it with proDD parameters etc. 
get_dp_standata <- function(sel_uprotein_id,
	matchedms_mat_rscl,ribomat,ribomatrscl,
	params=proddparam
	){
	stopifnot(all(sel_uprotein_id%in%rownames(matchedms_mat_rscl)))
	# msids <- ms_id2protein_id%>%filter(uprotein_id %in% sel_uprotein_id)%>%.$ms_id
	msdata <- matchedms_mat_rscl[sel_uprotein_id,,drop=F]
	# msdata <- sizefactnorm(msdata)
	# MSmed = 0
	# MSmed = median(msdata,na.rm=T)
	# msdata = msdata - MSmed
	msdata %<>% replace_na(-Inf)
	n_missing_ms <- sum(!is.finite(msdata))
	#
	tpvect <- colnames(msdata)%>%str_extract('[^_]+')%>%as.factor%>%as.numeric
	#clear unecessary columns from this object
	ribocols <- colnames(ribomatrscl)%>%str_subset('ribo')
	ribomatrscl <- ribomatrscl[,ribocols]#the ribomat had extra columns for other data that needed to be gotten rid of
	ribotpvect <- colnames(ribomatrscl)%>%str_extract('[^_]+')%>%as.factor%>%as.numeric
	ribomat <- ribomatrscl[sel_uprotein_id,ribocols,drop=F]
	ribomed = 0
	stopifnot(ncol(ribomat)==length(ribotpvect))
	# ribomed = median(ribomat,na.rm=T)
	# ribomat <- ribomat - ribomed
	#
	prodd_params <- params
	#
	require(splines2)
	time=1:n_distinct(tpvect)
	timeknots <- time[c(-1,-length(time))]
	mybs <- cbind(1,ibs(time, knots = timeknots,degree = 1, intercept = TRUE))
	mydbs = bs(time, knots = timeknots,degree = 1, intercept = TRUE)
	#
	standata = list(
	  nsamples=ncol(msdata),#number of samples
	  nribosamples=ncol(ribomat),#number of samples
	  G=nrow(msdata),#number of proteins
	  T=n_distinct(tpvect),#info on the number of conditions
	  totalmissing=n_missing_ms,#info on th total amount of missing data
	  lMS=msdata,#data
	  lribo=ribomat,#data
	  experimental_design=tpvect, #indicator variable matching the samples to conditions
	  experimental_design_r=ribotpvect,
	  zeta=prodd_params$hyper_params$zeta,#the spread of the dropout point for a library, gets combined with the varianc per protein
	  rho=prodd_params$hyper_params$rho,#rho the location of the dropout point for a given library
	  mu0=prodd_params$hyper_params$mu0,#fit by proDD the mean of means
	  sigma20=prodd_params$hyper_params$sigma20,#fit by proDD - the variance in means
	  eta=prodd_params$hyper_params$eta,#fit by proDD - th evariance in protein variances
	  nu=prodd_params$hyper_params$nu,#fit by proDD - the mean of protein variances
	  voom_sigma=allvoom_sigma[sel_uprotein_id,ribocols,drop=F],
	  mybs = mybs,
	  mydbs = mydbs
	)
	standata
	invisible(standata)
}
get_dp_standata <- partial(get_dp_standata,
	matchedms_mat_rscl=matchedms_mat_rscl,ribomat=ribomat,ribomatrscl=ribomatrscl,
	params=proddparam)
set.seed(0)
# uids4stan%<>%sample(100)
#Make sure these interesting genes are in it
uids4stan <- union(metainfo%>%filter(gene_name%in%c('Flna','Satb2'))%>%filter(isbest)%>%.$uprotein_id%>%unique,uids4stan)
#
#make sure we can pull all the data
allstandata <- get_dp_standata(uids4stan)
#First, try fitting Satb2 - well behaved, gets us our param vector
satb2fits<-replicate(simplify=F,20,{
	pre = Sys.time()
	# ponlyoptres<-rstan::optimizing(dp_model_protonly,data=get_dp_standata(best_satb2_uid),init=list(prot0=array(0,dim=1)))
	optres<-rstan::optimizing(dp_model,data=get_dp_standata('ENSMUSP00000110057_4528'),algorithm='Newton',as_vector=FALSE,hessian=TRUE)
	Sys.time() - pre
	optres
})


#get best fit results
optres <- satb2fits[order(satb2fits%>%map_dbl('value'))%>%rev]%>%.[[1]]
#get initial values for single genes - set all to zero except the variance parameter
initvals <- optres$par
for(p in names(initvals))initvals[[p]] <- initvals[[p]] - initvals[[p]]
initvals[['sigma2']] = initvals[['sigma2']]+0.2

# initvals <- optres$par%>%{.['sigma2[1]<-0.2;.}%>%stanpars_to_list
set.seed(0)

names(optres$par)%<>%str_replace_all('[,\\[\\]]','.')%>%str_replace('\\.$','')
optres$par[colnames(optres$hessian)] - 1.96*sqrt(diag(solve(-optres$hessian)))
optres$par[colnames(optres$hessian)] + 1.96*sqrt(diag(solve(-optres$hessian)))


optres$hessian


require(R.utils)
require(rstan)

################################################################################
########Run optimization over all individual genes
################################################################################
	
# 'data/geneoptfits.rds'%>%file.remove
# selgene='ENSMUSP00000110057_4528'

if(!file.exists('data/geneoptfits.rds')){
	#
	require('R.utils')
	genefits <- list()
	# 			
	genefits <- mclapply(mc.cores=20,uids4stan,(function(selgene){
		cat(paste0('.',which(uids4stan==selgene),'..'))
			#
			selgenefits<-lapply(rep(selgene,5),safely(function(selgene){
				withTimeout(timeout=15,{
					# message(selgene)
					optdat <- get_dp_standata(selgene)
					optdat$lMS
					if(sum(is.finite(optdat$lMS))<2) return(NULL)
					#apparently we need actual hessian - Netwon succeeds for individual genes where
					#L-BFGS does not. 
					optimizing(dp_model,data=optdat,init=initvals,algorithm='Newton',as_vector=F,hessian=TRUE) 
				})
			}))
		selgenefits
	}))
	saveRDS(genefits,file='data/geneoptfits.rds')
}else{
	genefits<- readRDS((file='data/geneoptfits.rds'))
}
genefits%>%length
#see how many crashed for each gene
genefits%>%map_dbl(.%>%map('result')%>%map_lgl(is.null)%>%sum)%>%table
#None.
genefitsres<- genefits%>%map(.%>%map('result'))
genefitsres%<>%setNames(uids4stan)

#see how many converge
convnums <- genefitsres%>%map_dbl(.%>%map_dbl('return_code')%>%`==`(0)%>%sum)
genefitsres <- genefitsres[which(convnums!=0)]
conv_uprotids <- genefitsres%>%names
stopifnot(all(convnums==5))

#get the best fit object for each gene
bestfits <- lapply(genefitsres,function(selgenefits){
	bestfit<-selgenefits%>%map_dbl('value')%>%which.max
	selgenefits[[bestfit]]
})
#and optimized parameters
bestfitinits <- bestfits%>%map('par')
#name these
bestfits%<>%setNames(conv_uprotids)
bestfitinits%<>%setNames(conv_uprotids)

diag(solve(-genefitsres[[1]][[1]]$hessian))

genefitsres[[1]]%>%map('par')%>%map('l_st')
genefitsres[[1]]%>%map('par')%>%map('l_pihalf')

#All are now converging, so probably not ncessary to do this subsetting
finaltestgenes <- conv_uprotids

#now select the best opt fit and optimize again witht that
#this should succeed without warnings
# options(warn=1)
# bestfitinits%>%names%>%setdiff(rownames(matchedms_mat_rscl))
# reopts <- imap(bestfitinits,safely(function(bestfitinit,selgene_uid){
	# optres<-rstan::optimizing(dp_model,data=get_dp_standata(selgene_uid),
	# init=bestfitinit,algorithm='Newton')
# }))
# reopts %<>% map('result')
# reopts%>%map('return_code')%>%`!=`(0)%>%which%>%names
# finaltestgenes <-reopts%>%map('return_code')%>%`==`(0)%>%which%>%names

##we now have the modes for each gene, let's combine these into a list of parameters for the joint model
get_comb_initvals <- function(bestfitinits){
	combinitvals <- lapply(names(bestfitinits[[1]])%>%setNames(.,.),function(argind){
		bestfitinits%>%
		map(argind)%>%
		setNames(.,seq_along(.))%>%
		do.call(what=partial(abind::abind,along=1))
	})
	combinitvals	
}
combinitvals <- get_comb_initvals(bestfitinits)

#Now what if we try to fit normalization factors for each timepoint, fitting over all
#genes at once?
ribonorm_combinitvals <- c(combinitvals,list('protnorm'=array(rep(0,5),dim=list(5)),'ribnorm'=array(rep(0,5),dim=list(5))))

stop('Now to try optimizing norm factors over all genes')
message(paste0('optimizing simple model on : ',length(finaltestgenes)),' genes')

# LBFGS doesn't seem to work at all
alloptres<-rstan::optimizing(dp_model,data=get_dp_standata(finaltestgenes),
		init=combinitvals,as_vector=F,verbose=T,algorithm='LBFGS')

#now try to optimize the Ribo-seq libnorm factors as well
set.seed(0)
#We can fit these, but they change minimally
ribonormopt <- optimizing(dp_model_ribnorm,
	data=get_dp_standata(finaltestgenes),
	init=ribonorm_combinitvals,
	as_vector=F,save_iterations=FALSE,algorithm='LBFGS',verbose=TRUE,iter=2e3,hessian=FALSE)

#So - the avlues for these don't seem to be changing very much - and are converging on pretty mininiscule values
#Which makes me think I can ignore them.
plot(ribonormopt$par$ribnorm)
plot(ribonormopt$par$protnorm)

stop()

mschangegenes <- finaltestgenes%>%interect(metainfo%>%filter(sig_MS_change))

ribonormopt <- optimizing(dp_model_ribnorm,
	data=get_dp_standata(finaltestgenes),
	init=ribonorm_combinitvals,
	as_vector=F,algorithm='LBFGS',verbose=TRUE,iter=1,hessian=FALSE)


ribonormopt$par$ribnorm%>%txtplot
ribonormopt$par$protnorm%>%txtplot


#Now let's also try our model that allows protein deviations - this will be 
#underdetermined, unless we have multiple genes in play
finaltestgenes%>%head
devtestgenes <- finaltestgenes%>%head(100)

#Okay so this goes nuts and fits massive deviations
devopt <- optimizing(dp_model_dev,
	data=get_dp_standata(devtestgenes),
	init=get_comb_initvals(bestfitinits[devtestgenes]),
	as_vector=F,algorithm='LBFGS',verbose=TRUE,iter=1e3,hessian=FALSE)

# BiocManager::install(c('tidybayes'))
protmeanvect <- allstandata$lMS[1,]%>%split(ceiling(seq_along(.)/3))%>%map_dbl(mean)
devopt$par$p_dev[1,]
devopt$par$p_dev[1,]

#so the deviancem model fits better
sum(abs(devopt$par$prot[1,] - protmeanvect))
sum(abs(optres$par$prot[1,] - protmeanvect))

#Okay but surely those deviances effect protein a LOT
teval <- optimizing(dp_model_dev,
	data=get_dp_standata(devtestgenes),
	init=devopt$par%>%{.$p_dev[,] <- 0; . },
	as_vector=F,algorithm='LBFGS',verbose=TRUE,iter=50e3,hessian=FALSE)

teval$par$prot[1,]
alloptres$par$prot[1,]
teval$par$p_dev[1,]

################################################################################
########Now, let's try actually sampling from the posterior 
################################################################################
#for one gene
i=1

#This doesn't work that well....
samplings <- imap(bestfitinits[finaltestgenes[1]],safely(function(bestfitinit,selgene_uid){
	optres<-rstan::sampling(dp_model,data=get_dp_standata(selgene_uid),chains=4,
	init=function(){bestfitinit},control=list(adapt_delta=.9),iter=1e3)
}))
samplings%<>%map('result')

#try sampling from a model with the deviations
pdevsamplings  <- sampling(dp_model_dev,
	chains=4,
	data=get_dp_standata(devtestgenes),
	init=function(){devopt$par},
	verbose=TRUE,control=list(adapt_delta=.9),iter=400)

#look at the n_eff for these genes
summary(ribonormopt)$summary%>%
	{x=.;rownames(x)%>%parse_stan_pars(c('gene','time'))%>%bind_rows%>%cbind(.,x)}%>%
	filter(parameter%>%str_detect('pi|Ks'))%>%
	group_by(parameter)%>%
	slice(1:3)%>%as.data.frame

#This emphatically does not help...
#now plot
plotfile<- here(paste0('plots/','tmp','.pdf'))
pdf(plotfile)

dev.off()
normalizePath(plotfile)



#now plot
plotfile<- here(paste0('plots/','tmp','.pdf'))
pdf(plotfile)
gname = filter(metainfo,uprotein_id==names(samplings)[[1]])%>%.$gene_name
pairplot <- mcmc_pairs(samplings[[i]],pars=vars(matches('lKs|sigma|st|pi|(prot\\[1,2\\])|(prot\\[1,3\\])|lp__')),save_gg_objects=TRUE)
# str(pairplot)
print(pairplot)
dev.off()
normalizePath(plotfile)




