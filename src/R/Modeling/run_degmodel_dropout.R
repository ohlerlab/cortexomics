# source('src/R/Rprofile.R')
try({library(colorout)})
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
library(proDD)
# # library(conflicted)
# # ({map(lsf.str("package:BiocGenerics")%>%as.character%>%setdiff('basename'),.f=conflict_prefer,'IRanges')})
# # suppressMessages({map(lsf.str("package:BiocGenerics"),.f=conflict_prefer,'BiocGenerics')})
# # load.image('data/intergrete')

#set up stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#
if(!exists('best_protein_ids')){
	load(here('data/integrate_exprdata2.Rdata'))
}

dp_stanfile = here('src/Stan/mod_proDD.stan')%T>%{stopifnot(file.exists(.))}
dp_model = rstan::stan_model(dp_stanfile)

dp_stanfile_protonly = here('src/Stan/dropout_protonly.stan')%T>%{stopifnot(file.exists(.))}
dp_model_protonly = rstan::stan_model(dp_stanfile_protonly)



################################################################################
########get our data, scale it, get prior parameters use proDD
################################################################################
	
if(!file.exists('data/proddparams.Rdata')){
	library(proDD)
	#get matrices of riboseq and ms data, median norm them, then subtract again so the values center on zero
	#(the former is necessary, the second is simply for numeric reasons - stan will freak out if values are high)
	best_ms_ids <- ms_id2protein_id$ms_id[match(best_uprotein_ids,ms_id2protein_id$uprotein_id)]
	matchedms_mat_rscl <- matchedms_mat[best_ms_ids,]
	matchedms_mat_rscl <- matchedms_mat_rscl%>%{proDD::median_normalization(.)}
	msmed <- matchedms_mat_rscl%>%apply(2,median,na.rm=T)%>%median(.,na.rm=T)
	matchedms_mat_rscl %<>% subtract(msmed)
	matchedms_mat_rscl%<>% set_rownames(best_uprotein_ids)
	length(best_protein_ids)
	length(best_uprotein_ids)

	best_protein_ids<-best_uprotein_ids%>%str_extract('[^_]+')
	ribomed <- mscountvoom$E[]%>%median(.,na.rm=T)
	ribomatrscl <- mscountvoom$E[best_uprotein_ids,]%>%set_rownames(best_uprotein_ids)
	ribomatrscl <- ribomatrscl%>%{proDD::median_normalization(.)}
	ribomed <- ribomatrscl%>%apply(2,median,na.rm=T)%>%median(.,na.rm=T)
	ribomatrscl %<>% subtract(ribomed)

	ribocols <- colnames(allcountmat)%>%str_subset('ribo')
	allvoom_sigma <- countvoom$weights%>%{1/.}%>%
		set_rownames(rownames(countvoom$E))%>%
		set_colnames(colnames(countvoom$E))%>%
		.[best_protein_ids,ribocols,drop=F]%>%
		set_rownames(best_uprotein_ids)
	tpvect<-matchedms_mat_rscl%>%colnames%>%str_extract('[^_]+')

	#
	proddparam <- proDD::fit_parameters(matchedms_mat_rscl,tpvect)
	save(matchedms_mat_rscl,ribomatrscl,ribomed,msmed,proddparam,tpvect,allvoom_sigma,file='data/proddparams.Rdata')
}else{
	library(zeallot)
	# c(matchedms_mat_rscl,ribomatrscl,ribomed,msmed,proddparam,tpvect) %<-% get()
	# file.remove(file='data/proddparams.Rdata')
	load(file='data/proddparams.Rdata')
}

#function that pulls out the stan data we need in the right shape, and pairs it with proDD parameters etc. 
get_dp_standata <- function(sel_uprotein_id,
	matchedms_mat_rscl,allcountmat,ribomatrscl,
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
	ribocols <- colnames(allcountmat)%>%str_subset('ribo')
	ribotpvect <- ribocols%>%str_extract('[^_]+')%>%as.factor%>%as.numeric
	#
	ribomat <- ribomatrscl[sel_uprotein_id,ribocols,drop=F]
	ribomed = 0
	# ribomed = median(ribomat,na.rm=T)
	# ribomat <- ribomat - ribomed
	#
	#
	prodd_params <- params
	#
	#
	#
	#
	require(splines2)
	#
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
	  voom_sigma=allvoom_sigma[sel_uprotein_id,,drop=F],
	  mybs = mybs,
	  mydbs = mydbs
	)
	standata
	invisible(standata)
}
get_dp_standata <- partial(get_dp_standata,
	matchedms_mat_rscl=matchedms_mat_rscl,allcountmat=allcountmat,ribomatrscl=ribomatrscl,
	params=proddparam)
#
#load the metadata on all our genes and the many IDs involved
metainfo<-read_tsv(here('pipeline/exprdata/limma_genemetadata.tsv'))
#Pull out the gene names we want to analyze
uids4stan <- metainfo%>%
	filter(isbest)%>%#these are thee final pairs of gene/cds/mass spec ids that we use for modeling MS
	filter(sig_MS_change)%>%
	filter(n_stagemissing<=2)%>%#now filter them further - those with too many missing values will crash Rstan
	.$uprotein_id
##Should match the input objects used above
stopifnot(all(uids4stan%in%rownames(matchedms_mat_rscl)))
stopifnot(all(uids4stan%in%rownames(ribomatrscl)))
stopifnot(all(uids4stan%in%rownames(allvoom_sigma)))

set.seed(0)
# uids4stan%<>%sample(100)
#Make sure these interesting genes are in it
uids4stan <- union(metainfo%>%filter(gene_name%in%c('Flna','Satb2'))%>%filter(isbest)%>%.$uprotein_id%>%unique,uids4stan)

#temporary - this can go when we recalculated the data scaling
# stopifnot(uids4stan%>%is_in(rownames(matchedms_mat_rscl)))

#make sure we can pull all the data
allstandata <- get_dp_standata(uids4stan)
#
#First, try fitting Satb2 - well behaved, gets us our param vector
satb2fits<-replicate(simplify=F,20,{
	pre = Sys.time()
	# ponlyoptres<-rstan::optimizing(dp_model_protonly,data=get_dp_standata(best_satb2_uid),init=list(prot0=array(0,dim=1)))
	optres<-rstan::optimizing(dp_model,data=get_dp_standata(best_satb2_uid),algorithm='Newton')
	Sys.time() - pre
	optres
})

#get best fit results
optres <- satb2fits[order(satb2fits%>%map_dbl('value'))%>%rev]%>%.[[1]]
#get initial values for single genes - 
initvals <- optres$par%>%{.[]<-0;.['sigma2[1]']<-0.2;.}%>%stanpars_to_list


require(R.utils)
require(rstan)

# 'data/geneoptfits.rds'%>%file.remove
if(!file.exists('data/geneoptfits.rds')){
	#
	require('R.utils')
	genefits <- list()

		# selgene=uids4stan[666]
		# selgene=uids4stan[4947]
		# selgene=best_satb2_uid
# 			
	genefits <- mclapply(mc.cores=20,uids4stan,(function(selgene){
		cat(paste0('.',which(uids4stan==selgene),'..'))
			#
			selgenefits<-lapply(rep(selgene,10),safely(function(selgene){
				withTimeout(timeout=15,{
					# message(selgene)
					optdat <- get_dp_standata(selgene)
					optdat$lMS
					if(sum(is.finite(optdat$lMS))<2) return(NULL)
					# optdat$lMS<-optdat$lMS[1,,drop=F]%>%{.[.== -Inf]= -1.5;.}
					# optdat$lMS<-optdat$lMS[1,,drop=F]%>%{.[.== -Inf]= -1.5;.}
					# optimizing(dp_model,data=optdat,iter=100,verbose=1) 
					#for a given gene, try 20 values
					# suppressWarnings({optres<-lapply(1:20,function(sd) rstan::optimizing(dp_model,data=optdat,iter=1,seed=sd,as_vector=F))})

					# initsworked <- optres%>%map_dbl('return_code')%>%`==`(0)
					# bestinitind <- optres[initsworked]%>%map_dbl('value')%>%which.max
					# initvals <- optres[initsworked][[bestinitind]]$par%>%{.}
					#apparently we need the adaptive step size - Netwon succeeds for individual genes where
					#L-BFGS does not. 
					optimizing(dp_model,data=optdat,init=initvals,algorithm='Newton',as_vector=F) 
				})
			}))
		selgenefits
	}))
	saveRDS(genefits,file='data/geneoptfits.rds')
}else{
	genefits<- readRDS((file='data/geneoptfits.rds'))
}
stop()
#see how many crashed for each gene
genefits%>%map_dbl(.%>%map('result')%>%map_lgl(is.null)%>%sum)%>%table
#None.
genefitsres<- genefits%>%map(.%>%map('result'))
genefitsres%<>%setNames(uids4stan)

#see how many converge
convnums <- genefitsres%>%map_dbl(.%>%map_dbl('return_code')%>%`==`(0)%>%sum)
genefitsres <- genefitsres[which(convnums!=0)]
conv_uprotids <- genefitsres%>%names
stopifnot(all(convnums==10))

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
finaltestgenes <- conv_uprotids%>%head(100)
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
combinitvals <- lapply(names(bestfitinits[[1]])%>%setNames(.,.),function(argind){
	bestfitinits[finaltestgenes]%>%map(argind)%>%setNames(.,seq_along(.))%>%do.call(what=partial(abind,along=1))
})

#And use this to optimize them all at once (they are independent yes, but) does this
#crash stan?
optres<-rstan::optimizing(dp_model,data=get_dp_standata(finaltestgenes),
		init=combinitvals,as_vector=F,verbose=T,algorithm='LBFGS')
#No, but it converges in one step... which I guess it should
perturbcombinitvals<-combinitvals
perturbcombinitvals$prot0[1]<-3

allopt <- rstan::optimizing(dp_model,data=get_dp_standata(finaltestgenes),
		init=perturbcombinitvals,as_vector=F,verbose=T,algorithm='Newton',hessian=TRUE,iter=10)



samplings <- imap(bestfitinits[finaltestgenes],safely(function(bestfitinit,selgene_uid){
	optres<-rstan::sampling(dp_model,data=get_dp_standata(selgene_uid),
	init=function(){bestfitinit},control=list(adapt_delta=.9),iter=200)
}))

#get combination initvals with that
stanpars_to_list<-function(stanparvect){
  #get the dimensions of each of the parameters
  iparnames <- stanparvect%>%names
  parnames <- iparnames%>%str_remove('(\\[|\\.)[0-9\\,]+(\\]|\\.)$')
  parinds <- iparnames%>%str_extract('(\\[|\\.)[0-9\\,]+(\\]|\\.)$')%>%str_extract_all(regex('\\d+'))%>%map(as.numeric)
  #
  parinds <- parinds%>%split(parnames)%>%map(.%>%simplify2array%>%t)
  stopifnot(all(unique(parnames)%in%c(c("sigma2", "l_st", "lcM", "l_pihalf", "prot0", "cv", "prot",
		"mRNA", "Kd", "lKd", "lKs", "zetastar", "cM", "lp__"))))
  if('zetastar' %in% names(parinds) )parinds[['zetastar']]%<>%t
  #
  stanlist <- parinds%>%map(.%>%apply(2,max)%>%replace_na(1)%>%array(NA,dim=.))
  #
  parname='zetastar'
  for(parname in unique(parnames)){
    stanlist[[parname]][parinds[[parname]]] <- stanparvect[parnames==parname]
  }
  stanlist
}
samplingmedsaslist<-samplings[finaltestgenes]%>%map('result')%>%map(.%>%summary%>%.[[1]]%>%{setNames(.[,'50%'],rownames(.))})%>%map(.%>%stanpars_to_list%>%.[names(.)!='lp__'])

samplingcombinitvals <- lapply(names(bestfitinits[[1]])%>%setNames(.,.),function(argind){
	samplingmedsaslist%>%map(argind)%>%setNames(.,seq_along(.))%>%do.call(what=partial(abind,along=1))
})

#Now let's try with a 
dp_stanfile_ribnorm = here('src/Stan/mod_proDD_ribnorm.stan')%T>%{stopifnot(file.exists(.))}
dp_model_ribnorm = rstan::stan_model(dp_stanfile_ribnorm)

#
ribonorm_combinitvals <- c(combinitvals,list('ribonorm'=array(rep(0,5),dim=list(5))))

#now try to optimize the Ribo-seq libnorm factors as well
ribonormopt <- optimizing(dp_model_ribnorm,
	data=get_dp_standata(finaltestgenes),
	init=ribonorm_combinitvals,
	as_vector=F,hessian=F,algorithm='Newton',verbose=TRUE,iter=10)

ribonormopt$par$ribnorm%>%txtplot
ribonormopt$par%>%.$sigma2%>%log10%>%txtplot


#and try to sample from it.
ribonormsamp <- sampling(
	dp_model_ribnorm,
	data=get_dp_standata(finaltestgenes),
	init=function(){ribonorm_combinitvals},
	iter=100,
	control=list(adapt_delta=0.95)
)

parse_stan_pars<-function(stanpars,indnames=c('time','gene')){
  if(any(str_detect(stanpars,'\\['))){
    stopifnot(max(str_count(stanpars,','))==1)
    parsedpars<-stanpars%>%str_match('([^\\[]+)\\[?(\\d*),?(\\d*)\\]?')%>%as.data.frame%>%
    .[,-1]
  }else{
    # browser()
    stopifnot(any(str_detect(stanpars,'\\.')))
    stopifnot(max(str_count(stanpars,'\\.'))==2)
    parsedpars<-stanpars%>%str_match('([^\\.]+)\\.?(\\d*)\\.?(\\d*)\\.?')%>%as.data.frame%>%
    .[,-1]
  }
  #
  n_inds <- length(colnames(parsedpars))-1
  #
    parsedpars%>%.[,-1]%>%.[,]%>%
            set_colnames(indnames)%>%
            as.data.frame%>%
          mutate(parameter=parsedpars[,1])%>%
      select(parameter,!!!indnames)%>%
      mutate_at(vars(!!!indnames),list(~as.numeric(as.character(.))))%>%
      split(.,seq_len(nrow(.)))
}


tidy_opt_prots<-function(pars,gnames){
	protdf<-pars$prot%>%set_rownames(gnames)%>%set_colnames(tps)%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%
		gather(time,estimate,-uprotein_id)%>%
		filter(!str_detect(uprotein_id,'\\.\\d$'))%>%
		safe_left_join(metainfo%>%distinct(uprotein_id,gene_name))
	protdf%<>%mutate(
		CI.L = estimate - 0.0001,
		CI.R = estimate - 0.0001,
		ntime = match(time,tps),
		assay='MS'
	)
	protdf
}

parse_sampling_pars<-function(sgenesampling,uids){
	out = sgenesampling%>%rstan::summary(.)%>%.[[1]]%>%as.data.frame(stringsAsFactors=FALSE)%>%select(estimate = `50%`, CI.L = `2.5%`,CI.R = `97.5%`)%>%
	rownames_to_column('par')%>%mutate(params=par%>%parse_stan_pars(c('gene','time')))%>%unnest%>%
	select(parameter,ntime=time,gene,estimate,CI.L,CI.R)
		singlegene = out%>%filter(parameter=='prot')%>%.$gene%>%is_in(1)%>%all
	if(singlegene){	out$gene%<>%replace_na(1)}
	out%<>%mutate(uprotein_id=uids[gene])%>%
	safe_left_join(allow_missing=TRUE,metainfo%>%filter(!is.na(uprotein_id))%>%distinct(uprotein_id,gene_name))%>%
	# filter(parameter=='prot')%>%
	mutate(assay='MS')
	out%>%as.data.frame
	out
}

data.frame(uprotein_id=finaltestgenes)%>%left_join(metainfo)%>%filter(gene_name=='Satb2')
data.frame(uprotein_id=finaltestgenes)%>%left_join(metainfo)%>%filter(gene_name=='Flna')
finaltestgenes%>%head(2)

ribonorm_combinitvals$prot[which(finaltestgenes=='ENSMUSP00000098997_239'),]
get_dp_standata(finaltestgenes)$lMS[which(finaltestgenes==best_satb2_uid),]%>%txtplot
ribonorm_combinitvals$prot[which(finaltestgenes==best_satb2_uid),]


#parse the samplings made on individual genes
get_sampling_par_cis<-lapply(seq_along(finaltestgenes),function(i){
	parse_sampling_pars(sgenesampling=samplings[[i]]$result,uid=finaltestgenes[i])})%>%
	setNames(.,seq_along(.))%>%bind_rows(.id='genenum')%>%mutate(gene=gene+as.numeric(genenum)-1)
get_sampling_prot_cis<-get_sampling_par_cis%>%filter(parameter=='prot')
get_sampling_prot_lps<-get_sampling_par_cis%>%filter(parameter=='lp__')

ribonormopt_protdf<-tidy_opt_prots(ribonormopt$par,finaltestgenes)

summary<-rstan::summary

get_sampling_par_cis%>%filter(parameter=='l_pihalf')

ribnorm_sampling_protdf <- ribonormsamp%>%parse_sampling_prot(finaltestgenes)%>%
	filter(parameter=='prot')
	# filter(parameter=='ribnorm')
ribnorm_sampling_protdf%>%filter(parameter=='prot')
ribonormsamp%>%parse_sampling_prot(finaltestgenes)%>%filter(parameter=='ribnorm')
ribonormsamp%>%parse_sampling_prot(finaltestgenes)

stanmodels<-list(
	indivopts=combinitvals%>%tidy_opt_prots(finaltestgenes),
	indivsamplings=get_sampling_prot_cis,
	ribonormopt=ribonormopt_protdf,
	ribonormsamp=ribnorm_sampling_protdf
)
correctionmodels<-c(names(stanmodels))


best_ms_ids= metainfo$ms_id[match(best_uprotein_ids,metainfo$uprotein_id)]
for(corectionmodel in correctionmodels){
   #correct from scaled down vals up to matched_ms_mat vals
   # corval = (postmeanmat[best_uprotein_ids,]%>%apply(2,median,na.rm=T)%>%median) - matchedms_mat_rscl%>%apply(2,median,na.rm=T) 
  corval=msmed
  stopifnot(!any(stanmodels[[corectionmodel]]$gene_name%>%is.na))
  stanmodels[[corectionmodel]]$estimate <- stanmodels[[corectionmodel]]$estimate + median(corval)
  stanmodels[[corectionmodel]]$CI.L <- stanmodels[[corectionmodel]]$CI.L + median(corval)
  stanmodels[[corectionmodel]]$CI.R <- stanmodels[[corectionmodel]]$CI.R + median(corval)  
}

warning('I should first at least get the rescaling right')


for(i in seq_along(stanmodels)){ 
  message(names(stanmodels)[i])
  stopifnot(
  	c('estimate','CI.R','CI.L','gene_name','assay','ntime')%in%colnames(stanmodels[[i]])
  )
  stopifnot('Flna'%in%stanmodels[[i]]$gene_name)
}
conflict_prefer("intersect", "BiocGenerics")
lowlikgenes = get_sampling_prot_lps%>%arrange(estimate)%>%.$gene_name%>%head(2)

source('src/R/Shiny/model_app.R')


# save.image('data/run_degmodel_dropout.Rdata')
# load('data/run_degmodel_dropout.Rdata')