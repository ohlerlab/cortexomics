library(purrr)
library(here)
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
#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(assertthat))
suppressMessages(library(limma))
message('...done')

setwd(here())
# BiocManager::install('doMC')
source(here('src/R/cortexomics_myfunctions.R'))

# message('temp commented out')

transformdexprfile=here('pipeline/exprdata/transformed_data.txt')
designmatrixfile=here('pipeline/exprdata/designmatrix.txt')
registerDoMC(32)


#get the data
source(here('src/R/Modeling/stan_predict_impute.R'))


hierarchstanfile <- file.path(root,'src/Stan/degmodel_hierarch.stan') %T>%{stopifnot(file.exists(.))}
indivfitstanfile <- file.path(root,'src/Stan/degmodel_simple.stan') %T>%{stopifnot(file.exists(.))}
linearstanfile <- file.path(root,'src/Stan/degmodel_simple_linear.stan')%T>%{stopifnot(file.exists(.))}

#Funciton to fit linear and nonlinear model with stan
get_stanfit <- function(standata,stanfile,modelsamplefile,pars=NA,sampledir){
  require(tools)

  ##First copy the model file to a new folder
  tmpdir <- tempdir() 
  stopifnot(file.exists(stanfile))
  #delete the rds file
  stanfile%>%str_replace('.stan','.rds')%>%{suppressWarnings(file.remove(.))}
  # stanfile<-file.path(root,'src/Stan/degmodel_nonhierach.stan')
  tdir <-   tempdir()
  dir.create(tdir)
  file.copy(stanfile,tdir)
  modelcopy=file.path(tdir,basename(stanfile))
  stopifnot(file.exists(modelcopy))
  modelcopy


  # modelnm=basename(stanfile)%>%file_path_sans_ext%>%paste0('_n',length(standata$genes),basename(tdir))

  message(paste0('fitting model ',basename(tdir),' on ',standata$G,' genes'))
  # modelsamplefile<-paste0('stansamples_',basename(tdir),'_',basename(stanfile),'_modsamples')
  #
  n_chains <- 4

  setwd(sampledir)

  stanfit <- rstan::stan(file=modelcopy,
          model_name=basename(tdir),seed=1,
          # model_name=modelnm,seed=1,
          data=standata,
              control=list(adapt_delta=0.95,max_treedepth=15),save_dso=FALSE,
              pars=c(pars,'lrTE'),
              # init = lapply(seq_len(n_chains),function(id)list('lrTE'=array(rep(20,length(g2fit))))),
              chains=n_chains,iter=1e3,cores=n_chains,verbose=TRUE,save_warmup=FALSE,
              sample_file=modelsamplefile,
              # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
            )

    
  chainfiles <- Sys.glob(file.path(sampledir,'*_[0-9].csv'))%>%str_subset(file_path_sans_ext(basename(modelsamplefile)))
  stopifnot(length(chainfiles)==n_chains)
  

  list(fit=stanfit,chainfiles= chainfiles)

}

library(testthat)

test_that("the model Im using reliably converges for easy genes - those that are increasing",{
  #pick such a set of genes
  msincreasing <- exprdata%>%group_by(gene_name,assay,time)%>%summarise(signal=median(na.omit(signal)))%>%filter(assay=='MS')%>%group_by(gene_name)%>%summarise(is_increasing=signal[time=='E13'] < (signal[time=='P0']-2) )%>%{setNames(.$is_increasing,.$gene_name)}
  riboincreasing <- exprdata%>%group_by(gene_name,assay,time)%>%summarise(signal=median(na.omit(signal)))%>%filter(assay=='ribo')%>%group_by(gene_name)%>%summarise(is_increasing=signal[time=='E13'] < (signal[time=='P0']-2) )%>%{setNames(.$is_increasing,.$gene_name)}
  low_variance <-   exprdata%>%group_by(gene_name,assay,time)%>%filter(assay=='MS')%>%summarise(coeffvar = sd(signal)/mean(signal))%>%ungroup%>%mutate(coeffvar = rank(coeffvar)/length(coeffvar))%>%group_by(gene_name,assay)%>%summarise(lowvar = all(coeffvar<0.8))%>%{setNames(.$lowvar,.$gene_name)}
  not_missing <- exprdata%>%group_by(gene_name,assay)%>%filter(assay=='MS')%>%summarise(missing = any(is.na(signal)))%>%{setNames(!.$missing,.$gene_name)}

  test_genes <- msincreasing & riboincreasing & low_variance &   not_missing
  
  map(list(msincreasing,riboincreasing,low_variance, not_missing,test_genes),table)


  test_genes[test_genes]

  #run the model on them


  #see if they converge

  expect_true(all(is_convergent))
})




# stop()


# allhierarchfit <- get_stanfit(allstandata,hierarchstanfile,pars=hierarchpars,modelsamplefile='allhierarch.csv',sampledir=here('pipeline/stansamples'))
# allhierarchfit%>%saveRDS(here('data/allhierarchfit.rds'))

# alllinfit <- get_stanfit(allstandata,linearstanfile,pars='lp__',modelsamplefile='alllinear.csv',sampledir=here('pipeline/stansamples'))
# alllinfit%>%saveRDS(here('data/alllinfit.rds'))

allhierarchfit_rna <- get_stanfit(allstandata_rna,hierarchstanfile,pars=hierarchpars,modelsamplefile='allhierarch_rna.csv',sampledir=here('pipeline/stansamples'))
allhierarchfit_rna%>%saveRDS(here('data/allhierarchfit_rna.rds'))

alllinfit_rna <- get_stanfit(allstandata_rna,linearstanfile,pars='lp__',modelsamplefile='alllinear_rna.csv',sampledir=here('pipeline/stansamples'))
alllinfit_rna%>%saveRDS(here('data/alllinfit_rna.rds'))


# testhierarchfit <- get_stanfit(singlesetstandata,hierarchstanfile,pars=hierarchpars,modelsamplefile='stansamples/testhierarch.csv')
testhierarchfit <- get_stanfit(testsetstandata,hierarchstanfile,pars=hierarchpars,modelsamplefile='stansamples/testhierarch.csv')
testhierarchfit%>%saveRDS('../data/testhierarchfit.rds')
testhierarchfit<-readRDS('../data/testhierarchfit.rds')

testlinfit <- get_stanfit(testsetstandata,linearstanfile,pars=NULL,modelsamplefile='stansamples/testlin.csv',sampledir=here('pipeline/stansamples'))
testlinfit%>%saveRDS(here('data/testlinfit.rds'))






# allstanfits <- mclapply(genes2fit,safely(getstanfits),genelengths,exprdata,lengthnorm=FALSE)


dir.create(paste0('stanfits/'),showWarnings=F)


# message('testing fit')
# # modob<-safely(getstanfits)('Satb2',genelengths,exprdata,pars=c('lrTE','ldeg','prot'),lengthnorm=FALSE)
# stopifnot(is.null(modob$error))
# saveRDS(modob, paste0('stanfits/','Satb2','nonh_vs_lin.stanfit.rds'))
# message( paste0('stanfits/','Satb2','nonh_vs_lin.stanfit.rds'))
# modob$result$genes

#Now fit all models together
allstanfitssim <- safely(getstanfits)(genes2fit,genelengths,exprdata,lengthnorm=TRUE,)

#Use do par to write the inidividual model fits to disc
parrun<-foreach(g=genes2fit) %dopar%{
  suppressMessages({modob<-safely(getstanfits)(g,genelengths,exprdata,lengthnorm=TRUE)})
  saveRDS(modob, paste0('stanfits/',g,'nonh_vs_lin.stanfit.rds'))
  message(paste0(which(genes2fit==g), '/',length(genes2fit)))
}


source('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/R/modeling/degredation_trajplot.R')

stop()

allstanfitssim[[1]]$kinetic%>%
  summary%>%
  .$summary%>%
  as.data.frame%>%
  rownames_to_column('par')%>%
  mutate(ppar=parse_stan_pars(par))%>%unnest%>%filter(is.na(time),is.na(gene))%>%
  mutate(parameter,time,gene,mean,`2.5%`,`97.5%`)

stop()


# allstanfitssim%>%saveRDS('../data/allstanfitssim.rds')




###########PRocessing results

fitfiles <- Sys.glob(paste0('stanfits/*nonh_vs_lin.stanfit.rds'))
fitfiles%<>%setNames(.,basename(.)%>%str_replace('nonh_vs.*',''))

####Now read in the individual fits
allstanfits <- mclapply(fitfiles,readRDS)
allstanfits%<>%setNames(names(fitfiles))

allstanfitsln <- mclapply(genes2fit,safely(getstanfits),genelengths,exprdata,lengthnorm=TRUE)
allstanfitsln%<>%setNames(names(fitfiles))

#save larger objects
allstanfits%>%saveRDS('../data/allstanfits.rds')
allstanfitsln%>%saveRDS('../data/allstanfitsln.rds')


allstanfits%>%map('result')%>%
  map('kinetic')%>%
  map(summary)%>%map(1)%>%object.size%>%divide_by(1e6)

stop()


#Now, don't store all in memory, just save the summary
allstanfitsumm<-fitfiles%>%mclapply(safely(.%>%readRDS(.)%>%.$result%>%.$kinetic%>%summary%>%.[[1]]%>%as.data.frame))
stopifnot(allstanfitsumm%>%map('error')%>%map_lgl(is.null)%>%all)#they all worked
allstanfitsumm%<>%map('result')

####Now read in the list of indiivdual fit summaries
allstanfitsumm%>%saveRDS('../data/allstanfitsumm.rds')
allstanfitsumm<-readRDS('../data/allstanfitsumm.rds')

# stopifnot(allstanfitsumm%>%map_lgl(is,'matrix')%>%all)

allsummary <- allstanfitsumm%>%
  map(as.data.frame)%>%
  map(rownames_to_column,'par')%>%
  bind_rows(.id='gene_name')%>%
  mutate(ppars=parse_stan_pars(par))%>%
  unnest%>%
  select(-gene)

allsummarysim <- allstanfitssim$result[[1]]%>%summary%>%.[[1]]%>%as.data.frame%>%
  rownames_to_column('par')%>%
  mutate(ppars=parse_stan_pars(par))%>%
  unnest%>%
  mutate(gene_name=genes2fit[gene])


hallsummary <- allstanfitssim$result$kinetic%>%
  summary%>%.[[1]]%>%as.data.frame%>%
  rownames_to_column('par')%>%
  mutate(ppars=parse_stan_pars(par))%>%
  unnest%>%
  mutate(gene_name=allstanfitssim$result$genes[gene])



#Plot distributipon of rTEs

#about a fifth of them had convergenet TE
stopifnot(allsummary%>%filter(parameter=='lrTE')%>%
  .$sd %>%`>`(2.5)%>%mean%>%between(0.1,0.3))





stop()
allstanfits%>%map


allstanfits$lrTE


realdata

fitfiles <- Sys.glob(paste0('stanfits/*nonh_vs_lin.stanfit.rds'))


gnamei<-gnamesi

#Getting the log likelihood 
getloglikratios <- function(fitfile){
  fit<-fitfile%>%readRDS
  realstanfit <- fit%>%.$result%>%.$kinetic
  realstanfit_lin <- fit%>%.$result%>%.$linear
  fitdata<-fit%>%.$result%>%.$data
  kin=extract(realstanfit,'lp__')%>%.[[1]]%>%min
  lin=extract(realstanfit_lin,'lp__')%>%.[[1]]%>%min
  kin-lin
}

loglikratios <- fitfiles%>%str_subset('Flna')%>%mclapply(mc.cores=20,getloglikratios)


##Now with fixed values
# realstanfit_test <- rstan::stan(file='src/Stan/degmodel_simple.stan',data=realdata,
#                                control=list(adapt_delta=0.95,max_treedepth=20),
#                                init = list(list(lrTE=array(0,dim=c(1)))),
#                                chains=1,iter=1,warmup = 0,
#                                # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
#                                verbose=TRUE)


stop()

# 
# stanfit <- rstan::stan(file='src/degmodel.stan',data=list(G=n_genes,T=tps,K=K,
#                         MS=ms_array,ribo=ribo_mat),
#                        control=list(adapt_delta=0.90,max_treedepth=10),
#                        chains=4,iter=2e3,
#                        # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
#                        verbose=TRUE)


##########Examining stanfits

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

hist(degs,50)





library(rstan)
funnel <- stan_demo("funnel", seed = 12345)   # has 5 divergent transitions
pairs(funnel, pars = c("y", "x[1]", "lp__"), las = 1) # below the diagonal
funnel_reparam <- stan_demo("funnel_reparam") # has no divergent transitions
pairs(funnel_reparam, pars = c("y", "x[1]", "lp__"), las = 1) # below the diagonal


#ldeg_rTE parameter did NOT help, it just blows up to minus infinity


#' Notes on nodel fitting.... 
#' Looks like scaling of the standard deviation makes things work just fine
#' Now I"ll try putting things on a log scale
#' Okay so applying a prior like that to transformed parameters is bad (warning about needing to add in the jacobian of the transform)
#' so I now have rTE and tau getting initilialized on a log scale.
#' Combining the two genes Isoc1 seems to make the second one kind of divergent... (low Rhat, neff)
#' 
#' 
#' Open questions - how many things have flatly contradictory evidence for the rTE? - now with the new data have to recheck
#' 
#' I need to show that the similtaneous model doesn't give different fits to the single ones, or at least not too badly.
#' 
#' I need to identify a couple of test genes - e.g. the one with contradictory TE, a few with vauge TE, and a few with certain TE. Then I should put them in a big model together
#' First without contradictory one, then with, then with plus a mixture model component if necessary
#' Then I should introduce some library specific scaling factors.
#' 
#' I should bear in mind there could be some gene name - id conflicts....
#' 
#' Looking at the fits - a kind of consistent pattern seems to be emerging where the linear fit looks okay just.... understated. This would not I think be fixed by simple coefficients for the timepoints,
#' Because It's in opposite direction depending on time.
#' The things I"m wondering about are a) Should I average between timepoints (seems like obviously yes, actually) and should I try restricting the rTE to the reasonable value space?
#' Let's try tp restriction with a
#' 
#' NOt much discernable difference from the tp averaging. Didn't hurt, at least
#' 
#' Now to try the blunt prior on rTE - 
#' Programming this is a nightmare, as it turns out. Getting initialization to work ont he stan object was hard, but getting it to actually recompile teh model whenever I run was harder... still haven't.
#'
#' Now the hierarchical model, which... seems to fit faster than the all the independent ones??? 
#' Okay so the model object fits at least.....
#'
#' Need to compare the hierarch and non-hierarch models
#' So need to be able to run the hierarch models 
#' 
#' It would be good if I could compare the RNA and the Riboseq based model - particularly for TE based cases
#' 
#' Also look at our top stalling candidates.