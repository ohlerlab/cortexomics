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
library(testthat)
library(parallel)
#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(assertthat))
suppressMessages(library(limma))
message('...done')

conflict_prefer("accumulate", "purrr")
conflict_prefer("setdiff", "BiocGenerics")
conflict_prefer("between", "dplyr")
conflict_prefer("clusterApply", "BiocGenerics")

setwd(here())

source(here('src/R/Functions/cortexomics_myfunctions.R'))
source(here('src/R/Rprofile.R'))
# message('temp commented out')

#processedexprdatafile=here('pipeline/exprdata/transformed_data.txt')
#designmatrixfile=here('pipeline/exprdata/designmatrix.txt')
registerDoMC(32)


hierarchCIstanfile <- here('src/Stan/degmodel_dataconfint_hierarch.stan') %T>%{stopifnot(file.exists(.))}
indivCIstanfile <- here('src/Stan/degmodel_dataconfint.stan') %T>%{stopifnot(file.exists(.))}
indivCI_lin_stanfile <- here('src/Stan/degmodel_dataconfint_linear.stan') %T>%{stopifnot(file.exists(.))}
st_indivCI_stanfile <- 'src/Stan/st_model_confint.stan'



processedtbl<-read_tsv(here('pipeline/exprdata/limma_proDD_CIs.tsv'))
metainfo<-read_tsv(here('pipeline/exprdata/limma_genemetadata.tsv'))


dir.create('exprdata',showWarnings = FALSE)
# exprtbl <- read_tsv(processedexprdatafile) 

# processedtbl<-exprtbl%>%gather(dataset,signal,-gene_name)%>%separate(dataset,into=c('time','assay','rep'))%>%
#   group_by(gene_name,time,assay)%>%
#   mutate(sd=sd(na.omit(signal)),signal=mean(na.omit(signal)))%>%
#   group_by(gene_name)%>%
#   filter(!any(sd==0))

# processedtbl %<>% select(gene_name, everything())
# processedtbl

test_that("we have datathat looksl ike the old but with confidence intervals",{
  expect_gt(nrow(processedtbl),10)
  data_has_confints = 'sd' %in% colnames(processedtbl)
  points_per_gene <- processedtbl%>%group_by(gene_name)%>%tally%>%.$n
  expect_true(all(points_per_gene>0))
  expect_equal(n_distinct(points_per_gene),1)
  expect_true(data_has_confints) 

})

g2fit="Satb2"

get_standata_confint<-function(g2fit,genelengths,exprdata,lengthnorm,assay2get='ribo'){
  require(splines2)
  stopifnot(g2fit%in%exprdata$gene_name)

  MS_array <- exprdata%>%
    filter(gene_name%in%g2fit)%>%
    filter(assay=='MS')%>%
    arrange(match(gene_name,g2fit),time)%>%
    .$signal%>%
    array(dim=c(5,length(g2fit)))%>%
    aperm(c(1,2))%>%
    identity
    # {2^.}

  MS_tau <- exprdata%>%
    filter(gene_name%in%g2fit)%>%
    filter(assay=='MS')%>%
    arrange(match(gene_name,g2fit),time)%>%
    .$sd%>%
    array(dim=c(5,length(g2fit)))%>%
    aperm(c(1,2))%>%
    identity
    # {2^.}

  ribo_array <- exprdata%>%
    filter(gene_name%in%g2fit)%>%
    filter(assay== !!assay2get)%>%
    arrange(match(gene_name,g2fit),time)%>%
    .$signal%>%
    array(dim=c(5,length(g2fit)))%>%
    aperm(c(1,2))%>%
    identity
    # {2^.}

  ribo_tau <- exprdata%>%
    filter(gene_name%in%g2fit)%>%
    filter(assay== !!assay2get)%>%
    arrange(match(gene_name,g2fit),time)%>%
    .$sd%>%
    array(dim=c(5,length(g2fit)))%>%
    aperm(c(1,2))%>%
    identity
    # {2^.}

  #Use gene lengths to get DENSITY - we care about this, not the total counts
  if(lengthnorm) ribo_array%<>%sweep(2,STATS=log2(genelengths[g2fit]/1e3),FUN='-')
  # if(lengthnorm) ribo_tau%<>%sweep(2,STATS=log2(genelengths[g2fit]/1e3),FUN='-')

  #timepoint averaging - we'll leave in the first, unused timepoint...
  ribo_array_old<-ribo_array

  for(i in 2:nrow(ribo_array)){
    ribo_array[i,] <- (ribo_array_old[i,]+ribo_array_old[i-1,])/2
  }

  standata <- list(
     G=length(g2fit),T=length(tps),K=3,
     # gene_name=gnamei,
     MS=MS_array,
     MS_tau=MS_tau,
     # ribo =exprdata%>%filter(gene_name==gnamei)%>%filter(assay=='ribo',rep==1)%>%.$predicted_signal%>%{2^.}%>%rollmean(k=2)%>%matrix
     ribo =ribo_array,
     ribo_tau =ribo_tau 
  )
  dimnames(standata$MS)=list(tps,unique(g2fit))
  dimnames(standata$MS_tau)=list(tps,unique(g2fit))
  dimnames(standata$ribo)=list(tps,unique(g2fit))
  dimnames(standata$ribo_tau)=list(tps,unique(g2fit))
  
  time=1:length(tps)
  timeknots <- time[c(-1,-length(time))]
  mybs <- cbind(1,ibs(time, knots = timeknots,degree = 1, intercept = TRUE))
  mydbs = bs(time, knots = timeknots,degree = 1, intercept = TRUE)

  gsetname <- if(length(g2fit)==1)g2fit[1] else paste0('gset',substr(digest::digest(letters,len=10),1,10))
  
  standata<-list(
    gsetname=gsetname,
    G=length(g2fit),          #  genes
    T=length(tps),          #  timepoints
    lMS=MS_array,  # mass spec data mean
    lMS_tau=MS_tau,  # mass spec data precision
    lribo=ribo_array, # lriboseq (synthesis) data mean
    lribo_tau=ribo_tau, # lriboseq (synthesis) data sd
    mybs=mybs,
    mydbs=mydbs
  )
  standata
}
# test_that("the model Im using reliably converges for easy genes - those that are increasing",{
  #pick such a set of genes
  # msincreasing <- processedtbl%>%group_by(gene_name,assay,time)%>%summarise(signal=median(na.omit(signal)))%>%filter(assay=='MS')%>%group_by(gene_name)%>%summarise(is_increasing=signal[time=='E13'] < (signal[time=='P0']-2) )%>%{setNames(.$is_increasing,.$gene_name)}
  # riboincreasing <- processedtbl%>%group_by(gene_name,assay,time)%>%summarise(signal=median(na.omit(signal)))%>%filter(assay=='ribo')%>%group_by(gene_name)%>%summarise(is_increasing=signal[time=='E13'] < (signal[time=='P0']-2) )%>%{setNames(.$is_increasing,.$gene_name)}
  # low_variance <-   processedtbl%>%group_by(gene_name,assay,time)%>%filter(assay=='MS')%>%summarise(coeffvar = sd(signal)/mean(signal))%>%ungroup%>%mutate(coeffvar = rank(coeffvar)/length(coeffvar))%>%group_by(gene_name,assay)%>%summarise(lowvar = all(coeffvar<0.8))%>%{setNames(.$lowvar,.$gene_name)}
  # not_missing <- processedtbl%>%group_by(gene_name,assay)%>%filter(assay=='MS')%>%summarise(missing = any(is.na(signal)))%>%{setNames(!.$missing,.$gene_name)}

  # test_genes <- msincreasing & riboincreasing & low_variance &   not_missing
  
  # map(list(msincreasing,riboincreasing,low_variance, not_missing,test_genes),table)


  # expect_true(all(is_convergent))
# })



#Funciton to fit linear and nonlinear model with stan
get_stanfit <- function(standata,stanfile,pars=NA,sampledir,iter=200, n_chains=4,stanpars_to_listi=stanpars_to_list){

stanpars_to_listi(paramvals[[1]]%>%.[1,]%>%names)

  require(tools)
  ##First copy the model file to a new folder
  tmpdir <- file.path(tempdir(),tempfile())
  stopifnot(file.exists(stanfile))
  #delete the rds file
  stanfile%>%str_replace('.stan','.rds')%>%{suppressWarnings(file.remove(.))}
  # stanfile<-file.path(root,'src/Stan/degmodel_nonhierach.stan')
  system(str_interp('rm -rf ${tmpdir}'))
  dir.create(rec=TRUE,tmpdir)
  #
  stopifnot(file.copy(stanfile,tmpdir))
  modelcopy=file.path(tmpdir,basename(stanfile))
  stopifnot(file.exists(modelcopy))
  message(normalizePath(stanfile))
  message(normalizePath(modelcopy))
  # modelnm=basename(stanfile)%>%file_path_sans_ext%>%paste0('_n',length(standata$genes),basename(tmpdir))
  message(paste0('fitting model ',basename(tmpdir),' on ',standata$G,' genes'))
  # modelsamplefile<-paste0('stansamples_',basename(tmpdir),'_',basename(stanfile),'_modsamples')
  ngenes <- standata$G
  #
  message('optimizing with newton to get initial values')
  optpars<-optimizing(stan_model(stanfile),data=standata,algorithm = c("Newton"),hessian=FALSE)
  #
  stanparvect <- optpars$theta_tilde
  # initvals = list(
  #     cM = pmax(standata$lribo,0.01),
  #     Kd = array( rep(4,ngenes) ,dim=ngenes),#i.e. a 3 fold reduction due to degredation per tp
  #     prot0 = array(standata$lMS[1,],dim=ngenes)
  # )
   #
  # initvals$lKs = (2^initvals$prot0)*initvals$Kd / (2^standata$lribo[1,])
  #
  initvals=list(optpars$theta_tilde[1,]%>%stanpars_to_list)
  #
  initvals = rep(initvals,  n_chains)
  #
  #model name
  imodel_name <- basename(tmpdir)
  #path of sample files and model
  modelfilename<-tools::file_path_sans_ext(basename(stanfile))
  samplefilepath <- file.path(sampledir,paste0(modelfilename,'_',standata$gsetname))
  #get chain files, clear old ones
  chainfiles <- paste0(samplefilepath,'_',1:n_chains,'.csv')
  try({file.remove(chainfiles)})
  message(paste0('fitting full model ',basename(tmpdir),' on ',standata$G,' genes'))
  #fit model
  stanfit <- stan(file=modelcopy,
          model_name=imodel_name,
          # model_name=modelnm,seed=1,
          data=standata,
              control=list(adapt_delta=0.99,max_treedepth=20),save_dso=FALSE,
              # pars=c(pars,'lrTE'),
              init =   initvals,
              # chains=n_chains,
              chains=  n_chains,
              iter=iter,
              cores=n_chains,
              verbose=FALSE,
              save_warmup=FALSE,
              sample_file=samplefilepath
            )
  stopifnot(all(file.exists(chainfiles)))
  list(fit=stanfit,chainfiles= chainfiles)
}

{
  require(BiocParallel)
  # BiocManager::install('batchtools')
  param<- BatchtoolsParam(workers=4, cluster="sge", resources=list(queue='all'))
  sge_template <- '~/tmp.tmpl' 
  param$template%>%readLines%>%str_replace(regex('\\#\\$ -q .*'),'')%>%cat(file=sge_template,sep='\n')
  param<- BatchtoolsParam(workers=50, cluster="sge", resources=list(queue='all'),template=sge_template)
  #bplapply(BPPARAM=param
}

test_gene_names <-  metainfo%>%filter(sig_MS_change)%>%.$gene_name
genelengths = metainfo$width%>%setNames(metainfo$gene_name)
mygenes2fit <- c('Satb2','Flna','Gamt',test_gene_names)

splitdata <- mclapply(mc.cores=20,mygenes2fit,mymemoise(get_standata_confint),genelengths,processedtbl,lengthnorm=T,assay2get='ribo')

length(splitdata)
#First lets get the optimal fit for 
indivCIstanmodel<-mymemoise(stan_model)(indivCIstanfile)
indivCI_lin_stanmodel<-mymemoise(stan_model)(indivCI_lin_stanfile)
st_indivCIstanmodel<-mymemoise(stan_model)(st_indivCI_stanfile)


################################################################################
########Messing with optimization
################################################################################
  
# nlopt<-vb(indivCIstanmodel,data=splitdata[[2]],algorithm = c("meanfield"))
# linopt<-vb(indivCI_lin_stanmodel,data=splitdata[[2]],algorithm = c("fullrank"))
# save.image('data/degredation_confint_data.R')

nlopt<-optimizing(indivCIstanmodel,data=splitdata[[2]],algorithm = c("Newton"),hessian=TRUE)
linopt<-optimizing(indivCI_lin_stanmodel,data=splitdata[[2]],algorithm = c("Newton"),hessian=TRUE)

st_opt<-optimizing(st_indivCIstanmodel,data=splitdata[[2]],algorithm = c("Newton"),hessian=TRUE)


goptsetnames <- splitdata%>%map_chr('gsetname')
splitdata%<>%setNames(goptsetnames)

st_opts <- tryCatch({
    readRDS(here('data/st_opts.rds'))
  },error=function(e){
     st_opts <- mclapply(splitdata[goptsetnames], safely(function(x){cat('.');optimizing(st_indivCIstanmodel,data=x,algorithm = c("Newton"),hessian=TRUE)}))
     st_opts%>%saveRDS(here('data/st_opts.rds'))
     st_opts
  })
# st_opts <- mclapply(splitdata[goptsetnames], safely(function(x){cat('.');optimizing(st_indivCIstanmodel,data=x,algorithm = c("Newton"),hessian=TRUE)}))
# st_opts%>%saveRDS(here('data/st_opts.rds'))

paramstderrs <- st_opts%>%map('result')%>%map('hessian')%>%map(safely(.%>%multiply_by(-1)%>%solve%>%diag%>%sqrt))%>%map('result')%>%setNames(goptsetnames)
paramvals <- st_opts%>%map('result')%>%map('theta_tilde')%>%setNames(goptsetnames)

paramvals[[1]]%>%.[1,]%>%stanpars_to_list

parnames <- iparnames%>%str_remove('(\\[|\\.)[0-9\\,]+(\\]|\\.)$')
parinds <- iparnames%>%str_extract('(\\[|\\.)[0-9\\,]+(\\]|\\.)$')%>%str_extract_all(regex('\\d+'))%>%map(as.numeric)
innummax <- parinds%>%map_dbl(length)%>%max

matrix(NA,nrow=length(parinds),ncol = )

stanmodel_hes_stder_df<-paramvals%>%map_df(.id='gene_name',~enframe(.[1,],'parameter','value'))%>%
  mutate(parameter = str_replace_all(parameter,'[\\[\\],](?=.)','.'))%>%
  mutate(parameter = str_replace_all(parameter,']$',''))%>%
  inner_join(
    paramstderrs%>%map_df(.id='gene_name',enframe,'parameter','stderr')
  )

stanmodel_hes_stder_df%<>%  mutate(myhgrp = case_when(
    value < -5 ~ 'low',
    value > 5 ~ 'high',
    TRUE ~ 'middle' 
  ))

pdf('tmp.pdf')
stanmodel_hes_stder_df%>%filter(parameter=='l_st.1')
stanmodel_hes_stder_df%>%filter(parameter=='l_pihalf.1')
dev.off()

stanmodel_hes_stder_df%>%filter(parameter=='l_pihalf.1')%>%filter(value>0)%>%arrange(desc(value))


mcshanedf<-fread('ext_data/mcshane_etal_2016_S1.csv')
#
mcshanethalfs<-mcshanedf%>%select(2,38,41)%>%set_colnames(c('gene_name','half_life','McShane_deg_cat'))
#
mcshanethalfs$half_life%<>%str_replace('> 300','300')%>%as.numeric
mcshanethalfs$half_life %<>% {./24}

mcshanethalfs%>%filter(McShane_deg_cat=='ED')%>%.$half_life%>%.[.<12]%>%txtdensity()
mcshanethalfs%>%filter(McShane_deg_cat=='ED')%>%.$half_life%>%log10%>%.[.<24]%>%txtdensity()

#
library(txtplot)
#
pdf('tmp.pdf')
stanmodel_hes_stder_df%>%filter(parameter=='l_pihalf.1')%>%
  mutate(value = value)%>%
  left_join(mcshanethalfs)%>%
  mutate(half_life=half_life)%>%
  {qplot(data=.,color=McShane_deg_cat,.$value,log(.$half_life),ylab = 'McShane estimated Thalf',xlab='Model_estimated_thalf')+theme_bw()}
dev.off()
normalizePath('tmp.pdf')
#
mcshanetestdf<-stanmodel_hes_stder_df%>%filter(parameter=='l_pihalf.1')%>%
  mutate(value = value)%>%
  left_join(mcshanethalfs)%>%
  mutate(half_life=half_life)%>%
  mutate(mcshnhigh = ifelse(half_life==max(half_life,na.rm=T),'high','nothigh'))

mcshanetestdf%>%
  {table(.$mcshnhigh,.$myhgrp)}%>%
  chisq.test

mcshanetestdf%>%
  {table(.$McShane_deg_cat,.$myhgrp)}%>%
  chisq.test




stanmodel_hes_stder_df

mcshanedf%>%colnames

# pdf('tmp.pdf')
# plot(st_opt)
# dev.off()


# nlopt$value
# linopt$value
# paramint<-c('lKs.1','Kd.1')
# nlopt$hessian[paramint,paramint]


# library(txtplot)

# nlprotvals<-nlopt$theta_tilde[1,]%>%enframe('par','value')%>%  mutate(ppar=parse_stan_pars(par))%>%unnest%>%
#   filter(parameter%in%c('prot'))%>%.$value%T>%txtplot

# linprotvals <- linopt$theta_tilde[1,]%>%enframe('par','value')%>%  mutate(ppar=parse_stan_pars(par))%>%unnest%>%
#   filter(parameter%in%c('prot'))%>%.$value%T>%txtplot

# stprotvals <- st_opt$theta_tilde[1,]%>%enframe('par','value')%>%  mutate(ppar=parse_stan_pars(par))%>%unnest%>%
#   filter(parameter%in%c('prot'))%>%.$value%T>%txtplot


        




indivclustfits<-bplapply(BPPARAM=param,splitdata,get_stanfit=get_stanfit,stanfile=st_indivCI_stanfile,
# indivclustfits<-lapply(splitdata[,get_stanfitfun=get_stanfit,stanfile_i=st_indivCI_stanfile,
    function(gsetstandata,get_stanfitfun,stanfile_i){
  #stindivclustfits<-lapply(splitdata[1:2],get_stanfit=get_stanfit,st_indivCI_stanfile=st_indivCI_stanfile,function(gsetstandata,get_stanfit,st_indivCI_stanfile){
  MSmed = median(gsetstandata$lMS) 
  gsetstandata$lMS <- gsetstandata$lMS - MSmed
  #
  Ribomed = median(gsetstandata$lribo) 
  gsetstandata$lribo <- gsetstandata$lribo - Ribomed
  #
  #
  require(rstan)
  require(tools)
  require(magrittr)
  require(tidyverse)
  require(here)
  #get teh data  #
  testindivfit<-get_stanfitfun(
    gsetstandata,
    stanfile_i,
    pars=NULL,iter=1000,n_chains=4,
    sampledir=here('pipeline/stansamples')
  )
})


indivclustfits

stop()


sampledir%>%list.files(full=TRUE)
# sampledir%>%list.files(full=TRUE)%>%file.remove

#function to extract parameters from list of stanfit objects and ploti
plot_sample_density<-function(indivclustfits,par='lKs'){
  indivclustfits%<>%setNames(mygenes2fit[1:length(indivclustfits)])
  vals = indivclustfits %>%map_df(.id='gene',function(x)x[[1]]@sim$samples[[1]][['lKs[1]']]%>%data.frame(val=.))
  qplot(data=vals,color=gene,fill=gene,x=val,geom='blank')+geom_histogram()+theme_bw()
}

#now plot
pdf('tmp.pdf')
# plot_sample_density(indivclustfits,'lKs')
plot_sample_density(indivclustfits,'')
dev.off()
normalizePath('tmp.pdf')

indivclustfits[[5]]

stop()

indivclustfits[[2]][[1]]@sim$samples[[2]][['lKs[1]']]%>%max

#see if they converge
testhierarchfit <- get_stanfit(
  testsetstandata,
  pars=NULL,iter=400,n_chains=4,
  modelsamplefile='stansamples/testhierarch.csv',
  sampledir=here('pipeline/stansamples')
)




testhierarchfit[[1]]%>%summary%>%.[[1]]%>%head

parsum%>%filter(is.na(gene))

parsum<-testhierarchfit[[1]]%>%
  summary%>%
  .$summary%>%
  as.data.frame%>%
  rownames_to_column('par')%>%
  mutate(ppar=parse_stan_pars(par))%>%unnest
  # filter(is.na(time),is.na(gene))%>%
  # filter(par%>%str_detect('rTE'))
  parsum$par%>%unique


parsum%>%filter(is.na(gene))
parsum%>%filter(parameter=='lKs')
parsum%>%filter(parameter=='Kd')

parsum%>%filter(gene==1)%>%filter(parameter=='mRNA')
  standata$ribo[,1]
  parsum%>%filter(gene==2)%>%filter(parameter=='mRNA')
  standata$ribo[,2]
  parsum%>%filter(gene==1)%>%filter(parameter=='prot')
  standata$MS[,1]
  parsum%>%filter(gene==2)%>%filter(parameter=='prot')
  standata$MS[,2]

  #Okay so if I don't use prot lik, the mRNA variable converges correctly
#  parsum%>%filter(is.na(time))

  standata
  standata

  parsum%>%head
  Rhat_lrTE<-parsum%>%filter(par=='hmu_lrTE')%>%.$Rhat
  
  expect_lt(  Rhat_lrTE,1.1)
  expect_gt(  Rhat_lrTE,0.9)

  parsum%>%filter(par%>%str_detect('lrTE|ldeg'))
  parsum%>%filter(par%>%str_detect('rTE'))
#  parsum%>%filter(par%>%str_detect('prot'))
  


# })
#' The first four genes seem to converge okay
#' totally converge if indep and I do 1e3 
#' lognormal dist...
#' Asdding in others of our test genes do not
#' I need to look at trajectories.


stop()
test_that("our various improvements to the expression data do in fact improve it...")
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


source('src/R/modeling/degredation_trajplot.R')

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
