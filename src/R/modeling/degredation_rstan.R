library(purrr)
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


root <- '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/'

get_prot_samples<-function(fit) fit %>%as.data.frame%>%select(matches('prot'))%>%mutate(sample=1:nrow(.))%>%gather(par,value,-sample)%>%mutate(ppars=parse_stan_pars(par))%>%unnest%>%
  filter(parameter=='prot')%>%select(time,value,sample,gene)

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

# message('temp commented out')

transformdexprfile=file.path('exprdata/transformed_data.txt')
designmatrixfile=file.path('exprdata/designmatrix.txt')



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


########Now, let's compare models of different complexity levels
designmatrix$ribo <- designmatrix$assay %in% c('ribo','MS')
designmatrix$MS <- designmatrix$assay %in% c('MS')

design = model.matrix( ~ time + assay + time:assay , designmatrix, xlev = list(assay = c('total','ribo','MS')) )

exprdata <- exprmatrix%>%
  set_colnames(designmatrix$dataset)%>%
  as.data.frame%>%
  rownames_to_column('gene_name')%>%
  gather(dataset,signal,-gene_name)%>%
  left_join(designmatrix)

increasinggenes <- exprdata%>%
  group_by(gene_name)%>%filter(assay=='MS')%>%
  mutate(increasingMS = (mean(signal[time=='P0']) - mean(signal[time=='E13'])) > 2 )%>%
  filter(increasingMS)%>%.$gene_name

decreasinggenes <- exprdata%>%
  group_by(gene_name)%>%filter(assay=='MS')%>%
  mutate(decreasingMS = (mean(signal[time=='P0']) - mean(signal[time=='E13'])) < -2 )%>%
  filter(decreasingMS)%>%.$gene_name


#Check we aren't using hte wrongly scaled data.
stopifnot(exprdata%>%filter(gene_name=='Satb2',assay=='ribo',time=='E13')%>%.$signal%>%`<`(10))

limmafits <- list()
designmatrix$time %<>% as_factor

#fit the full model
limmafits[['full']] = limma::lmFit(exprmatrix,
                                   design=model.matrix( ~ time*(ribo+MS), designmatrix)
)

#with 4 splines
limmafits[['spline_4']] = limma::lmFit(exprmatrix,
                                        design=model.matrix( ~ ns(as.numeric(time),4)*(ribo+MS), designmatrix)
)

#with fewer splines
limmafits[['spline_3']] = limma::lmFit(exprmatrix,
                                       design=model.matrix( ~ ns(as.numeric(time),3)*(ribo+MS), designmatrix)
)

preddf <- exprdata %>%select(gene_name,dataset,signal)%>%
  left_join(getlimmapredictions(limmafits[['full']],'full',designmatrix)) %>%
  left_join(getlimmapredictions(limmafits[['spline_4']],'spline_4',designmatrix))%>%
  left_join(getlimmapredictions(limmafits[['spline_3']],'spline_3',designmatrix))


####Quick and dirty imputation with splines
exprdata%<>%left_join(preddf%>%select(gene_name,time,assay,predicted_signal_spline_3))%>%
  mutate(signal = ifelse(is.na(signal),predicted_signal_spline_3,signal))%>%
  select(-predicted_signal_spline_3)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)






#sketch
#we need a 
library(zoo)
# 
tps <- c('E13','E145','E16','E175','P0')
#
gnamei <- 'Satb2'
#
#increasing genes
gnamesi <- c('Satb2')
#decreasing genes
'Isoc1'
# gnamesi <- c('Satb2')
gnamesi <- c('Isoc1')
#
gnamesi <- exprdata%>%.$gene_name%>%sample(7)
#
gnamesi <- 'Satb2'
g2fit <- gnamesi



#named vector containg gene legnths
genelengths <- fread('feature_counts/data/E13_ribo_1/feature_counts')%>%
  select(Geneid,Length)%>%
  left_join(fread('ids.txt'),by=c(Geneid='gene_id'))%>%
  arrange(sample(1:n()))%>%
  group_by(gene_name)%>%
  slice(1)%>%
  {setNames(.$Length,.$gene_name)}


max_cores <- 32
g2fit<-exprdata$gene_name%>%unique%>%head(3)
lengthnorm=TRUE



g2fit <- c("Acadvl", "Ak1", "Asna1", "Cald1", "Cst3", "Ctnnd2", "Cul2",
"Dclk1", "Dpysl3", "Epb41l1", "Flna", "Hspa12a", "Igsf21", "Nos1",
"Orc3", "Phactr4", "Rap1b", "Rasa3", "Rps6ka5", "Satb2", "Srcin1",
"Stx1b", "Tbc1d24", "Trmt1l", "Zc2hc1a", "Zmym3")
pars=NA
g2fit<-exprdata$gene_name%>%unique
# g2fit<-'Satb2'

#Funciton to fit linear and nonlinear model with stan
getstanfits <- function(g2fit,genelengths,exprdata,pars=NA,lengthnorm=TRUE,...){

  MS_array <- exprdata%>%
    filter(gene_name%in%g2fit)%>%
    filter(assay=='MS')%>%
    arrange(match(gene_name,g2fit),time,rep)%>%
    .$signal%>%
    array(dim=c(3,5,length(g2fit)))%>%
    aperm(c(1,2,3))%>%
    {2^.}

  riboarray<-preddf%>%
    filter(gene_name%in%g2fit)%>%
    filter(assay=='ribo')%>%
    arrange(match(gene_name,g2fit),time)%>%
    .$predicted_signal_full%>%
    {2^.}%>%
    matrix(ncol=length(g2fit))

  #Use gene lengths to get DENSITY - we care about this, not the total counts
  if(lengthnorm) riboarray%<>%sweep(2,STATS=genelengths[g2fit],FUN='/')

  #timepoint averaging - we'll leave in the first, unused timepoint...
  riboarray_old<-riboarray
  for(i in 2:nrow(riboarray)){
    riboarray[i,] <- (riboarray_old[i,]+riboarray_old[i-1,])/2
  }

  standata <- list(
     G=length(g2fit),T=length(tps),K=3,
     # gene_name=gnamei,
     MS=MS_array,
     # ribo =exprdata%>%filter(gene_name==gnamei)%>%filter(assay=='ribo',rep==1)%>%.$predicted_signal%>%{2^.}%>%rollmean(k=2)%>%matrix
     ribo =riboarray 
  )


  #Fit non-hierarchical model 
  modelnm=paste0(g2fit[1],'_',length(g2fit))
  file.path(root,'src/Stan/degmodel_nonhierach.stan')%>%str_replace('.stan','.rds')%>%file.remove
  # origmodelfile<-file.path(root,'src/Stan/degmodel_nonhierach.stan')
  origmodelfile<-file.path(root,'src/Stan/degmodel_hierarch.stan')
  tdir <-   tempdir()
  dir.create(tdir)
  file.exists(origmodelfile)
  file.copy(origmodelfile,tdir)
  modelcopy=file.path(tdir,basename(origmodelfile))
  stopifnot(file.exists(modelcopy))
  modelcopy
  message(paste0('fitting model ',basename(origmodelfile),' on ',length(g2fit),' genes'))
  modelsamplefile<-paste0('stansamples_',basename(tdir),'_',basename(origmodelfile),'_modsamples')

  #
  n_chains <- 4
  realstanfit <- rstan::stan(file=modelcopy,
          model_name=basename(tdir),seed=1,
          # model_name=modelnm,seed=1,
          data=standata,
              control=list(adapt_delta=0.95,max_treedepth=15),save_dso=FALSE,
              pars=pars,
              # init = lapply(seq_len(n_chains),function(id)list('lrTE'=array(rep(20,length(g2fit))))),
              chains=n_chains,iter=1e3,cores=n_chains,verbose=TRUE,save_warmup=FALSE,
              sample_file=modelsamplefile,
              # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
            )

  #now fit linear model
  modelnm=paste0(g2fit[1],'_linear_',length(g2fit))
  file.path(root,'src/Stan/degmodel_nonhierach.stan')%>%str_replace('.stan','.rds')%>%file.remove
 
  linmodelfile <- 'src/Stan/degmodel_simple_linear.stan'
  linmodelsamplefile<-paste0('stansamples/',basename(tdir),'_',basename(linmodelfile),'_modsamples')

  realstanfit_lin <- rstan::stan(file=file.path(root,linmodelfile),seed=1,model_name=modelnm,data=standata,
                             control=list(adapt_delta=0.95,max_treedepth=20),sample_file='tmpstansamples.txt',
                             chains=4,iter=1000,cores=4,verbose=F,
                              pars=pars%>%setdiff(c('deg','ldeg')),
                              linmodelsamplefile,
                             # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
                            )

  list(kinetic=realstanfit,linear=realstanfit_lin,data=standata,genes=g2fit)

}

genes2fit <- exprdata$gene_name%>%c('Satb2','Orc3')%>%unique%>%setNames(.,.)
testgeneset <- c("Acadvl", "Ak1", "Asna1", "Cald1", "Cst3", "Ctnnd2", "Cul2",
"Dclk1", "Dpysl3", "Epb41l1", "Flna", "Hspa12a", "Igsf21", "Nos1",
"Orc3", "Phactr4", "Rap1b", "Rasa3", "Rps6ka5", "Satb2", "Srcin1",
"Stx1b", "Tbc1d24", "Trmt1l", "Zc2hc1a", "Zmym3")
genes2fit<-testgeneset
which(genes2fit == 'Satb2')

# allstanfits <- mclapply(genes2fit,safely(getstanfits),genelengths,exprdata,lengthnorm=FALSE)

registerDoMC(32)

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