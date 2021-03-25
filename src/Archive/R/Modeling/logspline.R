library(here)
# install.packages('here')


spl_getprotdf<-function(
splfit,
splinedf =  4,
deg = 0.1,
rTE = 10,
orthns_imat,
ms0ratio = 1,
logfit=TRUE
 ) {
  ntps <- length(predict(splfit,orthns_imat))

   trans = if(logfit) exp else identity

  dPdt <- function(itime, state, parameters){
   trans = if(logfit) exp else identity
   dPdt <- with(as.list(c(state,parameters)) ,{
     mRNA <- trans(c(orthns_imat[itime,]) %*% scoefs)
     # mRNA <- max(0,mRNA)
     (rTE * mRNA) - (deg*P)
   })
   list(dPdt)
  }

   
  stopifnot(splinedf %in% 1:100)
  splineribo <- trans(predict(splfit,orthns_imat))


  #set up the actual ode
  state<-c('P'= min(
    ((rTE*splineribo[1])/deg),
    # ((rTE*splineribo[1])/0.1)
    Inf
  ))

  state = state*ms0ratio
  parameters = list(rTE=rTE,bounds=c(1,ntps),deg=deg,df=splinedf,scoefs = coef(splfit),orthns_imat=orthns_imat)%>%as.list
  Pdf = ode(y = state, times = 1:ntps, func = dPdt, parms = parameters)
  data.frame(ribo = splineribo,P=Pdf[,2],time=1:length(splineribo))


}



################################################################################
########Now test the above
################################################################################
library(MASS) 
library(deSolve) 
# install.packages('txtplot')
library(splines2)
library(magrittr)
library(txtplot)
library(nnls)
# BiocManager('nnls')



flankn=1
ribo <-c(rep(100,flankn),200,400,800,400,200,rep(100,flankn))
time = 1:length(ribo)
# mydbs = cbind(1, bs(time,df=(length(time)-1)))
#What if I use m spline?




################################################################################
########Define bases
################################################################################
  
mybs =  cbind(1,time-1,iSpline(time,df=(length(time)-1)))
mydbs = cbind(1,iSpline(time,df=(length(time)-1),deriv=1))
dim(mybs);dim(mydbs)
mybs%>%{ginv(.)%*%.}%>%round(10)%>%{.-diag(ncol(.))}%>%abs%>%max
mydbs%>%{ginv(.)%*%.}%>%round(10)%>%{.-diag(ncol(.))}%>%abs%>%max
stopifnot(!any(mydbs<0))
stopifnot((!all(mydbs[1,]==0)))

invmydbs = ginv(mydbs)
mydbs_lims  = cbind(mydbs) %>% {.[c(1,nrow(.)),]}
myddbs_lims  = cbind(0,0,iSpline(time,length(time)-1,derivs=2))%>%{.[c(1,nrow(.)),]}
mydddbs_lims  = cbind(0,0,iSpline(time,length(time)-1,derivs=3))%>%{.[c(1,nrow(.)),]}
# mydbs_ddlims = cbind(mSpline(time,length(time)-1,derivs=2))%>%{.[c(1,nrow(.)),]}
dim(myddbs_lims)

iSpline(1:5,intercept=TRUE)
iSpline(1:5,intercept=TRUE,deriv=1)

# dim(mybs);dim(mydbs);length(zv)


# basis_ribolinlim = rbind(
#   mybs,
#   myddbs_lims,
#   mydddbs_lims
# )

# dbasis_ribolinlim = rbind(
#   mydbs,
#   myddbs_lims,
#   mydddbs_lims
# )




################################################################################
########Linear?
################################################################################
library(splines)
library(splines2)
library(splines)
library(tidyverse)
library(magrittr)
timeknots <- time[c(-1,-length(time))]
mybs <- cbind(1,ibs(time, knots = timeknots,degree = 1, intercept = TRUE))
mydbs = bs(time, knots = timeknots,degree = 1, intercept = TRUE)


dim(mybs);dim(mydbs)
mybs%>%{ginv(.)%*%.}%>%round(10)%>%{.-diag(ncol(.))}%>%abs%>%max
mydbs%>%{ginv(.)%*%.}%>%round(10)%>%{.-diag(ncol(.))}%>%abs%>%max
stopifnot(!any(mydbs<0))
stopifnot((!all(mydbs[1,]==0)))

invmydbs = ginv(mydbs)
mydbs_lims  = cbind(mydbs) %>% {.[c(1,nrow(.)),]}
#this is all unncessary with piecewise linear but meh
myddbs_lims  = cbind(0,0,iSpline(time,length(time)-1,derivs=2)*0)%>%{.[c(1,nrow(.)),]}
mydddbs_lims  = cbind(0,0,iSpline(time,length(time)-1,derivs=3)*0)%>%{.[c(1,nrow(.)),]}
# mydbs_ddlims = cbind(mSpline(time,length(time)-1,derivs=2))%>%{.[c(1,nrow(.)),]}
dim(myddbs_lims)
library(tidyverse)
dzv=lm( rep(1,nrow(mydbs)) ~ 0+mydbs)%>%.$coef%>%round(10)%>%replace_na(0)



dim(mybs);
dim(mydbs);
dim(mybs);
dim(myddbs_lims);
length(mydddbs_lims)


basis_ribolinlim = rbind(
  mybs,
  myddbs_lims,
  mydddbs_lims
)

dbasis_ribolinlim = rbind(
  mydbs,
  myddbs_lims[,-1],
  mydddbs_lims[,-1]
)



###

dim(mydbs)
dim(mybs)
ribo <- ribo%>%rpois(length(.),.)
splfit <- (nnls(mydbs,log(ribo)))
# splfit <- lm(ribo ~ 0+mydbs)
txtplot(ribo)
txtplot(mydbs %*% coef(splfit))
# splfit <- (nnls(dbasis_ribolinlim,c(log(ribo),0,0)))
# txtplot(dbasis_ribolinlim %*% coef(splfit))

stopifnot(!any(is.na(splfit$coef)))

#looks right
# num_ode_res$ribo%>%txtplot
# num_ode_res$P%>%txtplot

###Get the mspline fit for our data that also achieves linearity 

#So can we reverse this, given parameeters/
#Why can't we go from

basis_with_divs = rbind(
  mybs,
  cbind(myddbs_lims),
  cbind(mydddbs_lims)
)





################################################################################
########stan - oriignally tested on R code below
################################################################################
mydeg = 0.8
myrTE  = 20
splfit = lm(log(ribo) ~ 0+ mydbs)
num_ode_res<-spl_getprotdf(splfit,deg = mydeg,rTE=myrTE,orthns_imat=mydbs,logfit=TRUE)
c

prot = num_ode_res$P

n_genes<-20
num_ode_res
ntps<-length(prot)
lMS <- matrix(log(rep(prot,n_genes)),nrow=ntps)
lMS_tau <- matrix(rep(0.1,ntps*n_genes),nrow=ntps)
lribo <- matrix(log(rep(ribo,n_genes)),nrow=ntps)
lribo_tau <- matrix(rep(0.1,ntps*n_genes),nrow=ntps)
# ntps <- nrow(lribo)
mybs
mydbs
n_chains=4
#arrange
standata<-list(
  G=n_genes,          #  genes
  T=ntps,          #  timepoints
  lMS=lMS,  # mass spec data mean
  lMS_tau=lMS_tau,  # mass spec data precision
  lribo=lribo, # lriboseq (synthesis) data mean
  lribo_tau=lribo_tau, # lriboseq (synthesis) data sd
  mybs=mybs,
  mydbs=mydbs
)
standata%>%map(dim)


stanmodel <- here('src/Stan/degmodel_dataconfint.stan')%T>%{stopifnot(file.exists(.))}
stanmodel <- here('src/Stan/degmodel_dataconfint_hierarch.stan')%T>%{stopifnot(file.exists(.))}

teststandir<-tempdir()
teststandir%>%dir.create

modelsamplefile <- here('pipeline/stan/degmodel_dataconfint')

modelsamplefile <- here('pipeline/stan/degmodel_dataconfint_hierarch')
library(rstan)

stanmodelob<-stan_model(file=stanmodel,model_name='foom')

opt <- optimizing(stanmodelob,data=standata)

opt$par%>%enframe('par','value')%>%filter(par%>%str_detect('hmu|hsig'))

stop()
stanfit <-   stan(file=stanmodel,
        model_name='food',seed=1,
        # model_name=modelnm,seed=1,
        data=standata,
            control=list(adapt_delta=0.95,max_treedepth=15),save_dso=FALSE,
            # pars=c(pars,'lrTE'),
            # init = lapply(seq_len(n_chains),function(id) list('lrTE'=array(rep(20,length(g2fit))))),
            chains=n_chains,iter=600,cores=n_chains,verbose=TRUE,save_warmup=FALSE,
            sample_file=modelsamplefile,
            # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
          )


stanfitop <-   rstan::optimizing(file=stanmodel,
        model_name='foodb',seed=1,
        # model_name=modelnm,seed=1,
        data=standata,
            control=list(adapt_delta=0.95,max_treedepth=15),save_dso=FALSE,
            # pars=c(pars,'lrTE'),
            # init = lapply(seq_len(n_chains),function(id)list('lrTE'=array(rep(20,length(g2fit))))),
            chains=n_chains,iter=600,cores=n_chains,verbose=TRUE,save_warmup=FALSE,
            sample_file=modelsamplefile,
            # init=function(z) list(rTE=array(c(10),dim=c(n_genes)),MS0=array(ribo_mat[1,]*rTEs,dim=c(n_genes))),
          )
#so, multiple genes in parallel works!
#how about a hierarchical model with three instances of our
#test data? Yes this makes it run like a god damned pig.

# runmodel()
#summary(stanfit)

#so the hierarchical parameters don't really converge
#This maybe the prior. Switch to gamma priors on them, rather than exponentiating them?

summary(stanfit)[[1]]%>%as.data.frame%>%rownames_to_column('par')%>%filter(par%>%str_detect('hmu|hsig'))

stop()

###############################################################################
#######Now the oopot
###############################################################################

# library(brms)

# dP(t) = Ks R(t) - Kd (P(t))
# R(t) = dP(t) + Kd (P(t)) / Ks
# log(R(t)) = log(dP(t)+Kd P(t)) - log(Ks)
# log(R(t)) = log(dP(t)+Kd P(t)) - log(Ks)
# log(R(t)) = log(P(t) * dlP(t) +Kd P(t)) - log(Ks)

# lR(t) = log(P(t)) + log(dlP(t) + Kd(t)) - log(Ks)

# lR(t) = lP(t) + log(dlP + Kd) - log(Ks)

# dlP + Kd or the fold synthesis at any given time, must be positive.

# We call it M

# dlP+Kd = dbs %*% cM = dbs %*% (zv+cv)

# => cM - zv = cv

# mybs %*% [a,cv] + log(mydbs %*% M) - log(Ks)
# BiocManager::install(c('tryCatchLog','brms'))


{

loglinearsplineoptfunc <- function(parvect,invmydbs,mybs,zv,ribo,
  prot,ribosd,protsd,returnprot=FALSE,returnribo=FALSE){
    cat('.')
    oparvect = parvect   
    #ddon't exponentiate th intercept for the protein
    parvect[-3] <- exp(parvect[-3])
    if(any(!is.finite(parvect))) browser()
    stopifnot(!any(is.na(parvect)))
    stopifnot(!any(!is.finite(parvect)))

  with(as.list(parvect),{
    Ks = parvect[1]
    Kd = parvect[2]
    m0 = parvect[3]
    cM = parvect[-c(1:3)]

    degzv = dzv*Kd
    #so cv are the actual possibly negative, fold change
    #cM aree the positivee fold changes after degreedatio is removed
    cv = cM - degzv

    # lR(t) = log(P(t)) + log(dlP(t) + Kd(t)) - log(Ks)

    est_log_ribo = (mybs %*% c(m0,cv)) + log(mydbs %*% cM) - log(Ks)
    est_log_ribo
    (mybs %*% (c(m0,cM) - c(0,degzv))) + log(mydbs %*% cM) - log(Ks)
    est_log_protein = mybs %*% c(m0,cv)


    if(returnribo) return(exp(est_log_ribo))
    if(returnprot) return(exp(est_log_protein))

    # browser()

    ribo_LL = dnorm(log(ribo),est_log_ribo,ribosd,log=TRUE)
    prot_LL = dnorm(log(prot),est_log_protein,protsd,log=TRUE)
    alphapriorLL = dnorm(cM[-1],0,sd=4,log=TRUE)
    # alphapriorLL = dnorm(cM[-1],0,sd=10000,log=TRUE)
    alphapriorLL = 0
    # alphapriorLL = dnorm(px[-1],0,sd=20,log=TRUE)
    # nsconstraintLL = 
    #   dnorm(myddbs_lims[,-1] %*% cv,0,sd=0.1,log=TRUE)+
    #   dnorm(mydddbs_lims[,-1] %*% cv,0,sd=0.1,log=TRUE)
    nsconstraintLL = 0

    LL = sum(c(ribo_LL ,prot_LL,alphapriorLL,nsconstraintLL))
    # browser()
    if(!is.finite(LL)) browser()
    if(length(LL)==0) browser()
    # cat(LL);cat('...')

    - LL
  })
}
}
mydeg = 0.8
myrTE  = 20
splfit = lm(log(ribo) ~ 0+ mydbs)
num_ode_res<-spl_getprotdf(splfit,deg = mydeg,rTE=myrTE,orthns_imat=mydbs,logfit=TRUE)
c

prot = num_ode_res$P
parvect = c(Ks=log(10),Kd=log(0.5),ms0=log(prot[1]),mypxrest=rep(1,ncol(mybs)-1))
# parvect[-c(1:2)]=1

# parvect%<>%

lowervect = rep(-100,length(parvect))%>%setNames(names(parvect))%>%
  {.['Kd']=log(0.001);.['ms0']=1;.}
stopifnot(all(is.finite(lowervect)))
uppervect = rep(100,length(parvect))%>%setNames(names(parvect))%>%
{.['Kd']=log(1);.['ms0']=1e12;.}
stopifnot(all(is.finite(uppervect)))

stopifnot(!any(is.na(parvect)))
stopifnot(!any(!is.finite(parvect)))

stopifnot(
  loglinearsplineoptfunc(parvect,invmydbs=invmydbs,zv=dzv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot) >0
)


par = optim(parvect,
  fn=loglinearsplineoptfunc,
  method='L-BFGS-B',
  lower = lowervect,
  upper = uppervect,
  # hessian = TRUE,
  invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  # control=c(pgtol=0.00000001)
  )

txtplot(num_ode_res$ribo)
txtplot(num_ode_res$P)
fitparvect=par$par
exp(fitparvect[-3])



fitll=loglinearsplineoptfunc(fitparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot)
initll=loglinearsplineoptfunc(parvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot)
fitll;initll
fitribo=loglinearsplineoptfunc(fitparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  returnribo=TRUE)%T>%txtplot
fitprot=loglinearsplineoptfunc(fitparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  returnprot=TRUE)%T>%txtplot

txtplot(ribo,fitribo)
txtplot(prot,fitprot)

myrTE;exp(fitparvect['Ks'])
mydeg;exp(fitparvect['Kd'])


#what if I fix these values?

forcepar = optim(parvect,
  fn=loglinearsplineoptfunc,
  method='L-BFGS-B',
  lower = lowervect%>%{.['Ks']=log(myrTE-1);.['Kd']=log(mydeg*0.99);.},
  upper = uppervect%>%{.['Ks']=log(myrTE+1);.['Kd']=log(mydeg*1.01);.},
  # hessian = TRUE,
  invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  # control=c(pgtol=0.00000001)
)
forceparvect=forcepar$par
forcell=loglinearsplineoptfunc(forceparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot)
initll=loglinearsplineoptfunc(parvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot)
fitll
forcell;initll
forceribo=loglinearsplineoptfunc(forceparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  returnribo=TRUE)%T>%txtplot
forceprot=loglinearsplineoptfunc(forceparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  returnprot=TRUE)%T>%txtplot
txtplot(ribo,forceribo)
txtplot(prot,forceprot)


# if I reinput these values od I recover them

reinputpar = optim(parvect,
  fn=loglinearsplineoptfunc,
  method='L-BFGS-B',
  lower = lowervect,
  upper = uppervect,
  # hessian = TRUE,
  invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,
  ribo=fitribo,prot=fitprot,
  # control=c(pgtol=0.00000001)
)
reinputparparvect<-reinputpar$par
reinputparribo=loglinearsplineoptfunc(reinputparparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  returnribo=TRUE)%T>%txtplot
reinputparprot=loglinearsplineoptfunc(reinputparparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  returnprot=TRUE)%T>%txtplot
txtplot(reinputparprot,fitprot)
txtplot(reinputparprot,prot)
txtplot(reinputparribo,fitribo)

#Yes!






################################################################################
########MLE - for confiddence intervals
################################################################################
  
nargs <-length(parvect)

mle_loglinearsplineoptfunc<-function(Ks, Kd, ms0, mypxrest1, mypxrest2, mypxrest3, mypxrest4,
  mypxrest5, mypxrest6, mypxrest7,invmydbs,mybs,zv,ribo,prot,ribosd,protsd,returnprot=FALSE,returnribo=FALSE){
  loglinearsplineoptfunc(
      c(Ks, Kd, ms0, mypxrest1, mypxrest2, mypxrest3, mypxrest4,
      mypxrest5, mypxrest6, mypxrest7),
    invmydbs,mybs,zv,ribo,prot,ribosd,protsd,returnprot=FALSE,returnribo=FALSE)
}
mle_loglinearsplineoptfunc<-function(Ks, Kd, ms0, mypxrest1, mypxrest2, mypxrest3, mypxrest4,
  mypxrest5, mypxrest6, mypxrest7){
  loglinearsplineoptfunc(
      c(Ks, Kd, ms0, mypxrest1, mypxrest2, mypxrest3, mypxrest4,
      mypxrest5, mypxrest6, mypxrest7),
    invmydbs,mybs,zv,ribo,prot,ribosd,protsd,returnprot=FALSE,returnribo=FALSE)
}


ribosd=0.1
protsd=0.5
library(bbmle)

out<-mle2(mle_loglinearsplineoptfunc,
  start=parvect%>%as.list,
  method='L-BFGS-B',
  # fixed=list(invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot),
  lower=lowervect,
  upper=uppervect,
  control=list(maxit=20000)
)

fitparvect=coef(out)

fitribo=loglinearsplineoptfunc(fitparvect,
  invmydbs=invmydbs,zv=zv,mybs=mybs,
  ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  returnribo=TRUE)%T>%txtplot

fitprot=loglinearsplineoptfunc(fitparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  returnprot=TRUE)%T>%txtplot

confints<-confint(out,names(parvect)[1])

exp(confints)

exp(coef(out))

exp(-17)
exp(-1.9)
warnings()
# loglinearsplineoptfunc(parvect%>%{.['Kd']=log(0.01);.},invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,returnribo=TRUE)%>%txtplot

# # linearsplineoptfunc <- function(parvect,invmydbs,mybs,zv,ribo,prot,ribosd,protsd){

#   with(as.list(parvect),{
#     Ks = parvect[1]
#     Kd = parvect[2]
#     pxopt = parvect[-c(1:2)]
#     degzv = zv*Kd

#     yv = pxopt[-1] + invmydbs %*% ((mybs %*% pxopt)  * (mybs %*% degzv)) 
#     est_synthesis = mydbs %*% yv
#     est_protein = mybs %*% pxopt

#     ribo_LL = dnorm(log(ribo),log(est_synthesis/Ks),ribosd,log=TRUE)
#     prot_LL = dnorm(log(prot),log(est_protein),protsd,log=TRUE)
#     LL = sum(c(ribo_LL ,prot_LL))
#     if(!is.finite(LL))browser()
#     LL
#   })
# }

