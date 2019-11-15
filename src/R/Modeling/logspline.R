

################################################################################
########Now test the above
################################################################################
library(MASS) 
library(deSolve) 
# install.packages('txtplot')
library(splines2)
library(txtplot)

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

dim(mybs);dim(mydbs);length(zv)


basis_ribolinlim = rbind(
  mybs,
  myddbs_lims,
  mydddbs_lims
)

dbasis_ribolinlim = rbind(
  mydbs,
  myddbs_lims,
  mydddbs_lims
)

################################################################################
########Now this fitting but in log space?
################################################################################
  
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
     mRNA <- max(0,mRNA)
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

predict.nnls <- function(nnls,A){
  A %*% coef(nnls)
}


dim(mydbs)
dim(mybs)
ribo <- ribo%>%rpois(length(.),.)
splfit <- (nnls(mydbs,log(ribo)))
txtplot(mydbs %*% coef(splfit))
splfit <- (nnls(dbasis_ribolinlim,c(log(ribo),0,0)))
txtplot(dbasis_ribolinlim %*% coef(splfit))

stopifnot(!any(is.na(splfit$coef)))
mydeg = 0.6
myrTE  = 20
num_ode_res<-spl_getprotdf(splfit,deg = mydeg,rTE=myrTE,orthns_imat=mydbs,logfit=TRUE)

#looks right
num_ode_res$ribo%>%txtplot
num_ode_res$P%>%txtplot

dim(mybs)
dim(cbind(myddbs_lims))
dim(cbind(mydddbs_lims))

###Get the mspline fit for our data that also achieves linearity 
basis_with_divs = rbind(
  mybs,
  cbind(myddbs_lims),
  cbind(mydddbs_lims)
)

coef(nnls(mydbs,log(ribo)))




################################################################################
########Now the oopot
################################################################################

# library(brms)
# brm(bf(ribo ~ mydbs),nl=TRUE,data=data.frame(ribo=ribo))

dP(t) = Ks R(t) - Kd (P(t))
R(t) = dP(t) + Kd (P(t)) / Ks
log(R(t)) = log(dP(t)+Kd P(t)) - log(Ks)
log(R(t)) = log(dP(t)+Kd P(t)) - log(Ks)
log(R(t)) = log(P(t) * dlP(t) +Kd P(t)) - log(Ks)

lR(t) = log(P(t)) + log(dlP(t) + Kd(t)) - log(Ks)

lR(t) = lP(t) + log(dlP + Kd) - log(Ks)

dlP + Kd or the fold synthesis at any given time, must be positive.
We call it M

dlP+Kd = dbs %*% cM = dbs %*% (zv+cv)

=> cM - zv = cv

mybs %*% [a,cv] + log(mydbs %*% M) - log(Ks)

dzv=lm( rep(1,nrow(mydbs)) ~ 0+mydbs)%>%.$coef%>%round(10)%>%replace_na(0)

{
loglinearsplineoptfunc <- function(parvect,invmydbs,mybs,zv,ribo,prot,ribosd,protsd,returnprot=FALSE,returnribo=FALSE){
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
    cv = cM - degzv




    est_log_ribo = (mybs %*% c(m0,cv)) + log(mydbs %*% cM) - log(Ks)
    est_log_protein = mybs %*% c(m0,cv)

    if(returnribo) return(exp(est_log_ribo))
    if(returnprot) return(exp(est_log_protein))

    # browser()

    ribo_LL = dnorm(log(ribo),est_log_ribo,ribosd,log=TRUE)
    prot_LL = dnorm(log(prot),est_log_protein,protsd,log=TRUE)
    alphapriorLL = dnorm(cv[-1],0,sd=4,log=TRUE)
    # alphapriorLL = dnorm(cM[-1],0,sd=10000,log=TRUE)
    # alphapriorLL = 0
    # alphapriorLL = dnorm(px[-1],0,sd=20,log=TRUE)
    nsconstraintLL = 
      dnorm(myddbs_lims[,-1] %*% cv,0,sd=0.1,log=TRUE)+
      dnorm(mydddbs_lims[,-1] %*% cv,0,sd=0.1,log=TRUE)
    # nsconstraintLL = 0

    LL = sum(c(ribo_LL ,prot_LL,alphapriorLL,nsconstraintLL))
    # browser()
    if(!is.finite(LL)) browser()
    if(length(LL)==0) browser()
    # cat(LL);cat('...')

    - LL
  })
}

parvect = c(Ks=log(10),Kd=log(0.5),ms0=log(5000),mypxrest=rep(1,ncol(mybs)-1))
# parvect[-c(1:2)]=1
prot = num_ode_res$P
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
}
txtplot(ribo)
txtplot(prot)
parvect
fitparvect=par$par
fitparvect
exp(fitparvect)

fitll=loglinearsplineoptfunc(fitparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot)
initll=loglinearsplineoptfunc(parvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot)
fitll;initll
fitribo=loglinearsplineoptfunc(fitparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  returnribo=TRUE)%T>%txtplot
fitprot=loglinearsplineoptfunc(fitparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  returnprot=TRUE)%T>%txtplot

txtplot(ribo,fitribo)
txtplot(prot,fitprot)

exp(fitparvect['Ks'])
exp(fitparvect['Kd'])
myrTE
mydeg

#what if I fix these values?

par = optim(parvect,
  fn=loglinearsplineoptfunc,
  method='L-BFGS-B',
  lower = lowervect%>%{.['Ks']=log(19);.['Kd']=log(0.59);.},
  upper = uppervect%>%{.['Ks']=log(21);.['Kd']=log(0.61);.},
  # hessian = TRUE,
  invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  # control=c(pgtol=0.00000001)
)

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


out<-mle2(mle_loglinearsplineoptfunc,
  start=parvect%>%as.list,
  method='L-BFGS-B',
  # fixed=list(invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot),
  lower=lowervect,
  upper=uppervect,
  control=list(maxit=20000)
)
fitparvect=coef(out)
fitribo=loglinearsplineoptfunc(fitparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  returnribo=TRUE)%T>%txtplot

fitprot=loglinearsplineoptfunc(fitparvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  returnprot=TRUE)%T>%txtplot

confint(out,names(parvect)[1])
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



lR = lP + log(dlP + Kd) - log(Ks)
