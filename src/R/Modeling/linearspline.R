library(deSolve)
library(magrittr)
library(txtplot)
# library(minpack.lm)
library(tidyverse)
library(splines)
# library(SplinesUtils)
library(data.table)
library(nlme)
# library(nlmeODE)
library(BiocParallel)
library(splines2)



################################################################################
########Math  - this time linear space splines
################################################################################
# So I think this is all fine if we just use linear space splines
# IF we can express a time constant function as a simple combination of splines.
# Now, our linear thing is a 


# time = 0:23 + 10
# sig = 0:23 + 10

# dbs%>%args

# fit=lm(sig ~ 0+ bs(time,intercept=TRUE))
# fit
# fit_x2=lm(sig*2 ~ bs(time,intercept=TRUE))
# fit_x2


# pdf('tmp.pdf')
# bs(time,intercept=TRUE)%>%matplot
# dev.off()
# normalizePath('tmp.pdf')

# dbs(time,intercept=TRUE) %*% (fit$coef)
# bs(time,intercept=TRUE) %*% (fit$coef)
# ibs(time,intercept=TRUE) %*% (fit$coef)

# #so if we want our RNA to be 1+4 dimensions for an arbitrary fit
# #then we can have 
# #and we want our prot to be 1+4 dimensions for arbitrary fit
# inter_dbasis = cbind(1,dbs(time,df=4))
# integ_inter_dbasis = cbind(1,1:(1*length(time)),bs(time,df=4))





# #So weith default df = 3
# #so note that with an intercept, bs has only 3
# princomp(bs(time))%>%summary
# princomp(bs(time,intercept=TRUE))%>%summary
# #And dbs only has two
# princomp(dbs(time,intercept=TRUE))%>%summary
# #so to fit e.g. 5 timepoint in our rate of change, we'd need 6 
# mydf = 6
# #we would like there to be a meaningful intercept at both levels actually
# princomp(dbs(time,df=mydf,intercept=TRUE))%>%summary

# let bs be the linear spline basis
# px = ([m0,a1,a2,a3])
# so that bs %*% pv, bs is a txD matrix
# dbs is now also txD
# let zv be teh coefficients of a linear fit in our spline basis

# dP(t) = Ks P(t) - Kd (P(t))



# zv = lm(time ~ 0 + inter_dbasis)$coef
# zv = lm(1:24 ~ 0 + integ_inter_dbasis[,,-5])$coef


# note that px1 is t+1 for a perfect  


# dbasis = inter_dbasis
# mybs = integ_inter_dbasis



# ###actually hadamard producct

# #Preeeety sure this works
# t,D    D,1    =   t,D     D,1 -  t,D     D,1  %.%  t,D1   D1,1
# dbasis %*% px = dbasis %*% yv  -  dbasis %*% zv %.%  mybs %*% px1

# dbasis %*% px = dbasis %*% yv -  dbasis %*% zv %.%  mybs %*% px1
# #add
# dbasis %*% px + dbasis %*% zv %.%  mybs %*% px = dbasis %*% yv
# #flilp
# dbasis %*% yv = dbasis %*% px + dbasis %*% zv %.%  mybs %*% px1 

# #Okay now can I factorize this???
# #note the hadamard (elementwise) product there
# #t,D      D,1 = t,D    %*% D,1 +     t,D   . D,1        t,D1 %*% D1,1
# dbasis %*% yv = dbasis %*% px +     dbasis %*% zv %.% (mybs %*% px1)

# dbasis %*% yv = dbasis %*% px +     dbasis %*% zv %.% (mybs %*% px1)

# #Okay so no, this isn't a linear operation, but it's lcose
# #t,D      D,1 =  t,D  %*% D,1 +  t,D   .    D,1  %*%    
# dbasis %*% yv = dbasis %*% px + diag(mybs %*% px1) %*% dbasis %*% t(zv) 
# #Now let's assume that there's a left inverse to dbasis - i.e. it has linearly independent columns
# dbasis %*% yv = dbasis %*% px + diag(mybs %*% px1) %*% dbasis %*% t(zv) 

# #now we end up with this.......
# yv = px + dbinv %*% diag(mybs %*% px1) %*% dbasis %*% t(zv) 


# #We COULD also just express our zv in terms of the protein level basis
# dbasis %*% yv = dbasis %*% px +  (mybs %*% zv) %.% (mybs %*% px1)

# dbasis %*% yv = dbasis %*% px +  mybs %*%  (zv %.% px1)

# #again assume linear independence in the columns of dbasis
# dbasis %*% yv = dbasis %*% px +  mybs %*% zv  %.% mybs %*% zv



# #what about ks...



# dbasis %*% yv = dbasis %*% yv +  mybs %*%  (zv %.% px1)
# D,1 =D,1 +       D,t   %*% t,D1  %*%   D1,1

# #this is WAY nicer

#                     t,D1    D1,1
# yv = px + invdbasis %*% mybs %*%  (zv %.% px1)

# #waht if we assume that zv[1] = 0?
# yv = px + invdbasis %*% mybs %*%  (zv %.% px)





#now in orthogonolaizing our basis... we need to think about how the coefficients are effected.

#hang on, if one of our basis terms is the intercept in dbs space, then every term but one in




################################################################################
########Numerically test the linear space spline solution
################################################################################
# install.packages('deSolve')  
# library(deSolve)
spl_getprotdf<-function(
splfit,
splinedf =  4,
deg = 0.1,
rTE = 10,
orthns_imat,
ms0ratio = 1,
logfit=TRUE
 ) {
  ntps <- length(predict(splfit))
  dPdt <- function(itime, state, parameters){
   trans = if(logfit) exp else identity
   dPdt <- with(as.list(c(state,parameters)) ,{
     mRNA <- trans(c(1,orthns_imat[itime,]) %*% scoefs)
     mRNA <- max(0,mRNA)
     (rTE * mRNA) - (deg*P)
   })
   list(dPdt)
  }

   trans = if(logfit) log else identity
   
  stopifnot(splinedf %in% 1:100)
  splineribo <- trans(predict(splfit))

  #set up the actual ode
  state<-c('P'= min(
    ((rTE*splineribo[1])/deg),
    # ((rTE*splineribo[1])/0.1)
    Inf
  ))


  state = state*ms0ratio
  parameters = list(rTE=rTE,bounds=c(1,ntps),deg=deg,df=splinedf,scoefs = splfit$coef,orthns_imat=orthns_imat)%>%as.list
  Pdf = ode(y = state, times = 1:ntps, func = dPdt, parms = parameters)
  data.frame(ribo = splineribo,P=Pdf[,2],time=1:length(splineribo))


}

##############

#now we end up with this.
# yv = px + dbinv %*% diag(mybs %*% px1) %*% dbasis %*% t(zv) 


ribo <-c(rep(100,20),400,800,400,rep(100,20))
time<-seq_along(ribo)
ztime = time -1
nz = last(ztime)
# nsbasis = evaluate(OrthogonalSplineBasis(orthogonalsplinebasis::expand.knots(c(0,nz))),ztime)
# dnsbasis = evaluate(deriv(OrthogonalSplineBasis(orthogonalsplinebasis::expand.knots(c(0,nz)))),ztime)

nsbasis <- bs(time,df=5)
dnsbasis <- dbs(time,df=5)
splfit <- (lm(log(ribo) ~ nsbasis))
exp(predict(splfit))%>%txtplot
stopifnot(!any(is.na(splfit$coef)))

num_ode_res<-spl_getprotdf(splfit,deg = 0.1,rTE=20,splinedf=4,orthns_imat=nsbasis)

num_ode_res$ribo%>%txtplot(ylim=c(0,max(.)))
num_ode_res$P%>%txtplot(ylim=c(0,max(.)))


#waht if we assume that zv[1] = 0?
# yv = px + invdbasis %*% mybs %*%  (zv %.% px)

# library(orthogonalsplinebasis)

#This gets us an orthogonal spline matrix
# inter_dbasis=evaluate(OBasis(expand.knots(seq(first(time),last(time),length.out=3),4)),time)
# round(ginv(inter_dbasis) %*% inter_dbasis,4)

#But what if we want one of the columns to be an interept? Then
# intb %*% diag(v) =  

# inter_dbasis = cbind(1,dbs(time,df=6,intercept=FALSE))
# inter_dbasis = princomp(inter_dbasis)
# inter_dbasis = inter_dbasis$scores 


# OBasis(expand.knots(0,))


# integ_inter_dbasis = cbind(1,1:(1*length(time)),bs(time,df=4))

# SplineBasis(expand.knots(seq(first(time),last(time),length.out=3),4))

#how can I arrive at an orthonormal basis that includes an intercept?
#What if I get an orthonormal basis for the differences.
#then I arrive at the orthonormal basis for the actual space
# inter_dbasis=evaluate(OBasis(expand.knots(seq(first(time),last(time),length.out=3),4)),time)
# BiocManager::install('pracma')

# library(pracma)
# A <- matrix(c(0,-4,2, 6,-3,-2, 8,1,-1), 3, 3, byrow=TRUE)
# gs <- gramSchmidt(A)
# (Q <- gs$Q); (R <- gs$R)
# Q %*% R  # = A

# mydbs = cbind(1, bs(time))
# mydbs%>%{ginv(.)%*%.}%>%round(10)
# invmydbs = ginv(mydbs)
# mybs = cbind(1,1:length(time),ibs(time))

#wait, fuck lol c(1,bs) without intercept option is already pretty much orthogonal





################################################################################
########Now test the above
################################################################################
library(MASS) 
library(deSolve) 
# install.packages('txtplot')
library(splines2)
library(txtplot)

flankn=20
ribo <-c(rep(100,flankn),200,400,800,400,200,rep(100,flankn))
time = 1:length(ribo)
# mydbs = cbind(1, bs(time,df=(length(time)-1)))
#What if I use m spline?




################################################################################
########Define bases
################################################################################
  
mybs =  cbind(1,time,iSpline(time,df=(length(time)-1)))
mydbs = cbind(1,mSpline(time,df=(length(time)-1)))
mybs%>%{ginv(.)%*%.}%>%round(10)%>%{.-diag(ncol(.))}%>%abs%>%max
mydbs%>%{ginv(.)%*%.}%>%round(10)%>%{.-diag(ncol(.))}%>%abs%>%max

invmydbs = ginv(mydbs)
mydbs_lims  = cbind(,mydbs) %>% {.[c(1,nrow(.)),]}
myddbs_lims  = cbind(0,mSpline(time,length(time)-1,derivs=1))%>%{.[c(1,nrow(.)),]}
# mydbs_ddlims = cbind(mSpline(time,length(time)-1,derivs=2))%>%{.[c(1,nrow(.)),]}

zv = lm(rep(1,length(ribo))~0+mybs)%>%.$coef%>%round(10)%>%replace_na(0)

dim(mybs);dim(mydbs);length(zv)



################################################################################
########Test data
################################################################################
  


dim(mydbs)
dim(mybs)
ribo <- ribo%>%rpois(length(.),.)
splfit <- (lm(ribo ~ 0+mydbs))
splfit
predict(splfit)%>%txtplot
stopifnot(!any(is.na(splfit$coef)))
mydeg = 0.04
myrTE  = 20
num_ode_res<-spl_getprotdf(splfit,deg = mydeg,rTE=myrTE,orthns_imat=mydbs[,-1],logfit=FALSE)
#looks right
num_ode_res$ribo%>%txtplot
num_ode_res$P%>%txtplot

###Get the mspline fit for our data that also achieves linearity 
basis_with_divs = rbind(
  mybs,
  mydbs_lims,
  myddbs_lims
)

#Fit a
mypx_lm = lm(c(num_ode_res$P,0,0,0,0) ~ 0 + basis_with_divs,weights = c(rep(1,length(num_ode_res$P)),rep(1e3,4)))

basis_with_divs = rbind(
  mybs,
  myddbs_lims
)

mypx_lm = lm(c(num_ode_res$P,0,0) ~ 0 + basis_with_divs,weights = c(rep(1,length(num_ode_res$P)),rep(1e3,2)))
mypx = mypx_lm$coef%>%replace_na(0)

predict(mypx_lm)%>%head(-4)%>%txtplot

txtplot(num_ode_res$P,predict(mypx_lm)%>%head(-4))
txtplot(log10(num_ode_res$ribo),log10(est))




#zv is just the c(deg,0,0,....)


#now use yv = px + invdbasis %*% mybs %*%  (zv %.% px)
length(ribo)
length(num_ode_res$P)
dim(invmydbs)
dim(mydbs)
dim(mybs)
length(mypx)

degzv = zv*mydeg

mybs
length(mypx)

dim(matrix(mypx[-1]))
dim(invmydbs)
dim(invmydbs)
dim(mybs)
dim((mybs %*% degzv) )
length(degzv)
dim((mybs %*% (matrix(degzv)) ))
dim(((mybs %*% mypx)  * (mybs %*% degzv)) )

yv = matrix(mypx[-1]) + invmydbs %*% ((mybs %*% mypx)  * (mybs %*% degzv)) 


#You would think thi sis the same - nope
# yv = tail(matrix(mypx),-1) + invmydbs %*% ((mybs %*% mypx)  * (mybs %*% zv)) 
#how about also no!
# yv = tail(matrix(mypx),-1) + ((invmydbs %*%mybs %*% mypx)  * (invmydbs %*% mybs %*% zv)) 
est = mydbs %*% yv
est/(ribo)
mean(est/(ribo))
Ks

est%>%txtplot

# dbasis %*% yv = dbasis %*% px +  mybs %*% zv  %.% mybs %*% zv
#estimate using spline of P



#so, armed with the above, assuming those 
#small innaccuracies aren't imporant (the mean was very accurate after all)
#we can then write a nice linear function to optimize
# 0
parvect = par$par
# invmydbs=invmydbs;zv=zv;mybs=mybs;ribosd=0.1;protsd=0.3;ribo=ribo;prot=prot

# myyv = pxopt[-1] + invmydbs %*% ((mybs %*% pxopt)  * (mybs %*% degzv)) 
{

linearsplineoptfunc <- function(parvect,invmydbs,mybs,zv,ribo,prot,ribosd,protsd,returnprot=FALSE,returnribo=FALSE){
    oparvect = parvect
    browser()
    parvect[1:2] = exp(parvect[1:2])
    stopifnot(!any(is.na(parvect)))
    stopifnot(!any(!is.finite(parvect)))
  with(as.list(parvect),{
    Ks = parvect[1]
    Kd = parvect[2]
    pxopt = parvect[-c(1:2)]
    degzv = zv*Kd
    
    yv = pxopt + invmydbs %*% ((mybs %*% pxopt)  * (mybs %*% degzv)) 
    est_ribo = (mydbs %*% yv) / Ks
    est_protein = mybs %*% pxopt

    if(returnprot) return(est_protein)
    if(returnribo) return(est_ribo)

    # browser()

    ribo_LL = dnorm(log(ribo),log(est_ribo),ribosd,log=TRUE)
    prot_LL = dnorm(log(prot),log(est_protein),protsd,log=TRUE)
    # browser()
    alphapriorLL = dnorm(pxopt[-1],0,sd=1,log=TRUE)
    alphapriorLL = dnorm(pxopt[-1],0,sd=10000,log=TRUE)
    alphapriorLL = 0
    # alphapriorLL = dnorm(px[-1],0,sd=20,log=TRUE)

    nsconstraintLL = dnorm(mydbs_dlims %*% pxopt[-1],0,sd=0.1,log=TRUE)+dnorm(mydbs_ddlims %*% pxopt[-1],0,sd=0.1,log=TRUE)
    nsconstraintLL = 0

    LL = sum(c(ribo_LL ,prot_LL,alphapriorLL,nsconstraintLL))
    # browser()
    if(!is.finite(LL))browser()
    # cat(LL);cat('...')
    LL
  })
}

parvect = c(Ks=log(10),Kd=log(0.5),ms0=log(5000),mypxrest=rep(0,ncol(mybs)-1))
parvect[-c(1:2)]=mypx$coef%>%replace_na(.,mean(na.omit(.)))
prot = num_ode_res$P
# parvect%<>%

lowervect = rep(-2e4,length(parvect))%>%setNames(names(parvect))%>%{.['Ks']=log(0.1);.['mypxrest1']=log(1);.}
stopifnot(all(is.finite(lowervect)))
uppervect = rep(2e4,length(parvect))%>%setNames(names(parvect))%>%{.['Kd']=log(1);.}
stopifnot(all(is.finite(uppervect)))
stopifnot(!any(is.na(parvect)))
stopifnot(!any(!is.finite(parvect)))


stopifnot(
  linearsplineoptfunc(parvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot) < 0
)
}
par = optim(parvect,
  fn=linearsplineoptfunc,
  method='L-BFGS-B',
  lower = lowervect,
  upper = uppervect,
  # hessian = TRUE,
  invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  # control=c(pgtol=0.00000001)
  )

par$par == uppervect
par$par == uppervect

linearsplineoptfunc(par$par,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,returnprot=TRUE)%>%txtplot
linearsplineoptfunc(par$par,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,returnribo=TRUE)%>%txtplot

linearsplineoptfunc(parvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,returnprot=TRUE)%>%txtplot
linearsplineoptfunc(parvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,returnribo=TRUE)%>%txtplot

prot%>%txtplot
ribo%>%txtplot
}


# install.packages('brm')
# library(brm)


# brm(ribo ~ ,data=data.frame(ribo=ribo))



# 0
# debugger()


# #FUCK
# linear splines end up negative. No way to 


# ribo = c(rep(100,10),200,400,800,400,200,rep(100,10))
# ##Okay so if we can use M splines for the above, what about goiing back to log level





# lR = lP + log(dlP + Kd) - log(Ks)

# now let F = dLP + Kd be an M spline

# lR = lP + log(F) - log(Ks)

# lR = int(M%*%(px-Z)) + log(F) - log(Ks)

# lR = I%*%(px-Z) + log(M %*% (px))) - log(Ks)

# lR = I%*%(px-Z) + log(M %*% (px))) - log(Ks)

# lR = I%*%(px-Z) + log(Ks %*% M %*% (px))) - log(Ks)


# yv = pxopt[-1] + invmydbs %*% ((mybs %*% pxopt)  * (mybs %*% degzv)) 

# yv = pxopt[-1] + invmydbs %*% ((mybs %*% pxopt)  * (mybs %*% degzv))

# library(brms)
# brm(bf(ribo ~ mydbs),nl=TRUE,data=data.frame(ribo=ribo))







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



