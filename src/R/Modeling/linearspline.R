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
library(conflicted)

conflict_prefer("last", "dplyr")

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
########Test data
################################################################################




nnls<-nnls::nnls
predict.nnls <- function(nnls,A){
  A %*% coef(nnls)
}
library(lsei)
dim(mydbs)
dim(mybs)
ribo <- ribo%>%rpois(length(.),.)
splfit <- (nnls(mydbs,ribo))
txtplot(mydbs %*% coef(splfit))


# splfit <- (nnls(dbasis_ribolinlim,c(ribo,0,0,0,0)))
# txtplot(dbasis_ribolinlim %*% coef(splfit))

stopifnot(!any(is.na(splfit$coef)))
mydeg = 0.6
myrTE  = 20
num_ode_res<-spl_getprotdf(splfit,deg = mydeg,rTE=myrTE,orthns_imat=mydbs,logfit=FALSE)
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

#Fit a
# mypx_lm = lm( ~ 0 + basis_with_divs,weights = c(rep(1,length(num_ode_res$P)),rep(1e3,4)))
# predict(mypx_lm)%>%head(-4)%>%txtplot

# mypx_nnls<-nnls(basis_with_divs,c(num_ode_res$P,0,0,0,0))
# mypx = coef(mypx_nnls)

# mypx_nnls<-nnls(mybs,c(num_ode_res$P))
# mypx = coef(mypx_nnls)

stop()

mypx_lm = lm(c(num_ode_res$P) ~ 0 + mybs,weights = c(rep(1,length(num_ode_res$P))))

mypx_lm = lm(c(num_ode_res$P)-min(num_ode_res$P) ~ 0 + mybs[,-1],weights = c(rep(1,length(num_ode_res$P))))


limbasis = basis_ribolinlim
prot = num_ode_res$P
extnum = nrow(limbasis) - length(prot)
mypx_lm = lm(c(prot,rep(0,extnum)) ~ 0 + limbasis,weights = c(rep(1,length(num_ode_res$P)),rep(1e3,extnum)))
mypx = mypx_lm$coef%>%replace_na(0)

txtplot((limbasis%*%mypx)%>%head(-extnum))
txtplot((prot))


limbasis = rbind(
  mydbs,
  myddbs_lims%>%t%>%head(-1)%>%t,
  mydddbs_lims%>%t%>%head(-1)%>%t
)
extnum = nrow(limbasis) - length(prot)
myrx_lm = lm(c(ribo,rep(0,extnum)) ~ 0 + limbasis,weights = c(rep(1,length(num_ode_res$P)),rep(1e3,extnum)))
myyv = myrx_lm$coef%>%replace_na(0)
txtplot((limbasis%*%myyv)%>%head(-extnum))
txtplot((prot))


#now use yv = px + invdbasis %*% mybs %*%  (zv %.% px)
length(ribo)
length(num_ode_res$P)
dim(invmydbs)
dim(mydbs)
dim(mybs)
length(mypx)

zv = lm( rep(1,nrow(mybs)) ~ 0+mybs)$coef%>%replace_na(0)
# zv = c(mydeg,rep(0,length(mypx)-1))
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

txtplot((mybs %*% mypx))

#yv given px
yv = matrix(mypx[-1]) + invmydbs %*% ((mybs %*% mypx)  * (mybs %*% degzv)) 
est = mydbs %*% yv
est/(ribo)
mean(est/(ribo))
#so what's px given yv AND m0?
# mydbs (yv) = mydbs mypx[-1] + mybs mypx * mybs (degzv)
#turn that haamardd prodduce in to a maaatrix mult
# mydbs (yv) = mydbs mypx[-1] + mybs Ideg mypx 

# #now make things conformable? things go wrong here...
# 0mydbs (0yv) = 0mydbs mypx + mybs Ideg mypx  
# 0mydbs (0yv) = (0mydbs + mybs Ideg)  mypx  
# #now flip
# mypv = inv(0mydbs + mybs Ideg) 0mydbs 0yv
# mypv[-1] = inv(0mydbs + mybs Ideg) 0mydbs 0yv


# cbind(mydbs) %*% c(myyv)
# Ideg = diag(mydeg,ncol(mybs))
# mypx = ginv( mydbs + (mybs %*% Ideg) ) cbind(1,mydbs) 0yv

# #but non conf.



#################################################################
########Linear?
################################################################################
library(MASS) 
library(deSolve) 
# install.packages('txtplot')
library(splines2)
library(txtplot)

flankn=1
ribo <-c(rep(100,flankn),200,400,800,400,200,rep(100,flankn))
time = 1:length(ribo)


library(splines)
library(splines2)
library(splines)
library(tidyverse)
library(magrittr)
timeknots <- time[c(-1,-length(time))]
mybsdf1 <- cbind(1,ibs(time, knots = timeknots,degree = 1, intercept = TRUE))
mydbsdf1 = bs(time, knots = timeknots,degree = 1, intercept = TRUE)


dim(mybsdf1);dim(mydbsdf1)
mybsdf1%>%{ginv(.)%*%.}%>%round(10)%>%{.-diag(ncol(.))}%>%abs%>%max
mydbsdf1%>%{ginv(.)%*%.}%>%round(10)%>%{.-diag(ncol(.))}%>%abs%>%max
stopifnot(!any(mydbsdf1<0))
stopifnot((!all(mydbsdf1[1,]==0)))

invmydbsdf1 = ginv(mydbsdf1)
mydbsdf1_lims  = cbind(mydbsdf1) %>% {.[c(1,nrow(.)),]}
#this is all unncessary with piecewise linear but meh
myddbs_lims  = cbind(0,0,iSpline(time,length(time)-1,derivs=2)*0)%>%{.[c(1,nrow(.)),]}
mydddbs_lims  = cbind(0,0,iSpline(time,length(time)-1,derivs=3)*0)%>%{.[c(1,nrow(.)),]}
# mydbsdf1_ddlims = cbind(mSpline(time,length(time)-1,derivs=2))%>%{.[c(1,nrow(.)),]}
dim(myddbs_lims)
library(tidyverse)

dzv=lm( rep(1,nrow(mydbsdf1)) ~ 0+mydbsdf1)%>%.$coef%>%round(10)%>%replace_na(0)
zv=lm( rep(1,nrow(mybsdf1)) ~ 0+mybsdf1)%>%.$coef%>%round(10)%>%replace_na(0)



dim(mybsdf1);dim(mydbsdf1);
dim(mybsdf1);dim(myddbs_lims);length(mydddbs_lims)


basis_ribolinlim = rbind(
  mybsdf1,
  myddbs_lims,
  mydddbs_lims
)

dbasis_ribolinlim = rbind(
  mydbsdf1,
  myddbs_lims,
  mydddbs_lims
)

extnum = nrow(limbasis) - length(prot)

txtplot(prot)
txtplot(ribo)



mypx_lm = lm(prot ~ 0+mybsdf1[,-1])




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





yv = matrix(mypx[-1]) + invmydbs %*% ((mybs %*% mypx))  * (mybs %*% degzv)

#how about also no!
# yv = tail(matrix(mypx),-1) + ((invmydbs %*%mybs %*% mypx)  * (invmydbs %*% mybs %*% zv)) 
est = mydbs %*% yv
est/(ribo)
mean(est/(ribo))

# dbasis %*% yv = dbasis %*% px +  mybs %*% zv  %.% mybs %*% zv
#estimate using spline of P


dbasis %*% yv = dbasis %*% px +  mybs %*% px  %.% mybs %*% zv

#we need to think about the protein fold change.
P = mybs %*% px
dP = mydbs %*% px

#fold change given 
dP/P = mydbs %*% px / mybs %*% px
mydbs %*% px / mybs %*% px = 

# 



################################################################################
########Optimization
################################################################################
  

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
    parvect[1:2] = exp(parvect[1:2])
    stopifnot(!any(is.na(parvect)))
    stopifnot(!any(!is.finite(parvect)))
  with(as.list(parvect),{
    Ks = parvect[1]
    Kd = parvect[2]
    pxopt = parvect[-c(1:2)]
    degzv = zv*Kd
    
    yv = pxopt[-1] + invmydbs %*% ((mybs %*% pxopt)  * (mybs %*% degzv)) 
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

    nsconstraintLL = 
      dnorm(myddbs_lims %*% pxopt,0,sd=0.1,log=TRUE)+
      dnorm(mydddbs_lims %*% pxopt,0,sd=0.1,log=TRUE)
    nsconstraintLL = 0

    LL = sum(c(ribo_LL ,prot_LL,alphapriorLL,nsconstraintLL))
    # browser()
    if(!is.finite(LL)) browser()
    if(length(LL)==0) browser()
    # cat(LL);cat('...')

    - LL
  })
}

parvect = c(Ks=log(10),Kd=log(0.5),ms0=log(5000),mypxrest=rep(0,ncol(mybs)-1))
parvect[-c(1:2)]=mypx%>%replace_na(.,mean(na.omit(.)))
prot = num_ode_res$P
# parvect%<>%

lowervect = rep(-2e4,length(parvect))%>%setNames(names(parvect))%>%{.['Kd']=log(0.001);.['mypxrest1']=log(1);.}
stopifnot(all(is.finite(lowervect)))
uppervect = rep(2e4,length(parvect))%>%setNames(names(parvect))%>%{.['Kd']=log(1);.}
stopifnot(all(is.finite(uppervect)))
stopifnot(!any(is.na(parvect)))
stopifnot(!any(!is.finite(parvect)))


stopifnot(
  linearsplineoptfunc(parvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot) >0
)

linearsplineoptfunc(parvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot)
linearsplineoptfunc(parvect%>%{.['Kd']=log(0.1);.},invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot)



par = optim(parvect,
  fn=linearsplineoptfunc,
  method='L-BFGS-B',
  lower = lowervect,
  upper = uppervect,
  # hessian = TRUE,
  invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,
  # control=c(pgtol=0.00000001)
  )
}

par$par == uppervect
par$par == uppervect

linearsplineoptfunc(par$par,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,returnprot=TRUE)%>%txtplot
linearsplineoptfunc(par$par,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,returnribo=TRUE)%>%txtplot

linearsplineoptfunc(parvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,returnprot=TRUE)%>%txtplot
linearsplineoptfunc(parvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot,returnribo=TRUE)%>%txtplot

linearsplineoptfunc(parvect,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot)
linearsplineoptfunc(par$par,invmydbs=invmydbs,zv=zv,mybs=mybs,ribosd=0.1,protsd=0.3,ribo=ribo,prot=prot)


prot%>%txtplot
ribo%>%txtplot
# }


# install.packages('brm')
# library(brm)


# brm(ribo ~ ,data=data.frame(ribo=ribo))



# 0
# debugger()


# #FUCK
# linear splines end up negative. No way to 


# ribo = c(rep(100,10),200,400,800,400,200,rep(100,10))
# ##Okay so if we can use M splines for the above, what about goiing back to log level





# lR = lP + log(dlP + Kd) - d

# now let F = dLP + Kd be an M spline

# lR = lP + log(F) - d

# lR = int(M%*%(px-Z)) + log(F) - d

# lR = I%*%(px-Z) + log(M %*% (px))) - d

# lR = I%*%(px-Z) + log(M %*% (px))) - d

# lR = I%*%(px-Z) + log(Ks %*% M %*% (px))) - d


# yv = pxopt[-1] + invmydbs %*% ((mybs %*% pxopt)  * (mybs %*% degzv)) 

# yv = pxopt[-1] + invmydbs %*% ((mybs %*% pxopt)  * (mybs %*% degzv))



