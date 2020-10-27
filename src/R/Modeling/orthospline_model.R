library(deSolve)
library(txtplot)
library(minpack.lm)
library(tidyverse)
library(splines)
library(SplinesUtils)
library(data.table)
library(nlme)
library(nlmeODE)
library(BiocParallel)
library(splines2)


################################################################################
########some math
################################################################################
 
#The basic euqtoins with splnes
dP(t) = Ks R(t) - Kd P(t)

#get this done with 
R(t) = dP(t) + Kd(t) P(t) / Ks

#
P(t) = exp(Pl(t))
dP(t) = dPl(t)exp(Pl(t))

log(R(t)) = log( dPl(t) exp(Pl(t)) + Kd exp(Pl(t))  ) - log(Ks)

log(R(t)) = Pl(t) + log( dPl(t) + Kd ) - log(Ks)

#This is a problem for parametrization. So what if instead of initializing the parameters of Pl(t)
#We instead initialize the parameters of a split f(t) which describes dPl(t) + Kd ? Then we'd
#Still have to avoid it getting negative sometimes.
#We could have a spline which defines the log value of dPl(t) + Kd. It's exponent would then be


log(R(t)) = int(exp(F((t)))-Kd) + log( exp(F(t)) ) - log(Ks)




log( dPl(t) + Kd ) - log(Ks)



Our solution is 

A) exp(-D t)(m0 + Integrate_0_t(T R(t) * exp(D t )))

#if we use a spline to model log R

exp(-D t)(m0 + Integrate_0_t(T exp(Rspl(t)) * exp(D t )))


#So  here's  what solfram tells us about a segment with  linear changes in  degredation and
#(log)  syntehsis
B) integral exp(a_0 + a_1 t + d_0 t) dt = e^(t (a_1 + d_0) + a_0)/(a_1 + d_0) + constant

Plug  this into A)

exp(- d_0 t)(m0 +  exp(t (a_1 + d_0) + a_0)/(a_1 + d_0)  )

m0 exp(-d_0 t) + { exp(a_1 t + d0 t - d0 t + a_ 0)  / (a_1 + d_0 ) }

m0 exp(-d_0 t) + { exp(a_1 t + a_0)  / ( a_1 + d_0 ) }


#So here we have a_0 and a_1 as the start and slope of snthesis respectively across a time period
P(0) = M0

P(t1) = m0 exp(- d_0 t) + { exp(a_1 t + a_0)  / ( a_1 + d_0 ) }
 
P(t2) = P(t1) exp(-d_0 (t_d2)) + { exp(a_1 t_d2 + a_0)  / ( a_1 + d_0 ) }  





#but now re arrange this - we want RNA parameters
#we can parametrize them    Rst(t), t in 2:5 =  Rsl(t-1), a1(t) = R(t) - R(t-1)

P(t1) = m0 exp(- d_0 t) + { exp(Rsl(t1) t + Rst(t1))  / ( Rsl(t) + d_0 ) }


#Now looking at just that integral
Integrate(exp(log(T) + Rspl(t) + Dt))

#We can easily integrate splines -  so can we make that thing above into another spline
Integrate(exp(log(T) + Rspl(t) + Dt))
Integrate(exp(log(T) + Rspl(t) + Dt))


# can get the integral of a spline basis easily


#The fold change needs to constrained...
#we are essentially breaking down the fold change into a constant component - the
#degredation, dP(t) must be greater or euqal to Kd
#but dP(t) is a linear combination of bunch of parameters


R(t) =  { dPl(t) exp(Pl(t)) + Kd exp(Pl(t))   }/  Ks


log(R(t)) = Pl(t) + log( dPl(t) + Kd ) - log(Ks)

R(t)) = Pl(t) + log( dPl(t) + Kd ) - log(Ks)


if (k_d >> dPl(t))

then log(R(t)) = Pl(t) + log(Kd) - log(Ks)





log(R(t)) = Pl(t) + log( dPl(t) + Kd ) - log(Ks)

log(R(t)) = Pl(t) + log( dPl(t) + Kd ) - log(Ks)







################################################################################
########Math  - this time linear space splines
################################################################################
# So I think this is all fine if we just use linear space splines
# IF we can express a time constant function as a simple combination of splines.
# Now, our linear thing is a 


time = 0:23 + 10
sig = 0:23 + 10

dbs%>%args

fit=lm(sig ~ 0+ bs(time,intercept=TRUE))
fit
fit_x2=lm(sig*2 ~ bs(time,intercept=TRUE))
fit_x2


pdf('tmp.pdf')
bs(time,intercept=TRUE)%>%matplot
dev.off()
normalizePath('tmp.pdf')

dbs(time,intercept=TRUE) %*% (fit$coef)
bs(time,intercept=TRUE) %*% (fit$coef)
ibs(time,intercept=TRUE) %*% (fit$coef)

#so if we want our RNA to be 1+4 dimensions for an arbitrary fit
#then we can have 
#and we want our prot to be 1+4 dimensions for arbitrary fit
inter_dbasis = cbind(1,dbs(time,df=4))
integ_inter_dbasis = cbind(1,1:(1*length(time)),bs(time,df=4))





#So weith default df = 3
#so note that with an intercept, bs has only 3
princomp(bs(time))%>%summary
princomp(bs(time,intercept=TRUE))%>%summary
#And dbs only has two
princomp(dbs(time,intercept=TRUE))%>%summary
#so to fit e.g. 5 timepoint in our rate of change, we'd need 6 
mydf = 6
#we would like there to be a meaningful intercept at both levels actually
princomp(dbs(time,df=mydf,intercept=TRUE))%>%summary

let bs be the linear spline basis
px = ([m0,a1,a2,a3])
so that bs %*% pv, bs is a txD matrix
dbs is now also txD
let zv be teh coefficients of a linear fit in our spline basis

dP(t) = Ks P(t) - Kd (P(t))



zv = lm(time ~ 0 + inter_dbasis)$coef
zv = lm(1:24 ~ 0 + integ_inter_dbasis[,,-5])$coef


note that px1 is t+1 for a perfect  


dbasis = inter_dbasis
mybs = integ_inter_dbasis



###actually hadamard producct

#Preeeety sure this works
t,D    D,1    =   t,D     D,1 -  t,D     D,1  %.%  t,D1   D1,1
dbasis %*% px = dbasis %*% yv -  dbasis %*% zv %.%  mybs %*% px1

dbasis %*% px = dbasis %*% yv -  dbasis %*% zv %.%  mybs %*% px1
#add
dbasis %*% px + dbasis %*% zv %.%  mybs %*% px = dbasis %*% yv
#flilp
dbasis %*% yv = dbasis %*% px + dbasis %*% zv %.%  mybs %*% px1 

#Okay now can I factorize this???
#note the hadamard (elementwise) product there
#t,D      D,1 = t,D    %*% D,1 +     t,D   . D,1        t,D1 %*% D1,1
dbasis %*% yv = dbasis %*% px +     dbasis %*% zv %.% (mybs %*% px1)

dbasis %*% yv = dbasis %*% px +     dbasis %*% zv %.% (mybs %*% px1)

#Okay so no, this isn't a linear operation, but it's lcose
#t,D      D,1 =  t,D  %*% D,1 +  t,D   .    D,1  %*%    
dbasis %*% yv = dbasis %*% px + diag(mybs %*% px1) %*% dbasis %*% t(zv) 
#Now let's assume that there's a left inverse to dbasis - i.e. it has linearly independent columns
dbasis %*% yv = dbasis %*% px + diag(mybs %*% px1) %*% dbasis %*% t(zv) 

#now we end up with this.......
yv = px + dbinv %*% diag(mybs %*% px1) %*% dbasis %*% t(zv) 


#We COULD also just express our zv in terms of the protein level basis
dbasis %*% yv = dbasis %*% px +  (mybs %*% zv) %.% (mybs %*% px1)

dbasis %*% yv = dbasis %*% px +  mybs %*%  (zv %.% px1)

#again assume linear independence in the columns of dbasis

dbasis %*% yv = dbasis %*% yv +  mybs %*%  (zv %.% px1)
D,1 =D,1 +       D,t   %*% t,D1  %*%   D1,1

#this is WAY nicer

                           t,D1       D1,1
yv = px + invdbasis %*% mybs %*%  (zv %.% px1)

#waht if we assume that zv[1] = 0?
yv = 


px + invdbasis %*% mybs %*%  (zv %.% px)


#now in orthogonolaizing our basis... we need to think about how the coefficients are effected.

#hang on, if one of our basis terms is the intercept in dbs space, then every term but one in




################################################################################
########Numerically test the linear space spline solution
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
yv = px + dbinv %*% diag(mybs %*% px1) %*% dbasis %*% t(zv) 


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
yv = px + invdbasis %*% mybs %*%  (zv %.% px)

library(orthogonalsplinebasis)

#This gets us an orthogonal spline matrix
inter_dbasis=evaluate(OBasis(expand.knots(seq(first(time),last(time),length.out=3),4)),time)
round(ginv(inter_dbasis) %*% inter_dbasis,4)

#But what if we want one of the columns to be an interept? Then
intb %*% diag(v) =  

inter_dbasis = cbind(1,dbs(time,df=6,intercept=FALSE))
inter_dbasis = princomp(inter_dbasis)
inter_dbasis = inter_dbasis$scores 


OBasis(expand.knots(0,))


integ_inter_dbasis = cbind(1,1:(1*length(time)),bs(time,df=4))

SplineBasis(expand.knots(seq(first(time),last(time),length.out=3),4))

#how can I arrive at an orthonormal basis that includes an intercept?
#What if I get an orthonormal basis for the differences.
#then I arrive at the orthonormal basis for the actual space
inter_dbasis=evaluate(OBasis(expand.knots(seq(first(time),last(time),length.out=3),4)),time)
BiocManager::install('pracma')

library(pracma)
A <- matrix(c(0,-4,2, 6,-3,-2, 8,1,-1), 3, 3, byrow=TRUE)
gs <- gramSchmidt(A)
(Q <- gs$Q); (R <- gs$R)
Q %*% R  # = A

mydbs = cbind(1, bs(time))
mydbs%>%{ginv(.)%*%.}%>%round(10)
invmydbs = ginv(mydbs)
mybs = cbind(1,1:length(time),ibs(time))

#wait, fuck lol c(1,bs) without intercept option is already pretty much orthogonal

################################################################################
########Now test the above
################################################################################
library(MASS) 
mydbs = cbind(1, bs(time,df=(length(time)-1)))
mydbs%>%{ginv(.)%*%.}%>%round(10)
invmydbs = ginv(mydbs)
mybs = cbind(1,1:length(time),ibs(time,df=(length(time)-1)))

splfit <- (lm(ribo ~ 0+mydbs))
predict(splfit)%>%txtplot
stopifnot(!any(is.na(splfit$coef)))
mydeg = 0.04
myrTE  = 20
num_ode_res<-spl_getprotdf(splfit,deg = mydeg,rTE=myrTE,orthns_imat=mydbs[,-1],logfit=FALSE)
#looks right
num_ode_res$ribo%>%txtplot
num_ode_res$P%>%txtplot
mypx = lm(num_ode_res$P ~ 0+mybs)
mypx = mypx$coef%>%replace_na(0)


#zv is just the c(deg,0,0,....)
zv = rep(0,length(mypx))
zv[1] = mydeg

#now use yv = px + invdbasis %*% mybs %*%  (zv %.% px)
length(ribo)
length(num_ode_res$P)
dim(invmydbs)
dim(mydbs)
dim(mybs)
length(mypx)
length(zv)

yv = mypx[-1] + invmydbs %*% mybs %*% (zv * mypx)

(mydbs%*%yv) / num_ode_res$ribo

txtplot(mydbs%*%yv,ribo)

mydbs %*% yv

num_ode_res$P

txtplot(mydbs %*% yv,ylim=c(0,max(num_ode_res$ribo)))

#this looks mostly correc tthough I wonder where those weird little innacuracies are coming from.
yv = px + invdbasis %*% mybs %*%  (zv %.% px1)




################################################################################
########Try minimizing it
################################################################################
ssq=function(parms,ribo,P,dPdt=dPdt){
  
  dPdt <- function(itime, state, parameters){
   dPdt <- with(as.list(c(state,parameters)) ,{
     mRNA <- c(1,bs(itime,Boundary.knots=bounds,degree=df)) %*% scoefs
     mRNA <- max(0,mRNA)
     (rTE * mRNA) - (deg*P)
   })
   list(dPdt)
  }


  ntps = length(ribo)
  splinedf =  4
  bsbasis =  bs(1:ntps,degree=splinedf)

  deg = exp(parms['ldeg'])
  rTE = exp(parms['lrTE'])
  time=1:ntps

  bsdata=data.frame(ribo,time)
  bsfit <- (lm(data=bsdata,ribo ~ bs(time,degree=splinedf)))

  conf_interval <- predict(bsfit,, interval="confidence",
                           level = 0.95)
  #qplot(x=1:ntps,y=ribo,geom='blank')+geom_ribbon(aes(ymin = conf_interval[,2],ymax = conf_interval[,3],fill=I('grey')))+geom_point()+theme_bw()


  state<-c('P'= min(
    ((rTE*ribo[1])/deg),
    ((rTE*ribo[1])/0.1)
    ))

  parameters = list(rTE=rTE,bounds=c(1,ntps),deg=deg,df=splinedf,scoefs = bsfit$coef)%>%as.list
  Pdf = ode(y = state, times = 1:ntps, func = dPdt, parms = parameters)

  Pdf[,2] - P


}

################################################################################
########Test the OLS fitting approach on some toy data
################################################################################
rTE_vect <- c(.1,1,10)
deg_vect <- c(0.1,0.5,0.9)
parmgrid <- expand.grid(rTE=rTE_vect,deg=deg_vect)
prow=20
ribovect<-c(100,200,1000,1200,1000)
testdata <- lapply(1:nrow(parmgrid),safely(function(prow){
  parms=c(deg=parmgrid[prow,'deg'],rTE=parmgrid[prow,'rTE'])
   
  out <- simdata<-lapply(1:10,function(i){
      ribo <- sample(ribovect,replace=TRUE)


    # parameter fitting using levenberg marquart algorithm
    # initial guess for parameters

    # fitting

    out<-getprotdf(deg = deg,rTE=rTE,ntps=length(ribo),ribo=ribo)
    out
  })

    attr(out,'deg') <-parms['deg']
    attr(out,'rTE') <-parms['rTE']
  out
}))
attr(testdata[[1]][['result']],'deg')
alltestdata <- testdata%>%map('result')%>%map(~ tibble(deg=attr(.,'deg'),rTE=attr(.,'rTE'),data=.))%>%bind_rows
# jtestdata <- testdatalist%>%bind_rows(.id='actual_deg')
# ptestdata<-testdata%>%gather(assay,signal,-actual_deg,-time)

# ptestdata%>% # {
#   qplot(data=.,color=actual_deg,x=time,y=signal,geom='blank')+facet_wrap(scale='free',ncol=1,~assay)+geom_line()+geom_point()+theme_bw()+scale_y_log10()
# }
library(parallel)

#
parmgrid <- expand.grid(rTE=rTE_vect,deg=deg_vect)
prow=20
# paramtries <- mclapply(1:nrow(parmgrid),safely(function(prow){
data = alltestdata$data[[1]]


################################################################################
########try this out (or on real data!)
################################################################################

{
  require(BiocParallel)
  BiocManager::install('batchtools')
  param<- BatchtoolsParam(workers=1, cluster="sge", resources=list(queue='all'))
  sge_template <- '~/tmp.tmpl' 
  param$template%>%readLines%>%str_replace(regex('\\#\\$ -q .*'),'')%>%cat(file=sge_template,sep='\n')
  param<- BatchtoolsParam(workers=100, cluster="sge", resources=list(queue='all'),template=sge_template)
  #bplapply(BPPARAM=param
}
library(GenomicRanges)
#result <- bplapply(BPPARAM=param, c('1:1','1:3'), function(x)GenomicRanges::GRanges(x) )

allexprdata <- fread('pipeline/exprdata/limma_proDD_CIs.tsv')
allexprdata_formed <- allexprdata%>%select(uprotein_id,gene_name,assay,signal,time)
library(magrittr)
allexprdata_formed$assay%<>%ifelse(.=='MS','P',.)
allexprdata_formed$signal%<>%{exp(log(2)*.)}
allexprdata_formed%<>%spread(assay,signal)%>%group_by(uprotein_id,gene_name)%>%nest
allexprdata_formed%<>%group_by(gene_name)%>%mutate(changeP=map_lgl(data,~abs(log(first(.$P)/last(.$P)))>1))
allexprdata_formed%<>%mutate(totest = changeP)
allexprdata_formed%<>%filter(totest)


data<-allexprdata_formed$data[[1]]

#estimatesbak<-estimates
estimates <- bplapply(BPPARAM=param,allexprdata_formed$data,ssq=ssq,parmgrid=parmgrid,function(data,ssq,parmgrid){
  library(deSolve)
  library(minpack.lm)
  library(tidyverse)
  library(splines)
  #using this data
  ribo <- data$ribo
  P <- data$P
  ssq_p = partial(ssq,ribo=ribo,P=P) 
  #search the parameter grid until we get some kind of reasonable answer
  for(prow in sample(1:nrow(parmgrid))){
    # estimate <- safely(function(prow){
    estimate <- safely(function(prow){
    # parameter fitting using levenberg marquart algorithm
    # initial guess for parameters
    parms=c(ldeg=log(parmgrid[prow,'deg']),lrTE=log(parmgrid[prow,'rTE']))
    # fitting
    fitval=nls.lm(par=parms,fn=ssq_p,upper=c(ldeg=log(0.999),lrTE=log(1e6)),lower=c(ldeg=log(0.0001),lrTE=log(1/1e6)))
    cinf=confint(fitval)
    list(fitval,cinf)
    }) (prow)
    cat('.')
    if(!is.null(estimate[['result']])) break
  }
  estimate
})

converged <- estimates%>%map('result')%>%map_lgl(is.null)%>%not
converged%>%table

test_that("the simulations make some kind of sense",{

  alltestdata$isinconfint<-map2( alltestdata$deg,estimates%>%map('result')%>%map(2)%>%map(~.['deg',]),
    possibly(~ (.y[1] < .x ) & (.x < .y[2]),otherwise=NA)
  )
  alltestdata$isinconfint[alltestdata$isinconfint%>%map_lgl(~length(.)==0)] <- NA
  alltestdata$isinconfint%>%flatten_lgl

  alltestdata%>%select(-data)%>%unnest%>%group_by(deg,rTE,isinconfint)%>%tally%>%spread(isinconfint,n)


  alltestdata%>%select(-data)%>%unnest%>%group_by(deg,isinconfint)%>%tally%>%spread(isinconfint,n)

  alltestdata%>%select(-data)%>%unnest%>%group_by(rTE,isinconfint)%>%tally%>%spread(isinconfint,n)

  alltestdata%>%filter(deg==0.9)%>%select(-data)%>%.$isinconfint

  paramtries%>%map('result')%>%keep(Negate(is.null))%>%sample(1)

  paramtries[[1]]
  debugger()

})




converged <- estimates%>%map('result')%>%map_lgl(is.null)%>%not
converged%>%table

estimates[converged]%>%map('result')%>%map(2)%>%map(~.['lrTE',TRUE]%>%unlist%>%range)%>%setNames(allexprdata_formed$gene_name[converged])%>%bind_rows%>%t
estimates[converged]%>%map('result')%>%map(2)%>%map(~.['ldeg',TRUE]%>%unlist%>%range)%>%setNames(allexprdata_formed$gene_name[converged])%>%bind_rows%>%t

estimates[converged]%>%map('result')%>%map(2)%>%map(~.['ldeg',TRUE]%>%unlist%>%range)%>%setNames(allexprdata_formed$gene_name[converged])%>%bind_rows%>%t%>%.[.[,2]<log(1),]%>%apply(1,function(x) abs(first(x)-last(x)))%>%log10%T>%txtdensity%>%identity%>%exp


estimates[converged]%>%map('result')%>%map(2)%>%map(~.['ldeg',TRUE]%>%unlist%>%range)%>%setNames(allexprdata_formed$gene_name[converged])%>%bind_rows%>%t%>%.[.[,2]<log(1),]%>%apply(1,range)%>%txtdensity




library(txtplot)
estimates%>%map('result')%>%map(2)%>%map(.)

allexprdata_formed


estimates[converged]%>%map('result')%>%map(2)%>%map(~.['ldeg',TRUE]%>%unlist%>%range)%>%setNames(allexprdata_formed$gene_name[converged])%>%bind_rows%>%t%>%set_colnames(c('lower','upper'))%>%
  as.data.frame%>%
  filter(upper < log(1),lower> log(0.01))

estimates[converged]+

estimates[converged]%>%map('result')%>%head(1)%>%s
estimates%<>%setNames(allexprdata_formed$gene_name)
ldegestimates <- estimates[converged]%>%map('result')%>%map(1)%>%map('par')%>%map_dbl('ldeg')

goodests <- ldegestimates[ldegestimates > -2]



mcshane_df<-'ext_data/mcshane_etal_2016_S1.csv'%>%fread

mcshane_df%>%select(gene_name=`Gene names`,mu =`Half-life (exponential) 1-state-model [h]`)%>%mutate(mu=as.numeric(mu))%>%inner_join(enframe(goodests,'gene_name','ldeg'))%>%{cor.test(.$mu,.$ldeg,use='complete')}




exp(-2)

# BiocManager::install(c('splines2'))


splinedf<-4
# poly(time,splinedf)
# poly%>%args

ntps <- 5
time <- seq_len(ntps)
ribo= c(100,200,500,1000,500)%>%{rep(.,each=ntps/length(.))}%>%rpois(length(.),.)
# Q
protodedatasub <- protodedata%>%filter(actual_deg==0.6)
ribo <- protodedatasub$ribo
pfit <- lm(data=protodedatasub, ribo ~ 1 + poly(time,splinedf,raw=TRUE))
plot(ribo)
points(predict(pfit),type='line')

# plot(ribo)
# points(map_dbl(protodedatasub$time, ~ pfit$coef[1] + (.)*pfit$coef[2] + (.^2)*pfit$coef[3] + (.^3)*pfit$coef[4] + (.^4)*pfit$coef[5]),type='l')


# map(1:35,~{
#   e<-new.env()
#   e$t<-.
#   eval(
#     formexpr,e)
# })
# polyformula <- as.formula(polystring)

# map(1:35,~{
#   e<-new.env()
#   e$t<-.
#   eval(
#     formexpr,e)
# })

# library(testthat)
# test_that("Formula works",{

# library(nlmeODE)


# mRNAstrng<- paste0(' + ',map_chr(1:(splinedf+1),~paste0( pfit$coef[1+.],'* (t ^ ',.,') ')) , colllapse='')%>%paste0(collapse='+ ')%>%paste0(pfit$coef[1],.)
# formexpr<-parse(text=mRNAstrng)[[1]][[2]]
# iformvals <- map(moddata$time[1:35],~{
#   e<-new.env()
#   e$t<-.
#   eval(
#     formexpr,e)
# })

# polystring <- c(
#    paste0('~ rTE * ( ',,
#    map(1:splinedf,~paste0( pfit$coef[1+.],'* (t ^ ',.,') '))
# )%>%
# paste0(collapse=' + ')%>%paste0(' ) - (deg * P )')
# polyformula <- as.formula(polystring)


# expect_gt( map_dbl(1:35,~{
#   e<-new.env()
#   e$t<-.
#   eval(
#     formexpr,e)
# }),0)

# })



#extract some
moddata <- testdata%>%filter(actual_deg%in%c(0.9))%>%select(P,time,actual_deg,ribo)

  pfit <- lm(data=moddata, ribo ~ 1 + poly(time,splinedf,raw=TRUE))

  polystring <- c(paste0('~ rTE * ( ',pfit$coef[1]),
    map(1:(splinedf+1),~paste0( pfit$coef[1+.],'* (t ^ ',.,') '))
  )%>%paste0(collapse=' + ')%>%paste0(' ) - (deg * P )')

  formexpr<-parse(text=polystring)[[1]][[2]]

protodedata <- moddata %>%list%>%rep(3)%>%bind_rows%>%mutate(P=exp(log(P)+rnorm(n(),sd=0.1)))%>%
  arrange(time)%>%
  groupedData( P ~ time | actual_deg, data =  .)

{

# setting 1st compartiment of diff. equation (1 of 1)
degsynth <- list( DiffEq = list(
                      # dPdt = ~ - (rTE * bSpline(t,Boundary.knots=bounds,degree=df) - (deg*P))
                      # dPdt = ~ - (rTE * t + t^2) - (deg*P)
                      dPdt = eval(polyformula)
                      # dTdt = ~ T / T
                  ),
                 ObsEq  = list( 
                   P ~ P),
                 Parms  = c("deg", "rTE"),
                 States = c("P"),
                 Init   = list(TRUE,TRUE))

# deriv <- function(...){
#   browser()
#   stats::deriv(...)
# }

dmodel = nlmeODE( degsynth, protodedata, 
                 LogParms = T, JAC = T, SEQ = F, rtol = 0.01, atol = 0.01)

grp_n <- protodedata$actual_deg%>%n_distinct

nlme_res <- degsynth <- nlme( P ~ dmodel( deg, rTE, time, actual_deg ),
                   data    = protodedata,
                   fixed   = list( deg+rTE ~ 1) ,
                   random  = pdDiag( ~ 1),
                   start   = c(deg = 0.5, rTE = 10),
                   control = list( returnObject = T, msVerbose = T ),
                   verbose = T)


intervals(nlme_res,which='fixed')[[1]]%>%exp


}


expression({
    .value <- -(rTE * T + T^2) - deg * P
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("rTE",
        "T", "deg", "P")))
    .grad[, "rTE"] <- -T
    .grad[, "T"] <- -(rTE + 2 * T)
    .grad[, "deg"] <- -P
    .grad[, "P"] <- -deg
    attr(.value, "gradient") <- .grad
    .value
})

library(Deriv)

Deriv%>%args

myexpr <- ~ rTE()

custom_deriv <- function(expr, namevec){
  Dout <- Deriv(expr,x=c('rTE','T','deg','P'))  
  expression({
    expr 
  })
}


custom_deriv( ~ rTE  + t^2)
Dout <- Deriv( ~ rTE  + t^2)

nlEnv%>%names


with(nlEnv,{
    pars <- getParsNlme(plist, fmap, rmapRel, bmap, groups, beta,
        bvec, b, level, N)
    # res <- eval(model, data.frame(data, pars))
    eval(model, data.frame(data, pars))
    # if (!length(grad <- attr(res, "gradient"))) {
        # grad <- finiteDiffGrad(model, data, pars)
    # }

    # for (nm in names(plist)) {
    #     gradnm <- grad[, nm]
    #     if (is.logical(f <- plist[[nm]]$fixed)) {
    #         if (f)
    #             X[, fmap[[nm]]] <- gradnm
    #     }
    #     else X[, fmap[[nm]]] <- gradnm * f
    #     for (i in 1:Q) {
    #         if (is.logical(r <- plist[[nm]]$random[[i]])) {
    #             if (r)
    #               Z[, rmap[[i]][[nm]]] <- gradnm
    #         }
    #         else {
    #             rm.i <- rmap[[i]][[nm]]
    #             if (data.class(rm.i) != "list") {
    #               Z[, rm.i] <- gradnm * r
    #             }
    #             else {
    #               for (j in seq_along(rm.i)) {
    #                 Z[, rm.i[[j]]] <- if (is.logical(rr <- r[[j]]))
    #                   gradnm
    #                 else gradnm * rr
    #               }
    #             }
    #         }
    #     }
    # }
    # result <- c(Z[naPat, ], X[naPat, ], res[naPat])
    # result[is.na(result)] <- 0
    # result
  }
)

library(SplinesUtils)
# Q


################################################################################
########The math for ouru integral of hte synthesis euqation
################################################################################
  
# library(SplinesUtils)

# library(magrittr)
# RegBsplineAsPiecePoly(bsfit,"bs(time, degree = splinedf)")
# bsfit2<-bsfit
# bsfit2$coef %<>% add(2)
# RegBsplineAsPiecePoly(bsfit2,"bs(time, degree = splinedf)")
# curve%>%args

#So - we can get an analytical solution if the degredation rate and synthesis are both polynomials
#exp(-(1+2t+3(t^2)))(mu +   integrate[   (( (a0)+(a1*(t))+(a2*(t^2))+(a3*(t^3)) ) * T * (exp( (1+2t+3(t^2)) ))  ) , t] )

#So - we can get an analytical solution if the degredation rate and synthesis are both polynomials
#exp(-(1+2t+3(t^2)))(mu +   integrate[   (( (a0)+(a1*(t))+(a2*(t^2))+(a3*(t^3)) ) * T * (exp( (1+2t+3(t^2)) ))  ) , t] )


#
e^(-3 t^2 - 2 t - 1) (1/324 e^(2/3) T (sqrt(3 π) erfi((3 t + 1)/sqrt(3)) (54 a0 - 18 a1 - 3 a2 + 7 a3) + 6 e^(1/3 (3 t + 1)^2) (9 a1 + a2 (9 t - 3) + a3 (9 t^2 - 3 t - 2))) + μ)

Protfun(t,T,){function(x) ( (exp(b1*T0) - 1) * exp(-b0 - (b1*T0) - (3*(b2^2))) * ( (1/2) *exp(2/3)*sqrt(π/3)*T*erfi( ((3*t) + 1)/sqrt(3)))*(a0 + x*(a1 + x*(a2 + a3*x))) + mu) /b1

ribovect <- c(100,100,1000,600,600)
testdata=data.frame(ribo=ribovect,time=1:5)
bsfit=lm(data=testdata,ribo~ 0+bs(time, 4))
ppoly<-RegBsplineAsPiecePoly(bsfit,"bs(time, 4)",shift=FALSE)
ppoly%>%str
ppoly

T_n=5
a_0=ppoly[[1]]$coef[1]
a_1=ppoly[[1]]$coef[2]
a_2=ppoly[[1]]$coef[3]
a_3=ppoly[[1]]$coef[4]
T=40
t=1
mu= (100*T/D)


predict(bsfit,newx=1,newdata=data.frame(time=1:5)) - predict(ppoly,newx=1:5)

a_0 + a_1*t + a_2*t^2 + a_3*t^3


Let Poly_n(x) be the polynomial, piecewise, defining the mRNA

Our integral is Int_t(Poly_n(t)*exp(exp(ldeg)) )  

Which is 



tmp<-CubicInterpSplineAsPiecePoly(1:5,ribo,'natural')
tmp%>%str
nsfit <- lm(data=data.frame(ribo,time=1:5),ribo ~ spline(time,4))
pfit <- lm(data=data.frame(ribo,time=1:5),ribo ~ poly(time,3))

RegBsplineAsPiecePoly(nsfit,'ns(time, 4)',shift=F)

predict(pfit)

predict(nsfit)

map_dbl(seq(1,2,len=20),function(t){
  a_0+a_1*t+a_2*t^2+a_3*t^3
})%>%txtplot(width=100)


map_dbl(seq(0,2,len=20),function(t){

exp(-D*T_n)*
(
  (T*(exp(D*t))* (a_3*D^3*t^3 + a_2*D^3*t^2 + a_0*D^3 - 3*a_3*D^2*t^2 + a_1*D^2*(D*t - 1) - 2*a_2*D^2*t + 6*a_3*D*t + 2*a_2*D - 6*a_3) )/
    (D^4) + mu
)
})%>%txtplot(width=100)


# undebug(SplinesUtils:::predict.PiecePoly)
library(zeallot)
library(tidyverse)


ribo = c(100,100,100,100,100,100,100,500,100,100,100,100,100,100,100,100)
pfit <- lm(data=data.frame(ribo,time=1:length(ribo)),ribo ~1+poly(time,3,raw=TRUE))


predict(pfit)
c(a_0,a_1,a_2,a_3) %<-% pfit$coef



#For cubics 

Poly <- function(t){ a_0 + a_1*t + a_2*t^2 + a_3*t^3}
dPoly <- function(t){      a_1 + 2*a_2*t + 3*a_3*t^2}
ddPoly <- function(t){           2* a_2 +   6*a_3*t}
dddPoly <- function(t){                      6*a_3}

library(txtplot)
plot(map_dbl(1:length(ribo),Poly))
dddPoly(1)

dl = log(0.5)
igrow = function(t) exp(exp(dl)*t)
igrow_int <- function(t) exp( (exp(dl)*t) )
ideg = function(t) exp(-exp(dl)*t)
ideg_int <- function(t) exp( -(exp(dl)*t)  )






map(1:5,ideg_int)
map(1:5,ideg_int)

#So Int(igrow is then) exp(d t )/d


#So for a cubicInt

Int(Poly(t)*igrow())
u = Poly(t)
v = igrow
u int(v) - int(du,int(v)) =>

Poly(t)*Int(igrow) - Int(dPoly(t)*int(igrow))

u = dPoly(t)
v = int(igrow)
u int(v) - int(du , int(v)) =>


u = ddPoly(t)
v = int(int(igrow))


Poly(t)*Int(igrow) - ( 
  dPoly(t)*int(Int(igrow)) - Int(ddPoly(t ) *int(int((igrow)))) 
)

u=dddPoly(t)
v=int(int(int(igrow)))

#So we arrive at a tractable integral
Poly(t)*Int(igrow) - ( 
  dPoly(t)*int(Int(igrow)) - (
      ddPoly(t) *int(int(int(igrow)))  - Int(dddPoly(t)*int(int(int(igrow))) 
)


Poly(t)*Int(igrow) - ( 
  dPoly(t)*int(Int(igrow)) - (
      ddPoly(t) *int(int(int(igrow)))  - dddPoly*Int(int(int(int(igrow)))) 
)


Poly(1:5)

dl=log(0.5)
dl
E=exp(1)
T1=5


int_igrow<-function(T1,d=dl)(-1 + E^(E^d*T1))/E^d
int2_igrow<-function(T1,d=dl)((-1 + E^(E^d*T1))*T1)/E^d
int3_igrow<-function(T1,d=dl)((-1 + E^(E^d*T1))*T1^2)/E^d
int4_igrow<-function(T1,d=dl)((-1 + E^(E^d*T1))*T1^3)/E^d

t=2



degsynth <- function(t) Poly(t)*int_igrow(t) - (
  dPoly(t)*int2_igrow(t) - (
      ddPoly(t)*int3_igrow(t)  - dddPoly(t)*int4_igrow(t) 
  )
)
ribo
degsynth(5)

t=length(ribo)
m0=100

Pfunc <- function(t) {ideg_int(t) * (m0 + degsynth(t))}

map_dbl(1:length(ribo),Pfunc)

library(txtplot)

#from ( integral_0^T0 exp(-(b0 + b1 t + 3 b2^2)) dt) (μ + integral(a0 + a1 x + a2 x^2 + a3 x^3) T exp(b0 + b t + b2 t^2) dt)
Integrate[Exp[-(b0 + b1 t + 3 b2^2)], {t, 0, T0}] (μ + Integrate[(a0 + a1 x + a2 x^2 + a3 x^3) T Exp[b0 + b t + b2 t^2], t])



ribo =  c(100,200,500,1000,500,500,500)%>%{rep(.,each=ntps/length(.))}%>%rpois(length(.),.)
ribo =  c(500,500,500,500,500,500,500)%>%{rep(.,each=ntps/length(.))}
splinedf=2
time=1:ntps
bsdata=data.frame(ribo,time)
bsfit <- (lm(data=bsdata,ribo ~ bs(time,degree=splinedf)))
#bsfit$coef[2]=1
#bsfit$coef[2]=1
ppoly<-RegBsplineAsPiecePoly(bsfit,"bs(time, degree = splinedf)")
plot(predict(bsfit))
plot(predict(ppoly,newx=time))
ppoly
bsfit$coef
curve({-7.88e-13 + 72.8 * (x - 1) + 1.76 * (x - 1) ^ 2} ,from=0,to=1)
plot(bsfit)






################################################################################
########exponential splines - don't work
################################################################################

####OOOOOKay problem the fold change can't be any lower than the degredation rate.
#So the intial spline parameterization needs to be in terms of the difference between the degredation
#rate and the fold change


orthns<-function(x,df,...){
  nsbasis <- ns(x,df,...)
  svd(nsbasis)$u
}
x=1:20

dorthns<-function(x,df,...){
  
  nsbasis <- ns(x,df=4)
  dnsbasis <- rstpm2::nsxD(1:20,df=4)
  svdns <- svd(nsbasis)
  nsbasis %*% t(svd(nsbasis)$v)
  svdns$u
}


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
# conf_interval <- predict(linfit,, interval="confidence",
                         # level = 0.95)
# qplot(x=1:ntps,y=ribo,geom='blank')+geom_ribbon(aes(ymin = conf_interval[,2],ymax = conf_interval[,3],fill=I('grey')))+geom_point()+theme_bw()
# conf_interval <- predict(splfit, interval="confidence",
                         # level = 0.95)
#qplot(x=1:ntps,y=ribo,geom='blank')+geom_ribbon(aes(ymin = conf_interval[,2],ymax = conf_interval[,3],fill=I('grey')))+geom_point()+theme_bw()


  dPdt <- function(itime, state, parameters){
   trans = if(logfit) log else identity
   dPdt <- with(as.list(c(state,parameters)) ,{
     mRNA <- exp(c(1,orthns_imat[itime,]) %*% scoefs)
     mRNA <- max(0,mRNA)
     (rTE * mRNA) - (deg*P)
   })
   list(dPdt)
  }


   trans = if(logfit) log else identity

  stopifnot(splinedf %in% 1:100)
  splineribo <- exp(predict(splfit))

  #set up the actual ode
  state<-c('P'= min(
    ((rTE*splineribo[1])/deg),
    # ((rTE*splineribo[1])/0.1)
    Inf
  ))

  # orthns_imat <- orthns((1:ntps)-1,splinedf)


  state = state*ms0ratio
  parameters = list(rTE=rTE,bounds=c(1,ntps),deg=deg,df=splinedf,scoefs = splfit$coef,orthns_imat=orthns_imat)%>%as.list
  Pdf = ode(y = state, times = 1:ntps, func = dPdt, parms = parameters)

  #discrete aprox
  # Pdf<-rep(NA,ntps) 
  # Pdf[0] <- state
  # for(i in 2:ntps){
  #   (rTE*((ribo[i-1]+ribo[i])/2))
  #   pdf[i] <- (rTE*((ribo[i-1]+ribo[i])/2))
  # }
  data.frame(ribo = splineribo,P=Pdf[,2],time=1:length(splineribo))


}

################################################################################
########Spline model basis tests
################################################################################
  
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


# log(R(t)) = Pl(t) + log(Kd) - log(Ks)

# orthns(time,)




splinemodel_E <- function(Ks,Kd,a_0,a_1,a_2,a_3,a_4,a_5,R,P,nsbasis,dnsbasis,rsig,psig){
  
  n = length(R)

  Pl_t = cbind(1,nsbasis) %*% c(a_0,a_1,a_2,a_3,a_4,a_5)
  dPl_t = dnsbasis %*% c(a_1,a_2,a_3,a_4,a_5)

  # Ks = max(c(abs(min(dPl_t)),Ks))+0.001

  #AHA ok so dPl_t cannot be any less than the degreddation rate
  E_lR = Pl_t + log(dPl_t + Kd) - log(Ks)

  E_lR

  # log(R(t)) = Pl_t + log(dPl(t) + Kd ) - log(Ks)

}

splinemodel_LL <- function(parvect,R,P,nsbasis,dnsbasis,rsig,psig){

  # stopifnot(c('Ks','Kd','a_0','a_1','a_2','a_3','a_4','a_5') %in% names(parvect) )
  stopifnot(c('Ks','Kd','a_0','a_1','a_2','a_3','a_4','a_5') %in% names(parvect) )

  with(as.list(parvect),{
  n = length(R)

  Pl_t = cbind(1,nsbasis) %*% c(a_0,a_1,a_2,a_3,a_4,a_5)
  dPl_t = dnsbasis %*% c(a_1,a_2,a_3,a_4,a_5)
  
  mindPl_t = min(dPl_t)
  Kd = mindPl_t + ex_Kd
  # Ks = max(c(abs(min(dPl_t)),Ks))+0.001

  #AHA ok so dPl_t cannot be any less than the degreddation rate
  E_lR = Pl_t + log(dPl_t + Kd) - log(Ks)
  LL_P = dnorm(P,Pl_t,psig,log=T)
  LL_R = dnorm(R,E_lR,rsig,log=T)
  -sum(LL_P,LL_R)
  })

}

splinemodel2_LL <- function(parvect,R,P,nsbasis,dnsbasis,rsig,psig){

  # stopifnot(c('Ks','Kd','a_0','a_1','a_2','a_3','a_4','a_5') %in% names(parvect) )
  stopifnot(c('Ks','ex_Kd','a_0','a_1','a_2','a_3','a_4','a_5') %in% names(parvect) )

  with(as.list(parvect),{
  n = length(R)

  Pl_t = cbind(1,nsbasis) %*% c(a_0,a_1,a_2,a_3,a_4,a_5)
  dPl_t = dnsbasis %*% c(a_1,a_2,a_3,a_4,a_5)
  
  mindPl_t = min(dPl_t)
  Kd = mindPl_t - ex_Kd
  # Ks = max(c(abs(min(dPl_t)),Ks))+0.001
  #AHA ok so dPl_t cannot be any less than the degreddation rate
  if(any(dPl_t+Kd < 0 ))browser()
  E_lR = Pl_t + log(dPl_t + Kd) - log(Ks)
  LL_P = dnorm(P,Pl_t,psig,log=T)
  LL_R = dnorm(R,E_lR,rsig,log=T)
  -sum(LL_P,LL_R)
  })

}

#ZIFA

splfit$coef
R = ribo
P = out$P

orthns(1:(last(ztime)*1000),4)[,]
seq(1,length(ztime)*1000)

P%>%txtplot
psplfit <- lm(log(P) ~ nsbasis)
psplfit%>%predict%>%exp%>%txtplot

eR <- splinemodel_E(Ks=20,Kd=0.01,
  a_0=psplfit$coef[1],a_1=psplfit$coef[2],a_2=psplfit$coef[3],a_3=psplfit$coef[4],a_4=psplfit$coef[5],
  a_5=psplfit$coef[6],
  R=ribo,
  P=out$P,
  nsbasis,dnsbasis,rsig,psig)

parvect=list(Ks=20,Kd=0.01,
  a_0=psplfit$coef[1],a_1=psplfit$coef[2],a_2=psplfit$coef[3],a_3=psplfit$coef[4],a_4=psplfit$coef[5],
  a_5=psplfit$coef[6])%>%{setNames(unlist(.,use.names=F),names(.))}

eR <- splinemodel_LL(parvect,
  R=log(ribo),
  P=log(out$P),
  nsbasis,dnsbasis,rsig=1,psig=20
)

rsig=1
psig=20

# parvect['Kd']=0.9

# optim(parvect,splinemodel_LL,R=log(ribo),P=log(out$P),nsbasis=nsbasis,dnsbasis=dnsbasis,rsig=rsig,psig=psig)

parvect2 <- parvect
names(parvect2)[2]<-'ex_Kd'
parvect2['ex_Kd'] <-0.2

results <- optim(parvect2,splinemodel2_LL,R=log(ribo),P=log(out$P),lower=list(ex_Kd=0.001),nsbasis=nsbasis,dnsbasis=dnsbasis,rsig=rsig,psig=psig)


with(as.list(parvect2),{min(dnsbasis %*% c(a_1,a_2,a_3,a_4,a_5))}) - results$par['ex_Kd']



log(R(t)) = int(exp(F((t)))-Kd) + log( exp(F(t)) ) - log(Ks)


################################################################################
########
################################################################################
log(R(t)) = Pl(t) + log( dPl(t) + Kd ) - log(Ks)

#This is a problem for parametrization. So what if instead of initializing the parameters of Pl(t)
#We instead initialize the parameters of a split f(t) which describes dPl(t) + Kd ? Then we'd
#Still have to avoid it getting negative sometimes.
#We could have a spline which defines the log value of dPl(t) + Kd. It's exponent would then be
log(R(t)) = Pl(t) + log( dPl(t) + Kd ) - log(Ks)
log(R(t)) = int(exp(F((t)))-Kd) + log( exp(F(t)) ) - log(Ks)
#This leaves us with this tricky exponential integral to deal with.....
log(R(t)) = int{exp(F((t)))-Kd}  + F(t) - log(Ks)
log(R(t)) = int{exp(F((t)))-Kd}  + F(t) - log(Ks)



exp(predict(splfit))
  