library(deSolve)
library(minpack.lm)
library(tidyverse)
library(splines)
library(SplinesUtils)
library(data.table)
# clear all variables
# rm(list=ls())

# necessary library for nlmeODE
library(nlme)
library(nlmeODE)

# # setting data
# data <- Theoph
# # Why are these needed??!?
# data$Dose[data$Time!=0] <- 0    # necessary to speed nlme and also to be able to plot(augPred(...))
# data$Cmt <- rep(1,dim(data)[1]) # without this optional data column, nlme fails!

# # create groupedData object that separates the data per Subject
# data.grp = groupedData( conc ~ Time | Subject, data = data )

# # get subset of Theoph that is equal to 1
# data.test = subset(data, Subject == 1 )

# # plot Time x conc
# plot(data.test$Time,data.test$conc)

# # ploting the data.grp and see all subjects
# labels = list(x = "Time since drug administration", y = "Theophylline serum concentration")
# units  = list(x = "(hr)", y = "(mg/L)")
# plot(data.grp, xlab = paste(labels$x, units$x), ylab = paste(labels$y, units$y) )

# # setting 1st compartiment of diff. equation (1 of 1)
# OneComp <- list( DiffEq = list( dy1dt = ~ -ka * y1,
#                                 dy2dt = ~ ka * y1 - ke * y2 ),
#                  #ObsEq  = list(c1 = ~ 0, c2 = ~ y2/CL*ke), #from pkmodels
#                  ObsEq  = list( SC ~ 0, Cp ~ y2 / CL * ke),
#                  Parms  = c("ka", "ke", "CL"),
#                  States = c("y1", "y2"),
#                  Init   = list(0,0)
#                  )
# data.grp%>%head(1)
# # create nlmeODE model
# model = nlmeODE( OneComp, data.grp, 
#                  LogParms = T, JAC = T, SEQ = F, rtol = 0.01, atol = 0.01)

# # apply nlme to model with random effects:
# #  ka + ke + CL (all variables)
data.nlme.ka_ke_CL <- nlme( conc ~ model( ka, ke, CL, Time, Subject ),
                   data    = data.grp,
                   fixed   = list( ka + ke + CL ~ 1) ,
                   random  = pdDiag(ka + ke + CL ~ 1),
                   start   = c( ka = log(1.65), ke = log(0.08), CL = log(0.05) ),
                   control = list( returnObject = T, msVerbose = T ),
                   verbose = T)

# # summary for the results of nlme
# summary( data.nlme.ka_ke_CL )
# # plot of the results, with population and subject fitting
# plot(augPred(data.nlme.ka_ke_CL, level=0:1))

# # updates nlme using random effects:
# #  ka + CL
# data.nlme.ka_CL = update(data.nlme.ka_ke_CL, random = pdDiag(ka + CL ~ 1))
# # summary for the results of nlme
# summary( data.nlme.ka_CL )
# # plot of the results, with population and subject fitting
# plot(augPred(data.nlme.ka_CL, level=0:1))

# # updates nlme using a more accurate error correction for this data
# #  it is observerd that data.nlme.ka_CL residuals increase over time
# #  as such weigths will reflect this
# data.nlme.ka_CL.error2 <- update(data.nlme.ka_CL, weights = varConstPower( const = 0.7, power = 0.3))
# # summary for the results of nlme
# summary( data.nlme.ka_CL.error2 )
# # plot of the results, with population and subject fitting
# plot(augPred(data.nlme.ka_CL, level=0:1))

library(splines)
dPdt <- function(itime, state, parameters){
 dPdt <- with(as.list(c(state,parameters)) ,{
   mRNA <- c(1,bs(itime,Boundary.knots=bounds,degree=df)) %*% scoefs
   mRNA <- max(0,mRNA)
   (rTE * mRNA) - (deg*P)
 })
 list(dPdt)
}

getprotdf<-function(
ntps = 35,
splinedf =  4,
deg = 0.1,
rTE = 10,
timeline =  1:ntps,
ms0ratio = 3,
ribo= c(100,200,500,1000,500,500,500) ) {
  ribo = ribo %>%{rep(.,each=ntps/length(.))}%>%rpois(length(.),.)

# conf_interval <- predict(linfit,, interval="confidence",
                         # level = 0.95)
# qplot(x=1:ntps,y=ribo,geom='blank')+geom_ribbon(aes(ymin = conf_interval[,2],ymax = conf_interval[,3],fill=I('grey')))+geom_point()+theme_bw()
time=1:ntps

bsdata=data.frame(ribo,time)

bsfit <- (lm(data=bsdata,ribo ~ ns(time,splinedf)))

# conf_interval <- predict(bsfit, interval="confidence",
                         # level = 0.95)
#qplot(x=1:ntps,y=ribo,geom='blank')+geom_ribbon(aes(ymin = conf_interval[,2],ymax = conf_interval[,3],fill=I('grey')))+geom_point()+theme_bw()



  state<-c('P'= min(
    ((rTE*ribo[1])/deg),
    ((rTE*ribo[1])/0.1)
  ))
  state = state*ms0ratio
  parameters = list(rTE=rTE,bounds=c(1,ntps),deg=deg,df=splinedf,scoefs = bsfit$coef)%>%as.list
  Pdf = ode(y = state, times = 1:ntps, func = dPdt, parms = parameters)
  # Pdf<-rep(NA,ntps) 
  # Pdf[0] <- state
  # for(i in 2:ntps){
  #   (rTE*((ribo[i-1]+ribo[i])/2))
  #   pdf[i] <- (rTE*((ribo[i-1]+ribo[i])/2))
  # }

  data.frame(ribo = ribo,P=Pdf[,2],time=1:ntps)


}

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
testdata <- lapply(1:nrow(parmgrid),safely(function(prow){
  parms=c(deg=parmgrid[prow,'deg'],rTE=parmgrid[prow,'rTE'])
   
  out <- simdata<-lapply(1:10,function(i){
      ribo <- sample(c(100,200,1000,1200,1000),replace=TRUE)


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

# ptestdata%>%
# {
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
  sge_template <- '~/tmp.tmpl' 
  param$template%>%readLines%>%str_replace(regex('\\#\\$ -q .*'),'')%>%cat(file=sge_template,sep='\n')
  paramfunc<- BatchtoolsParam(workers=1, cluster="sge", resources=list(queue='all'),template='tmp.tmpl')
}
library(GenomicRanges)
result <- bplapply(BPPARAM=param, c('1:1','1:3'), function(x)GenomicRanges::GRanges(x) )

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
estimates <- bplapply(BPPARAM=param,allexprdata_formed$data%>%head(1),ssq=ssq,parmgrid=parmgrid,function(data,ssq,parmgrid){
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
    fitval=nls.lm(par=parms,fn=ssq_p,upper=c(deg=log(0.999),rTE=log(1e6)),lower=c(deg=log(0.0001),rTE=log(1/1e6)))
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

estimates[converged]%>%map('result')%>%map(2)%>%map(~.['rTE',TRUE]%>%unlist%>%range)%>%setNames(allexprdata_formed$gene_name[converged])%>%bind_rows%>%t
estimates[converged]%>%map('result')%>%map(2)%>%map(~.['deg',TRUE]%>%unlist%>%range)%>%setNames(allexprdata_formed$gene_name[converged])%>%bind_rows%>%t




estimates%>%map('result')%>%map(2)%>%map(.)

allexprdata_formed





# BiocManager::install(c('splines2'))


# splinedf<-4
# # poly(time,splinedf)
# # poly%>%args

# ntps <- 5
# time <- seq_len(ntps)
# ribo= c(100,200,500,1000,500)%>%{rep(.,each=ntps/length(.))}%>%rpois(length(.),.)
# # Q
# protodedatasub <- protodedata%>%filter(actual_deg==0.6)
# ribo <- protodedatasub$ribo
# pfit <- lm(data=protodedatasub, ribo ~ 1 + poly(time,splinedf,raw=TRUE))
# plot(ribo)
# points(predict(pfit),type='line')

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

Protfun(t,T,){
  function(x) ( (exp(b1*T0) - 1) * exp(-b0 - (b1*T0) - (3*(b2^2))) * ( (1/2) *exp(2/3)*sqrt(π/3)*T*erfi( ((3*t) + 1)/sqrt(3)))*(a0 + x*(a1 + x*(a2 + a3*x))) + mu) /b1

bsfit=lm(c(100,100,1000,600,600)~bs(0:4,2))
ppoly<-RegBsplineAsPiecePoly(bsfit,"bs(0:4, 2)",shift=F)

D = 0.9999
T_n=5
a_0=ppoly[[1]]$coef[1]
a_1=ppoly[[1]]$coef[2]
a_2=ppoly[[1]]$coef[3]
a_3=ppoly[[1]]$coef[4]
T=40
t=3
mu=100


a_0 + a_1*t + a_2*t^2 + a_3*t^3

exp(-(T_n D)) (μ + integral(a_0 + a_1 t + a_2 t^2 + a_3 t^3) T exp(D t) dt)




exp(-D*T_n)*
(
  (T*(exp(D*t))* (a_3*D^3*t^3 + a_2*D^3*t^2 + a_0*D^3 - 3*a_3*D^2*t^2 + a_1*D^2*(D*t - 1) - 2*a_2*D^2*t + 6*a_3*D*t + 2*a_2*D - 6*a_3) )/
    (D^4) + mu
)



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

# lm_glm <- bsfit
# SplineTerm <- "bs(time, degree = splinedf)"

# RegBsplineAsPiecePoly <- function (lm_glm, SplineTerm, shift = FALSE){
#     if (!inherits(lm_glm, c("lm", "glm")))
#         stop("This function only works with models that inherits 'lm' or 'glm' class!")
#     tm <- terms(lm_glm)
#     tl <- attr(tm, "term.labels")
#     if (!(SplineTerm %in% tl))
#         stop(sprintf("Required SplineTerm not found! Available terms are:\n    %s",
#             toString(tl)))
#     is.bs <- grepl("bs(", SplineTerm, fixed = TRUE)
#     is.ns <- grepl("ns(", SplineTerm, fixed = TRUE)
#     if ((!is.bs) && (!is.ns))
#         stop("Provided SplineTerm is neither 'bs()' nor `ns()` from package 'splines'!")
#     pos <- match(SplineTerm, tl)
#     BsplineCoef <- unname(with(lm_glm, coefficients[assign ==
#         pos]))
#     na <- is.na(BsplineCoef)
#     if (any(na)) {
#         warning("NA coefficients found for SplineTerm; Replacing NA by 0")
#         BsplineCoef[na] <- 0
#     }
#     # allBsplineCoef <- BsplineCoef 
#     BsplineCoef[]<-0
#     lapply(1:length(BsplineCoef),function(coefn){
#       BsplineCoef[]<-0
#       BsplineCoef[coefn] <- 1
#       predvars <- attr(tm, "predvars")
#       BsplineCall <- predvars[[2L + pos]]
#       BsplineCall[[2]] <- quote(x)
#       degree <- if (is.bs)
#           BsplineCall[[3]]
#       else 3L
#       x <- with(as.list(BsplineCall[is.bs + 3:4]), c(Boundary.knots[1],
#           unname(knots), Boundary.knots[2]))
#       y <- base::c(eval(BsplineCall, data.frame(x = x)) %*% BsplineCoef)
#       n_pieces <- length(x) - 1L
#       PiecePolyCoef <- matrix(0, degree + 1L, n_pieces)
#       i <- 1L
#       while (i <= n_pieces) {
#           xg <- seq.int(x[i], x[i + 1L], length.out = degree +
#               1L)
#           yg <- base::c(eval(BsplineCall, data.frame(x = xg)) %*%
#               BsplineCoef)
#           Xg <- base::outer(xg - shift * x[i], 0:degree, "^")
#           U <- base::chol.default(base::crossprod(Xg))
#           PiecePolyCoef[, i] <- base::backsolve(U, base::forwardsolve(t.default(U),
#               base::crossprod(Xg, yg)))
#           i <- i + 1L
#       }
#       structure(list(PiecePoly = list(coef = PiecePolyCoef, shift = shift),
#           Bspline = list(call = BsplineCall, coef = BsplineCoef),
#           knots = x), class = c("PiecePoly", "RegBspline"))
#     })
# }
# ?SplinesUtils::RegBsplineAsPiecePoly
# Q
# SplinesUtils::RegBsplineAsPiecePoly(bsfit, "bs(time, degree = splinedf)")
# tmp<-RegBsplineAsPiecePoly(bsfit, "bs(time, degree = splinedf)")

# curve(-0.0606 + 0.0623 * x + 0.00173 * x ^ 2,1:35)
# curve(0.000865 - 0.00173 * x - 0.000865 * x ^ 2,1:35)

# plot(map_dbl(1:35,~ -0.0606 + 0.0623 * .x + 0.00173 * .x ^ 2))
# plot(map_dbl(1:35,~ 0.000865 - 0.00173 * .x - 0.000865 * .x ^ 2))

# plot()

# bs%>%args
# plot(bs(myx,3)[,3])
# lm(bs(1:35,3)[,3]~bs(myx,3))

# splinedf<-4

# lm(bs(myx,splinedf)[,splinedf]~bs(myx,splinedf))
# plot(lm)
# splinebasis <- bs(myx,splinedf)

# perffits <- map(1:splinedf,~ SplinesUtils::RegBsplineAsPiecePoly(lm(splinebasis[,.]~bs(myx,splinedf)),"bs(myx, splinedf)",shift=FALSE))

# #everything on row one 
# #(poly_deg,num_eq) 
# perffits[[1]][[1]]$coef

# perffits%>%length
# perffits[[1]][[1]]$coef%>%dim
# bsplinecoeffmat <- perffits%>%map(1)%>%map('coef')%>%bind_cols
# bsplinecoeffmat <- perffits%>%map(1)%>%map('coef')%>%map(rowSums)%>%bind_cols

# bsfit <- (lm(data=bsdata,ribo ~ bs(time,degree=splinedf)))


# fitpolycoefs <- as.matrix(bsplinecoeffmat) %*% matrix(bsfit$coef[-1])
# fitpolycoefs[1] <- fitpolycoefs[1] + bsfit$coef[1]

# plot(predict(bsfit))
# points(map_dbl(myx, ~fitpolycoefs[1] + (.)*fitpolycoefs[2] + (.^2)*fitpolycoefs[3] + (.^3)*fitpolycoefs[4]),type='l')



# bsfit$coef


# #BOOOM - we now have a formula for going from a polynomial described by splines to a polynomial described by the linear variable, I think
# points(map_dbl(myx, ~perffit[[1]]$coef[1] + (.)*perffit[[1]]$coef[2] + (.^2)*perffit[[1]]$coef[3] + (.^3)*perffit[[1]]$coef[4]),type='l')

# perffit[[1]]$coef


# ?SplinesUtils::RegBsplineAsPiecePoly
# Q
# SplinesUtils::RegBsplineAsPiecePoly(bsfit, "bs(time, degree = splinedf)")
# tmp<-RegBsplineAsPiecePoly(bsfit, "bs(time, degree = splinedf)")

# curve(-0.0606 + 0.0623 * x + 0.00173 * x ^ 2,1:35)
# curve(0.000865 - 0.00173 * x - 0.000865 * x ^ 2,1:35)

# plot(map_dbl(1:35,~ -0.0606 + 0.0623 * .x + 0.00173 * .x ^ 2))
# plot(map_dbl(1:35,~ 0.000865 - 0.00173 * .x - 0.000865 * .x ^ 2))

# plot()

# bs%>%args
# plot(bs(myx,3)[,3])
# lm(bs(1:35,3)[,3]~bs(myx,3))

# perffit <- SplinesUtils::RegBsplineAsPiecePoly(lm(bs(1:35,3)[,3]~bs(myx,3)),"bs(myx, 3)",shift=FALSE)
# perffit[[1]]$coef

# #BOOOM - we now have a formula for going from a polynomial described by splines to a polynomial described by the linear variable, I think
# points(map_dbl(myx, ~perffit[[1]]$coef[1] + (.)*perffit[[1]]$coef[2] + (.^2)*perffit[[1]]$coef[3] + (.^3)*perffit[[1]]$coef[4]),type='l')



# predict.PiecePoly <- function (object, newx, deriv = 0L, ...) {
#   browser()
#   ## change symbol
#   PiecePolyObject <- object
#   ## extract piecewise polynomial coefficients
#   PiecePolyCoef <- PiecePolyObject$PiecePoly$coef
#   shift <- PiecePolyObject$PiecePoly$shift
#   ## get degree
#   degree <- dim(PiecePolyCoef)[1L] - 1L
#   ## deriv validation
#   if (deriv > degree) return(numeric(length(newx)))
#   if (deriv > 2) stop("The function only computes up to 2nd derivatives!")
#   ## get power
#   power <- 0:(degree - deriv)
#   ## extract knots
#   x <- PiecePolyObject$knots
#   ## which piece?
#   piece_id <- base::findInterval(newx, x, TRUE)
#   ind <- base::split.default(seq_len(length(newx)), piece_id)
#   unique_piece_id <- as.integer(names(ind))
#   n_pieces <- length(unique_piece_id)
#   ## loop through pieces
#   y <- numeric(length(newx))
#   i <- 1L
#   while (i <= n_pieces) {
#     ii <- unique_piece_id[i]
#     xi <- newx[ind[[i]]] - shift * x[ii]
#     pc <- PiecePolyCoef[, ii]
#     if (deriv > 0) pc <- pc[-1] * seq_len(length(pc) - 1L)
#     if (deriv > 1) pc <- pc[-1] * seq_len(length(pc) - 1L)
#     y[ind[[i]]] <- base::outer(xi, power, "^") %*% pc
#     i <- i + 1L
#     }
#   y
#   }

# predict(perffits[[2]],newx=3.5)


