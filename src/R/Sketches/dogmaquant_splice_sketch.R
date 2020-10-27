# DogmaQuant © 200X, Broad Institute, Inc. All rights reserved.
# Script file 5 of 8

library(abind)

mydims <- function(x) {
  sapply(dimnames(x),function(xx) xx[1:5])
}

load("fits.Rdata")
if (names(pfits[[1]])[1]=="par") {
  pfits <- t(sapply(pfits,function(x) {
    if (is.null(names(x))) {
      rep(NA,8)
    } else {
      x$par
    }
  }))
} else {
  pfits <- do.call(rbind,pfits)
}
pfits[is.infinite(pfits)] <- NA
pfits <- apply(pfits,2,function(x) {x[is.na(x)] <- mean(x,na.rm=TRUE); x})
M <- apply(pfits,2,mean,na.rm=TRUE)
V <- apply(pfits,2,var,na.rm=TRUE)
M["DC"] <- M["DL"] <- mean(M[c("DC","DL")])
M["TC"] <- M["TL"] <- mean(M[c("TC","TL")])
V["DC"] <- V["DL"] <- mean(V[c("DC","DL")])
V["TC"] <- V["TL"] <- mean(V[c("TC","TL")])
S <- sqrt(V)
inits <- pfits

load("prot.Rdata")
load("rna3ptspline.Rdata")
dimnames(prot)[[1]] <- c("M","H")

tmp <- apply(prot,1:3,function(x) sum(!is.na(x)))
tmp <- apply(tmp,3,function(x) all(x>=4))
tmp <- names(tmp[tmp])
prot <- prot[,,tmp,]

rna.unstim <- apply(p0[,c("E0","CT","C1","C2")],1,function(p) {
  names(p) <- c("E0","ET","E1","E2")
  x <- c(0,0.1,p["ET"],12,12.1)
  y <- exp(c(p["E0"],p["E0"],p["E0"]+p["E1"],p["E0"]+p["E2"],p["E0"]+p["E2"]))
  splinefun(x,y,method="natural")
})
rna.stim <- apply(p0[,c("E0","LT","L1","L2")],1,function(p) {
  names(p) <- c("E0","ET","E1","E2")
  x <- c(0,0.1,p["ET"],12,12.1)
  y <- exp(c(p["E0"],p["E0"],p["E0"]+p["E1"],p["E0"]+p["E2"],p["E0"]+p["E2"]))
  splinefun(x,y,method="natural")
})

rm(p0)

###
### Declare fitting functions
### 

timepoints <- as.numeric(dimnames(prot)[[4]])
names(timepoints) <- dimnames(prot)[[4]]

load("Kfunc.Rdata")


#These rate functions are for dealing with changes in the rates
#over time, applicable to both 
# 12 is because those chnges in raates are over a 12 hr period
# e.g. i f initial rate is 1, and final is 12, then you can
# do ratefunc(c(log(1),log(24)),t) to get the rate at time t



ratefunc <- function(par,x) {
  Y0 <- exp(par[1])
  Y12 <- exp(par[1])*exp(par[2]) # so second arg is like percent change
  (Y12-Y0)/12*x+Y0
}
library(testthat)

#In the case where parameter is a constant, literally jsut
#exponentiates the answer

expect_true(0.5 == ratefunc(c(log(0.5),log(1)),1) ) 
expect_true(0.5 == ratefunc(c(log(0.5),log(1)),12) ) 




#same this is the integral of the above
Int.ratefunc <- function(par,x) {
  Y0 <- exp(par[1])
  Y12 <- exp(par[1])*exp(par[2])
  (Y12-Y0)/12*0.5*x^2+Y0*x # This is the integral
}
#So in the case wheree there's no change over the 12 tps
#This is literalaly just returning x * exp(par[1])
Int.ratefunc(c(log(1.5),log(1)),8)

Int.ratefunc()

#this looks okay in wolfram : 
#integrate( exp(BSplineBasis[4, t ] ) * T * exp(exp(d) * t )  ) dt

#The T factors out
#T integral e^(N_4(t) + e^d t) dt

#We get this definite integral
T ( (-1 + E^(T1 (E^d + Subscript[N, 4])))/(E^d + Subscript[N, 4]) )


#expand? No, this doesn't integrate:
integrate( exp(BSplineBasis[4,1,t ] + BSplineBasis[4,2, t ] ) * T * exp(exp(d) * t )  ) dt

#This does , but 
integrate( ((a1*BSplineBasis[4,0,t ]) + (a2*BSplineBasis[4,1, t ]) ) * T * exp(exp(d) * t ) dt ) from 0 to T1 

Integrate[
(a0 BSplineBasis[4, 0, t] + a1 BSplineBasis[4, 1, t] + a4 BSplineBasis[4, 4, t]) T Exp[Exp[d] t], {t}]


#But what if we want to first apply 
#The original complicated version - allows for changing raates
Mfunc <- function(tp,par,rfunc,kpar) {
  Di=par[1:2]
  Ti=par[3:4]
  M0 <- exp(par[5])
  grass <- exp(par[6])
  temp <- try({
  exp(-Int.ratefunc(Di,tp)) * (M0 +
  integrate(function(tt) {
      rfunc(tt)*ratefunc(Ti,tt)*Kfunc(tt,kpar)*exp(Int.ratefunc(Di,tt))
  },
    lower=0,upper=tp)$value)+grass
  },silent=TRUE)
  if (class(temp)=="try-error") NA else as.vector(temp)
}








#We'll do our simple version, no chang ein rates
#Rate func and Int.ratefunc now just become exp and x*exp
Mfunc <- function(tp,par,rfunc,kpar) {
  Di=par[1:2]
  Ti=par[3:4]
  M0 <- exp(par[5])
  grass <- exp(par[6])
  temp <- try({
  exp(-Int.ratefunc(Di,tp)) * (M0 +
  integrate(function(tt) {
      rfunc(tt)*T*exp(exp(logD)*t)
      rfunc(tt)*ratefunc(Ti,tt)*exp(Int.ratefunc(Di,tt))
  },
    lower=0,upper=tp)$value)+grass
  },silent=TRUE)
  if (class(temp)=="try-error") NA else as.vector(temp)
}

rfunc = function(tt) 1

par <-  c()
Mfunc(1)



Hfunc <- function(tp,par,rfunc,kpar) {
  Di=par[1:2]
  Ti=par[3:4]
  M0 <- exp(par[5])
  grass <- exp(par[6])
  temp <- try({
  exp(-Int.ratefunc(Di,tp)) * (0 +
  integrate(function(tt) {
    rfunc(tt)*ratefunc(Ti,tt)*(1-Kfunc(tt,kpar))*exp(Int.ratefunc(Di,tt))},
    lower=0,upper=tp)$value)+grass
  },silent=TRUE)
  if (class(temp)=="try-error") NA else as.vector(temp)
}

###
### Declare optimizing functions: sqe(), which calls getfitted()
###

getfitted <- function(tp,par,rfunc,kpar) {
  # par <- par[c("D0","DC","T0","TC","M0","grass")]
  # tp <- 6
  Hhat <- Hfunc(tp,par,rfunc,kpar)
  Mhat <- Mfunc(tp,par,rfunc,kpar)
  c("M"=as.vector(Mhat),"H"=as.vector(Hhat))
}

llike <- function(x,y,noise) if (any(is.na(x) | is.infinite(x))) NA else sum(dnorm(y,mean=x,sd=noise,log=TRUE),na.rm=TRUE)

priorpen <- function(par,M,S) {
  nn <- c("D0","DC","DL","T0","TC","TL","M0","grass")
  sum(dnorm(par[nn],mean=M[nn],sd=S[nn],log=TRUE))
}

sqe <- function(par,from,to,PROT.UNSTIM,PROT.STIM,RNA.UNSTIM,RNA.STIM,KPAR.UNSTIM,KPAR.STIM,M,S,noiseM,noiseH) {

  #cat(paste(par,collapse=" "),"\n")
  UNSTIM.fitted <- t(sapply(timepoints,getfitted,par[c("D0","DC","T0","TC","M0","grass")],RNA.UNSTIM,KPAR.UNSTIM))
  STIM.fitted <-  t(sapply(timepoints,getfitted,par[c("D0","DL","T0","TL","M0","grass")],RNA.STIM,KPAR.STIM))

  ym <- c(PROT.UNSTIM[,"M"],PROT.STIM[,"M"])
  xm <- c(UNSTIM.fitted[,"M"],STIM.fitted[,"M"])
  yh <- c(PROT.UNSTIM[,"H"],PROT.STIM[,"H"])
  xh <- c(UNSTIM.fitted[,"H"],STIM.fitted[,"H"])
  #plot(timepoints,PROT.STIM[,"H"])
  #lines(timepoints,STIM.fitted[,"H"])
  
  temp <- -llike(xm,ym,noiseM)-llike(xh,yh,noiseH)-priorpen(par,M,S)
  #temp <- -as.vector(logLik(lm(ym ~ xm))+logLik(lm(yh ~ xh)))
  #cat(signif(temp,4),"\n")
  temp
}

llratiotest <- function(LL.null,LL.alt,df) {
  D <- -2*LL.null+2*LL.alt
  pchisq(D,df=df,lower.tail=FALSE)
}

getparams <- function(PROT.UNSTIM,PROT.STIM,RNA.UNSTIM,RNA.STIM,KPAR.UNSTIM,KPAR.STIM,M,S,inits) {
  #PROT.UNSTIM <- DAT["UNSTIM","R1",3,,]
  #PROT.STIM <- DAT["STIM","R1",3,,]
  #RNA.UNSTIM <- RNA[3,"UNSTIM",]
  #RNA.STIM <- RNA[3,"UNSTIM",]
  #fits <- try(optim("par"=inits,sqe,"PROT.UNSTIM"=PROT.UNSTIM,"PROT.STIM"=PROT.STIM,"RNA.UNSTIM"=RNA.UNSTIM,"RNA.STIM"=RNA.STIM,"KPAR.UNSTIM"=KPAR.UNSTIM,"KPAR.STIM"=KPAR.STIM,"M"=M,"S"=S,method="BFGS"))
  tt <- timepoints
  MUNSTIM <- PROT.UNSTIM[,"M"]
  MSTIM  <-  PROT.STIM[,"M"]
  noiseM <- (summary(lm(MUNSTIM ~ tt + I(tt^2)))$sigma+summary(lm(MSTIM ~ tt + I(tt^2)))$sigma)/2
  HUNSTIM <- PROT.UNSTIM[,"H"]
  HSTIM  <-  PROT.STIM[,"H"]
  noiseH <- (summary(lm(HUNSTIM ~ tt + I(tt^2)))$sigma+summary(lm(HSTIM ~ tt + I(tt^2)))$sigma)/2
  fits <- try(optim("par"=inits,sqe,"PROT.UNSTIM"=PROT.UNSTIM,"PROT.STIM"=PROT.STIM,"RNA.UNSTIM"=RNA.UNSTIM,"RNA.STIM"=RNA.STIM,"KPAR.UNSTIM"=KPAR.UNSTIM,"KPAR.STIM"=KPAR.STIM,"M"=M,"S"=S,"noiseM"=noiseM,"noiseH"=noiseH,method="BFGS",hessian=TRUE))
  show(fits)
  fits
}

genes <- dimnames(prot)[[3]]
names(genes) <- genes

Sys.time()
pfits <- 
lapply(genes,function(gn) {
  cat("\n",gn,"\n")
  getparams("PROT.UNSTIM"=t(prot[,"UNSTIM",gn,]),
            "PROT.STIM"=t(prot[,"STIM",gn,]),
            "RNA.UNSTIM"=rna.unstim[[gn]],
            "RNA.STIM"=rna.stim[[gn]],
            "KPAR.UNSTIM"=Kpars[,"UNSTIM"],
            "KPAR.STIM"=Kpars[,"STIM"],
            "M"=M,
            "S"=S,
            "inits"=inits[gn,])
})
Sys.time()

save.image(file="fits.Rdata")
