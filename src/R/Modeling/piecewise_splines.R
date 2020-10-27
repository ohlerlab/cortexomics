
################################################################################
########Instead just use piecewise...
################################################################################
dPdt <- function(itime, state, parameters){
 dPdt <- with(as.list(c(state,parameters)) ,{
   mRNA <- exp(c(1,orthns_imat[itime,]) %*% scoefs)
   mRNA <- max(0,mRNA)
   (rTE * mRNA) - (deg*P)
 })
 list(dPdt)
}




spl_getprotdf<-function(
splfit,
splinedf =  4,
deg = 0.1,
rTE = 10,
orthns_imat,
ms0ratio = 1
 ) {
  ntps <- length(predict(splfit))

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

#So here we have a_0 and a_1 as the start and slope of snthesis respectively across a time period
P(0) = M0

P(t1) = m0 exp(- d_0 t) + { exp(a_1 t + a_0)  / ( a_1 + d_0 ) }
 
P(t2) = P(t1) exp(-d_0 (t_d2)) + { exp(a_1 t_d2 + a_0)  / ( a_1 + d_0 ) }  

#we can simplify the slope term 
if a_1 = (r[t]-r[t-1]))/t, a_0 = r[i-1]
P(t2) = P(t1) exp(-d_0 (t_d2)) + { exp( ((r[t]-r[t-1])/t_d2) t_d2 + a_0)  / ( ((r[t]-r[t-1])/t_d2) + d_0 ) }  
P(t2) = P(t1) exp(-d_0 (t_d2)) + { exp((r[t]-r[t-1] + r[i-1])  / ( () + d_0 ) }  
P(t2) = P(t1) exp(-d_0 (t_d2)) + { exp((r[t]))  / ( a_1 + d_0 ) }  



m0 = num_ode_res$P[1]
ribo = num_ode_res$ribo
Ks = 20
Kd = 0.1


# exp(- d_0 t) * integrate[  T exp(z + z_1 t) exp( d_0 t )   , {t,0,1}] 
#input that shit into wolfram alpha

#Given a vector of RNA, a startin gprotein, a Ks and a Kd, this funciton outputs the protein values
piecewise_P  <- function(ribo,m0,Ks,Kd){

  n = length(ribo)
  ribo = log(ribo)
  P = rep(0,n)
  P[1] = m0
  i=2

  for(i in 2:n){
      #assume equally spaced tps
      step_t = 1
      a_1 = (ribo[i] - ribo[i-1]) / step_t

      a_0 = (ribo[i-1])

      # produced = (Ks * exp( (a_1*step_t)+a_0 )) / (a_1 + Kd) 
      produced = (Ks * (exp(Kd+a_1) - 1) * exp(a_0 - (Kd * step_t)) ) / (a_1 + Kd)

      fromlast = (P[i-1]*exp(- (step_t * Kd)) )

      # produced + fromlast
      P[i] = produced+fromlast 
      
      # P[i] = ( (Ks * exp( (a_1*step_t)+a_0 )) / (a_1 - K_d) ) + (P[i-1] * exp())  
      
  }
  P
}



#as a vector
degredation_since 




exp(-d_0 *t) * (((exp(d_0)-1)*d_0*a_0)+(exp(d_0)*(d_0-1)+1)*a_1) / {d_0 ^ 2}

exp(a_0) / Kd



num_ode_res<-spl_getprotdf(splfit,deg = 0.1,rTE=20,splinedf=4,orthns_imat=nsbasis)

ribo <-c(rep(100,20),400,800,400,rep(100,20))
# ribo <-c(rep(100,20),100,100,100,rep(100,20))
time<-seq_along(ribo)
nsbasis <- bs(time,df=5)
dnsbasis <- dbs(time,df=5)
splfit <- (lm(log(ribo) ~ nsbasis))
exp(predict(splfit))%>%txtplot
stopifnot(!any(is.na(splfit$coef)))

microbenchmark::microbenchmark(
  spl_getprotdf(splfit,deg = 0.1,rTE=20,splinedf=4,orthns_imat=nsbasis),
  piecewise_P(num_ode_res$ribo,m0=num_ode_res$P[1],Kd=0.1,Ks=20)  
)

num_ode_res <- spl_getprotdf(splfit,deg = 0.1,rTE=1,splinedf=4,orthns_imat=nsbasis)
piecewisepres <- piecewise_P(num_ode_res$ribo,m0=num_ode_res$P[1],Kd=0.1,Ks=1)  
c

num_ode_res$ribo%>%txtplot(ylim=c(0,max(.)))
num_ode_res$P%>%txtplot(ylim=c(0,max(.)))
piecewisepres%>%txtplot(ylim=c(0,max(.)))





