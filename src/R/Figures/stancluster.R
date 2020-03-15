# conflict_prefer("accumulate", "foreach")
# conflict_prefer("anyMissing", "matrixStats")
# conflict_prefer("as_function", "rlang")
# conflict_prefer("clusterApply", "parallel")

library(rstan)
library(rstan)
source(here::here('src/R/Rprofile.R'))
source('src/R/cortexomics_myfunctions.R')

naivebayesstan <- rstan::stan_model(model_code = '
data {
  int<lower=1> K;               // num clusters
  int<lower=1> V;               // num features
  int<lower=0> M;               // num datapoints
  // int<lower=0> N;               // total word instances
  // int<lower=1,upper=K> z[M];    // topic for doc m
  // int<lower=1,upper=V> w[N];    // word n
  real feat[M,V]; // the normally distributed features
  // int<lower=1,upper=M> doc[N];  // doc ID for word n
  // hyperparameters
  vector<lower=0>[K] alpha;     // cluster prior
  // vector<lower=0>[V] beta;      // word prior
  real pmeanmu;
  real pmeansd;
  real psdscale;
}
parameters {
  simplex[K] theta;   // clust prevalence
  real means[K,V];//means in each cluster
  real<lower=0> sds[K,V]  ; //sd of each cluster
}
//
transformed parameters{
  real gamma[M,K];  
  for (m in 1:M){
    for (k in 1:K){
    		gamma[m,k] = categorical_lpmf(k | theta);
			gamma[m,k] = gamma[m,k] +  normal_lpdf(feat[m] | means[k],sds[k]);		
	}	
  }
}
//
model {
  theta ~ dirichlet(alpha); // cluster prior,  
  	// means[k] ~ normal(beta);
 	// means[k] ~ normal();
  for(k in 1:K){
	  means[k] ~ normal(pmeanmu,pmeansd);
	  sds[k] ~ exponential(1/psdscale);
	}
  for (m in 1:M){
	target += log_sum_exp(gamma[m]);
  }
}
')

	
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# simulate data (2 classes and 3 items)
# set.seed(876)
library(mvtnorm)
# 
dat_m <-rbind(
	mvtnorm::rmvnorm( 200, mean=c(-3,-2,1), sigma=diag(.3,3) ),
	mvtnorm::rmvnorm( 200, mean=c(3,2,-1), sigma=diag(.3,3) )
)

correctclust <- c(rep(1,200),rep(2,200))

kmeansclust<-kmeans(dat_m,2)
kmeansclust$cluster%>%length

table(correctclust,kmeansclust$cluster)


extract<-rstan::extract


nbstandata <- list(
  K=K,
  V=ncol(dat_m),
  M=nrow(dat_m),
  feat=dat_m,
  alpha = rep(1,K)%>%{./sum(.)},
  psdscale = apply(dat_m,2,sd)%>%max%>%multiply_by(10),
  pmeansd = apply(dat_m,2,sd)%>%max%>%multiply_by(10),
  pmeanmu =0
)


fit <- sampling(
	naivebayesstan,
	data=nbstandata,
	iter=1000,
	chains=1
)

extract(fit)[['gamma']]%>%apply(c(2,3),mean)

sampfit <- extract(fit)[['gamma']]%>%apply(c(2,3),mean)%>%apply(1,which.max)

table(sampfit,kmeansclust$cluster)

library(tidyverse)

summary(fit)[[1]]%>%as.data.frame%>%rownames_to_column('par')%>%filter(par%>%str_detect('means|sds'))



ofit <- optimizing(
	naivebayesstan,
	data=nbstandata,
)
library(magrittr)

ofitclusts<-ofit$par%>%enframe%>%mutate(parse_stan_pars(name,c('m','k')))%>%unnest%>%filter(parameter=='gamma')%>%
	select(value,m,k)%>%spread(k,value)%>%{as.matrix(.[,2:3])}%>%apply(1,which.max)
table(correctclust,ofitclusts)

ofit$par%>%enframe%>%mutate(parse_stan_pars(name))%>%unnest%>%filter(parameter%in%c('means','sds'))




# %>%apply(c(1,3),mean)%>%apply(1,which.max)




vbfit <- vb(
	naivebayesstan,
	data=nbstandata,
	iter=1000
)

library(rstan)
library(rstan)
source(here::here('src/R/Rprofile.R'))



naivebayesstan <- rstan::stan_model(model_code = '
data {
  int<lower=1> K;               // num clusters
  int<lower=1> V;               // num features
  int<lower=0> M;               // num datapoints
  // int<lower=0> N;               // total word instances
  // int<lower=1,upper=K> z[M];    // topic for doc m
  // int<lower=1,upper=V> w[N];    // word n
  real feat[M,V]; // the normally distributed features
  // int<lower=1,upper=M> doc[N];  // doc ID for word n
  // hyperparameters
  vector<lower=0>[K] alpha;     // cluster prior
  // vector<lower=0>[V] beta;      // word prior
  real pmeanmu;
  real pmeansd;
  real psdscale;
}
parameters {
  simplex[K] theta;   // clust prevalence
  real means[K,V];//means in each cluster
  real<lower=0> sds[K,V]  ; //sd of each cluster
}
//
transformed parameters{
  real gamma[M,K];  
  for (m in 1:M){
    for (k in 1:K){
    		gamma[m,k] = categorical_lpmf(k | theta) + normal_lpdf(feat[m] | means[k],sds[k]);		
	}	
  }
}
//
model {
  theta ~ dirichlet(alpha); // cluster prior,  
  	// means[k] ~ normal(beta);
 	// means[k] ~ normal();
  for(k in 1:K){
	  means[k] ~ normal(pmeanmu,pmeansd);
	  sds[k] ~ exponential(1/psdscale);
	}
  for (m in 1:M){
	target += log_sum_exp(gamma[m]);
  }
}
')

	
nt = 5
betas = sample(1:10,nt-1)
mus = sample(1:10,nt)
sds = sample( 1:3,nt,rep=T)
ndp = 1e3
dat=matrix(0,ncol=nt,nrow=ndp)

dat[,1]<-rnorm(ndp,mus[1],sds[1])
for(i in 2:nt){
	dat[,i] <- rnorm(ndp,mus[i]+betas[i-1]*mus[i-1],sds[i])
}

colMeans(dat),cov(dat)



#https://www.johndcook.com/blog/2012/10/29/product-of-normal-pdfs/
can I just do that for the sampling?



 data(learning.test)
 # this is an in-sample prediction with naive Bayes (parameter learning
 # is performed implicitly during the prediction).
 bn = naive.bayes(learning.test, "A")
 pred = predict(bn, learning.test)
 table(pred, learning.test[, "A"])

 # this is an in-sample prediction with TAN (parameter learning is
 # performed explicitly with bn.fit).
 tan = tree.bayes(learning.test, "A")
 fitted = bn.fit(tan, learning.test, method = "bayes")
 pred = predict(fitted, learning.test)
 table(pred, learning.test[, "A"])


dat




 # this is an out-of-sample prediction, from a training test to a separate
 # test set.
 training.set = learning.test[1:4000, ]
 test.set = learning.test[4001:5000, ]
 bn = naive.bayes(training.set, "A")
 fitted = bn.fit(bn, training.set)
 pred = predict(fitted, test.set)
 table(pred, test.set[, "A"])


##Get the distance between points in our limma object, accounting for uncertainty.
 #seee https://math.stackexchange.com/questions/917292/expected-distance-between-two-vectors-that-belong-to-two-different-gaussian-dist
fit<-mscountebayes
cols2use <- fit$coef%>%colnames%>%str_subset('time')
se.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled[,cols2use]
centdist <- dist(fit$coef[,cols2use],method='euclidian',diag=TRUE)
coef_sigma <- (se.coef^2) %>% apply(1,sum)
dim(coef_sigma)
dim(as.matrix(centdist))
dist_un <- t(t(centdist + coef_sigma) + coef_sigma)

test_hdbscan = hdbscan(dist_un,minPts = )



