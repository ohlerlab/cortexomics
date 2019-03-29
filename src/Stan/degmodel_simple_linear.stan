#Runs pretty acceptably
data {
  int<lower=0> G;          //  genes
  int<lower=0> T;          //  timepoints
  int<lower=0> K;          //   replicates
  vector<lower=0>[G] MS[K,T];  // mass spec data
  vector<lower=0>[G] ribo[T]; // riboseq (synthesis) data

}
parameters {
  real<lower=0> tau;
  vector[G] lrTE;
  
}

transformed parameters {

  vector[G] prot[T];
  vector<lower=0>[G] rTE;

  #rTE on log scale
  rTE = exp(lrTE);

  #defining deg in terms of logdeg
  
  prot[1] = ribo[1] .* rTE;
  
  for(t in 2:T){
    prot[t] = ((ribo[t] .* rTE));
  } 

}

model{
  #priors
  log(tau) ~ normal(0,1000);

  for(t in 1:T){//for each tp
    for(k in 1:K){//for each replicate
      log(MS[k,t]) ~ normal(log(prot[t]), tau ) ;//normally distributed MS signal
    }
  }
}

generated quantities{
   vector<lower=0>[G] sP[T-1];
   
   for(t in 2:T){
    sP[t-1] = ((ribo[t] .* rTE));
  } 
}
