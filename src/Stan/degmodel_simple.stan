#Runs pretty acceptably
data {
  int<lower=0> G;          //  genes
  int<lower=0> T;          //  timepoints
  int<lower=0> K;          //   replicates
  vector<lower=0>[G] MS[T,K];  // mass spec data
  vector<lower=0>[G] ribo[T]; // riboseq (synthesis) data

}
parameters {
  vector[G] ms0logratio;  // starting mass spec relative to rTE

  vector<lower=-10,upper=0>[G] ldeg;  //amount degraded each tp - log scale
  
  real<lower=0> tau;
  vector<lower=0>[G] rTE;
  
}

transformed parameters {

  vector<lower=0>[G] deg; 
  vector[G] prot[T];
  vector<lower=0>[G] MS0;
 
  #defining our starting parameter MS0 in terms of it's ratio to the production
  MS0 = rTE .* exp(ms0logratio);
  
  #defining deg in terms of logdeg
  deg = exp(ldeg);
  
  prot[1] = MS0;
  
  for(t in 2:T){
    prot[t] = ((ribo[t] .* rTE)) + ((prot[t-1,]) .* (1-deg));
  } 

}

model{
  #priors
  ms0logratio ~ normal(0,20);
  
  for(t in 1:T){//for each tp
    for(k in 1:K){//for each replicate
      MS[t,k] ~ normal(prot[t], tau ) ;//normally distributed MS signal
    }
  }
  
    // exp(tau) ~ lognormal(0.1,10);
}
