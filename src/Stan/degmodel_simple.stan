#Runs pretty acceptably
data {
  int<lower=0> G;          //  genes
  int<lower=0> T;          //  timepoints
  int<lower=0> K;          //   replicates
  vector<lower=0>[G] MS[K,T];  // mass spec data
  vector<lower=0>[G] ribo[T]; // riboseq (synthesis) data

}
parameters {
  vector[G] ms0logratio;  // starting mass spec relative to rTE

  vector<lower=-20,upper=0>[G] ldeg;  //amount degraded each tp - log scale
  
  real<lower=0> tau;
  vector[G] lrTE;
  
}

transformed parameters {

  vector<lower=0>[G] deg; 
  vector[G] prot[T];
  vector<lower=0>[G] MS0;
  vector<lower=0>[G] rTE;

  #rTE on log scale
  rTE = exp(lrTE);

  #defining our starting parameter MS0 in terms of it's ratio to the production
  MS0 = rTE .* exp(ms0logratio);
  
  #defining deg in terms of logdeg
  deg = exp(ldeg);
  
  prot[1] = MS0;
  
  for(t in 2:T){
    // prot[t] = ((((ribo[t]+ribo[t-1])/2) .* rTE)) + ((prot[t-1,]) .* (1-deg));
    // prot[t] = ((ribo[t] .* rTE)) + ((prot[t-1,]) .* (1-deg));
  } 

}

model{
  #priors
  ms0logratio ~ normal(0,20);
  log(tau) ~ normal(0,1000);

  for(t in 1:T){//for each tp
    for(k in 1:K){//for each replicate
      log(MS[k,t]) ~ normal(log(prot[t]), tau ) ;//normally distributed MS signal
    }
  }
}

generated quantities{
   vector<lower=0>[G] dP[T-1];
   vector<lower=0>[G] sP[T-1];
   vector[G] degfact;
   
   for(t in 2:T){
    sP[t-1] = ((ribo[t] .* rTE));
    dP[t-1] = ((prot[t-1,]) .* (1-deg));
    degfact = (1-deg);
  } 
}
