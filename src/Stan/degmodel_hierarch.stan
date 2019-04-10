#first attemps at a hierarch
data {
  int<lower=0> G;          //  genes
  int<lower=0> T;          //  timepoints
  int<lower=0> K;          //   replicates
  vector<lower=0>[G] MS[K,T];  // mass spec data
  vector<lower=0>[G] ribo[T]; // riboseq (synthesis) data

}
parameters {
  vector[G] ms0logratio;  // starting mass spec relative to rTE

  vector<lower= -10,upper=0>[G] ldeg;  //amount degraded each tp - log scale
  
  vector[G] ltau;
  vector[G] lrTE;

  real hmu_lrTE;
  real hsig_lrTE;
  real hmu_ldeg;
  real hsig_ldeg;
  
}

transformed parameters {

  vector<lower=0>[G] deg; 
  vector[G] prot[T];
  vector<lower=0>[G] MS0;
  vector<lower=0>[G] rTE;
  vector<lower=0>[G] tau;

  #rTE on log scale
  rTE = exp(lrTE+19);
  
  #tau on log scale
  tau = exp(ltau);

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

  ms0logratio ~ normal(0,40);
  ltau ~ normal(0,1000);

  hmu_lrTE ~ normal(0,1000);
  hsig_lrTE ~ normal(0,1000);
  hmu_ldeg ~ normal(0,1000);
  hsig_ldeg ~ normal(0,1000);

  lrTE ~ normal(hmu_lrTE,exp(hsig_lrTE));
  ldeg ~ normal(hmu_ldeg,exp(hsig_ldeg));
  
  for(t in 1:T){//for each tp
    for(k in 1:K){//for each replicate
      log(MS[k,t]) ~ normal(log(prot[t]), tau ) ;//normally distributed MS signal
    }
  }
}

// generated quantities{
//    vector<lower=0>[G] dP[T-1];
//    vector<lower=0>[G] sP[T-1];
//    vector[G] degfact;
//    
//    for(t in 2:T){
//     sP[t-1] = ((ribo[t] .* rTE));
//     dP[t-1] = ((prot[t-1,]) .* (1-deg));
//     degfact = (1-deg);
//   } 
// }
