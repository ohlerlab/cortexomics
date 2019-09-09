#first attemps at a hierarch
data {
  int<lower=0> G;          //  genes
  int<lower=0> T;          //  timepoints
  vector<lower=0>[G] MS[T];  // mass spec data mean
  vector<lower=0>[G] MS_tau[T];  // mass spec data precision
  vector<lower=0>[G] ribo[T]; // riboseq (synthesis) data mean
  //vector<lower=0>[G] ribo_tau[T]; // riboseq (synthesis) data sd
}

parameters {
  vector[G] ms0logratio;  // starting mass spec relative to rTE

  vector<lower= -10,upper=0>[G] ldeg;  //amount degraded each tp - log scale
  
  vector[G] lrTE;

  real hmu_lrTE;
  real<lower= 0,upper=100> hsig_lrTE;
  real hmu_ldeg;
  real<lower= 0,upper=100> hsig_ldeg;
  
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
    prot[t] = ((ribo[t] .* rTE)) + ((prot[t-1,]) .* (1-deg));
  } 

}

model{
  #priors

  ms0logratio ~ normal(0,40);

  hmu_lrTE ~ normal(0,1000);
  hsig_lrTE ~ normal(0,1000);
  hmu_ldeg ~ normal(0,1000);
  hsig_ldeg ~ normal(0,1000);

  //lrTE ~ normal(hmu_lrTE,hsig_lrTE);
  //ldeg ~ normal(hmu_ldeg,hsig_ldeg);
  
  for(t in 1:T){//for each tp
      log2(prot[t]) ~ normal(log2(MS[t]),  MS_tau[t] ) ;
      target += -log2(prot[t]);
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
