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
  
  
  real<lower=0> rTEmu; //
  real<lower=0> rTEsig; //
  
  real<lower=0> tau;
  vector<lower=0>[G] rTE;
  
  
  real<lower=0>ms0lr_sig;
  real<lower=0>ms0lr_mu;
}

transformed parameters {

  vector<lower=0>[G] deg; 
  vector[G] prot[T];
  vector<lower=0>[G] MS0;
 
  #defining our starting parameter MS0 in terms of it's ratio to the production
  MS0 = rTE .* exp(ms0logratio);
  
  #defining deg in terms of logdeg
  deg = exp(ldeg);
  
  prot[1,1:G] = MS0[1:G];
  
  for(t in 2:T){
    prot[t,1:G] = (ribo[t,1:G] .* rTE[1:G]);
    for(ti in 1:(t-1)){
      prot[t,1:G] += ((prot[ti,1:G]) .* exp(- ti * deg[1:G])) ;
    }
  }  
}

model{
  
  ms0lr_sig ~ gamma(1,1);// weak prior on the ms0lr variance
  
  deg ~ beta(0.25,0.175); //weak prior on amount degraded between timepoints.
  
  ms0logratio ~ normal(rep_vector(ms0lr_mu,G),rep_vector(ms0lr_sig,G)); # join distribution of the log ratios for m0
  rTE ~ normal(rep_vector(rTEmu,G),rep_vector(rTEsig,G)); # join distribution of the rTE valuess
  
  for(t in 1:T){//for each tp
    for(k in 1:K){//for each replicate
      MS[t,k,1:G] ~ normal(prot[t,1:G], tau ) ;//normally distributed MS signal
    }
  }
  
  
}