#this version does okay, but the lower diagnoanl is full of deviances
#bout to try reparamtrizing to remove colinearity between ldeg and rTE
data {
  int<lower=0> G;          //  genes
  int<lower=0> T;          //  timepoints
  int<lower=0> K;          //   replicates
  vector<lower=0>[G] MS[T,K];  // mass spec data
  vector<lower=0>[G] ribo[T]; // riboseq (synthesis) data

}
parameters {
  vector[G] ms0logratio;  // starting mass spec relative to rTE

  //real l10pi_mu; //
  //real l10pi_sig; //
  //vector[G] l10pi; 
  vector<lower=-10,upper=0>[G] ldeg; 
  
  
  real<lower=0> rTEmu;
  real<lower=0> rTEsig;
  
  real<lower=0> tau;
  vector<lower=0>[G] rTE;
  
  
  real<lower=0>ms0lr_sig;
  real<lower=0>ms0lr_mu;
}

transformed parameters {

  vector<lower=0>[G] deg; 
  vector[G] prot[T];
  vector<lower=0>[G] MS0;
 
  MS0 = rTE .* exp(ms0logratio);
  
  #defining deg in temrs of the l10half lives
  //for(g in 1:G){
  //  deg[g] = pow(10,l10pi[g]) ;
  //}
  //deg = exp( (log(0.5)./(deg*36) ));
  
  #defining deg in terms of logdeg
  deg = exp(ldeg);
  
  //MS0 as unobserved - seems to result in multimodal posterior
  //prot[1] = (1 - deg) .* MS0;
  //prot[1] = prot[1] + (ribo[1] .* rTE);
 
  prot[1,1:G] = MS0[1:G];
  
  for(t in 2:T){
    prot[t,1:G] = (ribo[t,1:G] .* rTE[1:G]);
    for(ti in 1:(t-1)){
      prot[t,1:G] += ((prot[ti,1:G]) .* exp(- ti * deg[1:G])) ;//convoluted but not sure how state works in rstan
    }
  }  
}

model{
  
  //l10pi_mu ~ normal(48,100); // weak prior on the log half life
  //l10pi_sig ~ normal(10,10); // weak prior on the log half life spread

  //l10pi ~ normal(l10pi_mu,l10pi_sig);// join distribution for the log10 half lives
  
  ms0lr_sig ~ gamma(1,1);// weak prior on the ms0lr variance
  
  ms0logratio ~ normal(rep_vector(ms0lr_mu,G),rep_vector(ms0lr_sig,G)); # join distribution of the log ratios for m0
  rTE ~ normal(rep_vector(rTEmu,G),rep_vector(rTEsig,G)); # join distribution of the rTE valuess
  
  for(t in 1:T){//for each tp
    for(k in 1:K){//for each replicate
      MS[t,k,1:G] ~ normal(prot[t,1:G], tau ) ;//normally distributed MS signal
    }
  }
  
  
}





