data {
  int<lower=0> G;          //  genes
  int<lower=0> T;          //  timepoints
  vector<lower=0>[G] MS[T];  // mass spec data mean
  vector<lower=0>[G] MS_tau[T];  // mass spec data precision
  vector<lower=0>[G] ribo[T]; // riboseq (synthesis) data mean
  vector<lower=0>[G] ribo_tau[T]; // riboseq (synthesis) data sd

}

parameters {
  vector[G] ms0logratio;  // starting mass spec relative to rTE
  vector<lower= -10,upper=0>[G] ldeg;  //amount degraded each tp - log scale
  vector [G] mRNA[T];
  vector<lower=-10,upper=10>[G] lrTE;

  //real hmu_lrTE;
  //real<lower= 0,upper=100> hsig_lrTE;
  //real hmu_ldeg;
  //real<lower= 0,upper=100> hsig_ldeg;
  
}

transformed parameters {


  vector<lower=0>[G] deg; 
  vector[G] prot[T];
  //vector<lower=0.00001>[G] lin_prot[T];
  vector [G] lin_prot[T];
  vector [G] MS0;
  vector<lower=0>[G] rTE;
  vector<lower=0.00000001>[G] lin_mRNA[T];

  
  //rTE on log scale
  rTE = exp(lrTE);
    
  //defining our starting parameter MS0 in terms of it's ratio to the production
  MS0 = lrTE + ms0logratio;
    
  //defining deg in terms of logdeg
  deg = exp(ldeg);
  

  //print(lrTE);
  //print(rTE);
  //print(ms0logratio);
  //print(ldeg);


  lin_mRNA[1]  = exp(log(2) * mRNA[1]) ;
  lin_prot[1] = exp(lrTE + ms0logratio) ; 
  prot[1] = exp(log(2) * mRNA[1]);
  
  for(t in 2:T){
    lin_mRNA[t]  = exp(log(2) * mRNA[t]) ;
    lin_prot[t] = ((lin_mRNA[t] .* rTE)) + ((lin_mRNA[t-1,]) .* (1-deg));

    prot[t] = log2(lin_prot[t]);
    // print(prot[t]);
  }
  //print(prot);
}

model{
  //priors

  ms0logratio ~ normal(0,10);

  //hmu_lrTE ~ normal(0,1000);
  //hsig_lrTE ~ normal(0,1000);
  //hmu_ldeg ~ normal(0,1000);
  //hsig_ldeg ~ normal(0,1000);

  //lrTE ~ normal(hmu_lrTE,hsig_lrTE);
  //ldeg ~ normal(hmu_ldeg,hsig_ldeg);
  

  //for(t in 1:T){//for each tp
  //    log2(mRNA[t]) ~ normal(ribo[t],  ribo_tau[t] ) ; 
  //    target += -log2(mRNA[t]);
  //    log2(prot[t]) ~ normal(MS[t],  MS_tau[t] ) ;
  //   target += -log2(prot[t]);
  //}
  for(t in 1:T){//for each tp
      mRNA[t] ~ normal(ribo[t],  ribo_tau[t] ) ; 
      prot[t] ~ normal(MS[t],  MS_tau[t] ) ;
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
