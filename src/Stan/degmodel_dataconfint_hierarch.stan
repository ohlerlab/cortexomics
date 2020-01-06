data {
  int<lower=0> G;          //  genes
  int<lower=0> T;          //  timepoints
  vector<lower=0>[G] lMS[T];  // mass spec data mean
  vector<lower=0>[G] lMS_tau[T];  // mass spec data precision
  vector<lower=0>[G] lribo[T]; // lriboseq (synthesis) data mean
  vector<lower=0>[G] lribo_tau[T]; // lriboseq (synthesis) data sd
  matrix[T,T+1] mybs;
  matrix[T,T] mydbs;
}


parameters {
  vector[G] lKs;
  matrix<lower=0>[T,G] cM;  // vector of fold changes due to synthesis
  vector<lower=0,upper=20>[G] Kd;  //amount degraded each tp - log scale
  vector<lower= 0>[G] prot0;

  real hmu_lrTE;
  real <lower=0>hsig_lrTE;
  real hmu_ldeg;
  real <lower=0>hsig_ldeg;
}

transformed parameters {
    matrix [T,G] cv; // vector of fold changes including deg
    matrix [T,G] prot; // amounts of protein
    matrix [T,G] mRNA; // amounts of protein

    for(g in 1:G){
      cv[,g] = cM[,g] - Kd[g];
    } // get fold changes of protein

    prot = mybs[,2:T+1] * cv + rep_matrix(prot0,T)'; // get the full protein trajectory

    mRNA = prot + log2(mydbs * cM) - rep_matrix(lKs,T)'; // get the mRNA trajectory this implies

}


model{
  //priors
  //the data - can probalby dispense with this loop

  hmu_lrTE ~ normal(0,100);
  // hsig_lrTE ~ normal(0,100);
  hsig_lrTE ~ gamma(2,1./20.);
  hmu_ldeg ~ normal(0,100);
  //hsig_ldeg ~ normal(0,100);
  hsig_lrTE ~ gamma(2,1./20.);

  lKs ~ normal(hmu_lrTE,hsig_lrTE);
  log2(Kd) ~ normal(hmu_ldeg,hsig_ldeg);
  target += -log2(Kd);

  for(t in 1:T){
      lribo[t] ~ normal(mRNA[t],lribo_tau[t]);
      lMS[t] ~ normal(prot[t],lMS_tau[t]);
  }
}
