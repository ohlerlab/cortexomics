data {
  int<lower=0> G;          //  genes
  int<lower=0> T;          //  timepoints
  vector[G] lMS[T];  // mass spec data mean
  vector[G] lMS_tau[T];  // mass spec data precision
  vector[G] lribo[T]; // lriboseq (synthesis) data mean
  vector[G] lribo_tau[T]; // lriboseq (synthesis) data sd
  matrix[T,T+1] mybs; // spline basis for the protein trajectory
  matrix[T,T] mydbs; // differentiated spline basis for the RNA (piecewise linear right now)
}


parameters {
  vector[G] l_st; // the steady state
  matrix<lower=0>[T,G] cM;  // vector of fold changes due to synthesis
  vector<lower=-20,upper=20>[G] l_pihalf;  //log half life
  vector[G] prot0; // initial LOG amount of protein
}

transformed parameters {
    matrix [T,G] cv; // vector of fold changes including deg
    matrix [T,G] prot; // amounts of protein
    matrix [T,G] mRNA; // amounts of mRNA
    vector[G] Kd; // the degred
    vector[G] lKd; // the log degrad
    vector[G] lKs; // the (log) synthesis constant

    lKd = log2(log(2)) -  l_pihalf;
    Kd = exp(log(2)*lKd);
    lKs = l_st -  lKd ;

    for(g in 1:G){
      cv[,g] = cM[,g] - Kd[g];
      // test
      // cv[,g] = cM[,g] - 20;
    } // get fold changes of protein

    prot = mybs[,2:T+1] * cv + rep_matrix(prot0,T)'; // get the full protein trajectory

    mRNA = prot + log2(mydbs * cM) - rep_matrix(lKs,T)'; // get the mRNA trajectory this implies
}

model{

  l_pihalf ~ normal(0,4);
  l_st ~ normal(0,10);

  for(t in 1:T){
      // cv[t] ~ normal(0,10);
      lribo[t] ~ normal(mRNA[t],lribo_tau[t]); // use confidence intervals to weight the observations
      lMS[t] ~ normal(prot[t],lMS_tau[t]); // use confidence intervals to weight the observations
  }
}

// our prior for steady state levels should be something like the spread of steady state levels we
// see in the protein, 


