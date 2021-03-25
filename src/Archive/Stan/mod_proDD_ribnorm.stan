data {
  int nsamples;// number of samples
  int nribosamples;// number of samples
  int G;// number of proteins
  int T;// info on the number of conditions
  int totalmissing;// info on th total amount of missing data
  matrix[G, nsamples] lMS;// data
  matrix[G, nribosamples] lribo;// data
  matrix[G, nribosamples] voom_sigma;  
  int experimental_design[nsamples];// indicator variable matching the samples to conditions
  int experimental_design_r[nribosamples];
  real zeta[nsamples];// the spread of the dropout point for a library, gets combined with the varianc per protein
  real rho[nsamples];// rho the location of the dropout point for a given library
  real mu0;// fit by proDD the mean of means
  real sigma20;// fit by proDD - the variance in means
  real eta;// fit by proDD - th evariance in protein variances
  real nu;// fit by proDD - the mean of protein variances
  matrix[T,T+1] mybs; // spline basis for the protein trajectory
  matrix[T,T] mydbs; // differentiated spline basis for the RNA (piecewise linear right now)
}

parameters {
  real<lower=0> sigma2[G];// the variance for that protein
  vector<lower=-10,upper=10>[G] l_st; // the ratio of steady state ratio to ribo
  matrix<lower=-10,upper=10>[G,T] lcM;  // log vector of fold changes due to synthesis
  vector<lower=-20,upper=20>[G] l_pihalf;  //log half life
  vector[G] prot0; // initial LOG amount of protein
  vector[T] ribnorm; // normalization factor for the riboseq
  vector[T] protnorm;  // normalization factor for the protein
}

transformed parameters{
    matrix [G,T] cv; // vector of fold changes including deg
    matrix [G,T] prot; // amounts of protein
    matrix [G,T] mRNA; // amounts of mRNA
    vector[G] Kd; // the degred
    vector[G] lKd; // the log degrad
    vector[G] lKs; // the (log) synthesis constant
    real zetastar[totalmissing];
    matrix [G,T] cM;  // vector of fold changes due to synthesis

    lKd = log2(log(2)) -  l_pihalf;
    Kd = exp(log(2)*lKd);
    lKs = l_st -  lKd ;

    cM = exp(log(2)*lcM);

    for(g in 1:G){
      cv[g,] = cM[g,] - Kd[g];
      // test
      // cv[,g] = cM[,g] - 20;
    } // get fold changes of protein
    #rep matrix adds columns
    prot = cv * (mybs[,2:T+1]') + rep_matrix(prot0,T) - (rep_matrix(protnorm,G)'); // get the full protein trajectory

    mRNA = prot + log2(cM * (mydbs)' ) - rep_matrix(lKs,T) - (rep_matrix(ribnorm,G)'); // get the mRNA trajectory this implies

  {
    int counter = 1;
    for(i in 1:G){
      for(j in 1:nsamples){
        if(is_inf(lMS[i, j])){
          zetastar[counter] = zeta[j] * sqrt(1 + sigma2[i]/zeta[j]^2);
          counter = counter + 1;
        }
      }
    }
  }
}

model {
  // for(c in 1:ncond){
    // l_st[,c] ~ normal(mu0, sqrt(sigma20));
  #prior distribution for the variance per gene
  sigma2 ~ scaled_inv_chi_square(nu, sqrt(eta));

  #prior distributions for the steady state levels, and the 
  l_st ~ normal(0,5);
  l_pihalf ~ normal(0,4);

  #prior distribution on the fold changes per gene
  for(t in 1:T) lcM[,t] ~ normal(0,5); // put a prior on fold changes
  #broad prior distribution of the protein starting points 
  prot0 ~ normal(mu0, sqrt(sigma20*10));

  #prior disribution on the differences in protein 
  ribnorm ~ normal(0,3);
  protnorm ~ normal(0,3);
  // }


  {
    int counter = 1;
    for(i in 1:G){
      for(jr in 1:nribosamples){
        lribo[i,jr] ~ normal(mRNA[i,experimental_design_r[jr]], voom_sigma[i,jr]); // use confidence intervals to weight the observations
      }
      for(j in 1:nsamples){
        if(is_inf(lMS[i, j])){
          target += normal_lccdf(prot[i, experimental_design[j]] | rho[j], fabs(zetastar[counter]));
          counter += 1;
        }else{
          lMS[i, j] ~ normal(prot[i, experimental_design[j]], sqrt(sigma2[i]));
        }
      }
    }
  }
}
