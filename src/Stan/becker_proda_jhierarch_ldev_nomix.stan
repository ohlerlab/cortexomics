data {
  int G;// number of proteins
  int T;// info on the number of conditions
  matrix[G,T] lMSmu;
  matrix[G,T] lSeqmu;
  matrix[G,T] lMSsigma;
  matrix[G,T] lSeqsigma;
  real l_st_priorsd;
  real l_ribo_priorsd;
  real l_pihalf_priormu;
  real l_pihalf_priorsd;
}

parameters {
  real  <lower=0,upper=10> var_l_phalf;
  real  <lower=-10,upper=10> mu_l_pihalf;
  real  <lower=-10,upper=10> lKs;
  matrix<lower=-10,upper=10>[G,T] lribo;  // log vector of ribo-seq levels
  vector<lower=-20,upper=20>[G] l_pihalf;  //log half life
  vector[G] lprot0; // initial LOG amount of protein
}

transformed parameters{
    matrix [G,T] ribo; // amounts of protein
    matrix [G,T] prot; // amounts of protein
    matrix [G,T] dprot; // amounts of protein
    vector[G] lKd; // the degred
    vector[G] Ks; // the synthesis constant
    vector[G] m; // the slope in ribo/mRNA
    real sd_l_phalf;
    sd_l_phalf = sqrt(var_l_phalf);
    lKd = log(log(2)) -  l_pihalf;
    Ks = rep_vector(exp(lKs),G);
    ribo = exp(lribo);
    prot[,1] = exp(lprot0);
    // l_st = lKs - lKd;
    for(i in 2:T){
      m = ribo[,i] - ribo[,i-1] ;
      prot[,i] = 
        (Ks .* ribo[,i-1])./exp(lKd) - 
        ((Ks .* m) ./ (exp(lKd*2))) + 
        ((Ks .* m)  ./ exp(lKd)) +
        ((prot[,i-1])-((Ks .*ribo[,i-1])./exp(lKd))+((Ks .*m)./(exp(lKd*2)))).*exp(-exp(lKd));
    }
    prot = prot ;
}
model {
  var_l_phalf ~ inv_gamma(1,1);
  l_pihalf ~ normal(mu_l_pihalf,sd_l_phalf);
  mu_l_pihalf ~ normal(l_pihalf_priormu,3);
  for(g in 1:G){
    for(t in 1:T){
      lSeqmu[g,t] ~ normal(lribo[g,t],lSeqsigma[g,t]);
      }
      lMSmu[g,t]  ~ normal(log(prot[g,t]),lMSsigma[g,t]);
    }
  }
}