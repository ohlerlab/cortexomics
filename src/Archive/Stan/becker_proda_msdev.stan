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
transformed data{
  vector[G] l_pihalf; // the ratio of steady state ratio to ribo
  for(g in 1:G) l_pihalf[g] = -20;
}

parameters {
  vector<lower=-10,upper=10>[G] l_st; // the ratio of steady state ratio to ribo
  matrix<lower=-10,upper=10>[G,T] lribo;  // log vector of ribo-seq levels
  vector[G] lprot0; // initial LOG amount of protein
  matrix[G,T] msdev;//msdev
}

transformed parameters{
    matrix [G,T] ribo; // amounts of protein
    matrix [G,T] prot; // amounts of protein
    vector[G] lKd; // the degred
    vector[G] Ks; // the synthesis constant
    vector[G] m; // the slope in ribo/mRNA
    //get Kd
    lKd = log(log(2)) -  l_pihalf;
    //get Ks
    Ks = exp(l_st + lKd);
    ribo = exp(lribo);
    prot[,1] = exp(lprot0);
    // print("Ks:");
    // print(Ks);
    // print("l_st:");
    // print(l_st);
    // print("lKd:");
    // print(lKd);
    for(i in 2:T){
      // Becker's code:
      //  c = y_model[idx] - a * b / l + a * m / l ** 2
      // y_model[idx + 1] = 
      //a * b / l - 
      //a * m / l ** 2 +
      //a * m * dt / l + 
      //c * np.exp(-l * dt)
      //dt is just 1 in our model
      // we also can't do vectorized exponentiation, so we worth with lKd
      m = ribo[,i] - ribo[,i-1] ;
      prot[,i] = 
        (Ks .* ribo[,i-1])./exp(lKd) - 
        ((Ks .* m) ./ (exp(lKd*2))) + 
        ((Ks .* m)  ./ exp(lKd)) +
        ((prot[,i-1])-((Ks .*ribo[,i-1])./exp(lKd))+((Ks .*m)./(exp(lKd*2)))).*exp(-exp(lKd));
        // print((Ks .* ribo[,i-1])./exp(lKd));
        // print(((Ks .* m) ./ (exp(lKd*2))) );
        // print(((Ks .* m)  ./ exp(lKd)) );
        // print(((prot[,i-1])-((Ks .*ribo[,i-1])./exp(lKd))+((Ks .*m)./(exp(lKd*2)))).*exp(-exp(lKd)));
    }
    prot = exp(log(prot) + msdev);
}

model {
  // l_st ~ normal(0,l_st_priorsd);
  // l_pihalf ~ normal(l_pihalf_priormu,l_pihalf_priorsd);
  for(g in 1:G){
    for(t in 1:T){
      lSeqmu[g,t] ~ normal(lribo[g,t],lSeqsigma[g,t]);
      lMSmu[g,t]  ~ normal(log(prot[g,t]),lMSsigma[g,t]);
    }
  }
}
