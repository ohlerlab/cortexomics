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
  real<lower=-10,upper=10> lKs; // the ratio of steady state ratio to ribo
  matrix<lower=-10,upper=10>[G,T] lribo;  // log vector of ribo-seq levels
  vector<lower=-20,upper=20>[G] l_pihalf;  //log half life
  vector<lower=-20,upper=20>[G] lprot0; // initial LOG amount of protein
  vector<lower=-3,upper=3>[T-1] ribooffset;
  vector<lower=-3,upper=3>[T-1] protoffset;
}

transformed parameters{
    matrix [G,T] ribo; // amounts of protein
    matrix [G,T] prot; // amounts of protein
    vector[G] lKd; // the degred
    vector[G] Ks; // the synthesis constant
    vector[G] m; // the slope in ribo/mRNA
    matrix<lower=-10,upper=10>[G,T] offsetlribo;  // log vector of ribo-seq levels
    vector<lower=-3,upper=3>[T] lng_protoffset;

    lng_protoffset[1] = 0;
    lng_protoffset[2:T] = protoffset;

    offsetlribo=lribo;
    for(t in 2:T){
      // offsetlribo[,t] = offsetlribo[,t]+ribooffset[t-1]
      offsetlribo[,t] = offsetlribo[,t]+rep_vector(ribooffset[t-1],G);
    }
    //get Kd
    lKd = log(log(2)) -  l_pihalf;
    //get Ks
    Ks = rep_vector(exp(lKs),G);
    ribo = exp(offsetlribo);
    prot[,1] = exp(lprot0);
    // print("lKs:");
    // print(lKs);
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
      // prot[,i] = prot[,i] + exp(-5);
       // for(g in 1:G){
       //    for(t in 1:T){
       //         if((prot[g,t]==0)||(log(prot[g,t])==-Inf)){
       //        print((Ks .* ribo[,i-1])./exp(lKd));
       //        print(((Ks .* m) ./ (exp(lKd*2))) );
       //        print(((Ks .* m)  ./ exp(lKd)) );
       //        print(((prot[,i-1])-((Ks .*ribo[,i-1])./exp(lKd))+((Ks .*m)./(exp(lKd*2)))).*exp(-exp(lKd)));
       //          }
       //    }
       //  }

    }
    // print("lprot:");
    // print(log(prot));
}

model {
  // l_st ~ normal(0,l_st_priorsd);
  l_pihalf ~ normal(l_pihalf_priormu,l_pihalf_priorsd);
  for(g in 1:G){
    for(t in 1:T){
      lSeqmu[g,t] ~ normal(offsetlribo[g,t],lSeqsigma[g,t]);
      lMSmu[g,t]  ~ normal(log(prot[g,t])+lng_protoffset[t],lMSsigma[g,t]);
    }
  }
}