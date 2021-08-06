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
  // real  <lower=-10,upper=10> mu_lks;
  // real  <lower=0.01,upper=10> sd_lks;
  real  <lower=0,upper=10> var_l_phalf;
  real  <lower=-10,upper=10> mu_l_pihalf;
  real<lower=-20,upper=20> lKs; // the ratio of steady state ratio to ribo
  //vector<lower=-10,upper=10>[G] l_st; // the ratio of steady state ratio to ribo
  matrix<lower=-10,upper=10>[G,T] lribo;  // log vector of ribo-seq levels
  vector<lower=-20,upper=20>[G] l_pihalf;  //log half life
  vector[G] lprot0; // initial LOG amount of protein
  matrix[G,T] msdev;//msdev
  vector<lower = 0, upper = 1>[G] theta; // prob it's deviant - responsibility parameter

}

transformed parameters{
//     matrix [G,T] ribo; // amounts of protein
//     matrix [G,T] prot; // amounts of protein
//     vector[G] lKd; // the degred
//     vector[G] Ks; // the synthesis constant
//     vector[G] m; // the slope in ribo/mRNA
//     // vector[G] lKs; // the slope in ribo/mRNA
    real sd_l_phalf;
//     //get Kd
//     lKd = log(log(2)) -  l_pihalf;
//     //get Ks
//     //lKs = l_st + lKd;
//     //Ks = rep_vector(exp(lKs),G);
//     //Ks = rep_vector(exp(lKs),G);
//     ribo = exp(lribo);
//     prot[,1] = exp(lprot0);
//     // print("Ks:");
//     // print(Ks);
//     // print("l_st:");
//     // print(l_st);
//     // print("lKd:");
//     // print(lKd);
//     for(i in 2:T){
//       // Becker's code:
//       //  c = y_model[idx] - a * b / l + a * m / l ** 2
//       // y_model[idx + 1] = 
//       //a * b / l - 
//       //a * m / l ** 2 +
//       //a * m * dt / l + 
//       //c * np.exp(-l * dt)
//       //dt is just 1 in our model
//       // we also can't do vectorized exponentiation, so we worth with lKd
//       m = ribo[,i] - ribo[,i-1] ;
//       prot[,i] = 
//         (Ks .* ribo[,i-1])./exp(lKd) - 
//         ((Ks .* m) ./ (exp(lKd*2))) + 
//         ((Ks .* m)  ./ exp(lKd)) +
//         ((prot[,i-1])-((Ks .*ribo[,i-1])./exp(lKd))+((Ks .*m)./(exp(lKd*2)))).*exp(-exp(lKd));
//         // print((Ks .* ribo[,i-1])./exp(lKd));
//         // print(((Ks .* m) ./ (exp(lKd*2))) );
//         // print(((Ks .* m)  ./ exp(lKd)) );
//         // print(((prot[,i-1])-((Ks .*ribo[,i-1])./exp(lKd))+((Ks .*m)./(exp(lKd*2)))).*exp(-exp(lKd)));
//     }
// }
    matrix [G,T] ribo; // amounts of protein
    matrix [G,T] prot; // amounts of protein
    matrix [G,T] dprot; // amounts of protein
    vector[G] lKd; // the degred
    vector[G] Ks; // the synthesis constant
    vector[G] m; // the slope in ribo/mRNA
    sd_l_phalf = sqrt(var_l_phalf);
    lKd = log(log(2)) -  l_pihalf;
    Ks = rep_vector(exp(lKs),G);
    ribo = exp(lribo);
    prot[,1] = exp(lprot0);
    for(i in 2:T){
      m = ribo[,i] - ribo[,i-1] ;
      prot[,i] = 
        (Ks .* ribo[,i-1])./exp(lKd) - 
        ((Ks .* m) ./ (exp(lKd*2))) + 
        ((Ks .* m)  ./ exp(lKd)) +
        ((prot[,i-1])-((Ks .*ribo[,i-1])./exp(lKd))+((Ks .*m)./(exp(lKd*2)))).*exp(-exp(lKd));
    }
    dprot = exp(log(prot) + msdev);
}
model {
  // l_st ~ normal(0,l_st_priorsd);
  // l_pihalf ~ normal(l_pihalf_priormu,l_pihalf_priorsd);
  // var_l_phalf ~ inv_gamma(0.001,0.001);
  // sd_l_phalf ~ normal(mu_l_pihalf,3)
  // lKs ~ normal(mu_lks,sd_lks);
  l_pihalf ~ normal(mu_l_pihalf,sd_l_phalf);
  // l_pihalf ~ normal(mu_l_pihalf,1);
  // mu_l_pihalf ~ normal(l_pihalf_priormu,3);
  // l_pihalf ~ normal(mu_l_pihalf,2);
  // theta ~ beta(1,5)
  // msdev ~ normal(0,1)
  for(g in 1:G){
    for(t in 1:T){
      lSeqmu[g,t] ~ normal(lribo[g,t],lSeqsigma[g,t]);
      // lMSmu[g,t]  ~ normal(log(prot[g,t]),lMSsigma[g,t]);
      if(dprot[g,t] <= 0){
        print("problem with gene:");
        print(g);
        print(" dprot is:");
        print(dprot[g,]);
        print("lprot is:");
        print(log(prot[g,]));
        print(" msdev is:");
        print(msdev[g,]);
        print(" lMSmu is:");
        print(lMSmu[g,]);
        print(" lribo is:");
        print(lribo[g,]);
        print(" l_pihalf is:");
        print(l_pihalf[g]);
        print(" theta is:");
        print(theta[g]);
      }
      target += log_sum_exp(
            log(1-theta[g]) + normal_lpdf(lMSmu[g,t] | log(prot[g,t]),lMSsigma[g,t]),
            log(theta[g]) + normal_lpdf(lMSmu[g,t] | log(dprot[g,t]),lMSsigma[g,t])
      );
    }
  }
}
// so with free mu_l_pihalf and forced variance of 1 this works very well
// try prio on mu_l_pihalf still works well

