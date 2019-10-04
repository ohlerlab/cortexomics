data {
  int<lower=0> G;          //  genes
  int<lower=0> T;          //  timepoints
  vector<lower=0>[G] MS[G,T];  // mass spec data mean
  matrix<lower=0>[G] MS_tau[G,T];  // mass spec data precision
  matrix<lower=0>[G] ribo[2] [G,T]; // riboseq (synthesis) data mean
  matrix<lower=0>[G] ribo_tau[G,T]; // riboseq (synthesis) data sd
}

parameters {
  matrix<lower=0>[T,G] R;
  vector<lower= 0>[G] lK_s;
  vector<lower= -10,upper=0>[G] lK_d;  //amount degraded each tp - log scale
  vector<lower= 0>[G] prot0;  
}

transformed parameters {
    vector<lower = 0, upper= 1>[G] Kd;
    vector<lower = 0, upper= 1>[G] Ks;
    vector [G] riboslope;
    matrix [G,T] prot;
    
    Kd = exp(lK_d);
    Kd = exp(Ks);
   #see  piecewise_splines.R 

    prot[,1] = prot0;
    for(i in 2:n){
        #assume equally spaced tps
        // step_t = 1
        riboslope = (R[,i] - R[,i-1])
        
        # produced = (Ks * exp( (a_1*step_t)+a_0 )) / (a_1 + Kd) 
        // produced = (Ks * (exp(Kd+a_1) - 1) * exp(a_0 - (Kd * step_t)) ) / (a_1 + Kd)
        produced = (Ks * (exp(Kd+riboslope) - 1) * exp(R[,i-1] - (Kd)) ) / (riboslope + Kd)

        // fromlast = (P[,i-1]*exp(- (step_t * Kd)) )
        fromlast = (P[,i-1]*exp(-(Kd)) )

        prot[,i] = (Ks * (exp(Kd+riboslope) - 1) * exp(R[,i-1] - (Kd)) ) / (riboslope + Kd) + (P[,i-1]*exp(-(Kd)) )

    }
}


model{
  ribo[1][,t] ~ normal(R,  ribo_tau ) ; 
  ribo[1][,t] ~ normal(R,  ribo_tau ) ; 
  MS[t] ~ normal(prot,  MS_tau ) ;

}

// generated quantities{
//    vector [G] l2mRNA;
//    vector [G] l2prot;
//    l2mRNA = log2(mRNA);
//    l2prot = log2(prot);
// }







