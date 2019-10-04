#I'd been planning to get the ode solved in there but it's not necessary with the right params

// functions {
//   real[] synth(real t,
//                   vector prot,
//                   vector rTE,
//                   vector deg,
//                   matrix lin_mRNA,
//                   int G
//                   ) {
//     vector[G] dprotdt;


//     dprotdt = (lin_mRNA[1:G,t] .* rTE[1:G] ) - ( prot .* deg ) ;
//     return dprotdt;
//   }
// }

data {
  int<lower=0> G;          //  genes
  int<lower=0> T;          //  timepoints
  vector<lower=0>[G] MS[T];  // mass spec data mean
  vector<lower=0>[G] MS_tau[T];  // mass spec data precision
  vector<lower=0>[G] ribo[T]; // riboseq (synthesis) data mean
  vector<lower=0>[G] ribo_tau[T]; // riboseq (synthesis) data sd

}

// #Given a vector of RNA, a startin gprotein, a Ks and a Kd, this funciton outputs the protein values
// piecewise_P  <- function(ribo,m0,Ks,Kd){

//   n = length(ribo)
//   ribo = log(ribo)
//   P = rep(0,n)
//   P[1] = m0
//   i=2

//   for(i in 2:n){
//       #assume equally spaced tps
//       step_t = 1
//       a_1 = (ribo[i] - ribo[i-1]) / step_t

//       a_0 = (ribo[i-1])

//       # produced = (Ks * exp( (a_1*step_t)+a_0 )) / (a_1 + Kd) 
//       produced = (Ks * (exp(Kd+a_1) - 1) * exp(a_0 - (Kd * step_t)) ) / (a_1 + Kd)

//       fromlast = (P[i-1]*exp(- (step_t * Kd)) )

//       # produced + fromlast
//       P[i] = produced+fromlast 
      
//       # P[i] = ( (Ks * exp( (a_1*step_t)+a_0 )) / (a_1 - K_d) ) + (P[i-1] * exp())  
      
//   }
//   P
// }

parameters {
  vector<lower= 0>[G] rTE;
  matrix<lower=0>[T,G] mRNA;
  vector<lower= -10,upper=0>[G] ldeg;  //amount degraded each tp - log scale
  vector<lower= 0>[G] prot0;
  //real hmu_lrTE;
  //real<lower= 0,upper=100> hsig_lrTE;
  //real hmu_ldeg;
  //real<lower= 0,upper=100> hsig_ldeg;
  
}

transformed parameters {
    vector<lower = 0, upper= 1>[G] deg;
    vector[G] = P;
    p[]
}


model{

  //priors
  real prot[T,G];
  prot = integrate_ode_rk45(synth, prot0, 1 , {1,2,3,4,5}, rTE, deg, mRNA, G);

  for(t in 1:T){//for each tp
      log2(mRNA[t]) ~ normal(ribo[t],  ribo_tau[t] ) ; 
      target += -log2(mRNA[t]);
      log2(prot[t]) ~ normal(MS[t],  MS_tau[t] ) ;
      target += -log2(prot[t]);

  }
}

generated quantities{
   vector [G] l2mRNA;
   vector [G] l2prot;
   l2mRNA = log2(mRNA);
   l2prot = log2(prot);
}





