#Runs pretty acceptably
data {
  int<lower=0> G;          //  genes
  int<lower=0> T;          //  timepoints
  int<lower=0> K;          //   replicates
  int<lower=0> Kr;          //   replicates
  vector<lower=0>[G] MS[T];  // mass spec data
  vector<lower=0>[G] MS[T];  // mass spec data
  vector<lower=0>[G] ribo[Kr,T]; // riboseq (synthesis) data
}
parameters {
  vector[G] ms0logratio;  // starting mass spec relative to rTE
  vector<lower= -10,upper=0>[G] ldeg;  //amount degraded each tp - log scale  
  vector[G] ltau;
  vector[G] lrTE;
  
}

transformed parameters {

  vector<lower=0>[G] deg; 
  vector[G] prot[T];
  vector<lower=0>[G] MS0;
  vector<lower=0>[G] rTE;
  vector<lower=0>[G] tau;

  #rTE on log scale
  rTE = exp(lrTE+19);
  
  #tau on log scale
  tau = exp(ltau);

  #defining our starting parameter MS0 in terms of it's ratio to the production
  MS0 = rTE .* exp(ms0logratio);
  
  #defining deg in terms of logdeg
  deg = exp(ldeg);
  
  prot[1] = MS0;
  
  for(t in 2:T){
    prot[t] = ((ribo[t] .* rTE)) + ((prot[t-1,]) .* (1-deg));
  } 

}

model{
  #priors
  ms0logratio ~ normal(0,40);
  ltau ~ normal(0,1000);
  lrTE ~ normal(0,1000);

  for(t in 1:T){//for each tp
    for(k in 1:K){//for each replicate
      log(MS[k,t]) ~ normal(log(prot[t]), tau ) ;//normally distributed MS signal
    }
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
#It can also be confusing whether to use matrices or arrays when writing your code. There are actually four different ways to specify a two-dimensional collection of real numbers! But which one should you pick? This largely depends on the operations you need to perform. If you need to do matrix computations, you should be using a matrix. However, if you frequently need to index into the rows of the matrix it is more efficient to use arrays. In this situation, it will save you a headache to declare an array of type row_vector than to work with matrices.
#Matrices and arrays should also be traversed in different order. Loops involving arrays should be traversed with the last index varying fastest, whereas the opposite is true for matrices. Additionally, traversing through matrices is slightly more efficient. If for example your code involves an I x J array of matrices each having dimension R x C, then the most efficient way to write a loop that traverses every element is: