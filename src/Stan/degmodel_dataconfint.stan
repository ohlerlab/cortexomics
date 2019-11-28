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
}

transformed parameters {
    matrix [T,G] cv; // vector of fold changes including deg
    matrix [T,G] prot; // amounts of protein
    matrix [T,G] mRNA; // amounts of protein

    for(g in 1:G){
      cv[,g] = cM[,g] - Kd[g];
    } // get fold changes of protein

    prot = mybs[,2:T+1] * cv + rep_matrix(prot0,T)'; // get the full protein trajectory

    mRNA = prot + log(mydbs * cM) - rep_matrix(lKs,T)'; // get the mRNA trajectory this implies

}


model{
  //priors
  //the data - can probalby dispense with this loop
  for(t in 1:T){
      lribo[t] ~ normal(mRNA[t],lribo_tau[t]);
      lMS[t] ~ normal(prot[t],lMS_tau[t]);
  }
}

// generated quantities{
//    vector [G] l2mRNA;
//    vector [G] l2prot;
//    l2mRNA = log2(mRNA);
//    l2prot = log2(prot);
// }


// loglinearsplineoptfunc <- function(parvect,invmydbs,mybs,zv,ribo,
//   prot,ribosd,protsd,returnprot=FALSE,returnribo=FALSE){
//     cat('.')
//     oparvect = parvect   
//     #ddon't exponentiate th intercept for the protein
//     parvect[-3] <- exp(parvect[-3])
//     if(any(!is.finite(parvect))) browser()
//     stopifnot(!any(is.na(parvect)))
//     stopifnot(!any(!is.finite(parvect)))

//   with(as.list(parvect),{
//     Ks = parvect[1]
//     Kd = parvect[2]
//     m0 = parvect[3]
//     cM = parvect[-c(1:3)]

//     degzv = dzv*Kd
//     #so cv are the actual possibly negative, fold change
//     #cM aree the positivee fold changes after degreedatio is removed
//     cv = cM - degzv

//     # lR(t) = log(P(t)) + log(dlP(t) + Kd(t)) - log(Ks)

//     est_log_ribo = (mybs %*% c(m0,cv)) + log(mydbs %*% cM) - log(Ks)
//     est_log_ribo
//     (mybs %*% (c(m0,cM) - c(0,degzv))) + log(mydbs %*% cM) - log(Ks)
//     est_log_protein = mybs %*% c(m0,cv)


//     if(returnribo) return(exp(est_log_ribo))
//     if(returnprot) return(exp(est_log_protein))

//     # browser()

//     ribo_LL = dnorm(log(ribo),est_log_ribo,ribosd,log=TRUE)
//     prot_LL = dnorm(log(prot),est_log_protein,protsd,log=TRUE)
//     alphapriorLL = dnorm(cM[-1],0,sd=4,log=TRUE)
//     # alphapriorLL = dnorm(cM[-1],0,sd=10000,log=TRUE)
//     alphapriorLL = 0
//     # alphapriorLL = dnorm(px[-1],0,sd=20,log=TRUE)
//     # nsconstraintLL = 
//     #   dnorm(myddbs_lims[,-1] %*% cv,0,sd=0.1,log=TRUE)+
//     #   dnorm(mydddbs_lims[,-1] %*% cv,0,sd=0.1,log=TRUE)
//     nsconstraintLL = 0

//     LL = sum(c(ribo_LL ,prot_LL,alphapriorLL,nsconstraintLL))
//     # browser()
//     if(!is.finite(LL)) browser()
//     if(length(LL)==0) browser()
//     # cat(LL);cat('...')

//     - LL
//   })
// }
// }



