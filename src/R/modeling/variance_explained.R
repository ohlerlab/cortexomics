#'
#'This script will explore how much of the variance in our mass spec amounts can be explained using different data and models...

#We exepct that the riboseq will explain the Mass Spec better than the RNAseq.
#We expect that the riboseq if corrected will be still better (or maybe correction necessary to get better explanatory power)
#we expect that maybe some genes exist where RNAseq is better, BUT These genes will be genes with weird, and/or low Riboseq
#We at least expect that for TE genes, the Riboseq will explain the mass spec better.


bigmod <- read_stan_csv('pipeline/stansamples/')

combine_samples = function(sflist){
  min_iter = min(unlist(lapply(sflist,function(x){attr(x,'sim')$iter})))
  chop = function(x,min_iter){
    attr(x,'sim')$iter = min_iter
    attr(x,'sim')$n_save = min_iter
    z = attr(x,'sim')$samples[[1]]
    z = z[1:min_iter,]
    attributes(z)$sampler_params = attributes(z)$sampler_params[1:min_iter,]
    attr(x,'sim')$samples[[1]] = z
    return(x)
  }
  sflist = lapply(sflist,chop,min_iter=min_iter)
  return(rstan::sflist2stanfit(sflist))
}
modelfiles<-Sys.glob('pipeline/stansamples/test_ribo_hierarchical_*.csv')
modelreads <- mclapply(modelfiles,read_stan_csv)
combined_reads <- combine_samples(modelreads)
summary <- combined_model%>%summary

"pipeline/stansamples/testhierarch_[0-9].csv"%>%Sys.glob



summary%>%hea