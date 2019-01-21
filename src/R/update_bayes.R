
#Let's assume our MS values have a fixed variance.


#Our ebayes process would be - start with a flat prior
#Now estimate the ML decay constant and rTE given the flat priors on them
#now r-estimate the prior given the ML estimates of these.
#etc.


Dhat(t) = deg * t

Prot(t)  = exp(-dhat(t)) * (ms0 + ribo(t) + rTE(t) + exp(dhat(t)))

log(Prot(t))  = log( exp(-dhat(t)) * (ms0 + (ribo(1) + rTE + exp(dhat(1))) +  (ribo(t) + rTE(t) + exp(dhat(t))) ) ) 

# if that shit is constant then we can say that the 
	ms0*exp(-deg*t) + ribo(1)*exp(-deg*(t-1)) ... ribo(t)*exp(-t - t)  


	ms0*(deg^t) + ribo(1)*(deg^t-1) ... ribo(t)*deg^(t - t)  


	ms0*(deg^t) + ribo(1)*(deg^t-1) ... ribo(t)*deg^(t - t)  




#nls





#simplifying assumption - what if we assume a univform rTE for all genes? Then we could optimize a 