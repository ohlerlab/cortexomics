

#We are going to play with using gaussian quadrature to approximate a function as per
#https://academic.oup.com/bioinformatics/article/32/17/i702/2450770
#Pham and Jimenez

#the likelihood of our data given a set of parameters and predictors
`p(y|x;theta)` <- function(x,theta){
    #can be parametrized as the likelihood of u given those parameters and data
    #times the likelihood of the data given u, for all u.
    #so u is the biological variable pointed to by model theta and predictors x
    #via which those two are linked tot y.
    #integrating over u is the tricky bit.
    for( u in u)  `p(y|u)`(y,u) `p(u|x;theta)`(u,x,theta)
}

#tehcnical variabilitygiven a concentration
# tvar = `p(y|u)`(y,u)

#predicted concentration given predictors and model
# bvar = `p(u|x;theta)`(u,x,theta)

#a closed form solution to this integral exists if both are normal.

#Here, we'll instead assume that we dont' have a closed form, and that bvar is normal,
#and use gaussian quadrature. to approximate the integral
# gauss<-read_tsv('~/Dropbox/projects/cortexomics/gaussian_abscissas.tsv',col_names = FALSE)
# gauss%>%arrange(rev(X1))%>%.[-1,]%>%mutate(X1 = seq(nrow(gauss)+1,by=1,len=n()))%>%
#   {bind_rows(gauss,.)}%>%as.data.frame%>%
  write_colnames(c('m','z_m','wstar_m'))
write_tsv('~/Dropbox/projects/cortexomics/gaussian_abscissas_sym.tsv')
gauss<-read_tsv('~/Dropbox/projects/cortexomics/gaussian_abscissas_sym.tsv',col_names = FALSE)

lquad <- function(y,x,theta){
  for(m in 1:M){
    z_m <-
    u_m <- u_hat+sigma_hat*z_m
    lambda_m <-tvar(y,u_hat + (sqrt(2)*sigma_hat*z_m) )
    bvar(u_m,x,theta)
  }
}

Lquad<-function(Y,X,theta){
  log(theta_m) - (0.5) (2*pi*sigsq) - (Um - x_/(2*sigsq))
}
#

#a closed form solution to 
i=1:3
for (i in i){
  cat(i)
}


