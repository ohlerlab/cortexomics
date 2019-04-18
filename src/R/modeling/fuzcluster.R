################################################################################
########fuzzy clustering with sirt fuzcluster
################################################################################
	
#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# simulate data (2 classes and 3 items)
set.seed(876)
library(mvtnorm)
Ntot <- 1000  # number of subjects
# define SDs for simulating uncertainty
sd_uncertain <- c( .2, 1, 2 )

dat_m <- NULL   # data frame containing mean of membership function
dat_s <- NULL   # data frame containing SD of membership function

# *** Class 1
pi_class <- .6
Nclass <- Ntot * pi_class
mu <- c(3,1,0)
Sigma <- diag(3)
# simulate data
dat_m1 <- mvtnorm::rmvnorm( Nclass, mean=mu, sigma=Sigma )
dat_s1 <- matrix( stats::runif( Nclass * 3 ), nrow=Nclass )
for ( ii in 1:3){ dat_s1[,ii] <- dat_s1[,ii] * sd_uncertain[ii] }
dat_m <- rbind( dat_m, dat_m1 )
dat_s <- rbind( dat_s, dat_s1 )

# *** Class 2
pi_class <- .4
Nclass <- Ntot * pi_class
mu <- c(0,-2,0.4)
Sigma <- diag(c(0.5, 2, 2 ) )
# simulate data
dat_m1 <- mvtnorm::rmvnorm( Nclass, mean=mu, sigma=Sigma )
dat_s1 <- matrix( stats::runif( Nclass * 3 ), nrow=Nclass )
for ( ii in 1:3){ dat_s1[,ii] <- dat_s1[,ii] * sd_uncertain[ii] }
dat_m <- rbind( dat_m, dat_m1 )
dat_s <- rbind( dat_s, dat_s1 )
colnames(dat_s) <- colnames(dat_m) <- paste0("I", 1:3 )


BiocManager::install(c('sirt'))

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# estimation

#*** Model 1: Clustering with 8 random starts
res1 <- sirt::fuzcluster(K=2,dat_m, dat_s, nstarts=8, maxiter=25)
summary(res1)



#Verify that fuzz 