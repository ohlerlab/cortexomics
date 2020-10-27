################################################################################
########fuzzy clustering with sirt fuzcluster - example
################################################################################
	
# #*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# # simulate data (2 classes and 3 items)
# set.seed(876)
# library(mvtnorm)
# Ntot <- 1000  # number of subjects
# # define SDs for simulating uncertainty
# sd_uncertain <- c( .2, 1, 2 )

# dat_m <- NULL   # data frame containing mean of membership function
# dat_s <- NULL   # data frame containing SD of membership function

# # *** Class 1
# pi_class <- .6
# Nclass <- Ntot * pi_class
# mu <- c(3,1,0)
# Sigma <- diag(3)
# # simulate data
# dat_m1 <- mvtnorm::rmvnorm( Nclass, mean=mu, sigma=Sigma )
# dat_s1 <- matrix( stats::runif( Nclass * 3 ), nrow=Nclass )
# for ( ii in 1:3){ dat_s1[,ii] <- dat_s1[,ii] * sd_uncertain[ii] }
# dat_m <- rbind( dat_m, dat_m1 )
# dat_s <- rbind( dat_s, dat_s1 )

# # *** Class 2
# pi_class <- .4
# Nclass <- Ntot * pi_class
# mu <- c(0,-2,0.4)
# Sigma <- diag(c(0.5, 2, 2 ) )
# # simulate data
# dat_m1 <- mvtnorm::rmvnorm( Nclass, mean=mu, sigma=Sigma )
# dat_s1 <- matrix( stats::runif( Nclass * 3 ), nrow=Nclass )
# for ( ii in 1:3){ dat_s1[,ii] <- dat_s1[,ii] * sd_uncertain[ii] }
# dat_m <- rbind( dat_m, dat_m1 )
# dat_s <- rbind( dat_s, dat_s1 )
# colnames(dat_s) <- colnames(dat_m) <- paste0("I", 1:3 )

# dat_s[666,] = dat_s[666,] + 100
# dat_s[66,] = dat_s[666,] + 100




timeMSeffectnm<-timeMSeffect%>%set_colnames(colnames(.)%>%str_replace('_MSdev(\\w+)','\\1_MSdev'))

effcontrasts<-cbind(alltimeeff,timeTEeffect,timeMSeffectnm)

timeff_ciddf <-	lapply(colnames(effcontrasts)%>%setNames(.,.),function(datagroup){
	message(datagroup)
	# prediction_ob$coef
	topTable(contrasts.fit(bestmscountebayes,effcontrasts[,datagroup,drop=F]),coef=1,number=Inf,confint=.95)%>%
	as.data.frame%>%rownames_to_column('uprotein_id')
})%>%bind_rows(.id='datagroup')
timeff_ciddf%<>%as_tibble
# timeff_ciddf%<>%separate(datagroup,c('assay','time'))

#time seperate the assays




fuzclustdata<-timeff_ciddf%>%filter(!str_detect(datagroup,'E13'))%>%select(uprotein_id,datagroup,logFC)%>%
	spread(datagroup,logFC)
fuzclustdata_sd<-timeff_ciddf%>%filter(!str_detect(datagroup,'E13'))%>%mutate(sd = abs(logFC - CI.L)/2 )%>%select(uprotein_id,datagroup,sd)%>%
	spread(datagroup,sd)

fuzclustdata%<>%dplyr::rename('gene'=uprotein_id)%>%{set_rownames(as.matrix(.[,-1]),.[[1]] ) }
fuzclustdata_sd%<>%dplyr::rename('gene'=uprotein_id)%>%{set_rownames(as.matrix(.[,-1]),.[[1]] ) }

# BiocManager::install(c('sirt'))

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
# estimation

#*** Model 1: Clustering with 8 random starts

res1 <- sirt::fuzcluster(K=20,fuzclustdata, fuzclustdata_sd, nstarts=8, maxiter=25)#Verify that fuzz 

res1$posterior%>%apply(1,which.max)%>%table%>%sort





###load and assemble the data
#1 limma data estimate
#2 estimates of the kinetic coordinates
#3 estimates of the deviance from those.


###Now, cluster at various cluster numbers.
res1 <- sirt::fuzcluster(K=2,dat_m, dat_s, nstarts=8, maxiter=25)#Verify that fuzz 


#Plot clusters at different numbers

#Run GO-term analysis on different clusters

#export clusters
library(tidyverse)
#data which simply decomposes to two dimensions
basicadata<-list(c(-1,-1),c(1,1),c(-.5,.5),c(.5,-.5))%>%simplify2array%>%t
basicadata_p<-list(c(-1,-1),c(1,1),c(-.5,.5),c(.5,-.5),c(1,0.7))%>%simplify2array%>%t

pcbasicadata<-princomp(basicadata)
pcbasicadata_p<-princomp(basicadata_p)

qplot(basicadata[,1],basicadata[,2])+geom_abline(slope=pcbasicadata$loadings[,1]/pcbasicadata$loadings[,2],intercept=0)
qplot(basicadata_p[,1],basicadata_p[,2])+geom_abline(slope=pcbasicadata_p$loadings[,2]/pcbasicadata_p$loadings[,1],intercept=0)

basicadataerror <- basicadata%>%{.[]<-0.000001;.}
basicadata_perror <- basicadata_p%>%{.[]<-0.25;.}

cov(basicadata_p)
var(basicadata_p[,1])

library(magrittr)
(basicadata_p[,1] - mean(basicadata_p[,1])) %>% {t(.)%*%.} %>% divide_by(nrow(basicadata_p))
var(basicadata_p[,1])
var(basicadata_p[,1]%>%rep(1000))
mean((basicadata_p[,1]%>%rep(1000) - mean(basicadata_p[,1]))^2)
#Here our 
# µt ← d-dimensional vector initialized to 0


get_ktt <- function(inputdata = basicadata,inputerror = inputdata%>%{.[]<-0;.}){
	d=ncol(inputdata)
	mu_t = rep(0,d)
	# 2 foreach t ∈ T do
	# 3 µt += t.mean()
	# 4 end
	#5 µt /= T.length()
	mu_t = colMeans(inputdata)
	#6 KTT ← d ×d matrix initialized to 0
	Ktt <- matrix(0,ncol=d,nrow=d)
	#7 foreach t ∈ T do
	i=1
	# browser()
	for(i in 1:nrow(inputdata)){
		#8 ~m ← t.mean()
		m <- inputdata[i,]
		# KTT += ~m~m T +s2·t.cov()− µtµt
		tcov = diag(inputerror[i,])
		Ktt  = Ktt +  m %*% t(m)  +   tcov  - (mu_t %*% t(mu_t)  )
	}
	#10 end
	#11 KTT /= T.length(
	Ktt = Ktt / nrow(inputdata-1)
	Ktt
}

cov(basicadata)

princomp(basicadata);princomp(basicadata)$loadings
eigen(get_ktt(basicadata))

princomp(basicadata_p%>%rep(100));princomp(basicadata_p%>%rep(100)$loadings

#resembles covariance function with zero uncertainty
eigen(get_ktt(basicadata_p%>%list%>%rep(100)%>%reduce(.f=rbind)))
eigen(cov(basicadata_p%>%list%>%rep(100)%>%reduce(.f=rbind)))


#Now can we make basicdata_p look like basicdata if we dial up error?
basicadata_perror <- basicadata_p%>%{.[]<-0.0001;.[5,]<-100;.}

eigen(get_ktt(basicadata_p))
eigen(get_ktt(basicadata_p,basicadata_perror))


#with no uncertainty, looks the same as pca
covbasdata <- cov(basicadata)
covbasdata%>%svd%>%{.$d<-sqrt(.$d);.}
svd()%>%{.$d<-sqrt(.$d);.}
cov(basicadata)
Q

cov(basicadata)
get_ktt(basicadata, basicadata%>%{.[]<-0.25;.})
#with some, we get biggger values
svd(get_ktt(basicadata, basicadata%>%{.[]<-0.5;.}))%>%{.$d<-sqrt(.$d);.}


basicadata_perror <- basicadata_p%>%{.[]<-0.0001;.[5,]<-0.0001;.}

#Now with other point
princomp(basicadata_p)
svd(get_ktt(basicadata_p, basicadata_p%>%{.[]<-1;.}))%>%{.$d<-sqrt(.$d);.}
svd(cov(basicadata_p))%>%{.$d<-sqrt(.$d);.}

princomp(basicadata_p)
princomp(basicadata_p)$loadings

eigen(get_ktt(basicadata_p))%>%{.$values<-(.$values^2);.}

eigen(get_ktt(basicadata_p,basicadata_p%>%{.[]<-0.0001;.[5,]<-100;.}))
eigen(get_ktt(basicadata_p,basicadata_p%>%{.[]<-0.0001;.[5,]<-0.0001;.}))


posctrl<-run_uapca(datawith_highcertpoint)

negctrl<-run_uapca(datawith_lowcertpain)


data.frame(d1=,d2=)