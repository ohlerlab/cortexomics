library(tidyverse)


#sim some data
tps = 5
riboreps = 2
msreps = 3

ribopatterns = list(
	c(100,100,300,300,200),
	c(200,200,200,200,200),
	c(300,300,100,100,200)
)

deg = 0.5

rTE = 1000

ms0 = 100

ribo = sample(ribopatterns,1)[[1]]
ms = rep(NA,tps)
ms[1] = (rTE*ribo[1])+ms0 - (ms0*deg)
for (i in 2:length(ms)){
	ms[i] = (rTE*ribo[i])+(ms[i-1]*deg)
}


simprot <- function(deg,rTE,ribo,ms0){
	ribo = sample(ribopatterns,1)[[1]]
	ms = rep(NA,tps)
	ms[1] = (rTE*ribo[1])+ms0 - (ms0*deg)
	for (i in 2:length(ms)){
		ms[i] = (rTE*ribo[i])+(ms[i-1]*deg)
	}
	ms
}


simdata <- data_frame(trans=ribo,prot=ms,time=1:tps,
	ms_1=rnorm(n=seq_along(ms),mean=ms,sd=ms*0.2),
	ms_2=rnorm(n=seq_along(ms),mean=ms,sd=ms*0.2),
	ms_3=rnorm(n=seq_along(ms),mean=ms,sd=ms*0.2),
	ribo_1=rpois(n=seq_along(ribo),lambda=rnorm(n=seq_along(ribo),mean=ribo,sd=1)),
	ribo_2=rpois(n=seq_along(ribo),lambda=rnorm(n=seq_along(ribo),mean=ribo,sd=1))
)

simdata%>%select(-trans,-prot)%>%gather(assay,value,-time)%>%separate(assay,into=c('assay','rep'))%>%
	ggplot(data=.,aes(y=value,x=time))+facet_grid(scale='free',assay~.)+geom_point()

simdata%>%select(-trans,-prot)%>%gather(assay,value,-time)%>%separate(assay,into=c('assay','rep'))%>%
	ggplot(data=.,aes(y=log2(value),x=time))+facet_grid(scale='free',assay~.)+geom_point()


simdata <-rbind(
		data_frame(value=ribo,time=seq_along(ribo))%>%mutate(assay='ribo'),
		data_frame(value=simprot(deg,rTE,ribo,ms0),time=seq_along(value))%>%mutate(assay='ms')
	)
	ggplot(data=.,aes(y=value,x=time))+facet_grid(scale='free',assay~.)+geom_point()


rbind(
		data_frame(value=ribo,time=seq_along(ribo))%>%mutate(assay='ribo'),
		data_frame(value=simprot(deg,1,ribo,1e6),time=seq_along(value))%>%mutate(assay='ms')
	)
	ggplot(data=.,aes(y=value,x=time))+facet_grid(scale='free',assay~.)+geom_point()+theme(text=element_text(size=20))



ribo = sample(ribopatterns,1)[[1]]
	
par= list(rTE=100,ms0=100)






simulate_data <- function(deg,rTE,ms0,ribo){
	deg=exp(deg)
	ms = rep(NA,tps)
	ms[1] = ((rTE*ms0)*ribo[1])+ms0 - (ms0*deg)
	for (i in 2:length(ms)){
		ms[i] = (rTE*ms0*ribo[i])+(ms[i-1]*deg)
	}
	simdata <- data_frame(trans=ribo,prot=ms,time=1:tps,
		ms_1=rnorm(n=seq_along(ms),mean=ms,sd=ms*0.2),
		ms_2=rnorm(n=seq_along(ms),mean=ms,sd=ms*0.2),
		ms_3=rnorm(n=seq_along(ms),mean=ms,sd=ms*0.2),
		ribo_1=rpois(n=seq_along(ribo),lambda=rnorm(n=seq_along(ribo),mean=ribo,sd=1)),
		ribo_2=rpois(n=seq_along(ribo),lambda=rnorm(n=seq_along(ribo),mean=ribo,sd=1))
	)
	simdata
}

#let's assume we know our ribosomal values...
log.lklh.prot <- function(data,par,prot0,tps=5){
	deg = exp(par[1])
	rTE = par[2]
	prot0 = par[3]

	prot = rep(NA,tps)

	prot[1] = (rTE*prot0*data$trans[1])+prot0 - (prot0*deg)
	for (i in 2:length(prot)){
		prot[i] = (rTE*prot0*data$trans[i])+(prot[i-1]*deg)
	}

	sum(c(
		dnorm(data$ms_1,mean=prot,sd=data$prot*0.2,log=T),
		dnorm(data$ms_2,mean=prot,sd=data$prot*0.2,log=T),
		dnorm(data$ms_3,mean=prot,sd=data$prot*0.2,log=T)
	))
}

simdata <- simulate_data(0.5,1,100,ribopatterns[[1]])

log.lklh.prot(simdata,par=c(0.75,1,100))
log.lklh.prot(simdata,par=c(0.50,1,100))
log.lklh.prot(simdata,par=c(0.25,1,100))

log.lklh.prot(simdata,par=c(0.50,2,100))
log.lklh.prot(simdata,par=c(0.50,1,100))
log.lklh.prot(simdata,par=c(0.50,.5,100))

log.lklh.prot(simdata,par=c(0.50,1,200))
log.lklh.prot(simdata,par=c(0.50,1,100))
log.lklh.prot(simdata,par=c(0.50,1,50))



simdata%>%select(-trans,-prot)%>%gather(assay,value,-time)%>%
	separate(assay,into=c('assay','rep'))%>%
	ggplot(data=.,aes(y=value,x=time))+facet_grid(scale='free',assay~.)+geom_point()

realdeg = 0.2
realrTE = 1
realprot0 = 100
ribo=ribopatterns[[1]]

# testconfints<-function(realdeg,realrTE,realprot0,ribo){

# 	simdata <- simulate_data(realdeg,realrTE,realprot0,ribo)
# 	fit = optim(par=c(0.5,1,100),method='L-BFGS-B',fn=log.lklh.prot,data=simdata,lower=c(0,0,0),upper=c(1,1e9,1e9),hessian=T)
# 	fit
# 	fisher_info<-solve(-fit$hessian)
# 	fisher_info
# 	prop_sigma<-sqrt(diag(fisher_info))
# 	prop_sigma
# 	prop_sigma<-diag(prop_sigma)
# 	upper<-fit$par+1.96*prop_sigma
# 	lower<-fit$par-1.96*prop_sigma

# 	interval<-data.frame(value=fit$par, upper=diag(upper), lower=diag(lower))
# 	c(
# 		between(realdeg,interval$lower[1],interval$upper[1]),
# 		between(realdeg,interval$lower[2],interval$upper[2]),
# 		between(realdeg,interval$lower[3],interval$upper[3])
# 	)
# }

# degrange = seq(0.01,0.99,by=0.05)

# replicate(100,testconfints(sample(degrange,1),realrTE,realprot0,ribo=ribopatterns[[3]]))%>%apply(1,mean)


test1dconfints<-function(logrealdeg,realrTE,realprot0,ribo,
	datasimfunc=simulate_data,
	likfunc =log.lklh.prot ) {
	simdata <- datasimfunc(logrealdeg,realrTE,realprot0,ribo)
	#method='L-BFGS-B'
	fit = optim(par=c(0.5,1,100),
		method='L-BFGS-B',
		control=list('fnscale'= -1),
		fn=log.lklh.prot,
		data=simdata,
		lower=c(log(0.00001),10+0.99,0),
		upper=c(log(1),10+1,1e9),
		hessian=T)

	fisher_info<-solve(-fit$hessian)
	prop_sigma<-sqrt(fisher_info[1,1])
	# if(is.nan(prop_sigma)){ browser()}
	upper<-fit$par[1]+1.96*prop_sigma
	lower<-fit$par[1]-1.96*prop_sigma


	logrealdeg
	interval<-data.frame(value=fit$par[1], upper=upper, lower=lower)

	# browser()
	data_frame(
		deginconf=between(logrealdeg,interval$lower[1],interval$upper[1]),
		degdist = sqrt(((logrealdeg - fit$par[1])/logrealdeg)^2),degup=exp(interval$upper[1]),deglow=exp(interval$lower[1]),
		logrealdeg=exp(logrealdeg),estdeg=exp(fit$par[1]),
		realrTE=realrTE,estrTE=fit$par[2],
		realprot0=realprot0,estprot0=fit$par[3],
		ribo=paste0(ribo,collapse=',')
	)

}

degrange = log(seq(0.01,0.99,by=0.05))
#test our model fitting
#currently with fixed TE - notice that prot0 just goes up to compensate
reps<-replicate(simp=F,20,test1dconfints(
	logrealdeg=sample(degrange,1),
	# realrTE=sample(seq(0.5,3,len=20),1),
	realrTE=1,
	realprot0=sample(seq(3e3,3e6,l=10),1),
	ribo=sample(ribopatterns[2],1)[[1]]
))
#show as a table
reps%>%setNames(.,seq_along(.))%>%bind_rows%>%mutate(estprot0/realprot0)%>%arrange(logrealdeg)