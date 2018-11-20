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
ms[1] = ms0
for (i in 2:length(ms)){
	ms[i] = (rTE*mean(ribo[i],ribo[i-1]))+(ms[i-1]*deg)
}


#function simulating ms given a degredation constant, a constant value for rTE, starting ms, the riboseq,
simulate_data <- function(deg,rTE,ribo,ms0){
	deg = exp(deg)
	ms = rep(NA,tps)
	ms[1] = ms0
	for (i in 2:length(ms)){
		ms[i] = ms[i-1]+ (rTE*mean(ribo[i],ribo[i-1]))-(ms[i-1]*deg) 
	}
	ms
	simdata <- data_frame(trans=ribo,time=1:tps,ms,
		ribo_1=rpois(n=seq_along(ribo),lambda=rnorm(n=seq_along(ribo),mean=ribo,sd=1)),
		ribo_2=rpois(n=seq_along(ribo),lambda=rnorm(n=seq_along(ribo),mean=ribo,sd=1)),
		ms_1=rnorm(n=seq_along(ms),mean=ms,sd=ms*0.2),
		ms_2=rnorm(n=seq_along(ms),mean=ms,sd=ms*0.2),
		ms_3=rnorm(n=seq_along(ms),mean=ms,sd=ms*0.2)
	)	
	simdata
}


#let's assume we know our ribosomal values...
log.lklh.prot <- function(data,par,tps=5){
	
	deg = exp(par[1])
	rTE = par[2]
	prot0 = par[3]

	prot = rep(NA,tps)

	prot[1] = prot0

	for (i in 2:length(prot)){
		prot[i] = prot[i-1] + (rTE*data$trans[i]) - (prot[i-1]*deg)
	}
	msmeans <- replace_na((data$ms_1+data$ms_2+data$ms_1),0) / 3

	sum(c(
		# dnorm(data$ribo_1,mean=prot,sd=data$prot*0.2,log=T),
		# dnorm(data$ribo_1,mean=prot,sd=data$prot*0.2,log=T),
		dnorm(data$ms_1,mean=prot,sd=msmeans*0.2,log=T),
		dnorm(data$ms_2,mean=prot,sd=msmeans*0.2,log=T),
		dnorm(data$ms_3,mean=prot,sd=msmeans*0.2,log=T)
	))
}



#choose a pattern
ribo = sample(ribopatterns,1)[[1]]
#begin near equilibrium
ms0 = ribo[1] * rTE * (1/deg)
#now simulate data
simulate_data(log(deg),rTE,ribo,ms0)

ggarrange(
	simdata%>%select(-trans,-ms)%>%gather(assay,value,-time)%>%separate(assay,into=c('assay','rep'))%>%
		mutate(assay=factor(assay,unique(assay)))%>%
		ggplot(data=.,aes(y=value,x=time))+facet_grid(scale='free',assay~.)+geom_point(),
	simdata%>%select(-trans,-ms)%>%gather(assay,value,-time)%>%separate(assay,into=c('assay','rep'))%>%
		mutate(assay=factor(assay,unique(assay)))%>%
		ggplot(data=.,aes(y=log2(value),x=time))+facet_grid(scale='free',assay~.)+geom_point()
)

deg=0.5
ribo = ribopatterns[[1]]
prot0 <- ribopatterns[[1]][1] * rTE * (1/deg)
simdata <- simulate_data(log(deg),rTE,ribo,prot0)

#this seems to indeed have a minimum there
log.lklh.prot(simdata,par=c(log(deg),rTE,prot0))


test1dconfints<-function(logrealdeg,realrTE,realprot0,ribo,
		datasimfunc=simulate_data,
		likfunc =log.lklh.prot 
	) {

	simdata <- datasimfunc(logrealdeg,realrTE,ribo,realprot0)
	
	#filter out things with no cor?
	ribomean <- simdata$ribo_1 +simdata$ribo_2
	msmean <- simdata$ms_1 +simdata$ms_2+simdata$ms_3 
	corcheck <- anova(lm(data=simdata,ribomean~msmean ))
	isposcor <- cor.test(ribomean,msmean)%>%{.$estimate > 0 & ( .$p.value<0.05)}

	if(!isposcor) return(NA)

	#method='L-BFGS-B'
	TEinitial <- 100

	fit = optim(par=c(log(0.5),TEinitial,ribo[1]*TEinitial),
		method='L-BFGS-B',
		control=list('fnscale'= -1),
		fn=log.lklh.prot,
		data=simdata,
		lower=c(log(0.00001),0,0),
		upper=c(log(1),1e12,1e12),
		hessian=T)

	# fit = optim(par=c(log(0.5),100,ribo[1]*100),
	# 	method='Nelder-Mead',
	# 	control=list('fnscale'= -1),
	# 	fn=log.lklh.prot,
	# 	data=simdata,
	# 	hessian=T)

	fit$par
	c(logrealdeg,realrTE,realprot0)

	log.lklh.prot(simdata,fit$par)
	log.lklh.prot(simdata,c(logrealdeg,realrTE,realprot0))

	fisher_info<-solve(-fit$hessian)

	prop_sigma<-try({ sqrt(fisher_info[1,1])})
	if(is.nan(prop_sigma)){ browser()}
	if(is(prop_sigma,'try-error')){ browser()}
	upper<-fit$par[1]+1.96*prop_sigma
	lower<-fit$par[1]-1.96*prop_sigma

	interval<-data.frame(value=fit$par[1], upper=upper, lower=lower)

	browser()
	data_frame(
		deginconf=between(logrealdeg,interval$lower[1],interval$upper[1]),
		degdist = sqrt(((logrealdeg - fit$par[1])/logrealdeg)^2),#distance actual estimate
		degup=exp(interval$upper[1]),deglow=exp(interval$lower[1]),
		logrealdeg=exp(logrealdeg),estdeg=exp(fit$par[1]),
		realrTE=realrTE,estrTE=fit$par[2],
		realprot0=realprot0,estprot0=fit$par[3],
		ribo=paste0(ribo,collapse=',')
	)
}

degrange = log(seq(0.01,0.99,by=0.05))
#test our model fitting
reps<-replicate(simp=F,100,{
		test1dconfints(
			logrealdeg=sample(degrange,1),
			# realrTE=sample(seq(0.5,3,len=20),1),
			realrTE=100,
			# realprot0=sample(seq(3e3,3e6,l=10),1),
			realprot0=100*100,
			ribo=sample(ribopatterns[1],1)[[1]])
}
)

#show as a table
reptable<-reps%>%setNames(.,seq_along(.))%>%bind_rows%>%mutate(estprot0/realprot0)%>%arrange(logrealdeg)

reptable$deginconf%>%mean


Simulating data with 

plotlist <- list()


logrealdeg <- log(0.5)
realprot0 <- 100*ribo[1]*100
realrTE <- 100
ribo<- ribopatterns[[3]]
ribo<- c(1600,800,400,200,100)


simdata <- simulate_data(logrealdeg,realrTE,ribo,realprot0)
#method='L-BFGS-B'
TEinitial <- 100

fit = optim(par=c(log(0.5),TEinitial,simdata$ms_1[1]),
	method='L-BFGS-B',
	control=list('fnscale'= -1),
	fn=log.lklh.prot,
	data=simdata,
	lower=c(log(0.00001),0,0),
	upper=c(log(1),1e12,1e12),
	hessian=T)

fitsimmeddata <- replicate(simp=F,100,simulate_data(fit$par[1],fit$par[2],ribo,fit$par[3]))
fitsimmeddata%<>%bind_rows
fitsimmeddata 

fitsimmeddata 

simdata2plot<-simdata%>%select(-trans,-ms)%>%gather(assay,value,-time)%>%separate(assay,into=c('assay','rep'))%>%
		mutate(assay=factor(assay,unique(assay)))
fitsimmeddata2plot<-fitsimmeddata%>%select(-trans,-ms)%>%gather(assay,value,-time)%>%separate(assay,into=c('assay','rep'))%>%
		mutate(assay=factor(assay,unique(assay)))

fit$par[2]
par<-fit$par%>%map_dbl(round,2)
par[1]%<>%exp%>%round(2)


fisher_info<-solve(-fit$hessian)

prop_sigma<-try({ sqrt(fisher_info[1,1])})

degupper <- exp(fit$par[1]+1.96*prop_sigma)%>%round(2)
deglower <- exp(fit$par[1]-1.96*prop_sigma)%>%round(2)

prop_sigma<-try({ sqrt(fisher_info[2,2])})
rTEupper <- (fit$par[2]+1.96*prop_sigma)%>%floor
rTElower <- (fit$par[2]-1.96*prop_sigma)%>%floor

prop_sigma<-try({ sqrt(fisher_info[3,3])})
ms0upper <- (fit$par[3]+1.96*prop_sigma)%>%floor
ms0lower <- (fit$par[3]-1.96*prop_sigma)%>%floor

simdataplottitle<-str_interp('Simulated Data: \n Actual Parameters: rTE=${realrTE}, Beta=${exp(logrealdeg)}, Ms0=${realprot0}
	Estimated rTE=${par[2]} (${rTElower} - ${rTEupper}) 
	Estimated Beta=${par[1]} (${deglower} - ${degupper})  
	Estimated Ms0=${par[3]} (${ms0lower} - ${ms0upper})')

#simulate our data log scale
svg(h=5,w=6,'../plots/simdata_degredationest_degdominantribo.svg'%T>%{normalizePath(.)%>%message})
ggpubr::ggarrange(
	simdata2plot%>%
		ggplot(data=.,aes(y=log2(value),x=time))+facet_grid(scale='free',assay~.)+geom_point(size=3)+theme_minimal()+
		geom_point(data=fitsimmeddata2plot%>%filter(assay=='ms'),position='jitter',alpha=I(0.1))+
		ggtitle(simdataplottitle)
	
)
dev.off()