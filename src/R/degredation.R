library(tidyverse)
library(zeallot)


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
log.lklh.prot <- function(data,par){
	
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


# Simulating data with 

plotlist <- list()


logrealdeg <- log(0.5)
realprot0 <- 100*ribo[1]*100
realrTE <- 100
ribo<- ribopatterns[[3]]
ribo<- c(1600,800,400,200,100)


simdata <- simulate_data(logrealdeg,realrTE,ribo,realprot0)

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


exprtable<-fread('exprdata/transformed_data.txt')
ribocols <- exprtable%>%colnames%>%str_subset('ribo')

sampdata<-exprtable%>%select(gene_name,one_of(ribocols))%>%
	gather(dataset,value,-gene_name)%>%
	separate(dataset,c('time','assay','rep'))%>%
	filter(gene_name==sample(gene_name,1))


	# %>%
	# ggplot(aes(x=time,y=value))+geom_point()


c(myvar1,myvar2,...myvar3) %<-% list(1,2,3,4)
myvar1==1
myvar2==2
identical(myvar3,list(3,4))



log.lklh.prot <- function(data,par){
	
	deg = exp(par[1])
	rTE = par[2]
	prot0 = par[3]
	tps = nrow(data)

	prot = rep(NA,tps)

	prot[1] = prot0

	trans <- (data$ribo1+data$ribo2)/2

	for (i in 2:length(prot)){
		prot[i] = prot[i-1] + (rTE*trans[i]) - (prot[i-1]*deg)
	}

	sum(c(
		# dnorm(data$ribo_1,mean=prot,sd=data$prot*0.2,log=T),
		# dnorm(data$ribo_1,mean=prot,sd=data$prot*0.2,log=T),
		dnorm(data$MS1,mean=prot,sd=0.1,log=T),
		dnorm(data$MS2,mean=prot,sd=0.1,log=T),
		dnorm(data$MS3,mean=prot,sd=0.1,log=T)
	))
}

prot[1] = ms0
prot[2] = ms0 * exp(deg*1) + rTE*ribo * exp(deg*0)
prot[3] = ms0 * exp(deg*2) + rTE*ribo[1,2]* exp(deg*1) + rTE*ribo[2,3* exp(deg*0)

#Now with real data....

#get our real data

#for each gene, pull out estimates of degredation and rTE

#update our prior distributions of these, weighting the observations by the inverse squared standard error (I think?)

#Now, update
exprdata<-fread('exprdata/transformed_data.txt')
mscols <- exprdata%>%colnames%>%str_subset('MS_')
ribocols <- exprdata%>%colnames%>%str_subset('ribo_')


ms_tall<-exprdata%>%
	select(one_of(mscols),gene_name)%>%
	select(1:3,gene_name)%>%
	gather(dataset,value,-gene_name)%>%
	separate(dataset,c('time','assay','rep'))%>%
	group_by(time,gene_name)

ms_meansd <- ms_tall%>%	summarise(mean=mean(na.omit(value)),sd=sd(na.omit(value)))

ms_meansd%>%filter(sd>0.5)%>%ungroup%>%arrange(gene_name)%>%as.data.frame%>%left_join(ms_tall)
#there is some higher var at the low end of the protein scores as well....
ms_meansd%>%ggplot(aes(x=sd))+geom_histogram()+coord_cartesian(xlim=c(0,0.5))
ms_meansd%>%ggplot(aes(x=mean,y=sd))+geom_point()+geom_smooth()

exprdatareshape<-exprdata%>%select(gene_name,one_of(c(ribocols,mscols)))%>%
	gather(dataset,value,-gene_name)%>%
	mutate(value = 2^value)%>%
	separate(dataset,c('time','assay','rep'))%>%
	identity%>%
	group_by(gene_name,time)%>%mutate(datname=paste0(assay,rep))%>%
	select(-assay,-rep)%>%
	spread(datname,value)


nloglikms <- function(ribo1,ribo2,MS1,MS2,MS3,ldeg=log(0.5),rTE = TEinitial,prot0 = sampdata$MS1[1] ){

	deg = exp(ldeg)	

	tps = length(ribo1)

	prot = rep(0,tps)

	prot[1] = prot0

	trans <- (ribo1+ribo2)/2


	for (i in 2:tps){
		prot[i] <- prot[i-1] + (rTE*trans[i]) - (prot[i-1]*deg)
	}
	prot <- log2(prot+1)
	-sum(c(
		dnorm(log2(MS1),mean=prot,sd=0.1,log=T),
		dnorm(log2(MS2),mean=prot,sd=0.1,log=T),
		dnorm(log2(MS3),mean=prot,sd=0.1,log=T)
	))
}

exprdatareshape$gene_name

gene_namei<-ugnames[1]

getcoefestimates<-function(gene_namei){
	sampdata <- exprdatareshape%>%ungroup%>%filter(gene_name==gene_namei)
	startpars <-list(ldeg=log(0.5),rTE=TEinitial,prot0=sampdata$MS1[1])
	fixedpars <- list(ribo1=sampdata$ribo1,ribo2=sampdata$ribo2,MS1=sampdata$MS1,MS2=sampdata$MS2,MS3=sampdata$MS3)%>%
		lapply(na.omit)

	do.call(nloglikms,c(startpars,fixedpars))
	
	out<-mle(nloglikms,
		start=startpars,
		# start=list(par=c(log(0.5),TEinitial,sampdata$MS1[1])),
		method='L-BFGS-B',
		fixed=fixedpars,
		lower=c(log(0.00001),0,0),
		upper=c(log(1),1e12,1e12))

	cint <- confint(out)
	cint[1,1:2]%<>%exp%>%round(3)
	rownames(cint)[1]%<>%str_replace('ldeg','deg')
	cint[2,1:2]%<>%round(5)
	cint[3,1:2]%<>%round(5)
	ests <- coef(out)[names(startpars)]
	ests[1] %<>% exp%>%round(5)
	cint%>%as.data.frame %>%mutate(estimate=round(ests,5))
}
library(parallel)

ugnames <- exprdatareshape$gene_name%>%unique%>%sample(1000)

tests<-mclapply(ugnames,function(i) safely(getcoefestimates)(i) )

tests%>%map('error')%>%map_lgl(is,'simpleError')
tests%>%map('error')%>%head
tests%>%map('result')%>%keep(Negate(is.null))



tests%>%map('result')%>%keep(Negate(is.null))%>%map(set_rownames,c('deg','rTE','prot0'))%>%map(rownames_to_column,'var')%>%bind_rows%>%
	filter(var=='rTE')%>%set_colnames(c('var','lower','upper','est'))%>%
	{ ggpubr::ggarrange(nrow=2,
		ggplot(.,aes(ymin=lower,ymax=upper,y=est,x=seq_along(est)))+geom_pointrange()+coord_cartesian(ylim=c(0,1)),
		ggplot(.,aes(est))+geom_histogram()+scale_x_continuous(limits=c(0,5)),
		)
	}

tests%<>%setNames(ugnames)


#plot the log distributions of rTE and rDeg.
tests%>%map('result')%>%keep(Negate(is.null))%>%map(set_rownames,c('deg','rTE','prot0'))%>%map(rownames_to_column,'var')%>%bind_rows(.id='gene_name')%>%
	filter(var!='prot0')%>%
	set_colnames(c('gene_name','var','lower','upper','est'))%>%
	select(var,est,gene_name)%>%
	spread(var,est)%>%
	filter(rTE < quantile(rTE,0.95))%>%
	{
	 # ggpubr::ggarrange(nrow=2,
		ggplot(.,aes(y=log(deg),x=log(rTE)))+geom_point()
		# )
	}

