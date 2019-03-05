
library(purrr)
library(ggExtra)


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

optimwrap <- function (minuslogl, start = formals(minuslogl), method = "BFGS",
    fixed = list(), nobs, ...)
{
    call <- match.call()
    n <- names(fixed)
    fullcoef <- formals(minuslogl)
    if (any(!n %in% names(fullcoef)))
        stop("some named arguments in 'fixed' are not arguments to the supplied log-likelihood")
    fullcoef[n] <- fixed
    if (!missing(start) && (!is.list(start) || is.null(names(start))))
        stop("'start' must be a named list")
    start[n] <- NULL
    start <- sapply(start, eval.parent)
    nm <- names(start)
    oo <- match(nm, names(fullcoef))
    if (anyNA(oo))
        stop("some named arguments in 'start' are not arguments to the supplied log-likelihood")
    start <- start[order(oo)]
    nm <- names(start)
    f <- function(p) {
        l <- as.list(p)
        names(l) <- nm
        l[n] <- fixed
        do.call("minuslogl", l)
    }
    oout <- if (length(start))
        optim(start, f, method = method, hessian = TRUE, ...)
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
	simdata <- data.frame(trans=ribo,time=1:tps,ms,
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

	data.frame(
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
reptable<-reps%>%setNames(.,seq_along(.))%>%keep(~ ! any(is.na(.)))%>%bind_rows%>%mutate(estprot0/realprot0)%>%arrange(logrealdeg)

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


exprtable<-fread('./exprdata/transformed_data.txt')
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

# prot[1] = ms0
# prot[2] = ms0 * exp(deg*1) + rTE*ribo * exp(deg*0)
# prot[3] = ms0 * exp(deg*2) + rTE*ribo[1,2]* exp(deg*1) + rTE*ribo[2,3* exp(deg*0)

#Now with real data....

#get our real data

#for each gene, pull out estimates of degredation and rTE

#update our prior distributions of these, weighting the observations by the inverse squared standard error (I think?)

#Now, update 
exprdata_all<-fread('./exprdata/transformed_data.txt')
mscols <- exprdata_all%>%colnames%>%str_subset('MS_')
ribocols <- exprdata_all%>%colnames%>%str_subset('ribo_')


ms_tall<-exprdata_all%>%
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

exprdatareshape<-exprdata_all%>%select(gene_name,one_of(c(ribocols,mscols)))%>%
	gather(dataset,value,-gene_name)%>%
	mutate(value = 2^value)%>% ####NOTE I'm de-logging it here... probably not good
	separate(dataset,c('time','assay','rep'))%>%
	identity%>%
	group_by(gene_name,time)%>%mutate(datname=paste0(assay,rep))%>%
	select(-assay,-rep)%>%
	spread(datname,value)

'Satb2' %in% exprdata_all$gene_name
'Satb2' %in% exprdatareshape$gene_name


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

ugnames <- exprdatareshape$gene_id

getcoefestimates<-function(gene_namei){
	sampdata <-  sampdata%>%ungroup%>%filter(gene_name==gene_namei)
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

ugnames <- exprdatareshape$gene_name%>%unique%>%sample(1000)%>%c('Satb2')%>%unique

tests<-mclapply(ugnames,function(i) safely(getcoefestimates)(i) )
tests%<>%setNames(ugnames)

tests%>%map('error')%>%map_lgl(is,'simpleError')%>%table
tests%>%map('result')%>%keep(Negate(is.null))

#Plot the distribution of parameter estimates
tests%>%map('result')%>%
	keep(Negate(is.null))%>%
	map(set_rownames,c('deg','rTE','prot0'))%>%
	map(rownames_to_column,'var')%>%
	bind_rows%>%
	filter(var=='rTE')%>%
	set_colnames(c('var','lower','upper','est'))%>%
	{ ggpubr::ggarrange(nrow=2,
		ggplot(.,aes(ymin=lower,ymax=upper,y=est,x=seq_along(est)))+geom_pointrange()+coord_cartesian(ylim=c(0,1)),
		ggplot(.,aes(est))+geom_histogram()+scale_x_continuous(limits=c(0,5)),
		)
	}



'Satb2' %in% ugnames
tests$Satb2

#plot the log distributions of rTE and rDeg.
indiv_ests <- tests%>%map('result')%>%
	keep(Negate(is.null))%>%
	map(set_rownames,c('deg','rTE','prot0'))%>%
	map(rownames_to_column,'var')%>%bind_rows(.id='gene_name')%>%
	filter(var!='prot0')%>%
	set_colnames(c('gene_name','var','lower','upper','est'))


 pdfexpr<-function(file,expr,...){
    dir.create(dirname(file))
    pdf(file,...)
     print(expr)
     dev.off()
     message(normalizePath(file))
 }
 

#  #plot distribution of estimated 
#  pdfexpr('../plots/modelling/deg_rTE.pdf',{

# 	indiv_ests%>%
# 	select(var,est,gene_name)%>%
# 	spread(var,est)%>%
# 	filter(rTE < quantile(rTE,0.95))%>%
# 	{
# 		ggMarginal(ggplot(.,aes(y=log(deg),x=log(rTE)))+geom_point(),type='histogram')
# 	}
# })



gene_name_i = 'Satb2'
'Satb2' %in% indiv_ests$gene_name

 #plot distribution of estimated 
 pdfexpr('../plots/modelling/deg_rTE.pdf',w=14,h=7,expr={
 	ggarrange(nrow=1,plotlist=list(
		indiv_ests%>%
			select(var,est,gene_name)%>%
			spread(var,est)%>%
			filter(rTE < quantile(rTE,0.95))%>%
			# {assert_that(gene_name_i %in% .$gene_name);.}%>%
			mutate(Gene_Name=fct_other(gene_name,gene_name_i))%>%
			{
				ggMarginal(ggplot(.,aes(y=log(deg),x=log(rTE),color=Gene_Name))+geom_point(),type='histogram')
			},

	exprdatareshape%>%
		set_colnames(c('gene_name','time','ms_1','ms_2','ms_3','ribo_1','ribo_2'))%>%
		filter(gene_name==gene_name_i)%>%
		gather(assay,value,-time,-gene_name)%>%
		separate(assay,into=c('assay','rep'))%>%
		mutate(assay=factor(assay,unique(assay)))%>%
		ggplot(data=.,aes(y=value,x=time))+facet_grid(scale='free',assay~.)+geom_point()
		)
 	)
})

pdf('../plots/modelling/deg_rTE.pdf',w=14,h=7)
 	ggarrange(nrow=1,ncol=2,plotlist=list(
		indiv_ests%>%
			select(var,est,gene_name)%>%
			spread(var,est)%>%
			filter(rTE < quantile(rTE,0.95))%>%
			# {assert_that(gene_name_i %in% .$gene_name);.}%>%
			mutate(Gene_Name=fct_other(gene_name,gene_name_i))%>%
			{
				ggMarginal(ggplot(.,aes(y=log(deg),x=log(rTE),color=Gene_Name))+geom_point(),type='histogram')
			},

	exprdatareshape%>%
		set_colnames(c('gene_name','time','ms_1','ms_2','ms_3','ribo_1','ribo_2'))%>%
		filter(gene_name==gene_name_i)%>%
		gather(assay,value,-time,-gene_name)%>%
		separate(assay,into=c('assay','rep'))%>%
		mutate(assay=factor(assay,unique(assay)))%>%
		ggplot(data=.,aes(y=value,x=time))+facet_grid(scale='free',assay~.)+geom_point()
		)
 	)
dev.off()




pdfexpr('../plots/modelling/mass_spec_sdplot.pdf',{

	exprdata%>%ungroup%>%
		select(one_of(mscols),gene_name,time)%>%
		gather(dset,signal,-time,-gene_name)%>%
		group_by(time)%>%

	log2(exprdata[,mscols])


})

# install.packages('ggalt')
library(ggalt)

##I need to work out what the proper SD is for the mass spec data, on the log 2 scale.....

#Load our LFQ values, ge the log2 mean and log2 sd for the complete cases
sdplotdf<-'ms_tables/ms_LFQ_total_ms_tall.tsv'%>%
	read_tsv%>%
	filter(fraction=='total')%>%
	group_by(time,fraction,Protein_IDs)%>%
	filter(length(signal)==3)%>%
	filter(!any(is.na(signal)))%>%
	summarise(mean_sig=log2(mean(signal)),sd_signal = log2(sd(signal)))


#now plot the relationship between log2 signal and sd of the log2 signal
pdfexpr('../plots/modelling/mass_spec_sdplot.pdf',{
{
	sdplotdf%>%	
	ggplot(aes(x=log2(mean_sig),y=log2(sd_signal)))+
	# geom_bkde2d()+
	facet_wrap(~time)+
	# geom_point(alpha=I(0.1))+
	geom_smooth(method = "lm", se = FALSE)+
	stat_density_2d(aes(fill = ..level..), geom = "polygon")+
	geom_abline(slope=1)+
	theme_bw()
}
})


#now plot the relationship between log2 signal and sd of the log2 signal
pdfexpr('../plots/modelling/mass_spec_sdplot.pdf',{
{
	sdplotdf%>%	
	ggplot(aes(x=(mean_sig),y=(sd_signal)))+
	# geom_bkde2d()+
	facet_wrap(~time)+
	# geom_point(alpha=I(0.1))+
	geom_smooth(method = "lm", se = FALSE,data=sdplotdf%>%filter(between(mean_sig,27,32) ))+
	stat_density_2d(aes(fill = ..level..), geom = "polygon")+
	geom_abline(slope=1)+
	theme_bw()
}
})



#now plot the relationship between log2 signal and sd of the log2 signal
pdfexpr('../plots/modelling/mass_spec_sdplot_all.pdf',{
{
	sdplotdf%>%	
	ggplot(aes(x=(mean_sig),y=(sd_signal)))+
	# geom_bkde2d()+
	geom_point(alpha=I(0.1))+
	geom_smooth(method = "lm", se = TRUE)+
	# stat_density_2d(aes(fill = ..level..), geom = "polygon")+
	geom_abline(slope=1)+
	theme_bw()
}
})


#seems to be a predictable enough linear function, let's get those fits

sd_mean_lmfits <- sdplotdf%>%	summarise(mean_sig=log2(mean(signal)),sd_signal = log2(sd(signal)))%>%
	group_by(time)%>%nest%>%mutate(fit = map(data, ~ lm (data=. , sd_signal ~ mean_sig)) )

sd_mean_lmfits[[3]]

meansiglopes <- sd_mean_lmfits%>%mutate(l2mean_sig_slp =  fit%>%map_dbl(~ .$coef['mean_sig']))%>%{setNames(.$l2mean_sig_slp,.$time)}
meansigint <- sd_mean_lmfits%>%mutate(l2mean_sig_slp =  fit%>%map_dbl(~ .$coef['Intercept']))%>%{setNames(.$l2mean_sig_slp,.$time)}

sd_mean_lmfits$fit[[1]]$coef


####Check math on this sd calc see https://en.wikipedia.org/wiki/Log-normal_distribution
#we want u and s^2, the relevant parameters for the values on a log scale
# we have m and and a means of getting v, the relevant parametrs for the values on the non-log scale
# u = log ((m)/(sqrt(1+ v/(m^2))))
# s^2 = log (1+v/(m^2))
#so....
# s = sqrt(log(1 + sd(m)^2 / (m^2)))

#see 'Simulated linear test applied to quantitativeproteomics' eq 16,17
#for the comments about coefficient of variation approach - but my lines aren't straight

cvslope = sd_mean_lmfits$fit[[1]]$coef['mean_sig']
cvint = sd_mean_lmfits$fit[[1]]$coef['(Intercept)']
l2means = (sdplotdf%>%filter(time=='E13')%>%.$mean_sig)

means = 2^(sdplotdf%>%filter(time=='E13')%>%.$mean_sig)
sds = 2^((l2means*cvslope)+cvint)
s = sqrt(log2(1 + (sds^2) / (means^2) ))
u = log2(means/(sqrt(1+(sds)/(means^2))))
simdata=replicate(3,rnorm(n=length(u),mean=u,sd=s))
# simdatamean = simdata%>%apply(1,mean)
simdatasd = 2^simdata%>%apply(1,sd)

get_l2_ms_params <- function(m,cvslope,cvint){
	l2m = log2(m)
	sds = 2^((l2m*cvslope)+cvint)
	s = sqrt(log2(1 + (sds^2) / (m^2) ))
	u = log2(m/(sqrt(1+(sds)/(m^2))))
	list(u=u,s=s)
}

ms_params <- 
 list(
 	cvslope=sd_mean_lmfits$fit%>%map_dbl(~.$coef['mean_sig'])%>%setNames(sd_mean_lmfits$time),
	cvint=sd_mean_lmfits$fit%>%map_dbl(~.$coef['(Intercept)'])%>%setNames(sd_mean_lmfits$time)
)

get_l2_ms_params(means[1:5],ms_params$cvslope,ms_params$cvint)

pdfexpr('../plots/modelling/mass_spec_sdplot_sim.pdf',{
	data.frame(mean_sig=l2means,sd_signal=log2(simdatasd))%>%	
		ggplot(aes(x=(mean_sig),y=(sd_signal)))+
		# geom_bkde2d()+
		geom_point(alpha=I(0.1))+
		geom_smooth(method = "lm", se = TRUE)+
		# stat_density_2d(aes(fill = ..level..), geom = "polygon")+
		geom_abline(slope=1)+
		theme_bw()
})







#this is the negative log likelihood functon (i.e., func to be minimized)
#ribo - a gene,timepoint matrix
#MS - a gene,timepoint,replicate matrix
#rTE a gene, vector - the 
#prot0 a gene, vector - nuisance parameter - the initial protein level
#ms_s a time, vector - used to get the SD for the log-MS values from the log-MS values
nLL_model <- function(ldeg,prot0,ribo,MS,rTE,ms_sd){
	#our degredation constant gets optimized in log space, but used in linear space
	deg = exp(ldeg)	
	#Our protein vector is shaped like a slice of the MS [gene,time,replicate] array
	prot = MS[,,1,drop=F]
	dim(prot)=dim(prot)[1:2]
	prot[] <- 0
	prot[,1] <- prot0
	#build up our protein array tp by tp
	i=2
	for (i in 2:ncol(ribo)){
		prot[,i] <- prot[,i-1,drop=F] + (rTE*ribo[,i,drop=F]) - (prot[,i-1,drop=F]*deg)
	}
	#transform our MS to log scale
	prot <- log2(prot+1)
	#the sd will be a linear func of the mean, as per plots
	ms_sd <- sweep(prot,2,ms_sd,FUN='*')
	#finally get our log likelihood
	out <- 
	-sum(
		dnorm(log2(MS+1),mean=prot,sd=prot * ms_sd,log=TRUE)
	,na.rm=T)
	if(!is.finite(out)) browser()
	out
}

nLL_model_deg <- function(ldeg_prot0,ribo,MS,rTE,ms_params){
	ldeg=ldeg_prot0[1]
	prot0=ldeg_prot0[2]
	#our degredation constant gets optimized in log space, but used in linear space
	deg = exp(ldeg)	
	#Our protein vector is shaped like a slice of the MS [gene,time,replicate] array
	prot = MS[,,1,drop=F]
	dim(prot)=dim(prot)[1:2]
	prot[] <- 0
	prot[,1] <- prot0
	#build up our protein array tp by tp
	i=2
	for (i in 2:ncol(ribo)){
		prot[,i] <- prot[,i-1,drop=F] + (rTE*ribo[,i,drop=F]) - (prot[,i-1,drop=F]*deg)
	}
	#transform our MS to log scale, using our parameters
	l2ms_params<-get_l2_ms_params(prot[1,],ms_params$cvslope,ms_params$cvint)

	#finally get our log likelihood
	out <- 
	-sum(
		dnorm(log2(MS+1),mean=l2ms_params$u,sd=l2ms_params$s,log=TRUE)
	,na.rm=T)
	# if(!is.finite(out)) browser()
	# browser()
	out
}


nLL_model_deg_sep <- function(ldeg,prot0,ribo,MS,rTE,ms_params){
	#our degredation constant gets optimized in log space, but used in linear space
	deg = exp(ldeg)	
	#Our protein vector is shaped like a slice of the MS [gene,time,replicate] array
	prot = MS[,,1,drop=F]
	dim(prot)=dim(prot)[1:2]
	prot[] <- 0
	prot[,1] <- prot0
	#build up our protein array tp by tp
	i=2
	for (i in 2:ncol(ribo)){
		prot[,i] <- prot[,i-1,drop=F] + (rTE*ribo[,i,drop=F]) - (prot[,i-1,drop=F]*deg)
	}
	#transform our MS to log scale, using our parameters
	l2ms_params<-get_l2_ms_params(prot[1,],ms_params$cvslope,ms_params$cvint)

	#finally get our log likelihood
	out <- 
	-sum(
		dnorm(log2(MS+1),mean=l2ms_params$u,sd=l2ms_params$s,log=TRUE)
	,na.rm=T)
	# if(!is.finite(out)) browser()
	# browser()
	out
}

#we are currently plotting the model's fit
nLL_model_deg_plot <- function(ldeg_prot0,ribo,MS,rTE,ms_params){
	ldeg=ldeg_prot0[1]
	prot0=ldeg_prot0[2]
	#our degredation constant gets optimized in log space, but used in linear space
	deg = exp(ldeg)	
	#Our protein vector is shaped like a slice of the MS [gene,time,replicate] array
	prot = MS[,,1,drop=F]
	dim(prot)=dim(prot)[1:2]
	prot[] <- 0
	prot[,1] <- prot0
	#build up our protein array tp by tp
	i=2
	for (i in 2:ncol(ribo)){
		prot[,i] <- prot[,i-1,drop=F] + (rTE*ribo[,i,drop=F]) - (prot[,i-1,drop=F]*deg)
	}
	#transform our MS to log scale
	# prot <- log2(prot+1)
	#the sd will be a linear func of the mean, as per plots
	# ms_sd <- sweep(prot,2,ms_sd,FUN='*')
	#finally get our log likelihood
	l2ms_params<-get_l2_ms_params(prot[1,],ms_params$cvslope,ms_params$cvint)
	
	melt(MS)%>%set_colnames(c('gene_id','time','rep','signal'))%>%
	mutate(var='MS',pred=F,signal=log2(signal+1))%>%
	{bind_rows(.,data.frame(gene_id=1,time=1:length(l2ms_params$u),rep=1,var='MS',pred=TRUE,signal=l2ms_params$u) )}%>%
	{bind_rows(.,data.frame(gene_id=1,time=1:length(ribo),pred=F,rep=1,var='ribo',signal=log2(ribo[1,]) ))}%>%
	{
		ggplot(.,aes(y=signal,x=time))+
		geom_point(data=filter(.,!pred))+
		#facet_grid(scale='free',var ~ .)
		geom_ribbon(data=filter(.,pred),aes(y=signal,ymax=signal+(2*l2ms_params$s),ymin=signal-(2*l2ms_params$s)),alpha=0.1)+
		geom_line(data=filter(.,pred),alpha=0.1)+
		facet_grid(scale='free',var ~. )
	}%>%
	identity
	# ggsave(file='../plots/modelling/example_ms_fit.pdf'%T>%{normalizePath(.)%>%message})
}

n_genes = dim(expr_array)[1]
# debug(nLL_model)

#test the likelihood function in single gene and group modes
igenevect <- c(1,188)

nLL_model(
	ribo=ribo_matrix[igenevect,,drop=F],
	MS=expr_array[igenevect,,,drop=F],
	ldeg=rep(log(0.5),length(igenevect)),
	rTE=100,
	prot0=expr_array[igenevect,1,1],
	ms_sd=meansiglopes
)

nLL_model(
	ribo=ribo_matrix[10,,drop=F],
	MS=expr_array[10,,,drop=F],
	ldeg=c(log(0.5),log(0.1))[1],
	rTE=100,
	prot0=rep(30986,2)[1],
	ms_sd=meansiglopes
)

nLL_model_deg(
	ribo=ribo_matrix[10,,drop=F],
	MS=expr_array[10,,,drop=F],
	ldeg_prot0=c(c(log(0.5),log(0.1))[1],rep(30986,2)[1]),
	rTE=100,
	ms_params=ms_params
)

nLL_model_deg_plot(
	ribo=ribo_matrix[5,,drop=F],
	MS=expr_array[5,,,drop=F],
	ldeg_prot0=c(c(log(0.5),log(0.1))[1],rep(30986,2)[1]),
	rTE=100,
	ms_params=ms_params
)%>%ggsave(file='../plots/modelling/example_ms_fit_2.pdf'%T>%{normalizePath(.)%>%message})

Q

#now you need to call as above in an mle loop to get the degredation constants for individual genes
estimate_degredation <- function(ribo,MS,ldeg,rTE,prot0,ms_params){
	#fixed vs start
	startpars <-list(ldeg=ldeg,prot0=prot0)
	fixedpars <- list(rTE=rTE,MS=MS,ribo=ribo,ms_params=ms_params)
	#checks
	stopifnot(is.array(MS))
	stopifnot(is.numeric(MS))
	stopifnot(is.numeric(ribo))
	stopifnot(between(exp(ldeg),0,1))
	#testrun	
	do.call(nLL_model_deg_sep,c(startpars,fixedpars))
	# out<-mle(nLL_model_deg_sep,
	# 	start=startpars,
	# 	# start=list(par=c(log(0.5),TEinitial,sampdata$MS1[1])),
	# 	# method='Brent',
	# 	method='L-BFGS-B',
	# 	fixed=fixedpars,
	# 	lower=c(log(0),0+(1/1e12)),
	# 	upper=c(log(1-0.1),1e20),
	# )
	out<-optimwrap(nLL_model_deg_sep,
		start=startpars,
		# start=list(par=c(log(0.5),TEinitial,sampdata$MS1[1])),
		# method='Brent',
		method='L-BFGS-B',
		fixed=fixedpars,
		lower=c(log(0),0+(1/1e12)),
		upper=c(log(1-0.1),1e20),
	)

	nLL_model_deg_sep


	if(out$convergence==0){
		return(out$par)
	}else{
		return(NA)
	}
}
# lapply(1:nrow(ribo_matrix),function(i){
lapply(1:10,function(i){
	estimate_degredation(
		ribo=ribo_matrix[i,,drop=F],
		MS=expr_array[i,,,drop=F],
		ldeg=c(log(0.5),log(0.1))[1],
		rTE=100,
		prot0=rep(30986,2)[1],
		ms_params=ms_params
	)	
})









pdf('tmp.pdf')
melt(MS)%>%set_colnames(c('gene_id','time','rep','signal'))%>%mutate(var='MS')%>%
	{bind_rows(.,data.frame(gene_id=1,time=1:length(ribo),rep=1,var='ribo',signal=ribo[1,]) )}%>%
	ggplot(aes(y=signal,x=time))+geom_point()+facet_grid(scale='free',var ~ .)
dev.off()

pdfexpr('tmp.pdf',{plot(1)
})







0
estimate_degredations <- function(exprdata,rTE){

	fits<-mclapply(exprdata$data,safely(estimate_degredation)(exprdata_g), rTE = rTE )

	worked <- fits%>%map('result')%>%Negate(is.null)

	stopifnot(all(worked))

	tests%>%
		map('result')%>%
		map(set_rownames,c('deg','prot0'))%>%
		map(rownames_to_column,'var')%>%
		bind_rows%>%
		filter(var=='rTE')%>%
		set_colnames(c('var','lower','upper','est'))%>%
#Plot the distribution of parameter estimates
tests%>%map('result')%>%
	keep(Negate(is.null))%>%
	map(set_rownames,c('deg','rTE','prot0'))%>%
	map(rownames_to_column,'var')%>%
	bind_rows%>%
	filter(var=='rTE')%>%
	set_colnames(c('var','lower','upper','est'))%>%
	{ ggpubr::ggarrange(nrow=2,
		ggplot(.,aes(ymin=lower,ymax=upper,y=est,x=seq_along(est)))+geom_pointrange()+coord_cartesian(ylim=c(0,1)),
		ggplot(.,aes(est))+geom_histogram()+scale_x_continuous(limits=c(0,5)),
		)
	}
}


#I need to look at the genes which are outliers and figure out what's wrong with them....
exprdata_g <- exprdatareshape%>%group_by(gene_name)%>%nest%>%filter(gene_name %in% 'Satb2')%>%unnest


library(assertthat)
assert_that( exprdatareshape%>%group_by(gene_name)%>%tally%>%pluck('n')%>%table%>%length%>%`==`(1) ,
	msg='all gene names should have the same number of assocaited time points')


degests <- exprdatareshape %>% group_by(gene_name)%>% nest %>% slice(1:10) %>% mutate(deg = map(data,~safely(estimate_degredations)(rTE,.)))
degests %<>% mutate(iserror = deg%>%map('result')%>%map_lgl(is.null))
degests[3,]$deg



exprdatareshape$gene_name

gene_namei<-ugnames[1]



iter = 0
exprdata <- exprdatareshape
message('using gene name for now')
exprdata$gene_id <- exprdata$gene_name

exprdata



#Function to get the ribo vector we use to infer the protein levels, can mod later with splines
get_trans_vect <- function(ribodata){
	ribodata%>%as.matrix%>%apply(1,mean,na.rm=T)
}

exprdata$ribo <- get_trans_vect(exprdata%>%ungroup%>%select(matches('ribo')))

mscols <- colnames(exprdatareshape)%>% str_subset('MS')

exprdata %>% 
	select(-matches('gene_name'))%>%
	select(gene_id,matches('ribo'),one_of(mscols))%>%
	group_by(gene_id)

#

#main loop of the optimization cycle
fit_rTE_degs <- function(exprdata,indiv_ests){

	rTE <- indiv_ests%>%spread(var,est)%>%pluck('rTE')%>%na.omit%>%{median(.[.!=min(.)])}
	start_deg_med <- indiv_ests%>%spread(var,est)%>%pluck('deg')%>%na.omit%>%{median(.[.!=min(.)])}
	degredations <- indiv_ests%>%spread(var,est)%>%pluck('deg')%>%replace_na(start_deg_med)


	#Step 2
	#Taking the median rTE from the last stage as the initial value, I now re estimate degredation rates, 
	#then reestimate rTE, etc, until convergence.....
	while(!convergence){
		
		degredations <- estimate_degredations(exprdata,rTE)

		rTE <- estimate_rTE(exprdata,degredations)

		conv <- test_convergence(degredations,rTE,old_degredations,old_rTE)

		if(conv > 5) convergence <- TRUE

		old_degredations <- degredations
		rTE <- old_rTE

		iter = iter + 1

		if(iter > max_iter) break
	}

	if(!convergence){
		if(iter > max_iter){
			message('maximum iterations reached without convergence')
		}else{
			stop('something went wrong, no convergence but iter < max')
		}
	}

	return(list(rTE,degredations))

}

expr_array <- exprdata%>%ungroup%>%select(one_of(mscols))%>%as.matrix%>%array(dim=c(length(.)/15,5,3))
ribo_matrix <- exprdata$ribo%>%matrix(ncol=5,)

debug(nLL_model)

#test the likelihood function in single gene and group modes
nLL_model(
	ribo_matrix[T,,drop=F],
	expr_array[T,,,drop=F],
	ldeg=c(log(0.5)),
	rTE=100,
	prot0=3098,
	ms_sd=meansiglopes,
)

nLL_model(
	ribo_matrix[10,,drop=F],
	expr_array[10,,,drop=F],
	ldeg=c(log(0.5),log(0.1))[1],
	rTE=100,
	prot0=rep(30986,2)[1],
	ms_sd=meansiglopes,
	tps=5
)

indiv_ests <- fit_indiv_rTE_degs(exprdata,indiv_ests)

params_fit <- fit_rTE_degs(exprdata,indiv_ests)

#Step 3
#Plot some of the weird looking ones....

###Not sure what the most efficient way to strucutre the esxpr data is..
####I'll use nest and unnest to start with - easy error handlng etc.


#Okay so this is working.... acceptably well.
#right now I can have my neg log likelihood function with all arguments sep, and I use optimwrap (just a truncated mle() ) to
#clobber the function into something optim can pass it's vectors into....
#I don't know if there's a speed down associated with the function overhead of having the optimized function wrapped like that in the
#optime loop
#My plotting function works okay, I should modify it to take seperate vairables again...
#I've also spent a while trying to get the mass spec variance right - it's too wide right now, from the looks of things.
#I should insepct that with my plotting function and then do something about it
#I could try just using the simpler CV = exp(-c) formula for the variance, although the slopes don'' reeeeeaaaaaalllly look even
#Or I could just leave it...
#Like, is there a constant variance within a single gene