
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

ugnames <- exprdatareshape$gene_
ugnames%>%intersect('Satb2')

getcoefestimates<-function(gene_namei){
	sampdata <-  %>%ungroup%>%filter(gene_name==gene_namei)
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
     expr
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
sdplotdf<-'ms_tables/ms_LFQ_total_ms_tall.tsv'%>%
	read_tsv%>%
	filter(fraction=='total')%>%
	group_by(time,fraction,Protein_IDs)%>%
	filter(length(signal)==3)%>%
	filter(!any(is.na(signal)))
sdplotdf%>%	summarise(mean_sig=mean(signal),sd_signal = sd(signal))%>%
	ggplot(aes(x=log2(mean_sig),y=log2(sd_signal)))+geom_point()+facet_wrap(~time)

pdfexpr('../plots/modelling/mass_spec_sdplot.pdf',{

	
{sdplotdf%>%	summarise(mean_sig=mean(signal),sd_signal = sd(signal))%>%
	ggplot(aes(x=log2(mean_sig),y=log2(sd_signal)))+
	# geom_bkde2d()+
	facet_wrap(~time)+
	# geom_point(alpha=I(0.1))+
	geom_smooth(method = "lm", se = FALSE)+
	stat_density_2d(aes(fill = ..level..), geom = "polygon")
}

})

sd_mean_lmfits <- sdplotdf%>%	summarise(mean_sig=log2(mean(signal)),sd_signal = log2(sd(signal)))%>%
	group_by(time)%>%nest%>%mutate(fit = map(data, ~ lm (data=. , sd_signal ~ mean_sig)) )

meansiglopes <- sd_mean_lmfits%>%mutate(l2mean_sig_slp =  fit%>%map_dbl(~ .$coef['mean_sig']))%>%{setNames(.$l2mean_sig_slp,.$time)}

#this is the negative log likelihood functon (i.e., func to be minimized)
#ribo - a gene,timepoint matrix
#MS - a gene,timepoint,replicate matrix
#rTE a gene, vector - the 
#prot0 a gene, vector - nuisance parameter - the initial protein level
#ms_s a time, vector - used to get the SD for the log-MS values from the log-MS values
nLL_model <- function(ribo,MS,ldeg,rTE,prot0,ms_sd){
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
	-sum(
		dnorm(log2(MS),mean=prot,sd=prot * ms_sd,log=TRUE)
	,na.rm=T)
}

n_genes = dim(expr_array)[1]
# debug(nLL_model)

#test the likelihood function in single gene and group modes
igenevect <- c(1,188)
nLL_model(
	ribo_matrix[igenevect,,drop=F],
	expr_array[igenevect,,,drop=F],
	ldeg=rep(log(0.5),length(igenevect)),
	rTE=100,
	prot0=expr_array[igenevect,1,1],
	ms_sd=meansiglopes
)

nLL_model(
	ribo_matrix[10,,drop=F],
	expr_array[10,,,drop=F],
	ldeg=c(log(0.5),log(0.1))[1],
	rTE=100,
	prot0=rep(30986,2)[1],
	ms_sd=meansiglopes
)

#this is similiar to the above but takes ALL gene's expression as the input
#and also a vector of degredation constants, 
nLL_rTE <- function(ribo,MS,ldeg,rTE,prot0,ms_sd){

	deg = exp(ldeg)	

	tps = nrow(ribo1)

	prot = MS1
	prot[] <- 0

	prot[1,] = prot0

	trans <- (ribo1+ribo2)/2

	for (i in 2:tps){
		prot[i,] <- prot[i-1,] + (rTE*trans[i,]) - (prot[i-1,]*deg)
	}

	prot <- log2(prot+1)

	-sum(c(
		dnorm(log2(MS1),mean=prot,sd=0.1,log=T),
		dnorm(log2(MS2),mean=prot,sd=0.1,log=T),
		dnorm(log2(MS3),mean=prot,sd=0.1,log=T)
	))
}


nLL_deg <- function(ribo,MS,ldeg,rTE,prot0,ms_sd,tps){

	deg = exp(ldeg)	

	

	prot[1] = prot0

	for (i in 2:tps){
		prot[i] <- prot[i-1] + (rTE*trans[i]) - (prot[i-1]*deg)
	}
	prot <- log2(prot+1)
	-sum(c(
		dnorm(log2(MS1),mean=prot,sd=ms_sd,log=T),
		dnorm(log2(MS2),mean=prot,sd=ms_sd,log=T),
		dnorm(log2(MS3),mean=prot,sd=ms_sd,log=T)
	))
}

estimate_degredation <- function(rTE,exprdata_g){
	startpars <-list(ldeg=log(0.5),prot0=exprdata_g$MS1[1])
	fixedpars <- list(rTE=rTE,ribo1=exprdata_g$ribo1,ribo2=exprdata_g$ribo2,MS1=exprdata_g$MS1,MS2=exprdata_g$MS2,MS3=exprdata_g$MS3)%>%
		lapply(na.omit)

	do.call(nloglikms,c(startpars,fixedpars))
	
	out<-mle(nloglikms,
		start=startpars,
		# start=list(par=c(log(0.5),TEinitial,sampdata$MS1[1])),
		method='L-BFGS-B',
		fixed=fixedpars,
		lower=c(log(0.00001),0,0),
		upper=c(log(1),1e20,1e20))

	# cint <- confint(out)
	# cint[1,1:2]%<>%exp%>%round(3)
	# rownames(cint)[1]%<>%str_replace('ldeg','deg')
	# cint[2,1:2]%<>%round(5)
	# cint[3,1:2]%<>%round(5)
	ests <- coef(out)[names(startpars)]
	# ests[1] %<>% exp%>%round(5)
	# cint%>%as.data.frame %>%mutate(estimate=round(ests,5))
	ests['ldeg']
}


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
fit_rTE_degs <- function(exprdata,indiv_ests)

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
	tps=5
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


