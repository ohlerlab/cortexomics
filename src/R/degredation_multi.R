library(tidyverse)
library(zeallot)
library(ggpubr)
library(stats4)

rawcounts <- fread('feature_counts/all_feature_counts')
rawlfq <- fread('ms_tables/ms_LFQ_total_ms_tall.tsv')
rawlfq$time%<>%str_replace('p','')
testlfq <- rawlfq%>%
	filter(gene_name%>%str_detect(regex('satb2|Raf1',ign=T)),!Protein_IDs%>%str_detect('H3BKH3'))

testlfq%<>%select(gene_name,time,replicate,signal)%>%
	mutate(datacol=paste0('MS',replicate))%>%
	select(-replicate)%>%
	spread(datacol,signal)

rawcounts%>%head

testcounts <- rawcounts%>%
	select(gene_id=feature_id,matches('ribo'))%>%
	left_join(fread('ids.txt'))%>%
	filter(gene_name%>%str_detect(regex('satb2|raf1',ign=T)))%>%
	gather(datacol,signal,-gene_id,-gene_name)%>%
	separate(datacol,c('time','dat','replicate'))%>%
	mutate(replicate=paste0('ribo',replicate))%>%
	select(-gene_id,-dat)%>%
	spread(replicate,signal)%>%
	identity

testdata<-left_join(testlfq,testcounts)


# nloglikms <- function(ribo1,ribo2,MS1,MS2,MS3,ldeg=log(0.5),rTE = TEinitial,prot0 = sampdata$MS1[1] ){

# 	deg = exp(ldeg)	

# 	tps = length(ribo1)

# 	prot = rep(0,tps)

# 	prot[1] = prot0

# 	trans <- (ribo1+ribo2)/2


# 	for (i in 2:tps){
# 		prot[i] <- prot[i-1] + (rTE*trans[i]) - (prot[i-1]*deg)
# 	}
# 	prot <- log2(prot+1)
# 	-sum(c(
# 		dnorm(log2(MS1),mean=prot,sd=0.1,log=T),
# 		dnorm(log2(MS2),mean=prot,sd=0.1,log=T),
# 		dnorm(log2(MS3),mean=prot,sd=0.1,log=T)
# 	))
# }


# startpars <-list(ldeg=log(0.5),
# 	rTE=testdata$MS1[1]/testdata$ribo1[1],
# 	prot0=testdata$MS1[1])
# fixedpars <- list(ribo1=testdata$ribo1,ribo2=testdata$ribo2,MS1=testdata$MS1,MS2=testdata$MS2,MS3=testdata$MS3)

# do.call(nloglikms,c(startpars,fixedpars))



# out<-mle(nloglikms,
# 	start=startpars,
# 	method='L-BFGS-B',
# 	fixed=fixedpars,
# 	lower=c(log(0.00001),0,0),
# 	upper=c(log(1),1e12,1e12)
# )

# testdata

# out@coef

# library(gganimate)


# testdata=expand.grid(time=1:5,rep=1:10)%>%
# 	data.frame%>%
# 	mutate(signal=time+rnorm(length(time),sd=0.1*time))

# p <- testdata%>%ggplot(aes(x=time,y=signal))+geom_point()
# p <- p + transition_states(time)
# rlang::last_error()



#scale up to multiple genes

nloglikms_matrix <- function(ribo1,ribo2,MS1,MS2,MS3,ldeg,rTE,prot0){

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

n_genes=2


startpars <-list(ldeg=log(0.5)%>%rep(n_genes),
	rTE=testdata%>%group_by(gene_name)%>%slice(1)%>%mutate(t=MS1/ribo1)%>%.$t,
	prot0=testdata%>%group_by(gene_name)%>%slice(1)%>%mutate(t=mean(MS1,MS2,MS3,na.rm=T))%>%.$t
)

fixedpars <- list(
	ribo1=matrix(testdata$ribo1,ncol=2),
	ribo2=matrix(testdata$ribo2,ncol=2),
	MS1=matrix(testdata$MS1,ncol=2),
	MS2=matrix(testdata$MS2,ncol=2),
	MS3=matrix(testdata$MS3,ncol=2)
)

do.call(nloglikms_matrix,c(startpars,fixedpars))

out<-mle(nloglikms_matrix,
	start=startpars,
	method='L-BFGS-B',
	fixed=fixedpars,
	lower=c(log(0.00001),0,0),
	upper=c(log(1),1e12,1e12)
)

unlist(startpars)

fit = optim(par=unlist(startpars),
	method='L-BFGS-B',
	# control=list('fnscale'= -1),
	fn=nloglikms_matrix,
	ribo1=matrix(testdata$ribo1,ncol=2),
	ribo2=matrix(testdata$ribo2,ncol=2),
	MS1=matrix(testdata$MS1,ncol=2),
	MS2=matrix(testdata$MS2,ncol=2),
	MS3=matrix(testdata$MS3,ncol=2),
	lower=unlist(rep(c(log(0.00001),0,0),n_genes)),
	upper=unlist(rep(c(log(1),1e12,1e12),n_genes)),
	hessian=T)

multiple genes

nloglikms_matrix <- function(ribo1,ribo2,MS1,MS2,MS3,ldeg,rTE,prot0){

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

Lorenz<-function(t, state, parameters) {
	with(as.list(c(state, parameters)),{
	# rate of change
	dP <- rTE * ribo - exp(ldeg) * prot   
	dY <- b * (Y-Z)
	dZ <- -X*Y +c*Y - Z
	# return the rate of change
	list(c(dX, dY, dZ))
	}) # end with(as.list ...
}






