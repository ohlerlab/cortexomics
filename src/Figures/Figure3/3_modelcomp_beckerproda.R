################################################################################
################################################################################
base::source(here::here('src/R/Rprofile.R'))
library(rstan)

base::source('src/R/Functions/rstan_functions.R')

if(!exists("cdsgrl")) {
	# base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
	load('data/1_integrate_countdata.R')
}


sel_prodpreds<-readRDS('data/sel_prodpreds.rds')

sel_prodpreds%<>%mutate(se = ifelse(is.na(se),(CI.R-CI.L)/3.92,se))

# sel_ms_mat<-readRDS('data/sel_ms_mat.rds')
countpred_df<-readRDS('data/countpred_df.rds')
countpred_df%<>%mutate(se = (CI.R-CI.L)/3.92)
# tx_countdata<-readRDS('data/tx_countdata.rds')


#Things we want ot sim - production slow half life, production high half life, 
bmodel <- rstan::stan_model('src/Stan/becker_proda.stan')
bmodel_stationary <- rstan::stan_model('src/Stan/becker_proda_stationary.stan')
bmodel_degonly <- rstan::stan_model('src/Stan/bmodel_degonly.stan')
bmodel_linear <- rstan::stan_model('src/Stan/becker_proda_linear.stan')
allTEchangedf<-read_tsv('tables/manuscript/go_all_highcount_updown.tsv')

timepoints = ribocols%>%str_replace('_ribo_\\d','')

downtegenes = allTEchangedf%>%filter(down==1)%>%.$gene_name
uptegenes = allTEchangedf%>%filter(up==1)%>%.$gene_name
# bmodel_stationary = fix_param(bmodel,vars2fix = c('l_st','l_pihalf'))%>%{f='src/Stan/bmodel_stationary.stan';cat(.,file=f);f}%>%stan_model
#bmodel_degonly = fix_param(bmodel,vars2fix = c('l_st'))%>%{f='src/Stan/bmodel_degonly.stan';cat(.,file=f);f}%>%stan_model

gid2gnm = mcols(cds)%>%as.data.frame%>%distinct(gene_id,gene_name)%>%{safe_hashmap(.[[1]],.[[2]])}
gnm2gid = mcols(cds)%>%as.data.frame%>%distinct(gene_id,gene_name)%>%{safe_hashmap(.[[2]],.[[1]])}
countpred_df$gene_name = gid2gnm[[countpred_df$gene_id]]
sel_prodpreds$gene_name = gid2gnm[[sel_prodpreds$gene_id]]

count_ests = countpred_df%>%distinct(contrast,gene_name,.keep_all=TRUE)%>%select(gene_name,contrast,logFC)%>%spread(contrast,logFC)%>%{set_rownames(as.matrix(.[,-1]),.$gene_name)}
count_ses = countpred_df%>%distinct(contrast,gene_name,.keep_all=TRUE)%>%select(gene_name,contrast,se)%>%spread(contrast,se)%>%{set_rownames(as.matrix(.[,-1]),.$gene_name)}
prot_ests = sel_prodpreds%>%distinct(time,gene_name,.keep_all=TRUE)%>%select(gene_name,time,diff)%>%spread(time,diff)%>%{set_rownames(as.matrix(.[,-1]),.$gene_name)}
prot_ses = sel_prodpreds%>%distinct(time,gene_name,.keep_all=TRUE)%>%select(gene_name,time,se)%>%spread(time,se)%>%{set_rownames(as.matrix(.[,-1]),.$gene_name)}
ribocols = count_ests%>%colnames%>%str_subset('ribo')
rnacols = count_ests%>%colnames%>%str_subset('total')

xtailfoldchange<-Sys.glob('pipeline/xtail/xtail_*')%>%setNames(timepoints[-1])%>%
	map_df(.id='time',fread)%>%
	mutate(changetype='translational_xtail')

teupgenes=xtailfoldchange%>%filter(time=='P0_ribo')%>%filter(adj_p_value<0.05,log2fc > 1)%>%.$gene_name
tedowngenes = xtailfoldchange%>%filter(time=='P0_ribo')%>%filter(adj_p_value<0.05,log2fc < -1)%>%.$gene_name

protscl = prot_ests%>%colMedians
prot_ests = sweep(prot_ests,2,'-',STAT=protscl)

countrscl = count_ests%>%colMedians
count_ests = sweep(count_ests,2,'-',STAT=countrscl)

make_standata <- function(g,prot_ests,c_ests,prot_ses,c_ses){
	l2e = log2(exp(1))
	nT=ncol(c_ests)
	sampdata=list()
	sampdata$G=1
	sampdata$T=nT
	sampdata$l_st_priorsd=3
	sampdata$l_ribo_priorsd=3
	sampdata$l_pihalf_priormu=0.5
	sampdata$l_pihalf_priorsd=3
	sampdata$lMSmu = prot_ests[g,,drop=FALSE]/l2e #- median(log(P))
	sampdata$lSeqmu = c_ests[g,,drop=FALSE]/l2e #- median(log(ribo))
	sampdata$lMSsigma = prot_ses[g,,drop=FALSE]/l2e
	sampdata$lSeqsigma = c_ses[g,,drop=FALSE]/l2e
	sampdata
}
make_ribodata = function(g)make_standata(g,prot_ests,count_ests[,ribocols],prot_ses,count_ses[,ribocols])
make_rnadata = function(g)make_standata(g,prot_ests,count_ests[,rnacols],prot_ses,count_ses[,rnacols])

models = list(
	production = bmodel,
	stationary = bmodel_stationary,
	degredation = bmodel_degonly,
	linear = bmodel_linear
)
cgenes = rownames(count_ests)

testgrpsize = Inf

dec_genes = (prot_ests[,1]-prot_ests[,5])%>%sort%>%names%>%intersect(cgenes)%>%intersect(c(downtegenes,uptegenes))%>%head(testgrpsize)
incr_genes = (prot_ests[,1]-prot_ests[,5])%>%sort%>%names%>%intersect(cgenes)%>%intersect(c(downtegenes,uptegenes))%>%tail(testgrpsize)
notedec_genes = (prot_ests[,1]-prot_ests[,5])%>%sort%>%names%>%intersect(cgenes)%>%setdiff(c(downtegenes,uptegenes))%>%head(testgrpsize)
noteincr_genes = (prot_ests[,1]-prot_ests[,5])%>%sort%>%names%>%intersect(cgenes)%>%setdiff(c(downtegenes,uptegenes))%>%tail(testgrpsize)


testgenes = c(dec_genes,incr_genes,notedec_genes,noteincr_genes)

datafuns = list(riboseq =  make_ribodata , rnaseq =  make_rnadata)
		
opts = lapply(testgenes%>%setNames(.,.),function(gene){
	lapply(datafuns,function(datafun){
		simdata = datafun(gene)
		lapply(models,function(model){
			initvals=list()
			initvals$l_pihalf <- array(log(0.5),1)
			initvals$lribo = array(simdata$lSeqmu,c(1,length(simdata$lSeqmu)))
			initvals$l_st = array(simdata$lMSmu[[1]]-simdata$lSeqmu,1)
			initvals$lprot0 = array(simdata$lMSmu[[1]],1)
			reopt <- rstan::optimizing(model,data=simdata,init=initvals,as_vector=F,save_iterations=TRUE)
			reopt
		})
	})
})

#now let's do a chi squared test
model_tests = lapply(names(opts)%>%setNames(.,.),function(gene){
	lapply(names(datafuns)%>%setNames(.,.),function(datafun){
	lapply(names(models)%>%setNames(.,.),function(modname){
		opt = opts[[gene]][[datafun]][[modname]]
		gdata = datafuns[[datafun]](gene)
		#
		perrors =(log(opt$par$prot) - gdata$lMSmu)
		w_perrors = perrors/gdata$lMSsigma
		rerrors = (opt$par$lribo - gdata$lSeqmu)
		w_rperrors = rerrors/gdata$lSeqsigma
		n_df = opt$par[get_stanpars(models[[modname]])]%>%unlist%>%length
		errorsum = sum(c(w_perrors,w_rperrors)^2)
		BIC = -2*opt$value+(log(n_df))*n_df
		pval = 1 - pchisq(errorsum,n_df)
		c(BIC=BIC,pval=pval,residuals = w_perrors)
	})
})
})

modeltestdf = model_tests%>%
	map_df(.id='gene',.%>%
		map_df(.id='data',.%>%
			bind_rows(.id='model')))%>%
	mutate(passtest = ! (pval < 0.01))%>%
	group_by(gene,data)%>%
	mutate(best=BIC==min(BIC))
residlistcol = modeltestdf%>%ungroup%>%select(matches('residuals'))%>%as.matrix%>%t%>%as.data.frame
modeltestdf$residuals = residlistcol%>%as.list
modeltestdf%<>%select(-matches('residuals\\d'))

modeltestdf

igene=modeltestdf$gene[1]




modeltestdf%>%filter(!any(passtest))


modeltestdf%>%filter(data=='rnaseq')%>%filter(best)%>%.$model%>%table
modeltestdf%>%filter(data=='riboseq')%>%filter(best)%>%.$model%>%table

modeltestdf%>%filter(data=='rnaseq')%>%summarise(reject=all(!passtest))%>%.$reject%>%table
modeltestdf%>%filter(data=='riboseq')%>%summarise(reject=all(!passtest))%>%.$reject%>%table

#okay so testing the residuals gives us nothing to go on, effectively al of the time
resid_tests = map_df(.id='gene',names(opts)%>%setNames(.,.),function(igene){
	gdf = modeltestdf%>%filter(gene==igene,best)
	tidy(t.test(
	gdf%>%filter(data=='riboseq')%>%.$residuals%>%.[[1]],
	gdf%>%filter(data=='rnaseq')%>%.$residuals%>%.[[1]]
))
})

modeltestdf%>%group_by(gene,data)%>%
	summarise(mbic = min(BIC))%>%summarise(bestdata = data[which.min(mbic)])%>%
	mutate(iste = gene%in%c(uptegenes,downtegenes))%>%group_by(bestdata)%>%tally

modeltestdf%>%filter(best,data=='riboseq')%>%.$model%>%table
modeltestdf%>%filter(best,data=='rnaseq')%>%.$model%>%table


rnagoodtest = modeltestdf%>%filter(data=='rnaseq',model=='production',best,passtest)%>%.$gene

Ksvals = opts[rnagoodtest]%>%map_dbl(function(opt){
	opt[['riboseq']][['production']]$par$Ks%>%log
})




meantes = count_ests[,colnames(count_ests)%>%str_subset('TE')]%>%rowMeans

enframe(Ksvals,'gene_name','rna_Ks')%>%left_join(meantes%>%enframe('gene_name','TE'))%>%
	{cor.test(.$rna_Ks,.$TE)}

enframe(Ksvals,'gene_name','rna_Ks')%>%left_join(meantes%>%enframe('gene_name','TE'))%>%
	{quicktest(.$rna_Ks,.$TE)}







# make_simdata
# #data looks like this
# list(
# 	G;// number of proteins
#   int T;// info on the number of conditions
#   matrix[G,T] lMSmu;
#   matrix[G,T] lSeqmu;
#   matrix[T,T] lMSsigma[G];
#   matrix[T,T] lSeqsigma[G];
#   real l_st_priorsd;
#   real l_ribo_priorsd;
#   real l_pihalf_priormu;
#   real l_pihalf_priorsd;

#goal - a chi-squared test to distinguish the different types
#generate the different types

################################################################################
########This mess
################################################################################
if(F){

	library(tidyverse)
library(here)
library(magrittr)
library(splines2)
library(splines)


################################################################################
########
################################################################################
library(tidyverse)	
library(here)	
sampdata <- here('data/sampdata.rds')%>%read_rds
ribo <- c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 100, 100, 30, 
          30, 30, 30, 30, 30, 30, 30, 30, 30)
nT = length(ribo)
Pe = ribo*3
P = Pe
stopifnot(length(Pe) == length(ribo))


sampdata$T=nT
sampdata$lMSmu = t(log(P)) #- median(log(P))
sampdata$lSeqmu = t(log(ribo)) #- median(log(ribo))
sampdata$lMSsigma = array(diag(rep(0.1,nT)),c(1,nT,nT))
sampdata$lSeqsigma = array(diag(rep(0.1,nT)),c(1,nT,nT))

time = 1:length(ribo)
timeknots <- time[c(-1,-length(time))]
sampdata$mybs =  cbind(1,ibs(time, knots = timeknots,degree = 1, intercept = TRUE))
sampdata$mydbs =bs(time, knots = timeknots,degree = 1, intercept = TRUE)

initvals=list()
initvals$prot0 <- array(sampdata$lMSmu[1],1)

d_0 = 0.6
Kd = d_0
lKd = log(Kd)
lpihalf = log(log(2)) - lKd
pihalf = log(2) / Kd
lpihalf = log(log(2)) - lKd
log(log(2)) - lpihalf
lKd
#log(log(2)) - lKd




sampdata[names(sampdata)%>%str_subset('sigma$')] %<>% map(multiply_by,0.001)
# 
# sampdata%>%map(dim)
# sampdata
# proDAsigmastan2

d_0 = 0.6;
Ks = 40;
ribo= c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 100, 100, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30);
P = rep(NA,length(ribo))
P[1] = 3* (ribo[1] * Ks)

for(i in 2:length(ribo)){
  P[i] = stepsynthdeg(d_0=d_0,Ks=Ks,T=1, s_0 = ribo[i-1],s_1 =  ribo[i],P_0 = P[i-1])
}

lribo = log(ribo)
l_st = log(Ks)-log(d_0)
sampdata$prot0 = array(initvals$l_st+initvals$lribo[1],1)
initvals$l_pihalf <- array(log(1),1)
initvals$lribo = array(log(ribo),c(1,length(ribo)))
initvals$l_st = array(log(2),1)
initvals$lprot0 = 1+array(initvals$l_st+initvals$lribo[1],1)
# initvals$lprot0 = NA
initvals$prot0=NULL
Kd = log(log(2)) -  initvals$l_pihalf;
Ks = exp(initvals$l_st - lKd);
ribo = exp(lribo);
exp(initvals$lprot0);
beckermodel <- rstan::stan_model('src/Stan/becker_proda.stan')
opt <- rstan::optimizing(beckermodel,data=sampdata,verbose=TRUE,init=initvals,as_vector=F,iter=0,save_iterations=TRUE)
library(txtplot)
txtplot(log(opt$par$prot))
txtplot(opt$par$lribo)
opt$par$l_st
opt$par$lprot0
opt$par$l_st
log(opt$par$prot)
sampdata$lMSmu
initvals$lprot0
initvals$l_st

sampdata$lMSmu
sampdata2 <- sampdata
sampdata2$lMSmu <- log(opt$par$prot)
opt2 <- rstan::optimizing(beckermodel,data=sampdata2,verbose=TRUE,as_vector=F,save_iterations=TRUE)
txtplot(log(opt$par$prot))
opt2$par$l_pihalf
opt$par$l_pihalf
initvals$l_pihalf

opt2$par$l_st
opt$par$l_st

opt2$par$prot/opt$par$prot
txtplot(log(opt2$par$prot),log(opt$par$prot))


}