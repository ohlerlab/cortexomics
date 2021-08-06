################################################################################
################################################################################
base::source(here::here('src/R/Rprofile.R'))
library(rstan)

base::source('src/Archive/R/Functions/rstan_functions.R')

if(!exists("tx_countdata")) {
	# base::source("src/Figures/Figure0/0_load_annotation.R")
	load('data/1_integrate_countdata.R')
}

#actually raw data doesn't really show differences in mean IBAQ much
# rawprotdata = readRDS('data/proteinmsdata.rds')
sel_prodpreds<-readRDS('data/sel_prodpreds.rds')

sel_prodpreds%<>%mutate(se = ifelse(is.na(se),(CI.R-CI.L)/3.92,se))

# sel_ms_mat<-readRDS('data/sel_ms_mat.rds')
countpred_df<-readRDS('data/countpred_df.rds')
countpred_df%<>%mutate(se = (CI.R-CI.L)/3.92)
# tx_countdata<-readRDS('data/tx_countdata.rds')

{
#Things we want ot sim - production slow half life, production high half life, 
bmodel <- rstan::stan_model('src/Archive/Stan/becker_proda.stan')
# bmodel_ribooffset <- rstan::stan_model('src/Archive/Stan/becker_proda_ribooff.stan')
bmodel_msdev <- rstan::stan_model('src/Archive/Stan/becker_proda_msdev.stan')
bmodel_stationary <- rstan::stan_model('src/Archive/Stan/becker_proda_stationary.stan')
bmodel_degonly <- rstan::stan_model('src/Archive/Stan/bmodel_degonly.stan')
bmodel_linear <- rstan::stan_model('src/Archive/Stan/becker_proda_linear.stan')
bmodel_accumulation <- rstan::stan_model('src/Archive/Stan/becker_proda_accumulation.stan')
models = list(
	production = bmodel,
	# production_offset = bmodel_ribooffset,
	stationary = bmodel_stationary,
	degredation = bmodel_degonly,
	# accumulation = bmodel_accumulation,
	linear = bmodel_linear,
	msdev = bmodel_msdev
)
}

allTEchangedf<-'tables/xtailTEchange.tsv'%>%read_tsv
downtegenes = allTEchangedf%>%filter(down==1)%>%.$gene_name
uptegenes = allTEchangedf%>%filter(up==1)%>%.$gene_name
# bmodel_stationary = fix_param(bmodel,vars2fix = c('l_st','l_pihalf'))%>%{f='src/Archive/Stan/bmodel_stationary.stan';cat(.,file=f);f}%>%stan_model
#bmodel_degonly = fix_param(bmodel,vars2fix = c('l_st'))%>%{f='src/Archive/Stan/bmodel_degonly.stan';cat(.,file=f);f}%>%stan_model

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
protscl = prot_ests%>%colMedians
# protscl = prot_ests%>%colMedians%>%median%>%rep(ncol(prot_ests))
prot_ests = sweep(prot_ests,2,'-',STAT=protscl)
countrscl = count_ests%>%colMedians
count_ests = sweep(count_ests,2,'-',STAT=countrscl)
# if(USE_RAW){
# }
make_standata <- function(g,prot_ests,c_ests,prot_ses,c_ses){
	l2e = log2(exp(1))
	nT=ncol(c_ests)
	sampdata=list()
	sampdata$G=1
	sampdata$T=nT
	sampdata$l_st_priorsd=3
	sampdata$l_ribo_priorsd=3
	sampdata$l_pihalf_priormu=0
	sampdata$l_pihalf_priorsd=5
	sampdata$lMSmu = prot_ests[g,,drop=FALSE]/l2e #- median(log(P))
	sampdata$lSeqmu = c_ests[g,,drop=FALSE]/l2e #- median(log(ribo))
	sampdata$lMSsigma = prot_ses[g,,drop=FALSE]/l2e
	sampdata$lSeqsigma = c_ses[g,,drop=FALSE]/l2e
	sampdata$ribooffset = c(0, -0.1169, -0.4907, -0.5889, -0.7206)
	sampdata
}
make_ribodata = function(g)make_standata(g,prot_ests,count_ests[,ribocols],prot_ses,count_ses[,ribocols])
make_rnadata = function(g)make_standata(g,prot_ests,count_ests[,rnacols],prot_ses,count_ses[,rnacols])
datafuns = list(riboseq =  make_ribodata , rnaseq =  make_rnadata)
#precompute the data
gdatas = lapply(datafuns,function(datafun){
	lapply(rownames(prot_ests)%>%setNames(.,.),function(gene){
		simdata = datafun(gene)
	})
})

###For selecting subsets of genes
contrdf<-readRDS('data/contrdf.rds')
siggenes = rownames(count_ests)%>%intersect(rownames(prot_ests))
siggenes = rownames(count_ests)%>%intersect(rownames(prot_ests))#%>%intersect(contrdf%>%filter(adj_pval<0.05)%>%.$gene_name)
#
testgrpsize = if(!getwd()%>%str_detect('Users/dharnet/')) Inf else 50
timepoints = ribocols%>%str_replace('_ribo_\\d','')
xtailfoldchange<-Sys.glob('pipeline/xtail/xtail_*')%>%setNames(timepoints[-1])%>%
	map_df(.id='time',fread)%>%
	mutate(changetype='translational_xtail')

teupgenes=xtailfoldchange%>%filter(time=='P0_ribo')%>%filter(adj_p_value<0.05,log2fc > 1)%>%.$gene_name
tedowngenes = xtailfoldchange%>%filter(time=='P0_ribo')%>%filter(adj_p_value<0.05,log2fc < -1)%>%.$gene_name
dec_genes = (prot_ests[,1]-prot_ests[,5])%>%sort%>%names%>%
	intersect(siggenes)%>%intersect(c(downtegenes,uptegenes))%>%head(testgrpsize)
incr_genes = (prot_ests[,1]-prot_ests[,5])%>%sort%>%names%>%
	intersect(siggenes)%>%intersect(c(downtegenes,uptegenes))%>%tail(testgrpsize)
notedec_genes = (prot_ests[,1]-prot_ests[,5])%>%sort%>%names%>%
	intersect(siggenes)%>%setdiff(c(downtegenes,uptegenes))%>%head(testgrpsize)
noteincr_genes = (prot_ests[,1]-prot_ests[,5])%>%sort%>%names%>%
	intersect(siggenes)%>%setdiff(c(downtegenes,uptegenes))%>%tail(testgrpsize)
testgenes = c(dec_genes,incr_genes,notedec_genes,noteincr_genes)%>%unique

##for comparing half life results
#get mcshane half life data too
mcshanedf<-fread('ext_data/mcshane_etal_2016_S1.csv')
#
mcshanethalfs<-mcshanedf%>%select(2,38,41)%>%set_colnames(c('gene','half_life','McShane_deg_cat'))
#
mcshanethalfs$half_life%<>%str_replace('> 300','300')%>%as.numeric
mcshanethalfs%<>%filter(half_life<299)
mcshanethalfs$half_life %<>% {./24}


