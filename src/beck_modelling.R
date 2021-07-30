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
	linear = bmodel_linear
	# msdev = bmodel_msdev
	
)


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

timepoints = ribocols%>%str_replace('_ribo_\\d','')
xtailfoldchange<-Sys.glob('pipeline/xtail/xtail_*')%>%setNames(timepoints[-1])%>%
	map_df(.id='time',fread)%>%
	mutate(changetype='translational_xtail')

teupgenes=xtailfoldchange%>%filter(time=='P0_ribo')%>%filter(adj_p_value<0.05,log2fc > 1)%>%.$gene_name
tedowngenes = xtailfoldchange%>%filter(time=='P0_ribo')%>%filter(adj_p_value<0.05,log2fc < -1)%>%.$gene_name

protscl = prot_ests%>%colMedians
prot_ests = sweep(prot_ests,2,'-',STAT=protscl)

countrscl = count_ests%>%colMedians
count_ests = sweep(count_ests,2,'-',STAT=countrscl)


#get mcshane half life data too
mcshanedf<-fread('ext_data/mcshane_etal_2016_S1.csv')
#
mcshanethalfs<-mcshanedf%>%select(2,38,41)%>%set_colnames(c('gene','half_life','McShane_deg_cat'))
#
mcshanethalfs$half_life%<>%str_replace('> 300','300')%>%as.numeric
mcshanethalfs%<>%filter(half_life<299)
mcshanethalfs$half_life %<>% {./24}


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
	sampdata$ribooffset = c(0, -0.1169, -0.4907, -0.5889, -0.7206)
	sampdata
}
make_ribodata = function(g)make_standata(g,prot_ests,count_ests[,ribocols],prot_ses,count_ses[,ribocols])
make_rnadata = function(g)make_standata(g,prot_ests,count_ests[,rnacols],prot_ses,count_ses[,rnacols])

contrdf<-readRDS('data/contrdf.rds')

siggenes = rownames(count_ests)%>%intersect(rownames(prot_ests))
siggenes = rownames(count_ests)%>%intersect(rownames(prot_ests))#%>%intersect(contrdf%>%filter(adj_pval<0.05)%>%.$gene_name)



testgrpsize = if(!getwd()%>%str_detect('Users/dharnet/')) Inf else 50

dec_genes = (prot_ests[,1]-prot_ests[,5])%>%sort%>%names%>%
	intersect(siggenes)%>%intersect(c(downtegenes,uptegenes))%>%head(testgrpsize)
incr_genes = (prot_ests[,1]-prot_ests[,5])%>%sort%>%names%>%
	intersect(siggenes)%>%intersect(c(downtegenes,uptegenes))%>%tail(testgrpsize)
notedec_genes = (prot_ests[,1]-prot_ests[,5])%>%sort%>%names%>%
	intersect(siggenes)%>%setdiff(c(downtegenes,uptegenes))%>%head(testgrpsize)
noteincr_genes = (prot_ests[,1]-prot_ests[,5])%>%sort%>%names%>%
	intersect(siggenes)%>%setdiff(c(downtegenes,uptegenes))%>%tail(testgrpsize)


testgenes = c(dec_genes,incr_genes,notedec_genes,noteincr_genes)%>%unique

datafuns = list(riboseq =  make_ribodata , rnaseq =  make_rnadata)
	

#precompute the data
gdatas = lapply(datafuns,function(datafun){
	lapply(testgenes%>%setNames(.,.),function(gene){
		simdata = datafun(gene)
	})
})


# model=models$msdev
model=models$production
datafun='riboseq'
gene=filteredgenes[1]
simdata  = gdatas[[datafun]][[gene]]

# safely(rstan::optimizing)(models$msdev,data=simdata,init=initvals,as_vector=F,hessian=TRUE)$result$par$prot
# safely(rstan::optimizing)(models$production,data=simdata,init=initvals,as_vector=F,hessian=TRUE)$result$par$prot
# reopt <- safely(rstan::optimizing)(model,data=simdata,init=initvals,as_vector=F,hessian=TRUE)

#file.remove('data/bmodelopts.rds')
# if(!file.exists(here('data/bmodelopts.rds'))){

	bmodelopts <- mclapply(mc.cores=10,
		# testgenes%>%setNames(.,.)%10.[filteredgenes],function(gene){
		testgenes%>%setNames(.,.)%>%.[TRUE],function(gene){
		cat('.')
		lapply(names(datafuns)%>%setNames(.,.),function(datafun){
			simdata  = gdatas[[datafun]][[gene]]
				reopt = lapply(models,function(model){
					for(i in 1:10){
						initvals=list()
						initvals$l_pihalf <- array(log(0.5),1)
						initvals$lribo = array(simdata$lSeqmu,c(1,length(simdata$lSeqmu)))
						initvals$l_st = array(simdata$lMSmu[[1]]-simdata$lSeqmu,1)
						initvals$lprot0 = array(simdata$lMSmu[[1]],1)
						initvals$msdev = array(0,c(1,5))
						reopt <- safely(rstan::optimizing)(model,data=simdata,init=initvals,as_vector=F,hessian=TRUE)
						reopt
						if(!is.null(reopt$result)) break
					}
					opt = reopt$result
					opt
				})
		})
	})

	# bmodelopts
	# bmodelopts = bmodelopts[unique(names(bmodelopts))]
	saveRDS(bmodelopts,here('data/bmodelopts.rds'))

# }else{
# 	bmodelopts<-readRDS(here('data/bmodelopts.rds'))
# }


modname=names(models)[1]
gene=names(bmodel)
#pre-compute the number of dfs in the model (note that the gene and datasource don't effect this, hence
#the two magic numbers)
modname = names(models)
n_dflist = lapply(names(models)%>%setNames(.,.),function(modname){
	 bmodelopts[[1]][[1]][[modname]]$par[get_stanpars(models[[modname]])]%>%unlist%>%length
})

# stopifnot('msdev' %in% names(bmodelopts[[1]][[1]]))
# opt = bmodelopts[[gene]][[datafun]][['msdev']]
# if(!file.exists(here('data/modeltestdf.rds'))){

	#now let's do a chi squared test
	# model_tests = mclapply(mc.cores=20,names(bmodelopts)%>%setNames(.,.),safely(function(gene){
	modeltestdf = mclapply(mc.cores=10,names(bmodelopts)%>%setNames(.,.),possibly(NULL,.f=function(gene){
	# modeltestdf = mclapply(mc.cores=1,names(bmodelopts)%>%setNames(.,.),identity(function(gene){
	# modeltestdf = map_df(.id='gene',names(bmodelopts)%>%setNames(.,.),function(gene){
		cat('.')
		map_df(.id='data',names(gdatas)%>%setNames(.,.),function(datafun){
			map_df(.id='model',names(models)%>%setNames(.,.),function(modname){
				opt = bmodelopts[[gene]][[datafun]][[modname]]
				gdata = gdatas[[datafun]][[gene]]
				#
				perrors =(log(opt$par$prot) - gdata$lMSmu)
				w_perrors = perrors/gdata$lMSsigma
				rerrors = (opt$par$lribo - gdata$lSeqmu)
				w_rperrors = rerrors/gdata$lSeqsigma
				n_df = n_dflist[[modname]]
				errorsum = sum(c(w_perrors,w_rperrors)^2)
				BIC = -2*opt$value+log(5)*n_df
				pval = 1 - pchisq(errorsum,n_df)
				sumstats = c(BIC=BIC,pval=pval,residuals = w_perrors,errorsum=errorsum)
				# sumstats -> tmpsumstatsp
				# tmpsumstatsp
				# browser()
				sumstats
			})
		})
	}))%>%bind_rows(.id='gene')

	modeltestdf <-modeltestdf%>%
		mutate(passtest = ! (pval < 0.05))%>%
		group_by(gene,data)%>%
		mutate(best=BIC==min(BIC))

	modeltestdf%<>%rowwise%>%mutate(sumresid = sum(residuals1^2+residuals2^2+residuals3^2+residuals4^2+residuals5^2))%>%
		group_by(gene,data)
	# modeltestdf%<>%distinct(gene,data,model,.keep_all=TRUE)
	# saveRDS(modeltestdf,here('data/modeltestdf.rds'))
# }else{
# 	modeltestdf<-readRDS(here('data/modeltestdf.rds'))
# }

modeltestdf <-modeltestdf%>%
	# filter(!model=='production_offset')%>%
	# filter(!model=='msdev')%>%
	# filter(!model=='production')%>%
	# mutate(model = ifelse(model=='production_offset','production',model))%>%
	mutate(passtest = ! (pval < 0.05))%>%
	group_by(gene,data)%>%
	mutate(best=BIC==min(BIC))

modeltestdf%>%filter(gene==filteredgenes[1])%>%filter(data=='riboseq')%>%filter(model=='production')%>%t
modeltestdf%>%filter(gene==filteredgenes[1])%>%filter(data=='riboseq')%>%filter(model=='msdev')%>%t
modeltestdf%>%filter(data=='riboseq')%>%filter(best)%>%filter(passtest)%>%.$model%>%table
modeltestdf%>%filter(passtest)%>%group_by(gene)%>%slice(which.min(BIC))%>%.$data%>%table

modeltestdf%>%filter(gene%in%filteredgenesold)%>%filter(data=='riboseq')%>%filter(best)%>%filter(passtest)%>%.$model%>%table
modeltestdf%>%filter(gene%in%filteredgenes)%>%filter(data=='riboseq')%>%filter(best)%>%filter(passtest)%>%.$model%>%table


diffgenes = filteredgenes%>%setdiff(filteredgenesold)
diffgenes[1]


modeltestdf%>%filter(gene%in%'Myl6b')%>%filter(data=='riboseq')%>%as.data.frame

# filteredgenesold <- filteredgenes
modfilteredgenes = modeltestdf%>%
	filter(data=='riboseq')%>%
	filter(
		passtest[model=='production'],
		!passtest[model=='degredation'],
		!passtest[model=='linear'],
		!passtest[model=='stationary'],
	# !passtest[model=='accumulation'],
		best[model=='production'])%>%
	.$gene%>%unique

modfilteredgenes%>%length
'Myl6b'%in% modfilteredgenes

deggenes = modeltestdf%>%
	filter(data=='riboseq')%>%
	filter(model=='degredation')%>%
	filter(passtest[model=='degredation'])%>%
	filter(best)%>%.$gene

# lowfiltergenes = modeltestdf%>%
# 	filter(data=='riboseq')%>%
# 	filter(passtest[model=='production'],best[model=='production'])%>%
# 	.$gene%>%
# 	unique

allfailgenes = modeltestdf%>%group_by(gene)%>%filter(all(!passtest))

estimate = bmodelopts[modfilteredgenes]%>%map_df(.id='gene',possibly(NULL,.f=.%>%.[['riboseq']]%>%.[['production']]%>%.$par%>%unlist))
estimate$l_pihalf%>%txtdensity
estimate$l_st%>%txtdensity

txtplot(estimate$l_pihalf,estimate$l_st)
estimate%>%filter(l_pihalf<5)%>%{txtplot(.$l_pihalf,.$l_st)}
estimate%>%filter(gene=='Myl6b')%>%as.data.frame

# highpigenes = estimate%>%filter(gene%in%modfilteredgenes)%>%filter(l_pihalf%>%`>`(5))%>%.$gene
# low_stgenes = estimate%>%filter(gene%in%modfilteredgenes)%>%filter(l_st%>%`<`(-5))%>%.$gene
filteredgenes = modfilteredgenes
# filteredgenes = setdiff(filteredgenes,highpigenes)%>%setdiff(low_stgenes)

'Myl6b'%in% highpigenes

estimate%>%filter(gene%in%filteredgenes)%>%.$l_pihalf%>%txtdensity
estimate%>%filter(gene%in%filteredgenes)%>%.$l_st%>%txtdensity

filteredgenes%>%length

pihalftest = estimate%>%filter(gene%in%filteredgenes)%>%
	inner_join(mcshanethalfs)%>%
	filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	{quicktest(.$l_pihalf,log(.$half_life));.}%>%
	{cor.test(.$l_pihalf,log(.$half_life))}%>%tidy
pihalftest

# filteredgenes=filteredgenesold

{

jointmodel1te = stan_model(here('src/Archive/Stan/becker_proda_oneKs.stan'))

stopifnot(length(filteredgenes)>100)

get_comb_initvals <- function(bestfitinits){
  combinitvals <- lapply(names(bestfitinits[[1]])%>%setNames(.,.),function(argind){
    bestfitinits%>%
      map(argind)%>%
      setNames(.,seq_along(.))%>%
      do.call(what=partial(abind::abind,along=1))
  })
  combinitvals	
}

combinitvals <- bmodelopts[filteredgenes]%>%map('riboseq')%>%map('production')%>%map('par')%>%get_comb_initvals
jointdata = datafuns$riboseq(filteredgenes)
jointdata$G = nrow(jointdata$lMSmu)
combinitvals$lKs = combinitvals$Ks%>%mean
#
combinitvals$l_pihalf%>%txtdensity
combinitvals$l_st%>%txtdensity

jopt = rstan::optimizing(jointmodel1te,data=jointdata,init=combinitvals,as_vector=F,save_iterations=TRUE,hessian=F,verbose=T)


pihalftest = jopt$par%>%.$l_pihalf%>%setNames(filteredgenes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	# filter(l_pihalf<3)%>%
	{quicktest(.$l_pihalf,log(.$half_life))}%>%tidy

message(paste0(sep='\n',capture.output(pihalftest)))

}





stop()
nedgenes <- (mcshanethalfs%>%filter(McShane_deg_cat=='NED')%>%.$gene)

modeltestdf%>%group_by(gene)%>%	
	filter(data=='riboseq')%>%
	filter(best)%>%
	mutate(isned = gene %in% nedgenes)%>%
	mutate(sumresid = sum(residuals1^2+residuals2^2+residuals3^2+residuals4^2+residuals5^2))%>%
	{split(.$sumresid,.$isned)}%>%
	{t.test(.[['TRUE']],.[['FALSE']])}
	# {table(.$isned,.$model==)}%>%fisher.test


codreltests = pihalftest%>%{paste0('rho = ',round(.$estimate,3),'\n','pval = ',ifelse(.$p.value > 0.001,round(.$p.value,2),format(.$p.value,format='e',digits=2)))}


#now plot
plotfile<- here(paste0('plots/jte_indiv_v_mcshane_pihalf','.pdf'))
dir.create(dirname(plotfile))
pdf(plotfile)
p = 
	jopt$par%>%.$l_pihalf%>%setNames(filteredgenes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	# filter(McShane_deg_cat!='NED')%>%
	# filter()
	# ggplot(.,aes(log(half_life),l_pihalf,color=McShane_deg_cat))+
	ggplot(.,aes(log(half_life),l_pihalf))+
	geom_point()+
	facet_grid(McShane_deg_cat~.)+
	scale_color_discrete(name='McShane_deg_cat')+
	scale_x_continuous(paste0('log2(Half Life) (McShane et al)'))+
	scale_y_continuous(paste0('Joint TE Model - estimated log(Half Life)'))+
	ggtitle(paste0('Measured Half Lives vs Estimated'),sub=codreltests)+
	geom_smooth(method='lm',color=I('black'))+
	theme_bw()
p
dev.off()
normalizePath(plotfile)





################################################################################
########Okay so the joint model with a point estimate works okay, what about
########If we optimize a hiearach model?
################################################################################

jointmodel_hierach = stan_model(here('src/Archive/Stan/becker_proda_jhierarch.stan'))

mu_lks = log(unique(jopt$par$Ks))
sd_lks = sd(combinitvals$l_st + (log(log(2))-combinitvals$l_pihalf))

jopth = rstan::optimizing(jointmodel_hierach,data=jointdata,init=c(mu_lks =mu_lks ,sd_lks = sd_lks,combinitvals),as_vector=F,save_iterations=TRUE,hessian=F,verbose=T)

jopth$par$mu_lks
jopth$par$sd_lks

 jopt$par%>%.$l_pihalf%>%setNames(filteredgenes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	# filter(l_pihalf<4)%>%
	{quicktest(.$l_pihalf,log(.$half_life));.}%>%
	{cor.test(.$l_pihalf,log(.$half_life))}%>%tidy

 jopth$par%>%.$l_pihalf%>%setNames(filteredgenes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	filter(l_pihalf<4)%>%
	{quicktest(.$l_pihalf,log(.$half_life));.}%>%
	{cor.test(.$l_pihalf,log(.$half_life))}%>%tidy

#hmmmmm, thi sdoesn't work that well...

lineargenes = modeltestdf%>%filter(data=='riboseq',model=='linear',best,passtest)%>%.$gene



pihalftest = estimate%>%filter(gene%in%filteredgenes)%>%filter(l_pihalf%>%between(-5,5),l_st%>%between(-5,5))%>%
	.$l_pihalf%>%setNames(filteredgenes)%>%enframe('gene','l_st')%>%
	mutate(gene=gene%>%tolower)%>%
	inner_join(mcshanethalfs%>%mutate(gene=gene%>%tolower))%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	{quicktest(.$l_st,log(.$half_life))}%>%tidy
pihalftest

'becker_proda_lKsfix.stan'


################################################################################
########If I optimize on the linear genes with my 
################################################################################
	

################################################################################
########
################################################################################
	
#that gives us suprisingly many!
#559so far....

hessianses = bmodelopts[filteredgenes]%>%map_df(.id='gene',possibly(NULL,.f=.%>%.[['riboseq']]%>%.[['production']]%>%.$hessian%>%{sqrt(diag(solve(-.)))}))
allhessianses = bmodelopts%>%map_df(.id='gene',possibly(NULL,.f=.%>%.[['riboseq']]%>%.[['production']]%>%.$hessian%>%{sqrt(diag(solve(-.)))}))
#These guys are weirdly multimodel...
#weird tripartate shape to these - seems like a lot of the time, one or the other of the parameters has it's se waaaaay down.
hessianses%>%filter(is.finite(l_st.1),is.finite(l_pihalf.1))%>%{txtplot(log10(.$l_pihalf.1),log10(.$l_st.1))}

allhes%>%filter(is.finite(l_st.1),is.finite(l_pihalf.1))%>%{txtplot(log10(.$l_pihalf.1),log10(.$l_st.1))}



allestimate = bmodelopts%>%map_df(.id='gene',possibly(NULL,.f=.%>%.[['riboseq']]%>%.[['production']]%>%.$par%>%unlist))

estimate_rnaseq = bmodelopts%>%map_df(.id='gene',possibly(NULL,.f=.%>%.[['rnaseq']]%>%.[['production']]%>%.$par%>%unlist))

allTEchangedf

##YAY - RNA estimate Ks correlates a bit, if we filter out the edge cases
estimate_rnaseq%>%
	filter(gene%in%lowfiltergenes)%>%
	filter(!(gene%in%teupgenes|gene%in%tedowngenes))%>%
	inner_join(tedf,by=c('gene'='gene_name'))%>%
	filter(abs(log(Ks))<5)%>%
	{quicktest(.$TE,log(.$Ks))}


#odd relationship of the estimates, pihalf is often going to something very small
txtplot(estimate%>%.$l_pihalf,estimate%>%.$l_st)
txtplot(allestimate%>%.$l_pihalf,allestimate%>%.$l_st)

#What about the ses and the estimates
#okay so very clear relationship where the low estimates have lower ses for the pihalf
estimate%>%left_join(hessianses,suffix=c('','_se'))%>%{txtplot(.$l_pihalf,log(.$l_pihalf.1))}
#also very clear relationship but two horned, very low and very high estimates have high se 
estimate%>%left_join(hessianses,suffix=c('','_se'))%>%{txtplot(.$l_st,log(.$l_st.1))}

estimate%>%filter(l_pihalf<5,abs(l_st)<5)

bmodelopts[filteredgenes]%>%map_df(.id='gene',possibly(NULL,.f=.%>%.[['riboseq']]%>%.[['production']]%>%.$hessian%>%{sqrt(diag(solve(-.)))}))%>%.$l_pihalf.1%>%na.omit%>%log%>%txtdensity



filteredgenes = modeltestdf%>%filter(data=='riboseq')%>%filter(passtest[model=='production'],!passtest[model=='degredation'],!passtest[model=='linear'],best[model=='production'])%>%.$gene%>%unique
lowfiltergenes = modeltestdf%>%filter(data=='riboseq')%>%filter(passtest[model=='production'],best[model=='production'])%>%.$gene%>%unique

estimate%>%inner_join(mcshanethalfs)%>%filter(l_pihalf<5,abs(l_st) < 5)%>%{quicktest(.$l_pihalf,log(.$half_life))}

allestimate%>%
	filter(gene %in% (lowfiltergenes))%>%
	# filter(gene %in% (filteredgenes))%>%
	inner_join(mcshanethalfs)%>%
	filter(abs(l_pihalf)<5,abs(l_pihalf) <  5)%>%
	# filter(McShane_deg_cat=='ED')%>%
	{quicktest(.$l_pihalf,log(.$half_life))}






stop()
residlistcol = modeltestdf%>%ungroup%>%select(matches('residuals'))%>%as.matrix%>%t%>%as.data.frame
modeltestdf$residuals = residlistcol%>%as.list
modeltestdf%<>%select(-matches('residuals\\d'))



modeltestdf%>%filter(data=='rnaseq')%>%filter(best)%>%.$model%>%table

modeltestdf%>%filter(data=='riboseq')%>%filter(best)%>%.$model%>%table

modeltestdf%>%filter(data=='rnaseq')%>%summarise(reject=all(!passtest))%>%.$reject%>%table
modeltestdf%>%filter(data=='riboseq')%>%summarise(reject=all(!passtest))%>%.$reject%>%table

#okay so testing the residuals gives us nothing to go on, effectively al of the time
resid_tests = map_df(.id='gene',names(bmodelopts)%>%setNames(.,.),function(igene){
	gdf = modeltestdf%>%filter(gene==igene,best)
	tidy(t.test(
	gdf%>%filter(data=='riboseq')%>%.$residuals%>%.[[1]],
	gdf%>%filter(data=='rnaseq')%>%.$residuals%>%.[[1]]
))
})
resid_tests%>%.$p.value%>%`<`(0.05)%>%table

modeltestdf%>%filter(best,data=='riboseq')%>%.$model%>%table
modeltestdf%>%filter(best,data=='rnaseq')%>%.$model%>%table


rnagoodtest = modeltestdf%>%filter(data=='rnaseq',model=='production',best,passtest)%>%.$gene

Ksvals = bmodelopts[rnagoodtest]%>%map_dbl(function(opt){
	opt[['riboseq']][['production']]$par$Ks%>%log
})




meantes = count_ests[,colnames(count_ests)%>%str_subset('TE')]%>%rowMeans

enframe(Ksvals,'gene_name','rna_Ks')%>%left_join(meantes%>%enframe('gene_name','TE'))%>%
	{cor.test(.$rna_Ks,.$TE)}

enframe(Ksvals,'gene_name','rna_Ks')%>%left_join(meantes%>%enframe('gene_name','TE'))%>%
	{quicktest(.$rna_Ks,.$TE)}


gene='Flna'

modeltestdf%<>%group_by(gene,data)%>%mutate(bestmodel = model[BIC==min(BIC)])


nonlgenes = modeltestdf%>%group_by(gene)%>%filter(data=='riboseq')%>%
	filter(bestmodel=='production')%>%
	summarise(bicdiff = BIC[model=='linear']-BIC[model=='production'])%>%
	arrange(bicdiff)%>%.$gene%>%tail(10)
modeltestdf%>%filter(gene%in%nonlgenes)


igene=nonlgenes[7]
modeltestdf%>%filter(gene%in%igene)%>%arrange(BIC)

modeltestdf%>%group_by(gene)%>%
	

# bigdiffgenes = modeltestdf%>%filter(data=='riboseq')%>%arrange(BIC)%>%summarise(BICdiff = BIC[2]-BIC[1],bestmodel = model[1],secondmodel=model[2])%>%arrange(-BICdiff)
# 
	# mutate(gbestmodel=sample(bestmodel[BIC==min(BIC)])
igene = filteredgenes%>%sample(1)
ntps = 1:5
modelcols<-c('data'='black','degredation'='green','production'='blue','stationary'='purple','linear'='lightblue')
models2plot<-c('degredation','production','linear')

{
datadf = map_df(.id='data' ,names(datafuns)%>%setNames(.,.),function(datafun){
	se = datafuns[[datafun]](igene)$lMSsigma%>%.[1,]
	mu=datafuns[[datafun]](igene)$lMSmu%>%.[1,]
	msdf=tibble(y=mu)%>%mutate(lower=y-1.96*se,upper=y+1.96*se)%>%mutate(assay='MS',model='data',time=ntps)
	se = datafuns[[datafun]](igene)$lSeqsigma%>%.[1,]
	mu=datafuns[[datafun]](igene)$lSeqmu%>%.[1,]
	seqdf = tibble(y=mu)%>%mutate(lower=y-1.96*se,upper=y+1.96*se)%>%mutate(assay='seq',model='data',time=ntps)
	datdf = rbind(msdf,seqdf)
	datdf
})

modelname='production'

modelms = bmodelopts[[igene]]%>%map_df(.id='data',.%>%map_df(.id='model',.%>%.$par%>%.$prot%>%.[1,]%>%log%>%{tibble(y=.)}%>%mutate(assay='MS',time=ntps)))
modelribo = bmodelopts[[igene]]%>%map_df(.id='data',.%>%map_df(.id='model',.%>%.$par%>%.$lribo%>%.[1,]%>%{tibble(y=.)}%>%mutate(assay='seq',time=ntps)))
modelms%<>%filter(model%in%models2plot)
ggdf = bind_rows(datadf,modelms)

ggplot(ggdf,aes(y=y,x=time,color=model,ymin=lower,ymax=upper))+geom_line(linetype=2)+
	geom_errorbar(width=I(0.2))+
	scale_color_manual(values=modelcols)+
	ggtitle(paste0('model fit',igene))+
	facet_grid(assay~data,scale='free')

}


# #now let's do a chi squared test
# model_tests = lapply(names(bmodelopts)%>%setNames(.,.),function(gene){
# 	lapply(names(datafuns)%>%setNames(.,.),function(datafun){
# 		lapply(names(models)%>%setNames(.,.),function(modname){
gene=testgenes[1]
datafun=names(datafuns)[1]

			opt = bmodelopts[[gene]][[datafun]][[modname]]
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
# 		})
# 	})
# })


estimate_rnaseq%>%filter(gene%in%filteredgenes)%>%inner_join(tedf,by=c('gene'='gene_name'))%>%filter(abs(log(Ks))<5)%>%{quicktest(.$TE,log(.$Ks))}



################################################################################
########Get real data for comparison
################################################################################
	
countlinearTEs <- get_contrast_cis(
	bestonlycountebayes,
	t(countonly_pred_te_design['TE',,drop=F])
)%>%select(protein_id = uprotein_id,logFC,CI.L,CI.R,adj.P.Val)

countlinearTEs%<>%safe_left_join(metainfo%>%distinct(protein_id,gene_name),by='protein_id')

mcshanedf<-fread('ext_data/mcshane_etal_2016_S1.csv')
#
mcshanethalfs<-mcshanedf%>%select(2,38,41)%>%set_colnames(c('gene_name','half_life','McShane_deg_cat'))
#
mcshanethalfs$half_life%<>%str_replace('> 300','300')%>%as.numeric
mcshanethalfs%<>%filter(half_life<299)
mcshanethalfs$half_life %<>% {./24}

mcshanethalfs%>%head

tedf = countpred_df%>%filter(contrast%>%str_detect('TE'))%>%group_by(gene_name)%>%summarise(TE=mean(logFC))





isource='ribo'
ipar='l_pihalf'

for(isource in names(fitlist)){
for(ipar in ipar ){

#now plot
plotfile<- here(paste0('plots/',source,'_',ipar,'_indiv_v_mcshane_pihalf','.pdf'))
dir.create(dirname(plotfile))
pdf(plotfile)
p = stansumdata%>%
	filter(source==isource)%>%
	filter(par%>%str_detect(ipar))%>%inner_join(metainfo%>%distinct(gene_name,uprotein_id))%>%inner_join(mcshanethalfs)%>%
	filter(McShane_deg_cat=='ED')%>%
	# filter(between(log2(half_life),-3,3))%>%
	# filter(uprotein_id%in%confuproteinids)%>%
	ggplot(.,aes(log2(half_life),estimate,color=McShane_deg_cat))+
	geom_point()+
	scale_color_discrete(name='McShane_deg_cat')+
	scale_x_continuous(paste0('log2(Half Life) (McShane et al)'))+
	scale_y_continuous(paste0('Estimated ',ipar,' individual model fits'))+
	ggtitle(paste0('Measured Half Lives vs Estimated'))+
	geom_smooth(method='lm')+
	theme_bw()
p
dev.off()
normalizePath(plotfile)
p$data%>%head
cor.test(p$data$half_life,p$data[,'estimate'],method='spearman')
cor.test(p$data$half_life,p$data[,'estimate'])
txtplot(p$data$half_life,p$data[,'estimate'],width=100)

}}

chrmcounts = Sys.glob('../cortexomics/pipeline/star/data/*/*.bam') %>% str_subset(neg=T,'transcript')%>%setNames(.,basename(.))%>% map_dbl(.%>%{x=.;system(str_interp('samtools view -c ${.} chrM'),intern=TRUE)}%>%as.numeric)

chrmcounts%>%enframe%>%mutate(name=str_replace(name,'\\.bam',''))%>%filter




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
beckermodel <- rstan::stan_model('src/Archive/Stan/becker_proda.stan')
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