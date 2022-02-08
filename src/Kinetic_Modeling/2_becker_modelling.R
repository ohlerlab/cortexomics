if(!exists('prot_ests'))source('src/Kinetic_Modeling/becker_modelling_prep.R')
model=models[['production']]
gene=rownames(prot_ests)[1]
datafun=names(datafuns)[[2]]
# file.remove('data/bmodelopts.rds')
if(!file.exists(here('data/bmodelopts.rds'))){

	bmodelopts <- mclapply(mc.cores=10,
		# testgenes%>%setNames(.,.)%10.[filteredgenes],function(gene){
		rownames(prot_ests)%>%setNames(.,.),function(gene){
		cat('.')
		lapply(names(datafuns)%>%setNames(.,.),function(datafun){
			simdata  = gdatas[[datafun]][[gene]]
				reopt = lapply(models,function(model){
					for(i in 1:20){
						# browser()
						getdraws = model@model_name=='becker_proda'
						idraws = ifelse(getdraws,500,0)
						# idraws=0
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

}else{
	bmodelopts<-readRDS(here('data/bmodelopts.rds'))
}

# bmodelopts%>%head(100)%>%map(2)%>%map('production')%>%map_dbl('return_code')%>%table
# bmodeloptsold<-bmodelopts

# modname=names(models)[1]
# gene=names(bmodel)
#pre-compute the number of dfs in the model (note that the gene and datasource don't effect this, hence
#the two magic numbers)
# modname = names(models)

# stopifnot('msdev' %in% names(bmodelopts[[1]][[1]]))
# opt = bmodelopts[[gene]][[datafun]][['msdev']]
# file.remove('data/modeltestdf.rds')
if(!file.exists(here('data/modeltestdf.rds'))){
	#now let's do a chi squared test
	# model_tests = mclapply(mc.cores=20,names(bmodelopts)%>%setNames(.,.),safely(function(gene){
	n_dflist = lapply(names(models)%>%setNames(.,.),function(modname){
	 		bmodelopts[[1]][[1]][[modname]]$par[get_stanpars(models[[modname]])]%>%unlist%>%length
		})
	modeltestdf = mclapply(mc.cores=1,names(bmodelopts)%>%setNames(.,.),possibly(NULL,.f=function(gene){
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
				sumstats = c(BIC=BIC,pval=pval,residuals = w_perrors,errorsum=errorsum,conv=opt$return_code==0)
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
	saveRDS(modeltestdf,here('data/modeltestdf.rds'))
	genebmodels = modeltestdf%>%filter(data=='riboseq')%>%filter(best)%>%filter(passtest)%>%
		{setNames(.$model,.$gene)}
	message(paste0(sep='\n',capture.output(genebmodels%>%table)))
}else{
	modeltestdf<-readRDS(here('data/modeltestdf.rds'))
	genebmodels = modeltestdf%>%filter(data=='riboseq')%>%filter(best)%>%filter(passtest)%>%
		{setNames(.$model,.$gene)}

}

gns_by_mod <- genebmodels%>%{split(names(.),.)}

notmsdevgenes = modeltestdf%>%
	filter(data=='riboseq')%>%
	group_by(gene)%>%filter(n()==5)%>%
	filter(BIC[model=='production']<BIC[model=='msdev'])%>%.$gene%>%unique

exgn='Satb2'
allgnms<-unique(names(bmodelopts))%>%setNames(.,.)
modelvals <- map_df(.id='gene',allgnms,function(exgn){
	model=genebmodels[exgn]
	optres = bmodelopts[[exgn]][['riboseq']][[model]]
	if(is.null(optres)){return(NULL)}
	modelests = tibble(assay='prot',estimate=optres$par[['prot']]%>%as.vector)%>%
		unnest(estimate)%>%
		mutate(estimate=log2(estimate))%>%
		mutate(estimate=estimate+protscl)%>%
		mutate(time=tps,model=model,gene_name=exgn)
	modelests
})

tps = c('E13','E145','E16','E175','P0')
bmodelvals <- modelvals
bmodelvals%>%write_tsv(here('data/bmodelvals.tsv'))