if(!exists('prot_ests'))source('src/becker_modelling_prep.R')
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
					# if(!is.null(opt) & !is.null(opt$theta_tilde)){
						# opt$cov <- cov(opt$theta_tilde)
						# opt$theta_tilde <- NULL
					# }
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
message(paste0(sep='\n',capture.output(genebmodels%>%table)))
bmodelopts%>%map('riboseq')%>%map('msdev')%>%map_lgl(is.null)%>%table
bmodelopts%>%map('riboseq')%>%map('production')%>%map_lgl(is.null)%>%table

notmsdevgenes = modeltestdf%>%filter(data=='riboseq')%>%group_by(gene)%>%filter(n()==5)%>%filter(2+BIC[model=='production']<BIC[model=='msdev'])%>%.$gene%>%unique

################################################################################
########DIagnostics
################################################################################
if(FALSE){	
	# modeltestdf%>%filter(gene==filteredgenes[1])%>%filter(data=='riboseq')%>%filter(model=='production')%>%t
	# modeltestdf%>%filter(gene==filteredgenes[1])%>%filter(data=='riboseq')%>%filter(model=='msdev')%>%t

	# modeltestdf%>%filter(passtest)%>%group_by(gene)%>%slice(which.min(BIC))%>%.$data%>%table

	# modeltestdf%>%filter(gene%in%filteredgenesold)%>%filter(data=='riboseq')%>%filter(best)%>%filter(passtest)%>%.$model%>%table
	# modeltestdf%>%filter(gene%in%filteredgenes)%>%filter(data=='riboseq')%>%filter(best)%>%filter(passtest)%>%.$model%>%table


	# modeltestdf%>%filter(gene%in%'Myl6b')%>%filter(data=='riboseq')%>%as.data.frame

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
		filter(conv==1)%>%
		.$gene%>%unique


	# lowfiltergenes = modeltestdf%>%
	# 	filter(data=='riboseq')%>%
	# 	filter(passtest[model=='production'],best[model=='production'])%>%
	# 	.$gene%>%
	# 	unique

estimate = bmodelopts[modfilteredgenes]%>%map_df(.id='gene',possibly(NULL,.f=.%>%.[['riboseq']]%>%.[['production']]%>%.$par%>%unlist))
estimate$l_pihalf%>%txtdensity
estimate$l_st%>%txtdensity

outlgenes = 	estimate%>%filter(!between(l_pihalf,-5,5))%>%.$gene
nonoutlgenes = 	estimate%>%filter(!gene%in%outlgenes)%>%.$gene
genenonconv = bmodelopts%>%map(2)%>%map('production')%>%discard(is.null)%>%map_dbl('return_code')%>%`==`(70)
convgenes = genenonconv%>%.[!.]%>%names
nonoutconvgenes = nonoutlgenes%>%intersect(convgenes)
estimate = bmodelopts[nonoutconvgenes]%>%map_df(.id='gene',possibly(NULL,.f=.%>%.[['riboseq']]%>%.[['production']]%>%.$par%>%unlist))
estimate$l_pihalf%>%txtdensity
estimate$l_st%>%txtdensity

# bmodelopts[convgenes]%>%map(2)%>%map('production')%>%discard(is.null)%>%map_dbl('return_code')
	allfailgenes = modeltestdf%>%group_by(gene)%>%filter(all(!passtest))

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
}