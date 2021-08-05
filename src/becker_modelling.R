model=models[['production']]
gene=rownames(prot_ests)[1]
datafun=names(datafuns)[[2]]
file.remove('data/bmodelopts.rds')
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
						reopt <- safely(rstan::optimizing)(model,data=simdata,init=initvals,as_vector=F,hessian=TRUE,draws=idraws)
						reopt
						if(!is.null(reopt$result)) break
					}
					opt = reopt$result
					if(!is.null(opt) & !is.null(opt$theta_tilde)){
						opt$cov <- cov(opt$theta_tilde)
						opt$theta_tilde <- NULL
					}
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
#now plot
modelcols <- c(
	msdev = '#c64730',
	production = '#fda27c',
	linear = '#ba5e47',
	degredation = '#1d4c9b',
	stationary = '#6289cc'
)
plotfile<- here(paste0('plots/','becker_trajectory_classes','.pdf'))
pdf(plotfile)
genebmodels%>%enframe%>%
	mutate(value = factor(value,levels=names(modelcols)))%>%
	ggplot(aes(x=factor(1),fill=value))+
	stat_count(geom='bar',width=I(0.3))+
	scale_fill_manual(values=modelcols)+
	theme_bw()
dev.off()
message(normalizePath(plotfile))

gclusters = readxl::read_xlsx(here('tables/S4.xlsx'))
gclusters$cluster%<>%as_factor
levels(gclusters$cluster)%<>%rev
plotfile<- here(paste0('plots/','becker_traj_classes_by_cluster','.pdf'))
pdf(plotfile)
genebmodels%>%enframe('gene_name')%>%inner_join(gclusters)%>%
	mutate(cluster=='NA')%>%
	mutate(value = factor(value,levels=names(modelcols)))%>%
	group_by(cluster,value)%>%tally%>%
	group_by(cluster)%>%mutate(prop = n/sum(n))%>%
	ggplot(aes(x=cluster,y=prop,fill=value))+
	stat_identity(geom='bar',width=I(0.3),position='stack')+
	scale_fill_manual(values=modelcols)+
	coord_flip()+
	theme_bw()
dev.off()
message(normalizePath(plotfile))

	
source('src/Figures/Figure4/1_go_term_funcs.R')
stopifnot(exists("GTOGO"))
gid_bmodels = setNames(genebmodels,gnm2gid[[names(genebmodels)]])

clustvect<-gid_bmodels
if(!'gene_id'%in%colnames(GTOGO))GTOGO%<>%mutate(gene_id=ensembl_gene_id)
get_cluster_gos = function(clustvect){
    out=map_df(.id='ontology',onts%>%setNames(.,.),function(ont){
      lapply(unique(clustvect)%>%setNames(.,.),safely(function(clustval){
            gids <- names(clustvect)
            stopifnot(mean(names(clustvect)%in%GTOGO$gene_id)>.9)
            filtGTOGO <- GTOGO %>%filter(gene_id %in%names(clustvect))
            # browser()
             rungo(
              gids[clustvect==clustval],
              filtGTOGO,
              ont,
              algo='classic'
            )
        }))%>%map('result')%>%bind_rows(.id='cluster')
    })
    out
}
trajclass_goterms <- get_cluster_gos(gid_bmodels)
trajclass_goterms%>%filter(elimFisher<0.05,Significant>10)%>%write_tsv('tables/becker_traj_goterms.tsv')

go_comparison_plot <-function(x){x%>%map_df(.id='cluster',.%>%arrange(elimFisher)%>%head(10)%>%select(Term,elimFisher))%>%
  mutate(Term = as_factor(Term))%>%
  mutate(cluster = as_factor(cluster))%>%
  ggplot(.,aes(x=cluster,y=Term,color=-log10(elimFisher),size=-log10(elimFisher)))+
  geom_point()+
  theme_bw()
}



goplotlist = trajclass_goterms%>%split(.,.$cluster)%>%lapply(.%>%
			arrange(-elimFisher)%>%tail%>%
		 	mutate(Term=as_factor(Term))%>%
		 	{ggplot(data=.,aes(x=Term,y=-log10(elimFisher),fill=Enrichment))+
		 	stat_identity(geom='bar')+
		 	geom_text(aes(label=GO.ID,y=4,color=I('white')))+
			ggplot2::scale_fill_continuous(
			      # limits = coltranslims,
			      # breaks = round(seq(from=min(coltranslims),to=max(coltranslims),len=n_colbreaks)),
			      # labels = col_label_format,
			      high = '#231f20',low='#d1d1d1'
			    )+
		 	coord_flip()+
		 	theme(axis.text.x=element_text(vjust=-1))+
		 	theme_bw()+
		 	ggtitle(.$cluster%>%unique)})
#now plot
plotfile<- here(paste0('plots/','beck_traj_class_goplots','.pdf'))
pdf(plotfile,h=7,w=7*5)
print(ggarrange(plotlist=goplotlist,ncol=length(goplotlist)))
dev.off()
message(normalizePath(plotfile))		 


trajclass_goterms%>%split(.,.$cluster)%>%.['msdev']%>%.[[1]]%>%
	filter(Significant>10)%>%filter(elimFisher<0.05)%>%arrange(-Enrichment)%>%
	head(12)

goplotlist = trajclass_goterms%>%split(.,.$cluster)%>%.['msdev']%>%lapply(.%>%
			arrange(-Enrichment)%>%
			filter(elimFisher<0.05)%>%
			filter(Significant>10)%>%
			head(12)%>%
		 	mutate(Term=as_factor(Term))%>%
		 	{ggplot(data=.,aes(x=Term,y=-log10(elimFisher),fill=Enrichment))+
		 	stat_identity(geom='bar')+
		 	geom_text(aes(label=GO.ID,y=4,color=I('white')))+
			ggplot2::scale_fill_continuous(
			      # limits = coltranslims,
			      # breaks = round(seq(from=min(coltranslims),to=max(coltranslims),len=n_colbreaks)),
			      # labels = col_label_format,
			      high = '#231f20',low='#d1d1d1'
			    )+
		 	coord_flip()+
		 	theme(axis.text.x=element_text(vjust=-1))+
		 	theme_bw()+
		 	ggtitle(.$cluster%>%unique)})
#now plot
plotfile<- here(paste0('plots/','msdev_goplot','.pdf'))
pdf(plotfile,h=7,w=7*1)
print(ggarrange(plotlist=goplotlist,ncol=length(goplotlist)))
dev.off()
message(normalizePath(plotfile))		 



################################################################################
########Now let's do some plots of the data
################################################################################
#let's get an example

c('Bcl11b','Flna','Satb2','Nes')
bicdiffdf = modeltestdf%>%filter(data=='riboseq')%>%
	group_by(gene)%>%
	# group_slice(1)%>%
	arrange(BIC)%>%
	summarise(bicdiff = BIC[n()] - BIC[n()-1],model=model[best])
exampgenes = bicdiffdf%>%group_by(model)%>%slice(which.max(bicdiff))%>%{setNames(.$gene,.$model)}
gns4plot=exampgenes%>%setNames(.,.)
modelnms=names(models)%>%setNames(.,.)
exgn<-gns4plot[[1]]
model <- modelnms[[1]]
modelvals <- map_df(.id='gene',gns4plot,function(exgn){
	map_df(.id='model',modelnms,function(model){
	optres = bmodelopts[[exgn]][['riboseq']][[model]]
	# parsds = 1.96*sqrt(diag(solve(-optres$hessian)))
	modelests = list(
		prot=log(optres$par[['prot']])%>%as.vector,
		ribo=optres$par$lribo%>%as.vector)%>%
		enframe('assay','value')%>%
		unnest(value)%>%
		mutate(time=c(1:5,1:5))
	modelests
	})
})
modeldata <- map_df(.id='gene',gns4plot,function(exgn){
	gdata <- gdatas[['riboseq']][[exgn]]
	datavals = list(prot = gdata$lMSmu%>%as.vector,ribo = gdata$lSeqmu%>%as.vector)%>%
	enframe('assay','value')%>%unnest(value)
	datasigs = list(prot = gdata$lMSsig%>%as.vector,ribo = gdata$lSeqsig%>%as.vector)%>%
	enframe('assay','value')%>%unnest(value)
	datavals$lower = datavals$value - datasigs$value*1.96
	datavals$upper = datavals$value + datasigs$value*1.96
	datavals%>%mutate(time=c(1:5,1:5))
})
modeldata$model='data'
gplotdf = bind_rows(modeldata,modelvals)

g2plotttile = paste0(exampgenes,'\n',names(exampgenes))%>%setNames(exampgenes)

modelcols2 <- c(modelcols,data='black')
modelcols2<-c('data'='black','degredation'='green','production'='blue','stationary'='purple','linear'='lightblue',msdev='red')

#now plot
plotfile<- here(paste0('plots/','trajectory_example_plots','.pdf'))
pdf(plotfile,w=5*3,h=1*4)
plotlist = gplotdf%>%
	split(.,.$gene)%>%
	lapply(function(x)x%>%
		filter(assay=='prot' | model=='data')%>%
		mutate(gene=g2plotttile[gene])%>%
		filter(gene%>%str_detect(model)|model=='data')%>%
		# group_slice(1)%>%
		ggplot(aes(y=value,x=time,color=model))+
		geom_line()+
		ylab('')+
		scale_color_manual(values=modelcols2)+
		facet_grid(scale='free',assay~gene)+
		geom_errorbar(width=(0.2),aes(ymin=lower,ymax=upper))+
		theme_bw()
	)
# plotlist[[5]] %<>% add(scale_color_manual(values=modelcols2))
print(ggarrange(ncol=length(plotlist),plotlist=plotlist,common.legend=TRUE))
dev.off()
message(normalizePath(plotfile))

##What if we classify into only 
devdf = modeltestdf%>%filter(data=='riboseq')%>%group_by(gene)%>%
	summarise(dev=BIC[model=='msdev']<BIC[model=='production'])


devdf%>%left_join(mcshanethalfs)%>%
	filter(!is.na(McShane_deg_cat))%>%
	glm(data=.,dev ~ McShane_deg_cat,family='binomial')%>%summary
	group_by(McShane_deg_cat)%>%summarise(sum(dev),sum(!dev))


residdf = modeltestdf%>%filter(data=='riboseq')%>%
	inner_join(devdf%>%filter(dev))%>%
	filter(model=='production')%>%
	{set_rownames(as.matrix(select(ungroup(.),matches('residuals'))),.[['gene']])}%>%head

princomp(residdf)%>%

################################################################################
########DIagnostics
################################################################################
if(TRUE){	
	modeltestdf%>%filter(gene==filteredgenes[1])%>%filter(data=='riboseq')%>%filter(model=='production')%>%t
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
		.$gene%>%unique


	# lowfiltergenes = modeltestdf%>%
	# 	filter(data=='riboseq')%>%
	# 	filter(passtest[model=='production'],best[model=='production'])%>%
	# 	.$gene%>%
	# 	unique

outlgenes = 	estimate%>%filter(!between(l_pihalf,-5,5))%>%.$gene
nonoutlgenes = 	estimate%>%filter(!gene%in%outlgenes)%>%.$gene
genenonconv = bmodelopts%>%map(2)%>%map('production')%>%discard(is.null)%>%map_dbl('return_code')%>%`==`(70)
convgenes = genenonconv%>%.[!.]%>%names
nonoutconvgenes = nonoutlgenes%>%intersect(convgenes)
# bmodelopts[convgenes]%>%map(2)%>%map('production')%>%discard(is.null)%>%map_dbl('return_code')
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
}