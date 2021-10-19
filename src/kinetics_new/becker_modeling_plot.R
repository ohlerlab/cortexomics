library(here)
if(!exists('genebmodels'))source('src/kinetics_new/becker_modelling.R')
################################################################################
########## 
################################################################################
#now plot
modelcols <- c(
	msdev = '#c64730',
	production = '#fda27c',
	linear = '#ba5e47',
	degredation = '#1d4c9b',
	stationary = '#6289cc'
)
plotfile<- here(paste0('plots/','kinetics_redo/becker_trajectory_classes','.pdf'))
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
	onts = c('BP','MF','')
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
plotfile<- here(paste0('plots/','kinetics_redo/beck_traj_class_goplots','.pdf'))
pdf(plotfile,h=7,w=7*5)
print(ggarrange(plotlist=goplotlist,ncol=length(goplotlist)))
dev.off()
message(normalizePath(plotfile))		 
#
trajclass_goterms%>%split(.,.$cluster)%>%.['msdev']%>%.[[1]]%>%
	filter(Significant>10)%>%filter(elimFisher<0.05)%>%arrange(-Enrichment)%>%
	head(12)

#
notmsdevgenes = modeltestdf%>%filter(data=='riboseq')%>%
	group_by(gene)%>%
	filter(n()==5)%>%
	filter(2+BIC[model=='production']<BIC[model=='msdev'])%>%
	.$gene%>%unique
#
msdevsplit = ifelse(names(genebmodels)%in%notmsdevgenes,'not_msdev','msdev')%>%
	setNames(names(genebmodels))
msdevsplitgoterms <- get_cluster_gos(msdevsplit)
goplotlist = msdevsplitgoterms%>%split(.,.$cluster)%>%.['msdev']%>%lapply(.%>%
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
plotfile<- here(paste0('plots/','kinetics_redo/msdev_goplot','.pdf'))
pdf(plotfile,h=7,w=7*1)
print(ggarrange(plotlist=goplotlist,ncol=length(goplotlist)))
dev.off()
message(normalizePath(plotfile))		 
msdevsplitgoterms%>%filter(elimFisher<0.05)%>%head


if(!exists('jhopth')) source('src/kinetics_redo/becker_modelling_hierarch.R')
#plotting absolute pihalfs
#
high_pihalf_genes <- jhopth$par$l_pihalf[match(notmsdevgenes,combinitvals$gene)]%>%
	setNames(notmsdevgenes)
high_pihalf_genes=high_pihalf_genes%>%unlist%>%.[!is.na(.)]%>%`>`(quantile(.,.75))
msdevsplit = ifelse(high_pihalf_genes,'high_pihalf','low_pihalf')%>%
	setNames(names(high_pihalf_genes))
msdevsplitgoterms <- get_cluster_gos(msdevsplit)
goplotlist = msdevsplitgoterms%>%filter(ontology=='BP')%>%split(.,.$cluster)%>%lapply(.%>%
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
plotfile<- here(paste0('plots/','high_pihalf_genes','.pdf'))
pdf(plotfile,h=7,w=7*2)
print(ggarrange(plotlist=goplotlist,ncol=length(goplotlist)))
dev.off()
message(normalizePath(plotfile))		 



#
#plotting change in pihalfs
#
high_pihalf_genes <- jhopth$par$l_pihalf[match(notmsdevgenes,combinitvals$gene)]%>%
	setNames(notmsdevgenes)
high_pihalf_genes=high_pihalf_genes%>%unlist%>%.[!is.na(.)]%>%`>`(quantile(.,.75))
msdevsplit = ifelse(high_pihalf_genes,'high_pihalf','low_pihalf')%>%
	setNames(names(high_pihalf_genes))
msdevsplitgoterms <- get_cluster_gos(msdevsplit)
goplotlist = msdevsplitgoterms%>%filter(ontology=='BP')%>%split(.,.$cluster)%>%lapply(.%>%
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
plotfile<- here(paste0('plots/','high_pihalf_genes','.pdf'))
pdf(plotfile,h=7,w=7*2)
print(ggarrange(plotlist=goplotlist,ncol=length(goplotlist)))
dev.off()
message(normalizePath(plotfile))		 



################################################################################
########Plot numbers of NED genes
################################################################################

msdevsplit = ifelse(names(genebmodels)%in%notmsdevgenes,'not_msdev','msdev')%>%
	setNames(names(genebmodels))
plotfile<- here(paste0('plots/kinetics_redo/','NED_vs_notmsdev_barplots','.pdf'))
pdf(plotfile,h=5,w=3)
ggdf=msdevsplit%>%enframe('gene','msdev')%>%inner_join(mcshanethalfs)%>%
	# filter(McShane_deg_cat!='UN')%>%
	group_by(McShane_deg_cat,msdev)%>%tally%>%
	group_by(msdev)%>%mutate(nn=n,tot=sum(n),num=paste0(n,'/',tot),n = n/sum(n))%>%
	filter(McShane_deg_cat=='NED')
pval = matrix(c(ggdf$nn,ggdf$tot),ncol=2)%>%fisher.test%>%tidy%>%{paste0('p = ',round(.$p.value,3))}
pval
ggdf%>%
	ggplot(aes(y=n,x=msdev,fill=McShane_deg_cat))+
	stat_identity(geom='bar',position='dodge',fill=I('black'))+
	geom_text(aes(label=num),vjust=2,color=I('white'))+
	ggtitle(pval)+
	xlab('')+
	scale_y_continuous('fraction in NED category',breaks=seq(0,0.125,by=0.025))+
	theme_bw()
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
exampgenes = bicdiffdf%>%group_by(model)%>%
	# filter(!gene%in%c('Ttr','Hist1h2bc'))%>%
	slice(which.max(bicdiff))%>%{setNames(.$gene,.$model)}
gns4plot=exampgenes%>%setNames(.,.)
modelnms=unique(bicdiffdf$model)%>%setNames(.,.)
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
# modelcols2 <- c(modelcols,data='black')
modelcols2<-c('data'='black','degredation'='green','production'='blue','stationary'='purple','linear'='lightblue',msdev='red')
#now plot
plotfile<- here(paste0('plots/kinetics_redo/','trajectory_example_plots','.pdf'))
pdf(plotfile,w=5*3,h=1*4)
# gplotdf$model%>%unique
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
devdf = modeltestdf%>%
	filter(data=='riboseq')%>%group_by(gene)%>%
	summarise(dev=BIC[model=='msdev']<BIC[model=='production'])

devdf%>%left_join(mcshanethalfs)%>%
	filter(!is.na(McShane_deg_cat))%>%
	glm(data=.,dev ~ McShane_deg_cat,family='binomial')%>%summary
	group_by(McShane_deg_cat)%>%summarise(sum(dev),sum(!dev))

ggdf=devdf%>%left_join(mcshanethalfs)%>%rename('msdev':=dev)%>%
	group_by(McShane_deg_cat,msdev)%>%
	tally%>%
	group_by(msdev)%>%mutate(nn=n,tot=sum(n),num=paste0(n,'/',tot),n = n/sum(n))%>%
	# filter(McShane_deg_cat!='UN')%>%
	filter(McShane_deg_cat=='NED')
pval = matrix(c(ggdf$nn,ggdf$tot),ncol=2)%>%fisher.test%>%tidy%>%{paste0('p = ',round(.$p.value,3))}
residdf = modeltestdf%>%filter(data=='riboseq')%>%
	inner_join(devdf%>%filter(dev))%>%
	filter(model=='production')%>%
	{set_rownames(as.matrix(select(ungroup(.),matches('residuals'))),.[['gene']])}%>%head



################################################################################
################################################################################
#now plot
plotfile<- here(paste0('plots/','pihalf_distributions','.pdf'))
pdf(plotfile,w=5,h=4)
jhopth$par$l_pihalf[jointdata$lMSmu%>%rownames%>%is_in(mcshanethalfs$gene)]%>%enframe%>%
	mutate(type='Indiv Estimates (McShane Genes)')%>%
	bind_rows(tibble(value = rnorm(10e3,jhopth$par$mu_l_pihalf,sqrt(jhopth$par$var_l_phalf)),type='hierarch_dist'))%>%
	bind_rows(mcshanethalfs%>%select(value=half_life)%>%mutate(type='mcshane'))%>%
	ggplot(aes(x=value,fill=type))+
		geom_density(alpha=I(0.5))+
		facet_grid(type~.)
dev.off()
message(normalizePath(plotfile))


################################################################################
########## Pihalf cors
################################################################################
#now plot
plotfile<- here(paste0('plots/kinetics_redo/hierarch_tefix_indiv_v_mcshane_pihalf','.pdf'))
dir.create(dirname(plotfile))
pdf(plotfile)
ggdf <- jhopth$par%>%.$l_pihalf%>%setNames(combinitvals$gene)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)
scatterlabls =ggdf%>%
	group_by(McShane_deg_cat)%>%
	nest%>%
	mutate(labl = paste0(
		'r = ',map_dbl(data,~cor(method='spearman',use='complete',log(.$half_life),.$l_pihalf))%>%round(3),'\n',
		'p = ',map_dbl(data,~cor.test(method='spearman',use='complete',log(.$half_life),.$l_pihalf)$p.value%>%round(3))
	))
p = 
	ggdf%>%
	ggplot(.,aes(log(half_life),l_pihalf))+
	geom_point(alpha=I(0.5),size=I(0.3))+
	facet_grid(McShane_deg_cat~.)+
	scale_color_discrete(name='McShane_deg_cat')+
	scale_x_continuous(paste0('log2(Half Life) (McShane et al)'))+
	scale_y_continuous(paste0('Joint TE Model - estimated log(Half Life)'))+
	geom_text(data=scatterlabls,hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))+
	ggtitle(paste0('Measured Half Lives vs Estimated'))+
	geom_smooth(method='lm',color=I('black'))+
	theme_bw()
p
dev.off()
normalizePath(plotfile)
ggdf%>%transmute(gene_name=gene,
		log2_est_half_life=l_pihalf,
		log2_mcshane_pihalf=log2(half_life),
		McShane_deg_cat)%>%
	write_tsv('tables/est_pihalf_v_mcshane.tsv')
