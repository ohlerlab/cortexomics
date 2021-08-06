get_comb_initvals <- function(bestfitinits){
  combinitvals <- lapply(names(bestfitinits[[1]])%>%setNames(.,.),function(argind){
    bestfitinits%>%
      map(argind)%>%
      setNames(.,seq_along(.))%>%
      do.call(what=partial(abind::abind,along=1))
  })
  combinitvals	
}

################################################################################
########Now do joint modelling
################################################################################
	

{
# filteredgenes=filteredgenesold
filteredgenes = modfilteredgenes%>%
	c(genebmodels%>%.[.=='production']%>%names)%>%
	c(genebmodels%>%.[.=='degredation']%>%names)%>%
	# c(genebmodels%>%.[.=='stationary']%>%names)%>%
	# c(genebmodels%>%.[.=='msdev']%>%names)%>%
	c(genebmodels%>%.[.=='linear']%>%names%>%{.})%>%
	unique%>%
	intersect(rownames(g_elong_mat))
filteredgenes=ext_elong_gns%>%intersect(filteredgenes)


jointmodel1te = stan_model(here('src/Archive/Stan/becker_proda_oneKs.stan'))
filteredgenes=nonoutconvgenes
combinitvals <- bmodelopts[filteredgenes]%>%map('riboseq')%>%map('production')%>%map('par')%>%get_comb_initvals
jointdata = datafuns$riboseq(filteredgenes)
jointdata$G = nrow(jointdata$lMSmu)
combinitvals$lKs = 0
#
combinitvals$l_pihalf%>%txtdensity
combinitvals$l_st%>%txtdensity
combinitvals$l_pihalf[]=0
combinitvals$l_st[]=0
jopt = NULL
while(is.null(jopt)){
	jopt <- rstan::optimizing(
	# jopt <- possibly(rstan::optimizing,NULL)(
		jointmodel1te,
		data=jointdata,
		init=combinitvals,
		as_vector=F,
		save_iterations=TRUE,
		hessian=F,verbose=T,)	
}
pihalftest = jopt$par%>%.$l_pihalf%>%setNames(filteredgenes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	# filter(l_pihalf<3)%>%
	{quicktest(.$l_pihalf,log(.$half_life))}%>%tidy

message(paste0(sep='\n',capture.output(pihalftest)))
message(paste0(sep='\n',capture.output(jopt$par$Ks[1])))

}


pihalftest = jopt$par%>%.$l_pihalf%>%setNames(filteredgenes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	# filter(l_pihalf<3)%>%
	{quicktest(.$l_pihalf,log(.$half_life))}%>%tidy
pihalftest
pihalftest = ixopt$par%>%.$l_pihalf%>%setNames(filteredgenes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	# filter(l_pihalf<3)%>%
	{quicktest(.$l_pihalf,log(.$half_life))}%>%tidy
pihalftest


jopt$value
ixopt$value
# nedgenes <- (mcshanethalfs%>%filter(McShane_deg_cat=='NED')%>%.$gene)

# modeltestdf%>%group_by(gene)%>%	
# 	filter(data=='riboseq')%>%
# 	filter(best)%>%
# 	mutate(isned = gene %in% nedgenes)%>%
# 	mutate(sumresid = sum(residuals1^2+residuals2^2+residuals3^2+residuals4^2+residuals5^2))%>%
# 	{split(.$sumresid,.$isned)}%>%
# 	{t.test(.[['TRUE']],.[['FALSE']])}
# 	# {table(.$isned,.$model==)}%>%fisher.test


# codreltests = pihalftest%>%{paste0('rho = ',round(.$estimate,3),'\n','pval = ',ifelse(.$p.value > 0.001,round(.$p.value,2),format(.$p.value,format='e',digits=2)))}


#now plot
plotfile<- here(paste0('plots/jte_indiv_v_mcshane_pihalf','.pdf'))
dir.create(dirname(plotfile))
pdf(plotfile)
ggdf <- jopt$par%>%.$l_pihalf%>%setNames(filteredgenes)%>%enframe('gene','l_pihalf')%>%
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
	# filter(abs(l_pihalf) <  5)%>%
	# filter(McShane_deg_cat!='NED')%>%
	# filter()
	# ggplot(.,aes(log(half_life),l_pihalf,color=McShane_deg_cat))+
	ggplot(.,aes(log(half_life),l_pihalf))+
	geom_point(alpha=I(0.5),size=I(0.3))+
	facet_grid(McShane_deg_cat~.)+
	scale_color_discrete(name='McShane_deg_cat')+
	scale_x_continuous(paste0('log2(Half Life) (McShane et al)'))+
	scale_y_continuous(paste0('Joint TE Model - estimated log(Half Life)'))+
	geom_text(data=scatterlabls,hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))+
	ggtitle(paste0('Measured Half Lives vs Estimated'),sub=codreltests)+
	geom_smooth(method='lm',color=I('black'))+
	theme_bw()
p
dev.off()
normalizePath(plotfile)


################################################################################
########Now let's use that fixed lKS to get half life estimates for eveyrthing
################################################################################
'becker_proda_lKsfix.stan'

nofailgenes = bmodelopts%>%map('riboseq')%>%map('production')%>%map_lgl(Negate(is.null))%>%.[.]%>%names

allinivals <- bmodelopts%>%.[nofailgenes]%>%map('riboseq')%>%map('production')%>%map('par')%>%get_comb_initvals
alldata = datafuns$riboseq(nofailgenes)
alldata$G = nrow(alldata$lMSmu)
#
alldata$lKs = jopt$par$lKs
#
allinivals$l_pihalf[]=0
allinivals$l_st[]=0
fopt <- NULL
while(is.null(fopt)){
	fopt <- possibly(rstan::optimizing,NULL)(
		stan_model('src/Archive/Stan/becker_proda_lKsfix.stan'),
		data=alldata,
		init=allinivals,
		as_vector=F,
		save_iterations=TRUE,
		hessian=F,verbose=T,)	
}

pihalftest = fopt$par%>%.$l_pihalf%>%setNames(nofailgenes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	# filter(l_pihalf<3)%>%
	{quicktest(.$l_pihalf,log(.$half_life))}%>%tidy

pihalftest

stop()

fopt$par%>%.$l_pihalf%>%setNames(nofailgenes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	summarise(mean(l_pihalf),mean(log(half_life)),sd(l_pihalf),sd(log(half_life)))


#now plot
plotfile<- here(paste0('plots/jte_indiv_v_mcshane_pihalf','.pdf'))
dir.create(dirname(plotfile))
pdf(plotfile,h=21/2,w=7/2)
ggdf <- fopt$par%>%.$l_pihalf%>%setNames(nofailgenes)%>%enframe('gene','l_pihalf')%>%
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
	# filter(abs(l_pihalf) <  5)%>%
	# filter(McShane_deg_cat!='NED')%>%
	# filter()
	# ggplot(.,aes(log(half_life),l_pihalf,color=McShane_deg_cat))+
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

################################################################################
########Now let's try with ribo offsets per timepoint?
################################################################################
	


{
# filteredgenes=filteredgenesold
filteredgenes = modfilteredgenes%>%
	c(genebmodels%>%.[.=='production']%>%names)%>%
	# c(genebmodels%>%.[.=='degredation']%>%names)%>%
# 	# c(genebmodels%>%.[.=='stationary']%>%names)%>%
# 	c(genebmodels%>%.[.=='msdev']%>%names)%>%
	# c(genebmodels%>%.[.=='linear']%>%names%>%{.})%>%
	unique
jointmodel5te = stan_model(here('src/Archive/Stan/becker_proda_fiveKs.stan'))
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
combinitvals$lKs = 0
#
combinitvals$l_pihalf%>%txtdensity
combinitvals$l_st%>%txtdensity
combinitvals$l_pihalf[]=0
combinitvals$l_st[]=0
combinitvals$ribooffset=array(0,4)
# combinitvals = combinitvals$ribooffset
jopt = NULL
while(is.null(jopt)){

	jopt <- possibly(rstan::optimizing,NULL)(
	# jopt <- rstan::optimizing(
		jointmodel5te,
		data=jointdata,
		init=combinitvals,
		as_vector=F,
		save_iterations=TRUE,
		hessian=F,verbose=T)

}


pihalftest = jopt$par%>%.$l_pihalf%>%setNames(filteredgenes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	# filter(l_pihalf<3)%>%
	{quicktest(.$l_pihalf,log(.$half_life))}%>%tidy

message(paste0(sep='\n',capture.output(pihalftest)))
message(paste0(sep='\n',capture.output(jopt$par$Ks[1])))



jopt$par$ribooffset%>%txtplot
jopt$par$protoffset%>%txtplot


}


################################################################################
########
################################################################################

trelongs = Sys.glob('../Ribotransformer/pipeline/ixnos_elong/*/*.elong.csv')%>%
	setNames(.,basename(dirname(.)))%>%
	map_df(.id='sample',read_csv)

iso_prop_df = iso_tx_countdata$abundance%>%as.data.frame%>%rownames_to_column('tr_id')%>%pivot_longer(-tr_id)%>%
	left_join(ids_nrgname%>%select(tr_id=transcript_id,gene_id,gene_name)%>%distinct)%>%
	group_by(gene_name,gene_id,name)%>%mutate(prop=value/sum(value+1e-12))%>%
	select(tr_id,gene_id,sample=name,prop)

table(filteredgenes%in%iso_prop_df$gene_name)

iso_prop_df%<>%filter(is.finite(prop))

g_elong_df = trelongs%>%mutate(tr_id = str_replace(tr_id,'\\.\\d+$',''))%>%
	left_join(ids_nrgname%>%select(tr_id=transcript_id,gene_id,gene_name)%>%distinct)%>%	
	left_join(iso_prop_df%>%ungroup%>%select(-gene_name),by=c('sample','gene_id'))%>%
	group_by(sample,gene_id,gene_name)%>%
	summarise(elong=weighted.mean(x=elong,w=prop))

g_elong_mat = g_elong_df%>%ungroup%>%select(gene_name,sample,elong)%>%
	pivot_wider(names_from='sample',values_from='elong')%>%
	{set_rownames(as.matrix(.[,-1]),.[['gene_name']])}
g_elong_mat%<>%replace_na(1)

g_elong_mat <- sapply(seq(1,10,by=2),function(i){
	g_elong_mat[,i]+g_elong_mat[,i+1]/2
})
ext_elong_gns =  g_elong_mat%>%rowMeans%>%sort%>%{c(head(.,100),tail(.,100))}%>%names

table(filteredgenes%in%rownames(g_elong_mat))



{
# filteredgenes=filteredgenesold
filteredgenes = modfilteredgenes%>%
	c(genebmodels%>%.[.=='production']%>%names)%>%
	c(genebmodels%>%.[.=='degredation']%>%names)%>%
	# c(genebmodels%>%.[.=='stationary']%>%names)%>%
	# c(genebmodels%>%.[.=='msdev']%>%names)%>%
	c(genebmodels%>%.[.=='linear']%>%names%>%{.})%>%
	unique%>%
	intersect(rownames(g_elong_mat))
filteredgenes=ext_elong_gns%>%intersect(filteredgenes)

ixnos_offsetmodel = stan_model(here('src/Archive/Stan/becker_proda_ixnos.stan'))
stopifnot(length(filteredgenes)>20)

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
jointdata$lSeq_offset = g_elong_mat[filteredgenes,]
combinitvals$lKs = 0
#
combinitvals$l_pihalf%>%txtdensity
combinitvals$l_st%>%txtdensity
combinitvals$l_pihalf[]=0
combinitvals$l_st[]=0
combinitvals$ribooffset=array(0,4)
# combinitvals = combinitvals$ribooffset
ixopt = NULL
while(is.null(ixopt)){

	ixopt <- possibly(rstan::optimizing,NULL)(
	# ixopt <- rstan::optimizing(
		ixnos_offsetmodel,
		data=jointdata,
		init=combinitvals,
		as_vector=F,
		save_iterations=TRUE,
		hessian=F,verbose=T)

}


pihalftest = ixopt$par%>%.$l_pihalf%>%setNames(filteredgenes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	# filter(l_pihalf<3)%>%
	{quicktest(.$l_pihalf,log(.$half_life))}%>%tidy

message(paste0(sep='\n',capture.output(pihalftest)))
message(paste0(sep='\n',capture.output(ixopt$par$Ks[1])))

}


################################################################################
########Okay so the joint model with a point estimate works okay, what about
########If we optimize a hiearach model?
################################################################################
{
jointmodel_hierach = stan_model(here('src/Archive/Stan/becker_proda_jhierarch_ldev_nomix.stan'))
gns4model = gns_by_mod$production
# gns4model = nonoutconvgenes
combinitvals <- bmodelopts[gns4model]%>%map('riboseq')%>%map('production')%>%map('par')%>%get_comb_initvals
jointdata = datafuns$riboseq(gns4model)
jointdata$G = nrow(jointdata$lMSmu)
combinitvals$lKs = 0
combinitvals$mu_lks=0
combinitvals$mu_l_pihalf=0
combinitvals$sd_lks=1
combinitvals$sd_l_phalf=1
combinitvals$theta = array(0.5,combinitvals$G)
#
jhopth = rstan::optimizing(jointmodel_hierach,data=jointdata,init=combinitvals,as_vector=F,save_iterations=TRUE,hessian=F,verbose=T)
#
jhopth$par%>%.$l_pihalf%>%setNames(gns4model)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	# filter(l_pihalf<4)%>%
	{quicktest(.$l_pihalf,log(.$half_life));.}%>%
	{cor.test(.$l_pihalf,log(.$half_life))}%>%tidy
}

jhopth$par%>%names
jhopth$par['var_l_phalf']
jhopth$par['mu_l_pihalf']
jhopth$par$l_pihalf%>%txtdensity
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