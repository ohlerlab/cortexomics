
################################################################################
########Code for plotting trajectories
################################################################################

samplingmedsaslist<-samplings[finaltestgenes]%>%map('result')%>%map(.%>%summary%>%.[[1]]%>%{setNames(.[,'50%'],rownames(.))})%>%map(.%>%stanpars_to_list%>%.[names(.)!='lp__'])
#Now get initial values from the samplings
samplingcombinitvals <- lapply(names(bestfitinits[[1]])%>%setNames(.,.),function(argind){
  samplingmedsaslist%>%map(argind)%>%setNames(.,seq_along(.))%>%do.call(what=partial(abind,along=1))
})	

data.frame(uprotein_id=finaltestgenes)%>%left_join(metainfo)%>%filter(gene_name=='Satb2')
data.frame(uprotein_id=finaltestgenes)%>%left_join(metainfo)%>%filter(gene_name=='Flna')
finaltestgenes%>%head(2)

ribonorm_combinitvals$prot[which(finaltestgenes=='ENSMUSP00000098997_239'),]
get_dp_standata(finaltestgenes)$lMS[which(finaltestgenes==best_satb2_uid),]%>%txtplot
ribonorm_combinitvals$prot[which(finaltestgenes==best_satb2_uid),]


#parse the samplings made on individual genes
get_sampling_par_cis<-lapply(seq_along(finaltestgenes),function(i){
	parse_sampling_pars(sgenesampling=samplings[[i]]$result,uid=finaltestgenes[i])})%>%
	setNames(.,seq_along(.))%>%bind_rows(.id='genenum')%>%mutate(gene=gene+as.numeric(genenum)-1)
get_sampling_prot_cis<-get_sampling_par_cis%>%filter(parameter=='prot')
get_sampling_prot_lps<-get_sampling_par_cis%>%filter(parameter=='lp__')

ribonormopt_protdf<-tidy_opt_prots(ribonormopt$par,finaltestgenes)

summary<-rstan::summary

get_sampling_par_cis%>%filter(parameter=='l_pihalf')

ribnorm_sampling_protdf <- ribonormsamp%>%parse_sampling_prot(finaltestgenes)%>%
	filter(parameter=='prot')
	# filter(parameter=='ribnorm')
ribnorm_sampling_protdf%>%filter(parameter=='prot')
ribonormsamp%>%parse_sampling_prot(finaltestgenes)%>%filter(parameter=='ribnorm')
ribonormsamp%>%parse_sampling_prot(finaltestgenes)

stanmodels<-list(
	indivopts=combinitvals%>%tidy_opt_prots(finaltestgenes),
	indivsamplings=get_sampling_prot_cis,
	ribonormopt=ribonormopt_protdf,
	ribonormsamp=ribnorm_sampling_protdf
)
correctionmodels<-c(names(stanmodels))


best_ms_ids= metainfo$ms_id[match(best_uprotein_ids,metainfo$uprotein_id)]
for(corectionmodel in correctionmodels){
   #correct from scaled down vals up to matched_ms_mat vals
   # corval = (postmeanmat[best_uprotein_ids,]%>%apply(2,median,na.rm=T)%>%median) - matchedms_mat_rscl%>%apply(2,median,na.rm=T) 
  corval=msmed
  stopifnot(!any(stanmodels[[corectionmodel]]$gene_name%>%is.na))
  stanmodels[[corectionmodel]]$estimate <- stanmodels[[corectionmodel]]$estimate + median(corval)
  stanmodels[[corectionmodel]]$CI.L <- stanmodels[[corectionmodel]]$CI.L + median(corval)
  stanmodels[[corectionmodel]]$CI.R <- stanmodels[[corectionmodel]]$CI.R + median(corval)  
}

warning('I should first at least get the rescaling right')


for(i in seq_along(stanmodels)){ 
  message(names(stanmodels)[i])
  stopifnot(
  	c('estimate','CI.R','CI.L','gene_name','assay','ntime')%in%colnames(stanmodels[[i]])
  )
  stopifnot('Flna'%in%stanmodels[[i]]$gene_name)
}
conflict_prefer("intersect", "BiocGenerics")
lowlikgenes = get_sampling_prot_lps%>%arrange(estimate)%>%.$gene_name%>%head(2)

# source('src/R/Shiny/model_app.R')
# save.image('data/run_degmodel_dropout.Rdata')
# load('data/run_degmodel_dropout.Rdata')