

################################################################################
########Okay so the joint model with a point estimate works okay, what about
########If we optimize a hiearach model?
################################################################################

# estimate%>%filter(gene%in%gns4model)%>%.$l_pihalf%>%txtdensity
# estimate%>%filter(gene%in%gns4model)%>%.$l_pihalf%>%min
get_comb_initvals <- function(bestfitinits){
  combinitvals <- lapply(names(bestfitinits[[1]])%>%setNames(.,.),function(argind){
    bestfitinits%>%
      map(argind)%>%
      setNames(.,seq_along(.))%>%
      do.call(what=partial(abind::abind,along=1))
  })
  combinitvals$genes = names(bestfitinits)
  combinitvals
}

{
jointmodel_hierach = stan_model(here('src/Archive/Stan/becker_proda_jhierarch_ldev.stan'))
# jointmodel_hierach = stan_model(here('src/Archive/Stan/becker_proda_jhierarch.stan'))
gns4model = gns_by_mod[c('production','linear','degredation')]%>%unlist
# gns4model = nonoutconvgenes
combinitvals <- bmodelopts[gns4model]%>%map('riboseq')%>%map('production')%>%map('par')%>%discard(is.null)%>%get_comb_initvals
jointdata = datafuns$riboseq(combinitvals$genes)
jointdata$G = nrow(jointdata$lMSmu)
combinitvals$lKs = 0
combinitvals$mu_lks=0
combinitvals$mu_l_pihalf=0
combinitvals$sd_lks=1
combinitvals$sd_l_phalf=1
combinitvals$theta = array(0.99,jointdata$G)
combinitvals$l_pihalf = array(0,jointdata$G)
combinitvals$l_st_dev = (jointdata$lMSmu - jointdata$lSeqmu)%>%rowMeans%>%array
combinitvals$msdev = jointdata$lMSmu - jointdata$lSeqmu - rep(combinitvals$l_st_dev,5)
#
jhopth = rstan::optimizing(jointmodel_hierach,data=jointdata,init=combinitvals,as_vector=F,save_iterations=TRUE,hessian=F,verbose=T)
#
jhopth$par%>%.$l_pihalf%>%setNames(combinitvals$genes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	# filter(l_pihalf<4)%>%
	{quicktest(.$l_pihalf,log(.$half_life));.}%>%
	{cor.test(.$l_pihalf,log(.$half_life))}%>%tidy
}
jhopth$par$l_pihalf%>%which.min