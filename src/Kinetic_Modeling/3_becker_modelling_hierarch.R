
#run prep then modeling
if(!exists('bmodelopts'))source('src/Kinetic_Modelling/becker_modelling.R')

mcshanetest <- function(x) x$l_pihalf%>%
	setNames(combinitvals$genes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	filter(McShane_deg_cat!='NED')%>%
	{quicktest(.$l_pihalf,log(.$half_life));.}%>%
	{cor.test(.$l_pihalf,log(.$half_life))}%>%tidy

################################################################################
########Okay so the joint model with a point estimate works okay, what about
########If we optimize a hiearach model?
################################################################################

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
jointmodel_hierach = stan_model(here('src/Stan/becker_proda_jhierarch_ldev_nomix.stan'))
gns4model = gns_by_mod[c('production','linear','degredation','msdev','stationary')]%>%unlist
combinitvals <- bmodelopts[gns4model]%>%
	map('riboseq')%>%map('production')%>%
	map('par')%>%discard(is.null)%>%get_comb_initvals
jointdata = datafuns$riboseq(combinitvals$genes)
jointdata$G = nrow(jointdata$lMSmu)
combinitvals$lKs = 0
combinitvals$mu_lks=0
combinitvals$mu_l_pihalf=0
combinitvals$sd_lks=1
combinitvals$sd_l_phalf=1

combinitvals$l_pihalf = array(0,jointdata$G)

jhopth = rstan::optimizing(jointmodel_hierach,
	data=jointdata,
	init=combinitvals,
	as_vector=F,save_iterations=TRUE,hessian=F,verbose=T,
)

}
