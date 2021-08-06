#run prep then modeling
if(!exists('bmodelopts'))source('src/becker_modelling.R')

mcshanetest <- function(x) x$l_pihalf%>%setNames(combinitvals$genes)%>%enframe('gene','l_pihalf')%>%
	inner_join(mcshanethalfs)%>%
	# filter(abs(l_pihalf) <  5)%>%
	filter(McShane_deg_cat!='NED')%>%
	# filter(l_pihalf<4)%>%
	{quicktest(.$l_pihalf,log(.$half_life));.}%>%
	{cor.test(.$l_pihalf,log(.$half_life))}%>%tidy

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
# jointmodel_hierach = stan_model(here('src/Archive/Stan/becker_proda_jhierarch_ldev_pmass.stan'))
# jointmodel_hierach = stan_model(here('src/Archive/Stan/becker_proda_jhierarch_ldev_5ks.stan'))
jointmodel_hierach = stan_model(here('src/Archive/Stan/becker_proda_jhierarch_ldev_nomix.stan'))
# jointmodel_hierach = stan_model(here('src/Archive/Stan/becker_proda_jhierarch.stan'))
gns4model = gns_by_mod[c('production','linear','degredation','msdev','stationary')]%>%unlist
# gns4model = notmsdevgenes
# gns4model = nonoutconvgenes
combinitvals <- bmodelopts[gns4model]%>%map('riboseq')%>%map('production')%>%map('par')%>%discard(is.null)%>%get_comb_initvals
jointdata = datafuns$riboseq(combinitvals$genes)
# txtdensity(jointdata$lMSmu - log(combinitvals$prot))
jointdata$G = nrow(jointdata$lMSmu)
combinitvals$lKs = 0
combinitvals$mu_lks=0
combinitvals$mu_l_pihalf=0
combinitvals$sd_lks=1
combinitvals$sd_l_phalf=1
# combinitvals$theta = array(0.5,jointdata$G)
# combinitvals$l_pihalf %<>% {.[!between(.,-5,5)]=0;.} 
# combinitvals$l_pihalf %<>% {.[!between(.,-5,5)]=0;}.
combinitvals$l_pihalf = array(0,jointdata$G)
# jointdata$lMSmu%<>% multiply_by(1/log2(exp(1)))
# jointdata$lSeqmu%<>% multiply_by(1/log2(exp(1)))
# jointdata$lMSsigma%<>% multiply_by(1/log2(exp(1)))
# jointdata$lSeqsigma%<>% multiply_by(1/log2(exp(1)))
# combinitvals$l_st_dev = (jointdata$lMSmu - jointdata$lSeqmu)%>%rowMeans%>%array
# combinitvals$msdev = jointdata$lMSmu - jointdata$lSeqmu - rep(combinitvals$l_st_dev,5)
#
jhopth = rstan::optimizing(jointmodel_hierach,
	data=jointdata,
	init=combinitvals,
	as_vector=F,save_iterations=TRUE,hessian=F,verbose=T,
)

}
jhopth$par$l_pihalf%>%which.min
# jhopth$par$theta%>%c(0,1)%>%txtdensity
# jhopth$par$theta
jhopth$par$Ks%>%.[1]
gns4model%>%length
jhopth$par%>%mcshanetest
jhopth$par$l_pihalf%>%txtdensity
jhopth$par$ribooffset%>%txtplot


#things we are tweaking, theta init at 0.5
#I think we probably need a prior on theta to stop the model just sticking everything into msdev
#maybe this gets better if we don't optimize though...
#we could also try a prior on msdev - penilize it for getting too large...
#anyway we have a fit with only production, linear and degredation, with theta at 0.5 to start,
#that's pretty fucking good already.

#will that fit with the hierarch parameters free?

#current state - production onnly or with deg and linear theta starts at 0.5, no init of the other parameters - we get convergence this way.

#prior on msdev normal 1, all genes - we end up with initialization errors again, and once again theta going to 1
#okay so what if weeeeeeee, tried to use pihalfs but onlyt the sane ones. No that still just goes right to theta=1
#I think this is an inevitable problem from optimization. Let's sample instead


#so with free parameters, can I at least opt without the msdev genes, hmmm seems to create problems?

#adding a minimum val to prot and dprot? Doesn't hurt, seems like... but also results in a basically all theta = 0 

#scaling the input data seems to cause failure, and yup, remove it and we are fine
# confirm we need the mixture? No we donät, if we leave out ms dev.
# do i need to go back and use target in the comparison scripts? - can' thurt. Done
# okay so let's figure out why theta is still in jopth....