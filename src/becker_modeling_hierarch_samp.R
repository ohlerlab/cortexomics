#run prep then modeling
if(!exists('jhopth'))source('src/becker_modelling_hierarch.R')

################################################################################
########Okay so the joint model with a point estimate works okay, what about
########If we optimize a hiearach model?
################################################################################
{
jhsamp = rstan::sampling(jointmodel_hierach,
    data=jointdata,
    init=function(){jhopth$par},
    control=list(adapt_delta=.98,max_treedepth=15),iter=2e3,
    chains=4
)
#
jhsamp$par%>%.$l_pihalf%>%setNames(combinitvals$genes)%>%enframe('gene','l_pihalf')%>%
    inner_join(mcshanethalfs)%>%
    # filter(abs(l_pihalf) <  5)%>%
    filter(McShane_deg_cat!='NED')%>%
    # filter(l_pihalf<4)%>%
    {quicktest(.$l_pihalf,log(.$half_life));.}%>%
    {cor.test(.$l_pihalf,log(.$half_life))}%>%tidy
}
jhsamp$par$l_pihalf%>%which.min
jhsamp$par$theta%>%txtdensity



mgsampling<-rstan::optimizing(dp_model_ribnorm,data=get_dp_standata_withpriors(c(mgids),ribomatrscl=countmats$ribo,ribo_sigma=sigmas$ribo),
                              init=mginits$ribo,iter=10)

mgsampling<-rstan::sampling(dp_model_ribnorm,data=get_dp_standata_withpriors(c(mgids),ribomatrscl=countmats$ribo,ribo_sigma=sigmas$ribo),chains=4,
                            init=function(){mginits$ribo},control=list(adapt_delta=.98,max_treedepth=15),iter=2e3)