library(multitaper)
library(ggpubr)

mmt<-function(vect){
	tapers=24
	bw=12
	if(length(vect)<25){slepians<-dpss(n=length(vect)+(50-length(vect)),k=tapers,nw=bw)}
	if(length(vect)>=25){slepians<-dpss(n=length(vect),k=tapers,nw=bw)}
	vals<-take_Fvals_spect(x = vect,n_tapers = tapers,time_bw = bw,slepians_values = slepians)
	pval<-pf(q=vals[1],df1=2,df2=(2*24)-2,lower.tail=F)
	return(list(vals[2],pval))
}          
pvect <- function(vlen,phase_props=c(5,2,1)){
	vect <- rep(0,vlen*3)

	phase_props <- phase_props%>%{./sum(.)}

	codonstrengths <- c(1:64)%>%sample(vlen,rep=T)%>%rep(each=3)

	vect_per <- phase_props%>%rep(vlen)%>%multiply_by(codonstrengths)%>%rpois(n=length(.),lambda=.)

	vect_per
}

# spec<-spec.mtm(ts(vect_per))

# spec$spec

counts = c(5,10,50,100,500)%>%setNames(.,.)

testreps<-map_df(.id='prob',counts,function(count){map_df(1:100,~{
	map_df(.id='periodic',c('periodic','aperiodic'),function(periodic){
	vect = pvect(300)
	if(!periodic=='periodic') vect = sample(vect)

	vect = rmultinom(n=1,size=count,prob=vect)
	

	poismod<-glm(vect~ rep_along(vect,0:2),family=poisson)
	nbmod<-glm.nb(vect~ rep_along(vect,0:2))
	mmt<-mmt(vect)

	data_frame(
		mtm=mmt[[1]],
		mtm_sig=mmt[[2]],
		fft=vect%>%fft%>%abs%>%tail(.,length(.)/(2/3))%>%max,
		pois=poismod%>%anova%>%as.list%>%.[[4]]%>%diff%>%abs,
		pois_sig=	poismod%>%summary%>%coef%>%.[2,4],
		negbin=nbmod%>%anova%>%as.list%>%.[[4]]%>%diff%>%abs,
		negbin_sig=nbmod%>%summary%>%coef%>%.[2,4],
		count=count
	)

})})})


testreps$pois%<>%pmax(1e-3)
testreps$negbin%<>%pmax(1e-3)
testreps$negbin_sig%<>%`<`(0.05)
testreps$mtm_sig%<>%`<`(0.05)

# testreps%>%split(.,.$prob)%>%lapply( .%>%qplot(data=.,mtm,fft,geom='point',log='xy'))%>%ggarrange(plotlist=.)
# testreps%>%split(.,.$prob)%>%lapply( .%>%qplot(data=.,mtm,negbin,color=periodic,geom='point',log='xy'))%>%ggarrange(plotlist=.)
# testreps%>%split(.,.$prob)%>%lapply( .%>%qplot(data=.,mtm,pois,color=periodic,geom='point',log='xy'))%>%ggarrange(plotlist=.)
X11(12,12)
ggarrange(
	testreps%>%qplot(data=.,x=mtm,y=negbin,color=negbin_sig,shape=periodic,geom='point',log='xy')+facet_grid(periodic~count,scale='free'),
	testreps%>%qplot(data=.,x=mtm,y=negbin,color=mtm_sig,shape=periodic,geom='point',log='xy')+facet_grid(periodic~count,scale='free'),
	nrow=2
)
dev.off()

with(testreps%>%filter(count==5),table(periodic,mtm_sig))
with(testreps%>%filter(count==5),table(periodic,negbin_sig))


with(testreps%>%filter(count==10),table(periodic,mtm_sig))
with(testreps%>%filter(count==10),table(periodic,negbin_sig))

with(testreps%>%filter(count==50),table(periodic,mtm_sig))
with(testreps%>%filter(count==50),table(periodic,negbin_sig))


table(testreps$periodic,testreps$negbin_sig)


testreps%>%qplot(data=.,mtm,pois,color=pois_sig,shape=periodic,geom='point',log='xy')+facet_grid(.~prob,scale='free')


testreps%>%qplot(data=.,mtm,pois,geom='point')

testreps%>%qplot(data=.,mtm,negbin,geom='point')+ggtitle('Multitaper Spectral Coefficient vs \nVariance Explained by Phase (glm.neg_bin)',sub=round(cor(testreps$mtm,testreps$negbin),3))
testreps%>%qplot(data=.,mtm,pois,geom='point')+ggtitle('Multitaper Spectral Coefficient vs \nVariance Explained by Phase (glm.poisson)',sub=round(cor(testreps$mtm,testreps$pois),3))





qplot(y=vect_per%>%fft%>%abs,x=seq_along(vect_per),geom='line')



vlens<-c(3e1,3e2,3e3,3e4+1,5e4,3e5)
vlens%>%setNames(.,.)

library(MASS)

bmarks<-map(vlens,function(vlen){
	vect = pvect(vlen)	
	microbenchmark::microbenchmark(times=5L,
		mtm=vect%>%spec.mtm(plot=F)%>%.$spec%>%tail(.,length(.)/(2/3))%>%max,
		fft=vect%>%fft%>%abs%>%tail(.,length(.)/(2/3))%>%max,
		pois=glm(vect~ rep_along(vect,0:2))%>%anova%>%as.list%>%.[[4]]%>%diff%>%abs,
		negbin=glm.nb(vect~ rep_along(vect,0:2))%>%anova%>%as.list%>%.[[4]]%>%diff%>%abs
	)%>%summary%>%as.data.frame%>%dplyr::select(expr,lq,median,uq)%>%mutate(N=length(vect))
})

bmarks%>%bind_rows%>%qplot(data=.,x=as.numeric(N),ymin=lq,ymax=uq,y=median,fill=expr,color=expr,geom=c('line','point'),log='xy')+theme_bw()+geom_ribbon(alpha=I(0.5))






