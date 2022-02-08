

metasignaldf_stgrp <- get_metasignaldf(bindatamats[mainsamps],longcdstrs) %>% 
	group_by(stage,section,start)%>%
	summarise(signal=mean(signal))
{
'plots/Positional_Analysis/'%>%dir.create(showWarn=F,rec=T)
library(rlang)
plotfile<-'plots/Positional_Analysis/fig1c_myribowaltz_allsec_stageov.pdf'%T>%pdf(h=6,w=12)
rwplot <- metasignaldf_stgrp%>%get_metaplot(c(0,0.012))
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message
}

# bindatamats%>%names
# metasignaldf_stgrp <- get_metasignaldf(bindatamats[nonmainsamps],longcdstrs) %>% 
# 	mutate(fraction = str_extract(sample,'80S|Poly'))%>%
# 	mutate(stage = sample%>%str_extract('(?<=80S|Poly).*?(?=_)'))%>%
# 	group_by(stage,section,start,fraction)%>%
# 	summarise(signal=mean(signal))
# {
# library(rlang)
# plotfile<-'plots/QC_plots/fig1c_myribowaltz_frac_metaplots.pdf'%T>%pdf(h=6,w=12)
# rwplot <- metasignaldf_stgrp%>%get_metaplot(c(0,0.006))
# print(rwplot)
# dev.off()
# normalizePath(plotfile)%>%message
# }

inf.omit = function(x) x[is.finite(x)]

