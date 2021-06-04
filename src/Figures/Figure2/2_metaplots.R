get_metasignaldf<-function(bindatamats,genelist=TRUE){
	metasignaldf = bindatamats%>%
		map_df(.id='sample',.%>%.[genelist,]%>%
			rownorm%>%
			# {./sum(.)}%>%
			colMeans(na.rm=T)%>%
			enframe('start','signal'))%>%
		group_by(sample)%>%
		mutate(section=c(rep('AUG',STARTWINDSIZE),'middle',rep('end',STOPWINDSIZE)))
	#assign stage
	metasignaldf$stage <- metasignaldf$sample%>%str_extract('[^_]+')
	# metasignaldf%<>%mutate(start = as.numeric(start) - ifelse(section=='AUG',STARTWINDSIZE+1,TOTBINS- STOPWINDSTART ))
	metasignaldf%<>%mutate(start = as.numeric(start) - ifelse(section=='AUG',1+FPEXT,1+FPEXT+STARTCDSSIZE+STOPCDSSIZE-2 ))
	metasignaldf
}

get_metaplot <- function(metasignaldf,ylims=NULL){
	if(!'fraction'%in%colnames(metasignaldf)) metasignaldf$fraction='total'
	if(!'sample'%in%colnames(metasignaldf)) metasignaldf$sample=metasignaldf$stage
	
	metasignaldf%>%
	filter(section!='middle')%>%
	split(.,.$section)%>%map( .%>%
		# slice_by(sample,8)%>%
		# filter(position%%3 == 0)%>%
		{
			isfirst = .$section[1]=='AUG'
			qplot(data=.,group=sample,color=stage,x=start,y=signal,geom='blank')+
			geom_line()+
			scale_x_continuous(name='position (bp)',
				limits=if(.$section[1]=='AUG') c(-FPEXT,STARTCDSSIZE) else c(-STOPCDSSIZE,TPUTREXT) ,
				minor_breaks=number_ticks,breaks=partial(number_ticks,n=12)
			)+
			facet_grid( fraction ~ section,scale='free_x')+
			scale_color_manual(values=stagecols)+
			scale_y_continuous(name='Mean Psite Count / CDS Total',limits=ylims)+
			theme_bw()+
			# theme_minimal()+
			theme(panel.grid = element_line(color=I('grey')))
		}
	)%>%
	ggarrange(plotlist=.,ncol=2,common.legend=T)
}

metasignaldf_stgrp <- get_metasignaldf(bindatamats[mainsamps],longcdstrs) %>% 
	group_by(stage,section,start)%>%
	summarise(signal=mean(signal))
{
'plots/Figures/Figure2/'%>%dir.create(showWarn=F,rec=T)
library(rlang)
plotfile<-'plots/Figures/Figure2/fig1c_myribowaltz_allsec_stageov.pdf'%T>%pdf(h=6,w=12)
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
# plotfile<-'plots/Figures/Figure2/fig1c_myribowaltz_frac_metaplots.pdf'%T>%pdf(h=6,w=12)
# rwplot <- metasignaldf_stgrp%>%get_metaplot(c(0,0.006))
# print(rwplot)
# dev.off()
# normalizePath(plotfile)%>%message
# }

#Let's define a subset of the data - those with the apparently biggest shift over time.
#Let's also define a 
inf.omit = function(x) x[is.finite(x)]

