mytpnames = c(E13 = "E12.5", E145 = "E14", E16 = "E15.5", E175 = "E17", P0 = "P0")
myassaynames = c('MS'='Mass Spec','total'='RNAseq','ribo'='Ribo-Seq','TE'='TE')
genenamelist =  c('Satb2')

{
make_traj_plot <- function(gene2plot,msdf,postmeanmat,postprecsdmat,exprdf,nonredgnames,cds,msrescale2lfq,best_uprotein_ids,ms_id2protein_id,prediction_df,tpnames,assaynames){
	plotlistlist = list()
	genenamelist = list()
	u = list('Bcl11b'=c(-2,4),'Flna'=c(-4,1),'Nes'=c(-4,1),'Satb2'=c(-1,5))
	

	assert_that(gene2plot %in% nonredgnames$gene_name)
	testpids <- cds%>%subset(gene_name==gene2plot)%>%.$protein_id%>%unique
	test_uids<-unique(exprdf$uprotein_id)%>%str_subset(testpids%>%paste0(collapse='|'))
	#get data for that test gene
	ggdf <- exprdf%>%filter(uprotein_id%in%test_uids)
	ggdf%<>%bind_rows(ggdf%>%filter(assay=='ribo')%>%mutate(assay='TE',signal = signal - (ggdf%>%filter(assay=='total')%>%.$signal)))
	postmeanmat_scaled <- postmeanmat - msrescale2lfq
	ggdf_msconf <- 
		safe_left_join(
			((postmeanmat_scaled[test_uids,,drop=F])-(1.96*postprecsdmat[test_uids,,drop=F]))%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,CI.L,-uprotein_id)%>%separate(dataset,c('time','assay')),
			((postmeanmat_scaled[test_uids,,drop=F])+(1.96*postprecsdmat[test_uids,,drop=F]))%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,CI.R,-uprotein_id)%>%separate(dataset,c('time','assay'))
		)%>%
		safe_left_join(
			((postmeanmat_scaled[test_uids,,drop=F]))%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,signal,-uprotein_id)%>%separate(dataset,c('time','assay')),
			)
	# ggdf_msconf$time%<>%factor%>%time
	#get ms-protein id pairs
	test_uids<-ggdf$uprotein_id%>%unique
	scaledata = function(x,assay2plot){
		if(assay2plot=='MS'){
			out = x%>%mutate_at(vars(one_of(c('signal','logFC','CI.R','CI.L'))),list(function(x)x+msrescale2lfq))
		}else{
			out = x	
		}
		rescale = if('signal' %in% colnames(x)){ median(x$signal[x$time=='E13']) } else {median(x$logFC[x$time=='E13'])}
		out = x %>%mutate_at(vars(one_of(c('signal','logFC','CI.R','CI.L'))),list(function(x)x-rescale))
		out
	}
	# plotdf%>%filter(assay==assay2plot)%>%scaledata(assay2plot='MS')
	# scaledata = identity
	#
	trajfile = './plots/tmp.pdf'
	for(testuid in test_uids){
		testuid %in% best_uprotein_ids
		uidfilt<-.%>%filter(uprotein_id==testuid)
		ms_id2protein_id%>%filter(uprotein_id==testuid)
		plotdf<-ggdf%>%uidfilt
		assays2plot <- unique(ggdf$assay)%>%sort%>%rev%>%.[order(.=='TE')]
		trajectoryplots<-lapply(assays2plot,function(assay2plot){
			if(assay2plot!='MS'){
				points = geom_point()
				linerange=NULL
			}else{
				points = geom_point(data=msdf%>%semi_join(ms_id2protein_id%>%filter(gene_name==gene2plot))%>%scaledata('dont unscale'))
				linerange=geom_linerange(color=I('blue'),data=ggdf_msconf%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),aes(y=signal,ymin=CI.L,ymax=CI.R))
			}
			ggplot(
			data = plotdf%>%filter(assay==assay2plot)%>%scaledata(assay2plot),
			aes(
				x=as.numeric(as_factor(time)),
				y=signal
			))+
			points+linerange+
			geom_ribbon(data=prediction_df%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),aes(x=as.numeric(as_factor(time)),y=logFC,ymin=CI.L,ymax=CI.R),fill='darkgreen',alpha=I(0.5))+	
			geom_line(data=prediction_df%>%uidfilt%>%filter(assay==assay2plot)%>%scaledata(assay2plot),linetype=2,aes(x=as.numeric(as_factor(time)),y=logFC))+	
			# geom_line(data=prediction_df%>%uidfilt,aes(x=time,y=logFC))+	
			scale_x_continuous(name='Stage',labels=tpnames)+
			theme_bw()+
			ggtitle(assaynames[assay2plot])
			# facet_wrap( ~ assay,scales='free')+
		})
		# trajectoryplots<-commonyspanplots(trajectoryplots,breakint=0.25,minrange=8)
		trajectoryplots[c(1:3)]%<>%lapply(function(plot) plot+ coord_cartesian(ylim=u[[gene2plot]]))
		trajectoryplots[c(4)]%<>%lapply(function(plot) plot+ coord_cartesian(ylim=c(-3,3)))

		plotlistlist = append(plotlistlist,list(trajectoryplots))
		genenamelist = append(genenamelist,gene2plot)
	}
	
	trajectoryplot<-ggarrange(plotlist=plotlistlist,ncol=4)
	trajectoryplot<-annotate_figure(trajectoryplot,top  = str_interp('Data vs Linear Model - ${testname}'))
	trajectoryplot
}
}

plotlistlist<-purely(make_traj_plot)(gene2plot='Satb2',
	nonredgnames=nonredgnames,
	cds=cds,
	msdf=msdf,
	postmeanmat=postmeanmat,
	postprecsdmat=postprecsdmat,
	exprdf=exprdf,
	tpnames=mytpnames,
	msrescale2lfq,
	best_uprotein_ids,
	ms_id2protein_id,
	prediction_df,
	assaynames=myassaynames
)

plotlistlist%<>%setNames(genenamelist)

for(testname in genenamelist){
	trajfile <- str_interp('plots/figures/figure2/traject_${testname}.pdf')
	# plotlistlist%<>%mapply(HARDCODELIMS,FUN=function(plotlist,lims) plotlist = lapply(plotlist,function(plot) plot+ coord_cartesian(ylim=lims)))
	pdf(trajfile,w=12,h=4)
	print(trajectoryplot)
	dev.off()
	system(str_interp('cp ${trajfile} plots/figures/figure2/traject_${testname}.pdf'))
	message(normalizePath(trajfile))
	message(normalizePath(str_interp('plots/figures/figure2/traject_${testname}.pdf')))
}

exprdf%<>%left_join(ms_id2protein_id%>%distinct(uprotein_id,gene_name))
exprdf%<>%mutate(protein_id=str_replace(uprotein_id,'_\\d+$',''))

save(
	cds,
	msdf,
	postmeanmat,
	postprecsdmat,
	exprdf,
	mytpnames,
	msrescale2lfq,
	best_uprotein_ids,
	prediction_df,
	myassaynames,file='data/trajplotobjects.Rdata'
)