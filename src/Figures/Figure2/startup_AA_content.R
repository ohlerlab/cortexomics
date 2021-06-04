
################################################################################
########AA content in the start up vs down genes?
################################################################################
poseffectdf<-read_tsv('tables/ribo_position_effect.tsv')
AAnum = 10
cdsstartseq = cdsgrl%>%sort_grl_st%>%resize_grl(3*AAnum,'start')%>%extractTranscriptSeqs(x=fafileob)
cdsstartaas = cdsstartseq%>%translate


##analysis of aa frqeuency.
startuptrs = poseffectdf%>%filter(str_pvalue<0.05)%>%filter(str_lfc>0)%>%.$tr_id
notuptrs = poseffectdf%>%filter(str_pvalue>0.05)%>%.$tr_id
#p
stupmat = cdsstartaas[startuptrs]%>%as.matrix%>%.[,-1]%>%apply(2,table)
stnotupmat = cdsstartaas[notuptrs]%>%as.matrix%>%.[,-1]%>%apply(2,table)
aas = stupmat%>%rownames%>%setNames(.,.)
codpos = setNames(1:9,paste0('codon_',2:10))
#
#neg results - no enrichment of any particular 
aaenrich=map_df(.id='codon',codpos,function(cod){
	map_df(.id='AA',aas,function(aa){
		nonaarows = aas!=aa
		tablemat = matrix(c(
			sum(stnotupmat[nonaarows,cod]),
			stnotupmat[aa,cod],
			sum(stupmat[nonaarows,cod]),
			stupmat[aa,cod]),ncol=2)
		tidy(fisher.test(tablemat))
	})
})
#
aaenrichallcods=
	map_df(.id='AA',aas,function(aa){
		nonaarows = aas!=aa
		tablemat = matrix(c(
			sum(stnotupmat[nonaarows,1:4]),
			sum(stnotupmat[aa,1:4]),
			sum(stupmat[nonaarows,1:4]),
			sum(stupmat[aa,1:4])),ncol=2)
		# print(tablemat)
		tidy(fisher.test(tablemat))
	})
aaenrichallcods%>%filter(p.value<0.05)

# aaenrich%>%filter(p.value<0.05)%>%
# 	# select(codon,AA,estimate)%>%
# 	arrange(AA)%>%
# 	select(-method,-alternative)%>%
# 	as.data.frame%>%
# 	mutate(is_acidic= AA %in% c('D','E'))%>%
# 	mutate(is_basic= AA %in% c('R','K','H'))

aaenrichdf = aaenrichallcods%>%
	mutate(is_acidic= AA %in% c('D','E'))%>%
	mutate(is_basic= AA %in% c('R','K','H'))
aaenrichdf%<>%left_join(readxl::read_xlsx('tables/S3.xlsx',2)%>%filter(time=='E12.5')%>%select(dwell_time,AA))
#
pdf<-grDevices::pdf
#now plot
plotfile<- here(paste0('plots/','Figures/Figure3/AAenrichbarplots','.pdf'))
pdf(plotfile)
enrichbarplot = aaenrichdf%>%
	arrange(dwell_time)%>%mutate(AA=as_factor(AA))%>%
	ggplot(.,aes(fill=AA,x=AA,y=log2(estimate),ymin=log2(conf.low),ymax=log2(conf.high)))+
	stat_identity(geom='bar')+
	geom_errorbar()+
	scale_fill_discrete(name='Amino Acid')+
	scale_x_discrete(paste0('Amino Acid'))+
	scale_y_continuous(paste0('Enrichment in Start-Up Group'))+
	ggtitle(paste0('AA enrichment in first 2-5 amino acids of start up ORFs'))+
	theme_bw()
dtscatter = aaenrichdf%>%
	arrange(dwell_time)%>%mutate(AA=as_factor(AA))%>%
	ggplot(.,aes(x=AA,y=dwell_time,color=AA))+
	geom_point()+
	scale_color_discrete(name='Amino Acid')+
	scale_x_discrete(paste0('Amino Acid'))+
	scale_y_continuous(paste0('Dwell time at E12.5'))+
	ggtitle(paste0('Dwell time at E12.5'))+
	theme_bw()
print(ggarrange(plotlist=list(enrichbarplot,dtscatter),nrow=2))
dev.off()
message(normalizePath(plotfile))

#dwell times
readxl::read_xlsx('tables/S3.xlsx',2)%>%filter(time=='E12.5')%>%select(dwell_time,AA)%>%
	left_join(aaenrichallcods)%>%{quicktest(log2(.$estimate),.$dwell_time)}

################################################################################
########codon-codon basis
################################################################################

ntpos = 1+(0:4)*3
ntpos %<>% setNames(paste('codon_',((.-1)/3)))
stupmat = sapply(ntpos[-1],function(i){ cdsstartseq[startuptrs]%>%substr(i,i+2)%>%table})
stnotupmat = sapply(ntpos[-1],function(i){ cdsstartseq[notuptrs]%>%substr(i,i+2)%>%table})
codons = rownames(stupmat)%>%setNames(.,.)

aas = stupmat%>%rownames%>%setNames(.,.)
codpos = setNames(1:4,paste0('codon_',2:5))

pos = codpos[1]
codon = codons[1]
#neg results - no enrichment of any particular 
codenrich=map_df(.id='pos',codpos,function(pos){
	map_df(.id='codon',codons,function(codon){
		# message(pos)
		# message(codon)
		nonaarows = codons!=codon
		tablemat = matrix(c(
			sum(stnotupmat[nonaarows,pos]),
			stnotupmat[codon,pos],
			sum(stupmat[nonaarows,pos]),
			stupmat[codon,pos]),ncol=2)
		tidy(fisher.test(tablemat))
	})
})

codenrichallcods=
	map_df(.id='codon',codons,function(codon){
		nonaarows = codons!=codon
		tablemat = matrix(c(
			sum(stnotupmat[nonaarows,pos]),
			sum(stnotupmat[codon,pos]),
			sum(stupmat[nonaarows,pos]),
			sum(stupmat[codon,pos])),ncol=2)
		tidy(fisher.test(tablemat))
	})

codenrichallcods%>%filter(p.value<0.05)

codenrichdf = codenrichallcods%>%
	left_join(readxl::read_xlsx('tables/S3.xlsx',2)%>%filter(time=='E12.5')%>%select(dwell_time,AA,codon))%>%
	mutate(is_acidic= AA %in% c('D','E'))%>%
	mutate(is_basic= AA %in% c('R','K','H'))


AAorder = c('NFMQKRVWILPCHASGYETD')
AAorder=AAorder%>%str_split('')%>%.[[1]]%>%setNames(seq_along(.),.)
pdf<-grDevices::pdf
#now plot
plotfile<- here(paste0('plots/','Figures/Figure3/start_codon_enrichbarplots','.pdf'))
pdf(plotfile,w=14,h=14)
enrichbarplot = codenrichdf%>%
	arrange(AAorder[AA],dwell_time)%>%
	mutate(codon=as_factor(codon))%>%
	mutate(AA=as_factor(AA))%>%
	filter(!is.na(dwell_time))%>%
	ggplot(.,aes(fill=AA,x=codon,y=log2(estimate),ymin=log2(conf.low),ymax=log2(conf.high)))+
	stat_identity(geom='bar')+
	geom_errorbar()+
	scale_fill_discrete(name='codon')+
	scale_x_discrete(paste0('Amino Acid'))+
	scale_y_continuous(paste0('Enrichment in Start-Up Group'))+
	facet_grid(.~AA,scale='free_x')+
	ggtitle(paste0('codon enrichment in amino acids 2-5 of start up ORFs'))+
	theme_bw()
dtscatter = codenrichdf%>%
	arrange(AAorder[AA],dwell_time)%>%
	mutate(AA=as_factor(AA))%>%
	mutate(codon=as_factor(codon))%>%
	filter(!is.na(dwell_time))%>%
	ggplot(.,aes(x=codon,y=dwell_time,color=AA))+
	geom_point()+
	scale_color_discrete(name='Amino Acid')+
	scale_x_discrete(paste0('Amino Acid'))+
	scale_y_continuous(paste0('Dwell time at E12.5'))+
	facet_grid(.~AA,scale='free')+
	ggtitle(paste0('Dwell time at E12.5'))+
	theme_bw()
print(ggarrange(plotlist=list(enrichbarplot,dtscatter),nrow=2))
dev.off()
message(normalizePath(plotfile))



# aaenrich%>%filter(p.value<0.05)%>%
# 	# select(codon,AA,estimate)%>%
# 	arrange(AA)%>%
# 	select(-method,-alternative)%>%
# 	as.data.frame%>%
# 	mutate(is_acidic= AA %in% c('D','E'))%>%
# 	mutate(is_basic= AA %in% c('R','K','H'))

# aaenrichdf = aaenrichallcods%>%
# 	mutate(is_acidic= AA %in% c('D','E'))%>%
# 	mutate(is_basic= AA %in% c('R','K','H'))
