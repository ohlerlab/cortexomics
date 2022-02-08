################################################################################
########Run tests for positional effects with Ribodipa
################################################################################

if(!exists('patterngenes')) patterngenes = longcdstrs	
# use main samples only
mainbindatamats = bindatamats[mainsamps]
#compress our bins into 4 codon units
data_binned <- mainbindatamats%>%lapply(function(mat){cbind(
	mat%>%.[,1:(FPEXT+STOPCDSSIZE)]%>%codmerge%>%codmerge(4),
	mat%>%.[,FPEXT+STOPCDSSIZE+1,drop=FALSE],
	mat%>%.[,(FPEXT+STOPCDSSIZE+2):TOTBINS,drop=FALSE]%>%codmerge%>%codmerge(4)
)})
data_binned<- abind(data_binned,along=0)

#name the compressed dimensions
startbinnum = (FPEXT+STOPCDSSIZE)/12
stopbinnum = (TPUTREXT+STOPCDSSIZE)/12
dimnames(data_binned)[[3]] <- c(paste0('start_',1:startbinnum),'mid',paste0('stop_',1:stopbinnum))
data_binned%<>%{lapply(matgenes%>%setNames(.,.),function(i).[,i,])}

#add in an extra row so that ribo dipa works, (it crashes otherwise...)
highc = data_binned[patterngenes]%>%map_dbl(sum)%>%which.max
highc = patterngenes[highc]
nullmat = matrix(0,nrow=length(mainbindatamats))
data_binned[[highc]] = cbind(data_binned[[highc]],nullmat)
data_binned[[highc]]%<>%{colnames(.)[ncol(.)]=ncol(.)+1;.}

# patterngenes = c(moststart_decr,moststart_incr)
ribodipadata = data_binned[patterngenes]
ribodipadata = data_binned%>%map(~.[mainsamps,])

#put some fake input in here
classlabel = data.frame(condition = mainsamps%>%str_extract('(E13|E145)|(E16|E175|P0)'))%>%
		mutate(comparison = ifelse(condition %in% c('E16','E175','P0'),2,1))

#now
if(!file.exists(here('data/result_ribodipa.rds'))){
	result_ribodipa  <- RiboDiPA::diffPatternTest(
		data = ribodipadata, 
		classlabel = classlabel, 
		method=c('gtxr', 'qvalue')
		)
	saveRDS(result_ribodipa,here('data/result_ribodipa.rds'))
}else{
	result_ribodipa<-readRDS(here('data/result_ribodipa.rds'))
}

ribodipares = result_ribodipa$gene%>%
	as.data.frame%>%rownames_to_column('tr_id')%>%
	# mutate(class = case_when(
		# tr_id%in%moststart_incr~'startincr',
		# tr_id%in%moststart_decr~'startdecr',
		# TRUE~'other')
	# )
	identity

ribodipares%<>%arrange(qvalue)
# ribodipares%>%head
# ribodipares%>%filter(qvalue<0.05)%>%.$tvalue%>%log10%>%txtdensity
# ribodipares%>%filter(qvalue>0.05)%>%.$tvalue%>%log10%>%inf.omit%>%txtdensity

dipabinres = result_ribodipa$bin
dipabinres%<>%map_df(.id='tr_id',.%>%as.data.frame%>%filter(!is.na(pvalue))%>%rownames_to_column('bin'))
dipabinres%<>%mutate(bin = bin%>%str_replace('\\.\\d+',''))
dipabinres%<>%arrange(bin%>%str_detect('end'),bin%>%str_detect('mid'))

#PLot the spatial distribution of the ribodipa fold changes
plotfile<- here(paste0('plots/','startlog2foldchange','.pdf'))
pdf(plotfile,w=10,h=5)
dipabinres%>%
	filter(bin%>%str_detect('start|stop'))%>%group_by(bin)%>%
	summarise(log2FoldChange=median(log2FoldChange))%>%
	mutate(section = str_extract(bin,'start|stop'))%>%
	mutate(position = ifelse(section=='start',
		((as.numeric(str_extract(bin,'\\d+$'))-1)*12)-FPEXT,
		((as.numeric(str_extract(bin,'\\d+$'))-1)*12)-STARTCDSSIZE
	))%>%
	{split(.,.$section)}%>%
	map( .%>%
	{ggplot(data=.,aes(x=position,y=log2FoldChange))+
	geom_point()+
	scale_x_continuous(name='position (bp)',
				limits=if(.$bin[1]%>%str_detect('start')) c(-FPEXT,STARTWINDSIZE) else c(-STOPCDSSIZE,TPUTREXT) ,
				# limits=if(isfirst) c(0,(STARTWINDSIZE)-1) else c(1-(STOPWINDSIZE),0) ,
				minor_breaks=number_ticks,breaks=partial(number_ticks,n=12)
			)+	
	scale_y_continuous('Average Binned log2FoldChange(relative Density)')+
	ggtitle('Average Ribodipa Fold Change, for Binned Positions')+
	theme_bw()})%>%
	ggarrange(plotlist=.,ncol=2,common.legend=T)
dev.off()
message(normalizePath(plotfile))
#
dipabinres$gene_id = trid2gid[[dipabinres$tr_id]]
dipabinres$is_sig <- dipabinres$pvalue<0.05

#get a data frame of transcripts with the start_up effect
poseffectdf = left_join(
	dipabinres%>%filter(bin==paste0('start_',1+(FPEXT/12)))%>%
		select(tr_id,str_pvalue=pvalue,str_lfc = log2FoldChange),
	dipabinres%>%filter(bin==paste0('stop_',STOPCDSSIZE/12))%>%
		select(tr_id,end_pvalue=pvalue,end_lfc = log2FoldChange),
	by='tr_id'
)
#
poseffectdf$gene_id<-trid2gid[[poseffectdf$tr_id]]
poseffectdf%>%write_tsv('tables/ribo_position_effect.tsv')

dirfuncs =list(up = function(x) x > 0,down=function(x) x < 0,both=TRUE)
dirname='up'
my_ontology <- 'BP'


run_go_ensembl <- function(geneList,my_ontology){
set.seed(1234)
require(org.Mm.eg.db)
require(DBI)
require(topGO)

table(geneList)

# Create topGO object

GOdata <-
  new(
    "topGOdata",
    ontology = my_ontology,
    allGenes = geneList,
    description = "Test",
    annot = annFUN.org,
    mapping = "org.Mm.eg.db",
    ID = "ENSEMBL"
  )
GOdata
}


if(!file.exists(here('data/poseffgores.rds'))){
		#Go term analysis of Ribodipa effects
	poseffgores <- imap(dirfuncs,function(dirfunc,dirname){
		godata<-poseffectdf%>%
			mutate(is_sig = (str_pvalue<0.05)&(dirfunc(str_lfc)))%>%
			{setNames(as.factor(as.numeric(.$is_sig)),.$gene_id)}

		stopifnot((c(0,1)%in%names(table(godata))))

		goob<-run_go_ensembl(godata,my_ontology)
		# goob$result
		results     <- runTest(goob, algorithm = 'elim', statistic = "fisher")
		results.tab <- GenTable(object = goob, elimFisher = results,topNodes = 100)
		results.tab%<>%mutate(Enrichment = Significant / Expected )
		results.tab%<>%mutate(elimFisher = elimFisher%>%str_replace("<","")%>%as.numeric )
		results.tab%<>%mutate(elimFisher = as.numeric(elimFisher) )
		results.tab%<>%mutate(direction=dirname)
		results.tab%<>%mutate(ont=my_ontology)
		results.tab	
	})
	saveRDS(poseffgores,here('data/poseffgores.rds'))
}else{
	poseffgores<-readRDS(here('data/poseffgores.rds'))
}


poseffectdf%<>%mutate(gene_name= gid2gnm[[gene_id]])
poseffectdf%>%filter(gene_name%>%is_in(genesofinterest))
poseffgores[[1]]%>%head(20)%>%write_tsv('tables/startup_poseffect_go.tsv')
poseffgores[[2]]%>%mutate(direction='down')%>%head(20)%>%write_tsv('tables/startdown_poseffect_go.tsv')

tubgenes = GTOGO%>%filter(go_id=='GO:1904528')%>%.$ensembl_gene_id

 p = poseffgores[[1]]%>%arrange(-elimFisher)%>%tail(10)%>%
		 	mutate(Term=as_factor(Term))%>%
		 	ggplot(aes(x=Term,y=-log10(elimFisher),fill=Enrichment))+
		 	stat_identity(geom='bar')+
		 	geom_text(aes(label=GO.ID,y=1.5,color=I('white')))+
		 	scale_fill_continuous(low='lightgrey',high='black')+
		 	# scale_fill_continuous()+
		 	scale_x_discrete('')+
		 	coord_flip()+
		 	theme(axis.text.x=element_text(vjust=-1))+
		 	theme_bw()
		stopifnot(length(my_catnm)==1)
		plotfile <- paste0('plots/Positional_Analysis/ggoplot_',as.character(my_ontology),'_',my_stat,'_','start_poseff','_','up','_','.pdf')
		pdf(plotfile,w=12,h=4)
		print(p)
		dev.off()
		message(normalizePath(plotfile))

 p = poseffgores[[2]]%>%arrange(-elimFisher)%>%tail(10)%>%
		 	mutate(Term=as_factor(Term))%>%
		 	ggplot(aes(x=Term,y=-log10(elimFisher),fill=Enrichment))+
		 	stat_identity(geom='bar')+
		 	geom_text(aes(label=GO.ID,y=1,color=I('white')))+
		 	scale_x_discrete('')+
		 	scale_fill_continuous(low='grey',high='darkgrey')+
		 	coord_flip()+
		 	theme(axis.text.x=element_text(vjust=-1))+
		 	theme_bw()
		stopifnot(length(my_catnm)==1)
		plotfile <- paste0('plots/Positional_Analysis/ggoplot_',as.character(my_ontology),'_',my_stat,'_','start_poseff','_','down','_','.pdf')
		pdf(plotfile,w=12,h=4)
		print(p)
		dev.off()
		message(normalizePath(plotfile))

poseffectdf<-read_tsv('tables/ribo_position_effect.tsv')

ribodipa_up = poseffectdf%>%filter((str_pvalue<0.05)&((str_lfc)>0))%>%.$gene_id
ribodipa_up = matgenes[trid2gid[[matgenes]]%>%is_in(ribodipa_up)]
metasignaldf_stgrp <- get_metasignaldf(mainbindatamats,ribodipa_up) %>% 
	group_by(stage,section,start)%>%
	summarise(signal=mean(signal)) 
{
library(rlang)
plotfile<-'plots/Positional_Analysis/metaplottesribodipaup.pdf'%T>%pdf(h=6,w=12)
rwplot <- metasignaldf_stgrp%>%get_metaplot(c(0,0.012))
print(rwplot+ggtitle('ribodipa - sig start up genes'))
dev.off()
normalizePath(plotfile)%>%message
}
ribodipa_notup = poseffectdf%>%filter((between(str_lfc,-0.3,0.3)))%>%.$gene_id
ribodipa_notup = matgenes[trid2gid[[matgenes]]%>%is_in(ribodipa_notup)]
metasignaldf_stgrp <- get_metasignaldf(mainbindatamats,ribodipa_notup) %>% 
	group_by(stage,section,start)%>%
	summarise(signal=mean(signal)) 
{
library(rlang)
plotfile<-'plots/Positional_Analysis/metaplottesribodipa_notup.pdf'%T>%pdf(h=6,w=12)
rwplot <- metasignaldf_stgrp%>%get_metaplot(c(0,0.012))
print(rwplot+ggtitle('ribodipa - start nochange genes'))
dev.off()
normalizePath(plotfile)%>%message
}


################################################################################
########Also do metaplots for the mono and poly
################################################################################

start_up_trs = poseffectdf%>%filter(str_lfc>0,str_pvalue<0.05)%>%.$tr_id
start_notup_trs = poseffectdf%>%filter((between(str_lfc,-0.3,0.3)))%>%.$tr_id
metasignaldf_stgrp<-bind_rows(.id='section',list(
		start_up_trs=get_metasignaldf(bindatamats[nonmainsamps],start_up_trs)%>%
			filter(section=='AUG')%>%select(-section),
		start_notup_trs =get_metasignaldf(bindatamats[nonmainsamps],start_notup_trs)%>%
			filter(section=='AUG')%>%select(-section)
	))
metasignaldf_stgrp <- metasignaldf_stgrp %>% 
	mutate(fraction = str_extract(sample,'80S|Poly'))%>%
	mutate(stage = sample%>%str_extract('(?<=80S|Poly).*?(?=_)'))%>%
	group_by(stage,section,start,sample,fraction)%>%
	summarise(signal=mean(signal))	

{
library(rlang)
plotfile<-'plots/Positional_Analysis/metaplottesribodipa_notup_polymono.pdf'%T>%pdf(h=6,w=12)
rwplot <- metasignaldf_stgrp%>%get_metaplot()
print(rwplot+ggtitle('ribodipa - start nochange genes'))
dev.off()
normalizePath(plotfile)%>%message
}

######################################################################
########And high/low transcriptional change
################################################################################

metasignaldf_stgrp <-  bind_rows(.id='fraction',list(
		dtxn_topq=get_metasignaldf(bindatamats[mainsamps],dtxn_topqtrs),
		dtxn_midqs =get_metasignaldf(bindatamats[mainsamps],dtxn_midqstrs),
		dtxn_lowqs=get_metasignaldf(bindatamats[mainsamps],dtxn_lowqstrs)
	))%>% 
	group_by(stage,section,start,fraction)%>%
	summarise(signal=mean(signal))
metasignaldf_stgrp$fraction%<>%factor(levels=c('dtxn_topq','dtxn_midqs','dtxn_lowqs'))
{
library(rlang)
plotfile<-'plots/dtxn_quart_metaplot.pdf'%T>%pdf(h=6,w=12)
rwplot <- metasignaldf_stgrp%>%get_metaplot()+facet_grid(scale='free')
print(rwplot)
dev.off()
normalizePath(plotfile)%>%message
}



################################################################################
########Test associations with ebp1 stuff, TE change
################################################################################
allxtail = Sys.glob('pipeline/xtail/*')%>%map_df(.id='time',fread)%>%group_by(gene_name)
allxtail$gene_id = gnm2gid[[ allxtail$gene_name]]
techange = allxtail%>%filter(time==3)%>%select(gene_id,log2fc,adj_p_value)

poseffectdf%>%filter(!is.na(str_pvalue),!is.na(end_pvalue))%>%{quicktest(.$str_lfc,.$end_lfc)}

#
#
library(ggrepel)
posvtedf = poseffectdf%>%mutate(gene_id=trid2gid[[tr_id]])%>%
	select(-gene_name)%>%
	filter((str_pvalue)<0.05)%>%
	left_join(
		techange%>%
		filter(adj_p_value<0.05)%>%
		identity
	,by='gene_id')%>%
	filter(log2fc>-5)
posvtedf%<>%mutate(glabl = ifelse(
	(str_lfc < 0) | str_lfc>2.5,
	gene_name,
	''
))
#now plot
plotfile<- here(paste0('plots/','start_up_vs_TEchange','.pdf'))
pdf(plotfile)
	# filter(padj<0.05)%>%
	ggplot(data=posvtedf,aes(str_lfc,log2fc,label=glabl))+
	geom_label_repel()+
	geom_point()+
	geom_smooth(method='lm')+
	scale_x_continuous(paste0('Relative Start Usage'))+
	scale_y_continuous(paste0('Overal TE Change'))+
	ggtitle(
		paste0('Start Usage vs. TE Change'),
		sub = posvtedf%>%{cor.test(.$str_lfc,.$log2fc)}%>%tidy%>%mutate_if(is.numeric,round,5)%>%
		# {str_interp('pearsons rho =\n ${.$conf.low} - ${.$conf.high};p =${.$p.value} ')}
		{paste0('rho = ',round(.$conf.low,3),' - ',round(.$conf.high,3),'\np = ',format(.$p.value,sci=T,digits=2))}%>%
		str_replace('0e\\+00',' < 1.0e-5')
	)+
	theme_bw()
dev.off()
message(normalizePath(plotfile))



##Now plot relationship with TE
{
tevals = iso_tx_countdata$abundance%>%as.data.frame%>%
	rownames_to_column('tr_id')%>%
	gather(dataset,val,-tr_id)%>%
	separate(dataset,into=c('time','assay','replicate'))%>%
	group_by(time,assay,tr_id)%>%summarise(val=mean(val))%>%
	spread(assay,val)%>%
	mutate(te = log2(ribo) - log2(total))%>%
	select(tr_id,time,ribo,total,te)

startwinds = trcds%>%resize(12,'start')%>%subset(seqnames%in%names(psitecovnorm[[1]]))
startpsitedens = psitecovnorm%>%
	map_df(.id='sample',~.[startwinds]%>%sum%>%enframe('tr_id','start_dens'))
startpsitedens%<>% separate(sample,into=c('time','assay','replicate'))
startpsitedens%<>%select(-assay)
startpsitedens%<>%group_by(tr_id,time)%>%summarise(start_dens=mean(start_dens))

#now plot
plotfile<- here(paste0('plots/','startprop_te_scatter','.pdf'))
pdf(plotfile)
tevals%>%
	inner_join(poseffectdf%>%filter(str_pvalue<0.05,str_lfc<0))%>%
	filter(time=='E13')%>%
	filter(between(te,-12,12))%>%
	left_join(startpsitedens)%>%
	filter(start_dens>0)%>%
	# filter(start_count>1)%>%
	group_by(time)%>%
	filter(ribo>quantile(ribo,0.25))%>%
	filter(total>quantile(total,0.25))%>%
	# filter(ribo>32)%>%
	# filter(log2(ribo)%>%between(-12,50))%>%
	# {quicktest(.$te,log2(.$start_sig))}%>%
	ggplot(.,aes(te,log2(start_dens)))+
	geom_point()+
	facet_grid(time~.)+
	geom_smooth(method='lm')+
	# scale_color_discrete(name='colorname',colorvals)+
	scale_x_continuous(paste0('TE'))+
	scale_y_continuous(paste0('Log2(start/CDS)'))+
	ggtitle(paste0('title'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))
}

#now plot
plotfile<- here(paste0('plots/Positional_Analysis/','te_startup_density','.pdf'))
pdf(plotfile)
tevals%>%
	# filter(time=='E13')%>%
	left_join(poseffectdf%>%filter(str_pvalue<0.05,str_lfc>0)%>%mutate(startup=TRUE))%>%
	mutate(startup = ifelse(is.na(startup),FALSE,TRUE))%>%
	filter(is.finite(te))%>%
	filter(between(te,-12,12))%>%
	ggplot(aes(x=te,fill=startup))+
	geom_density(,alpha=I(0.5))	+
	facet_grid(time~.)
	# {split(.$te,.$startup)}%>%
	# {t.test(.[[1]],.[[2]])}
dev.off()
message(normalizePath(plotfile))


#now plot
plotfile<- here(paste0('plots/Positional_Analysis/','te_startdown_density','.pdf'))
pdf(plotfile)
tevals%>%
	# filter(time=='E13')%>%
	left_join(poseffectdf%>%filter(str_pvalue<0.05,str_lfc<0)%>%mutate(startdown=TRUE))%>%
	mutate(startdown = ifelse(is.na(startdown),FALSE,TRUE))%>%
	filter(is.finite(te))%>%
	filter(between(te,-12,12))%>%
	ggplot(aes(x=te,fill=startdown))+
	geom_density(,alpha=I(0.5))	+
	facet_grid(time~.)
	# {split(.$te,.$startdown)}%>%
	# {t.test(.[[1]],.[[2]])}
dev.off()
message(normalizePath(plotfile))



##Okay so maybe the read length distributions at the starts differ over time?
uptrs = poseffectdf%>%filter(str_pvalue<0.05,str_lfc>0)%>%.$tr_id
startwinds = trcds%>%resize(12,'start')
startwinds = startwinds[uptrs]

if(!file.exists(here('data/startrlcountdf.rds'))){

	startrlcountdf = fpcovlist[mainsamps]%>%
		mclapply(mc.cores=4,.%>%
			map_df(.id='length',~.[startwinds]%>%
				sum%>%
				enframe('tr_id','count')
			)
		)%>%bind_rows(.id='sample')
	saveRDS(startrlcountdf,here('data/startrlcountdf.rds'))

}else{
	startrlcountdf<-readRDS(here('data/startrlcountdf.rds'))
}


#now plot
plotfile<- here(paste0('plots/Positional_Analysis/','startrldistplot','.pdf'))
pdf(plotfile,w=15,h=3)
startrlcountdf%>%group_by(length,sample)%>%summarise(count=sum(count))%>%
	group_by(sample)%>%mutate(freq = count/sum(count))%>%
	separate(sample,c('time','assay','rep'))%>%
	select(-assay)%>%	
	ggplot(.,aes(x=length,freq))+
	# stat_identity(geom='bar',position='dodge')+
	geom_point(geom='bar',position='dodge')+
	scale_x_discrete(paste0('read_length'))+
	scale_y_continuous(paste0('freq'))+
	facet_grid(.~time,scale='free_y')+
	ggtitle(paste0('Read length dists for samples at start window'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))




#
posvtedf = poseffectdf%>%
	select(-gene_name)%>%
	mutate(gene_id=trid2gid[[tr_id]])%>%
	filter((end_pvalue)<0.05)%>%
	left_join(
		techange%>%
		filter(adj_p_value<0.05)%>%
		identity
	,by='gene_id')%>%
	filter(log2fc>-5)
#now plot
plotfile<- here(paste0('plots/Positional_Analysis/','end_up_vs_TEchange','.pdf'))
pdf(plotfile)
	# filter(padj<0.05)%>%
	ggplot(data=posvtedf,aes(end_lfc,log2fc))+
	geom_point()+
	geom_smooth(method='lm')+
	scale_x_continuous(paste0('Relative Start Usage'))+
	scale_y_continuous(paste0('Overal TE Change'))+
	ggtitle(
		paste0('Start Usage vs. TE Change'),
		sub = posvtedf%>%{cor.test(.$str_lfc,.$log2fc)}%>%tidy%>%mutate_if(is.numeric,round,5)%>%
		# {str_interp('pearsons rho =\n ${.$conf.low} - ${.$conf.high};p =${.$p.value} ')}
		{paste0('rho = ',round(.$conf.low,3),' - ',round(.$conf.high,3),'\np = ',format(.$p.value,sci=T,digits=2))}%>%
		str_replace('0e\\+00',' < 1.0e-5')
	)+
	theme_bw()
dev.off()
message(normalizePath(plotfile))




poseffectdf%>%mutate(gene_id=trid2gid[[tr_id]])%>%
	filter(str_pvalue<0.05)%>%
	left_join(ebp1_ko_startstop,by='gene_id')%>%
	filter(log2fold_si_mock%>%between(-5,5))%>%
	filter(padj<0.05)%>%
	{quicktest(.$str_lfc,.$log2fold_si_mock)}

poseffectdf%>%mutate(gene_id=trid2gid[[tr_id]])%>%
	filter(str_pvalue<0.05)%>%
	left_join(ebp1_ip_dex%>%mutate(gene_id=gnm2gid[[gene_name]]),by='gene_id')%>%
	filter(log2fold_Total_IP%>%between(-5,5))%>%
	{quicktest(.$str_lfc,.$log2fold_Total_IP)}

poseffectdf%>%mutate(gene_id=trid2gid[[tr_id]])%>%
	left_join(ebp1_ko_startstop,by='gene_id')%>%
	filter(log2fold_si_mock%>%between(-5,5))%>%
	{quicktest(.$end_lfc,.$log2fold_si_mock)}

poseffectdf%>%mutate(gene_id=trid2gid[[tr_id]])%>%
	left_join(ebp1_ip_dex%>%mutate(gene_id=gnm2gid[[gene_name]]),by='gene_id')%>%
	filter(log2fold_Total_IP%>%between(-5,5))%>%
	{quicktest(.$end_lfc,.$log2fold_Total_IP)}

ebp1_ko_startstop = read_tsv("../Ebp1_ribo/tables/start_change_tbl.tsv")
ebp1_ip_dex = read_tsv('../Ebp1_ribo/data/ip_dexdf_s100.tsv')



#now plot
plotfile<- here(paste0('plots/Positional_Analysis/','ribodipares','.pdf'))
pdf(plotfile)
plotTest(result = result.pst, genes.list = ribodipares$tr_id%>%head(2), threshold = 0.05)[[1]]+scale_x_continuous(limits=c(1,50))
dev.off()
message(normalizePath(plotfile))

lapply(fpcovlist)
fpcovlist[[1]][[1]]%>%as("GRanges")%>%subset(score!=0)%>%shift(1)%>%coverage(weight='score')


fpcovtrs = fpcovlist[[1]][[1]]%>%names

ribocovtrs = ribocovtrs%>%intersect(fpcovtrs)

negcountdf = tibble(
	tr_id=names(cdsstartaas),
	negcount = cdsstartaas%>%as.character%>%str_count('R|K')
	)
#
startgrpdf = poseffectdf%>%inner_join(negcountdf)%>%
	mutate(startgrp = case_when(
		  (str_pvalue<0.05) & (str_lfc>0) ~ 'start_up',
		start_down = (str_pvalue<0.05) & (str_lfc<0) ~ 'start_down',
		TRUE ~ 'nochange'
	))
startgrpdf%<>%group_by(negcount,startgrp)%>%tally


#now plot
plotfile<- here(paste0('plots/Positional_Analysis/','startgrp_vsnegnum','.pdf'))
pdf(plotfile,w=14)
startgrpdf%>%
	filter(startgrp!='start_down')%>%
	group_by(negcount)%>%mutate(n=n/sum(n))%>%
	ggplot(.,aes(x=factor(negcount),fill=startgrp,y=n))+
	stat_identity(geom='bar',position='dodge')+
	scale_fill_discrete(name='Start Change Class')+
	scale_x_discrete(paste0('R/K/H residues in first ',AAnum,' AAs'))+
	scale_y_continuous(paste0('Freq'))+
	ggtitle(paste0('Start class vs basic residues'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))

poseffneuritedf = poseffectdf%>%
	inner_join(neurites%>%select(tr_id=transcript_id,RiboSeq_log2FC_Neurite_Soma))%>%
	mutate(startgrp = case_when(
		  (str_pvalue<0.05) & (str_lfc>0) ~ 'start_up',
		start_down = (str_pvalue<0.05) & (str_lfc<0) ~ 'start_down',
		TRUE ~ 'nochange'
	))

#now plot
plotfile<- here(paste0('plots/Positional_Analysis/','zapp_violin_startgrps','.pdf'))
pdf(plotfile)
poseffneuritedf%>%
	ggplot(.,aes(y=RiboSeq_log2FC_Neurite_Soma,x=startgrp,))+
	geom_violin()+
	# scale_color_discrete(name='colorname',colorvals)+
	scale_x_discrete(paste0('start group'))+
	scale_y_continuous(paste0('Zapp_RPF_Neur_v_Soma'))+
	# ggtitle(paste0('title'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))

startuptrs = poseffectdf%>%filter(str_pvalue<0.05,str_lfc>0)%>%.$tr_id
nonstartuptrs = poseffectdf%>%filter(!(str_pvalue<0.05&str_lfc>0))%>%.$tr_id
#
allaatable = cdsstartaas%>%as.matrix%>%.[,2:10]%>%table
#
aacontentmat = matrix(ncol=2,
	c(cdsstartaas[nonstartuptrs]%>%as.matrix%>%.[,2:4]%>%table%>%.[names(allaatable)],
	cdsstartaas[startuptrs]%>%as.matrix%>%.[,2:4]%>%table%>%.[names(allaatable)]
	))
aacontentmat%<>%set_rownames(names(allaatable))
aacontentmat[-1,]%>%sweep(2,STAT=colSums(.),F='/')%>%as.data.frame%>%mutate(l2ratio=log2(V2/V1))
aacontentmat%>%set_rownames()

