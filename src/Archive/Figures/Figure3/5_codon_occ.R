# # Setup

dir.create('plots/Figures/Figure3/',showWarn=F,rec=T)
#src/R/Figures/Figure2/codon_coverage.R
FLANKCODS<-15
{
if(!exists('here')) base::source(here::here('src/Rprofile.R'))
intersect <- BiocGenerics::intersect
# if(!exists('codonprofiles')) load(here('data/codon_coverage.Rdata'))
if((!exists('allcodsigmean_isomerge'))||(!'availability'%in%colnames(allcodsigmean_isomerge))){
	base::source(here('src/Figures/Figure3/3_tRNA_array_analysis.R'))
} 
library(rlang)
stopifnot('availability' %in% colnames(allcodsigmean_isomerge))
getwd()

stagecolsdisplay <- c(E12.5 = "#214098", E14 = "#2AA9DF", E15.5 = "#F17E22", E17 = "#D14E28", P0 = "#ED3124")
displaystageconv <- names(stagecolsdisplay)%>% setNames(c("E13", "E145", "E16", "E175", "P0"))
stagecols <- stagecolsdisplay %>% setNames(c("E13", "E145", "E16", "E175", "P0"))
stageconv <- names(stagecols) %>% setNames(c("E13", "E145", "E16", "E175", "P0"))
GENETIC_CODE<-Biostrings::GENETIC_CODE

# +
codonprofiles <- readRDS(here('data/codonprofiles.rds'))

codonprofiles$fraction = case_when(
	codonprofiles$sample%>%str_detect('ribo') ~ 'total',
	codonprofiles$sample%>%str_detect('Poly') ~ 'poly',
	TRUE ~ 'NA'
)
#codonprofiles%<>%rename('occupancy':=signal)

stopifnot(!any(codonprofiles$fraction=='NA'))

mainribosamples = read_csv('src/sample_parameter.csv')%>%filter(assay=='ribo',is.na(fraction))%>%.$sample_id
offsets <- read_tsv('ext_data/offsets_manual.tsv')
allTEchangedf<-fread('tables/xtailTEchange.tsv')
# -

# ## prepare occupancy data

# +


offsets %<>% select(offset, compartment, length, readlen)
codonprofiledat = codonprofiles %>% select(sample,codon,readlen,position,occ_nonorm,occupancy,one_of('fraction'))
trna_ab_df_samp = allcodsigmean_isomerge[c("fraction", "time", "sample", "anticodon", "abundance", "codon",
"weightedusage", "availability","rep")]
trna_ab_df_samp$fraction%<>%tolower
#
trna_ab_df_samp%<>%mutate(AA = as.character(translate(DNAStringSet(codon)))) 
offsets%<>%mutate(readlen=paste0('rl',length))
trna_ab_df = trna_ab_df_samp%>%group_by(fraction,time,codon,anticodon)%>%summarise_at(vars(one_of(c("abundance", "weightedusage", "availability"))),mean)
#
codonprofiledat$fraction%>%unique
codonoccs_samp<- codonprofiledat%>%
    filter(sample%in%mainribosamples)%>%
    safe_left_join(offsets%>%select(readlen,offset))%>%
    filter(position== -offset-3)%>%
    separate(sample,c('time','assay','rep'))%>%
    group_by(fraction,time,codon,rep)%>%
    summarise(occupancy = mean (occupancy))

codonoccs_samp%>%group_by(time,fraction,rep)%>%tally

codonoccs <- codonoccs_samp%>%
    group_by(fraction,time,codon)%>%
    summarise(occupancy = mean (occupancy))

codonoccs%>%write_tsv('tables/codonoccs_orig.tsv')

# #occpancy vs AA
codonproftppos<-codonprofiles%>%distinct(readlen,codon)%>%mutate(tppos = 1-as.numeric(str_replace(readlen,'rl','')))

samplestage <- unique(codonprofiles$sample)%>%{setNames(stageconv[str_replace(.,'_.*?_.*?$','')],.)}

}

# -

i=20
ir=100
options(repr.plot.width = i, repr.plot.height = i, repr.plot.res = ir)
	
	pdf('plots/Figures/Figure3/codonprof_rlindiv.pdf',w=24,h=14)
	codonprofiles%>%
		# mutate(time=stageconv[time])%>%
		# filter(codon%>%str_detect(c('GTC|AAC|ATG')))%>%
		mutate(stage=displaystageconv[samplestage[sample]])%>%
		filter(codon%>%str_detect(c('Glu-TTC|GTC|Val-CAC|AAC|ATG')))%>%
		ungroup%>%
		mutate(codon = as_factor(codon))%>%
		# filter(sample%>%str_detect(c('ribo')))%>%
		{print(ggplot(.,aes(position,occupancy,group=sample,
			color=stage))+
			scale_color_manual(values=stagecolsdisplay)+
			scale_x_continuous(minor_breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=3),breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=9) )+
			scale_y_continuous(limits=c(0,6))+
			facet_grid(codon~readlen)+
			geom_line(aes())+
			geom_vline(xintercept=0,linetype=2)+
			geom_vline(data=filter(codonproftppos,codon%in%.$codon),aes(xintercept=tppos),linetype=2)+
			coord_cartesian(xlim=c(-39,12))+
			geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
			geom_vline(data=offsets,aes(xintercept= -offset-5),color=I('green'),linetype=2)+
			# geom_rect(data=offsets,color=I('black'),alpha = I(0.1),aes(x=NULL,y=NULL,xmin= -8-offset, xmax = -offset, ymin = 0 ,ymax = Inf))+
			theme_bw())}
	dev.off()
	normalizePath('plots/Figures/Figure3/codonprof_rlindiv.pdf')

	pdf('plots/Figures/Figure3/offsetwindow_rlmerge_codprofs.pdf',w=14,h=7)
	codonprofiles%>%
		inner_join(offsets)%>%
		group_by(sample,codon,position)%>%summarise(occupancy=sum(occupancy))%>%
			filter(codon%>%str_detect(c('Glu-TTC|GTC|Val-CAC|AAC|ATG')))%>%
		# filter(codon%>%str_detect(c('GTC|AAC|ATG')))%>%
		# filter(codon%>%str_detect(c('Glu-TTC|GTC|Val-CAC|AAC|ATG')))%>%
		# filter(sample%>%str_detect(c('ribo')))%>%
		{print(ggplot(.,aes(position,occupancy,group=sample,
			color=displaystageconv[samplestage[sample]]))+
			# geom_rect(color=I('black'),alpha = I(0.3),aes(xmin= -6, xmax = 3, ymin = 0 ,ymax = Inf))+
			scale_color_manual(values=stagecolsdisplay)+
			scale_x_continuous(name = "Position of RPF 5' end relative to codon first nt",
				minor_breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=3),breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=9) )+
			facet_grid(codon~.)+
			geom_line(aes())+
			# geom_vline(xintercept=0,linetype=2)+
			# geom_vline(data=codonproftppos,aes(xintercept=tppos),linetype=2)+
			# coord_cartesian(xlim=c(-39,12))+
			coord_cartesian(xlim=c(-36,6))+
			# geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
			theme_bw())}
	dev.off()
	normalizePath('plots/Figures/Figure3/offsetwindow_rlmerge_codprofs.pdf')

	pdf('plots/Figures/Figure3/offsetwindow_rlmerge_codprofs.pdf',w=12,h=7)
	codonprofiles%>%
		inner_join(offsets)%>%
		mutate(time=samplestage[sample])%>%
		# mutate(position = position+offset)%>%
		filter(readlen=='rl29')%>%
		mutate(read_length=readlen%>%str_extract('\\d+'))%>%
		group_by(sample,read_length,codon,position)%>%summarise(occupancy=sum(occupancy))%>%
		filter(codon%>%str_detect(c('Glu-TTC|GTC|Val-CAC|AAC|ATG')))%>%
		# filter(occupancy%>%{.<quantile(.,0.1) | .>quantile(.,0.9)})%>%
		# filter(sample%in%c('E13_ribo_1','P0_ribo_2'))%>%
		mutate(stage = displaystageconv[samplestage[sample]])%>%
		filter(!is.na(stage))%>%
		{print(ggplot(.,aes(x=as.numeric(position),y=occupancy,
			color=stage,group=sample))+
			# geom_rect(color=I('black'),alpha = I(0.3),aes(xmin= -6, xmax = 3, ymin = 0 ,ymax = Inf))+
			# scale_color_manual(values=codonspeed[codon])+
			scale_color_manual(values=stagecolsdisplay)+
			scale_y_continuous("Normalized RPF-5'-edge Density ",limits=c(0,5))+
			# scale_color_gradient(name='Dwell Time',low='blue',high='red')+
			scale_x_continuous(name = "Position of RPF 5' end relative to codon first nt",
				minor_breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=3),breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=9) )+
			facet_grid(codon~read_length)+
			geom_line()+
			# geom_vline(xintercept=0,linetype=2)+
			# geom_vline(data=codonproftppos,aes(xintercept=tppos),linetype=2)+
			# coord_cartesian(xlim=c(-21,21))+
			coord_cartesian(xlim=c(-36,6))+
			geom_vline(data=offsets%>%filter(length=='29'),aes(xintercept= -offset-3),color=I('blue'),linetype=2)+
			# geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
			theme_bw())}
	dev.off()
	normalizePath('plots/Figures/Figure3/offsetwindow_rlmerge_codprofs.pdf')

	topbottomcods = codonoccs%>%filter(fraction=='total')%>%select(stage=time,codon,sumocc=occupancy)%>%
	group_by(codon)%>%summarise(sumocc=mean(sumocc))%>%
	mutate(dtrank=rank(sumocc))%>%
	filter(dtrank %in% c(1,2,3,61,60,59))%>%.$codon

	pdf('plots/Figures/Figure3/metametacodon_supp.pdf',w=12,h=7)
	codonprofiles%>%
		inner_join(offsets)%>%
		mutate(stage=samplestage[sample])%>%
		# mutate(position = position+offset)%>%
		filter(readlen=='rl29')%>%
		mutate(read_length=readlen%>%str_extract('\\d+'))%>%
		group_by(sample,stage,read_length,codon,position)%>%summarise(occupancy=sum(occupancy))%>%
		left_join(codonoccs%>%filter(fraction=='total')%>%select(stage=time,codon,sumocc=occupancy),by=c('stage','codon'))%>%
		# .$sumocc%>%unique%>%rank%>%sort
		filter(sample%in%c('E13_ribo_1','P0_ribo_2'))%>%
		filter(codon %in% topbottomcods)%>%
		# .$codon%>%n_distinct
		# group_by(codon,sample)%>%group_slice(1)
		{print(ggplot(.,aes(x=position,y=occupancy,
			color=codon,group=codon))+
			# geom_rect(color=I('black'),alpha = I(0.3),aes(xmin= -6, xmax = 3, ymin = 0 ,ymax = Inf))+
			# scale_color_manual(values=codonspeed[codon])+
			scale_y_continuous("Normalized RPF-5'-edge Density ")+
			# scale_color_gradient(name='Dwell Time',low='blue',high='red')+
			scale_color_discrete(name='codon')+
			scale_x_continuous(name = "Position of RPF 5' end relative to codon first nt",
				minor_breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=3),breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=9) )+
			facet_grid(sample~read_length)+
			geom_line()+
			# geom_vline(xintercept=,linetype=2)+
			# geom_vline(data=codonproftppos,aes(xintercept=tppos),linetype=2)+
			# coord_cartesian(xlim=c(-21,21))+
			coord_cartesian(xlim=c(-36,6))+
			geom_vline(data=offsets%>%filter(length=='29'),aes(xintercept= -offset-3),color=I('blue'),linetype=2)+
			theme_bw())}
	dev.off()
	normalizePath('plots/Figures/Figure3/metametacodon_supp.pdf')


	pdf('plots/Figures/Figure3/metametacodon_supp_grad.pdf',w=12,h=7)
	codonprofiles%>%
		inner_join(offsets)%>%
		mutate(stage=samplestage[sample])%>%
		# mutate(position = position+offset)%>%
		filter(readlen=='rl29')%>%
		mutate(read_length=readlen%>%str_extract('\\d+'))%>%
		group_by(sample,stage,read_length,codon,position)%>%summarise(occupancy=sum(occupancy))%>%
		left_join(codonoccs%>%filter(fraction=='total')%>%select(stage=time,codon,sumocc=occupancy),by=c('stage','codon'))%>%
		# .$sumocc%>%unique%>%rank%>%sort
		filter(sample%in%c('E13_ribo_1','P0_ribo_2'))%>%
		# filter(codon %in% topbottomcods)%>%
		# .$codon%>%n_distinct
		# group_by(codon,sample)%>%group_slice(1)
		{print(ggplot(.,aes(x=position,y=occupancy,
			color=sumocc,group=codon))+
			# geom_rect(color=I('black'),alpha = I(0.3),aes(xmin= -6, xmax = 3, ymin = 0 ,ymax = Inf))+
			# scale_color_manual(values=codonspeed[codon])+
			scale_y_continuous("Normalized RPF-5'-edge Density ")+
			scale_color_gradient(name='Dwell Time',low='blue',high='red')+
			# scale_color_discrete(name='codon')+
			scale_x_continuous(name = "Position of RPF 5' end relative to codon first nt",
				minor_breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=3),breaks = seq(0-(3*FLANKCODS),2+(3*FLANKCODS),by=9) )+
			facet_grid(sample~read_length)+
			geom_line(alpha=I(0.5))+
			# geom_vline(xintercept=,linetype=2)+
			# geom_vline(data=codonproftppos,aes(xintercept=tppos),linetype=2)+
			# coord_cartesian(xlim=c(-21,21))+
			coord_cartesian(xlim=c(-36,6))+
			geom_vline(data=offsets%>%filter(length=='29'),aes(xintercept= -offset-3),color=I('blue'),linetype=2)+
			theme_bw())}
	dev.off()
	normalizePath('plots/Figures/Figure3/metametacodon_supp_grad.pdf')



	#plot with the AAs as colors
	pdf('plots/Figures/Figure3/stripplot_aa_codon.pdf',w=12,h=5)
	codonoccs%>%
		mutate(AA = GENETIC_CODE[codon])%>%
		# filter(fraction=='Total')%>%
		mutate(codonname = paste0(AA,'-',codon))%>%
		group_by(codonname)%>%
		mutate(codmean=mean(occupancy))%>%group_by(AA)%>%mutate(aamean=mean(occupancy))%>%
		ungroup()%>%arrange(aamean,codmean)%>%
		mutate(codonname=as_factor(codonname),AA=as_factor(AA),codonname=as_factor(codonname))%>%
		{
		ggplot(filter(.,(codon%>%table(.)[.]) > 1),aes(y=occupancy,x=codonname,color=time))+
		geom_point(size=1)+
		# geom_line(color=I('grey'))+
		scale_color_manual(values=stagecols)+
		facet_grid(.~AA,scale='free_x')+
		# geom_rect(xmax='R-AGG',xmin='R-CGT',ymin=-Inf,ymax=Inf,alpha=0.2)+
		theme_bw()+
		theme(axis.text.x=element_text(angle=45,vjust=.5))+
		scale_y_continuous(name='Normalized Occupancy Over Codon A site')
	}
	dev.off()
	normalizePath('plots/Figures/Figure3/stripplot_aa_codon.pdf')

# # tRNA vs Occupancy Figures

	################################################################################
	########now both
	################################################################################
	{

	if('occupancy'%in%colnames(codonoccs))codonoccs%<>%dplyr::rename('dwell_time':=occupancy)

	#Now merge with the tRNA data
	tRNA_occ_df<-trna_ab_df%>%
		left_join(codonoccs)
	##Get the amino acid for each one
	tRNA_occ_df%<>%mutate(AA = GENETIC_CODE[codon]%>%qs('S+'))
	#
	tRNA_occ_df%<>%	left_join(codonfreqs%>%colSums%>%enframe('codon','freq'))%>%
		ungroup%>%
		mutate(common = freq>median(freq))
	# tRNA_occ_df%>%ungroup%>%filter(!is.na(abundance))	
	# tRNA_occ_df%<>%filter(fraction=='total')%>%group_by(codon)%>%filter(!all(is.na(abundance)))
    #also add tAI
    tRNA_occ_df%<>%group_by(AA)%>%mutate(tAI = freq/max(freq))%>%ungroup
    tRNA_occ_df%>%filter(time=='P0',fraction=='total')%>%mutate(isna=is.na(abundance))%>%lm(data=.,dwell_time~isna)%>%tidy%>%mutate(p.value<0.05)
    tRNA_occ_df%>%filter(time=='E13',fraction=='total')%>%mutate(isna=is.na(abundance))%>%lm(data=.,dwell_time~isna)%>%tidy%>%mutate(p.value<0.05)
    tRNA_occ_df%>%filter(time=='E13',fraction=='total')%>%ungroup%>%mutate(abundance=replace_na(abundance,min(abundance,na.rm=T)-1))%>%lm(data=.,dwell_time~AA+abundance)%>%tidy%>%
    	filter(term=='abundance')

    tRNA_occ_df%>%filter(time=='E13',fraction=='total')%>%ungroup%>%mutate(availability=replace_na(availability,min(availability,na.rm=T)-1))%>%lm(data=.,dwell_time~AA+availability)%>%tidy%>%
    	filter(term=='availability')
    tRNA_occ_df%>%filter(time=='P0',fraction=='total')%>%ungroup%>%mutate(abundance=replace_na(abundance,min(abundance,na.rm=T)-1))%>%lm(data=.,dwell_time~AA+abundance)%>%tidy%>%
    	filter(term=='abundance')
	tRNA_occ_df%>%filter(time=='P0',fraction=='total')%>%ungroup%>%mutate(abundance=replace_na(abundance,min(abundance,na.rm=T)-1))%>%lm(data=.,dwell_time~AA+abundance)%>%tidy%>%
    	filter(term=='abundance')
    tRNA_occ_df%>%filter(time=='P0',fraction=='total')%>%ungroup%>%mutate(abundance=replace_na(abundance,min(abundance,na.rm=T)-1))%>%lm(data=.,dwell_time~AA+abundance)%>%tidy%>%
    	filter(term=='abundance')

	tRNA_occ_df%<>%filter(is.finite(abundance))
	tRNA_occ_df$codon%>%n_distinct
    tRNA_occ_df_en <- tRNA_occ_df%>%left_join(tRNAenrichdf)

    tRNA_occ_df_en%>%write_tsv('tables/tRNA_stat_df.tsv')

	aatrna_occ_df<-tRNA_occ_df%>%
		group_by(time,AA)%>%
		summarise(dwell_time = mean(dwell_time),abundance_aa = log2(sum(2^(-abundance),na.rm=T)))
    #
 	tRNA_occ_df_en%<>%group_by(time,fraction,AA)%>%mutate(aacor_dwell_time = dwell_time - mean(dwell_time))
	codons4occ <- tRNA_occ_df$codon%>%unique%>%setdiff(c('TAG','TAA','TGA'))

	aa_corrected_occ_av_df <- tRNA_occ_df%>%
		# safe_filter(fraction=='Poly')%>%
		# filter(time=='E13')%>%
		# group_by(AA)%>%filter(n_distinct(codon)>1)%>%
		group_by(fraction,time,AA,codon)%>%summarise(dwell_time = mean(dwell_time),availability=mean(availability))%>%
		# group_by(AA)%>%mutate(availability = availability - mean(availability))%>%
		# group_by(AA)%>%mutate(availability = availability - mean(availability))%>%
		group_by(fraction,time,AA)%>%mutate(dwell_time_aacor = dwell_time - mean(dwell_time))


	# tmp1 <- tRNA_occ_df%>%
	# 	safe_filter(fraction=='poly')%>%
	# 	filter(time=='E13')%>%
	# 	# group_by(AA)%>%filter(n_distinct(codon)>1)%>%
	# 	group_by(AA,codon)%>%summarise(mr = mean(dwell_time),availability=mean(availability))%>%
	# 	# group_by(AA)%>%mutate(availability = availability - mean(availability))%>%
	# 	# group_by(AA)%>%mutate(availability = availability - mean(availability))%>%
	# 	group_by(AA)%>%mutate(mr = mr - mean(mr))

	# pdf('tmp.pdf')
	# aa_corrected_occ_av_df%>%
	# 	ggplot(data=.,aes(x=dwell_time_aacor,y=availability))+
	# 	geom_point(size=1)+
	# 	facet_grid(fraction ~ time)+
	# 	geom_text(show.legend=F,aes(label=codon,color=NULL),size=4)
	# dev.off()
	# normalizePath('tmp.pdf')



	}



	# #total occ, row abundance
	# testvals <- map_df(.id='fraction',unique(tRNA_occ_df_en$fraction)%>%setNames(.,.),function(fractioni){
	# 	map_df(.id='sigcol',sigcols,function(sigcol){
	# 		message(fractioni)
	# 		message(sigcol)
	# 		selvals<-tRNA_occ_df_en%>%
	# 		filter(fraction==fractioni)%>%
	# 			filter(time=='E13')%>%
	# 			select(AA,codon,dwell_time,!!sigcol)%>%
	# 			group_by(AA)%>%
	# 			filter(n_distinct(codon)>1)%>%
	# 			group_by(AA)%>%mutate(aacor_dwell_time = dwell_time - mean(dwell_time))%>%
	# 			select(aacor_dwell_time,!!sigcol)%>%
	# 			# group_by(AA)%>%summarise(dwell_time = mean(dwell_time),abundance=mean(!!sigcol))%>%
	# 			# group_by(AA)%>%mutate(abundance = abundance - mean(abundance))%>%
	# 			identity%>%
	# 			ungroup
	# 		if(all(is.na(selvals[[quo_text(sigcol)]]))) return(NULL)
	# 		selvals%>%{tidy(cor.test(use='complete',.$aacor_dwell_time,.[[quo_text(sigcol)]]))}
	# 	})
	# })
	# testvals%>%arrange(p.value)%>%filter()

	# #okay so availability and dwell time correlate - just
	# tRNA_occ_df_en%>%group_by(AA,codon)%>%filter(fraction=='total',time=='E13')%>%{tidy(cor.test(.$availability,.$dwell_time))}
	# tRNA_occ_df_en%>%group_by(AA,codon)%>%filter(fraction=='total',time=='E175')%>%{tidy(cor.test(.$availability,.$dwell_time))}
	# #abundance not really
	# tRNA_occ_df_en%>%group_by(AA,codon)%>%filter(fraction=='total',time=='E13')%>%{tidy(cor.test(.$abundance,.$dwell_time))}
 

	# sigcols <- c('availability','abundance','weightedusage','abundance_enrich','availability_enrich')
 #    sigcols <- syms(sigcols)%>%setNames(sigcols)
 #    sigcol=sigcols[1]
 #    fractioni='poly'
	tps=names(stagecols)
	makecorlabel = function(x) paste0('rho = ',round(x$estimate,3),'\n','pval = ',ifelse(x$p.value > 0.001,round(x$p.value,4),format(x$p.value,format='e',digits=4)))
	#make the grid of plots for the different codons
	# for(ifraction in c('poly','total')){
	for(ifraction in c('total')){
		for(timeset in list(all=tps)) {
		tRNA_occ_df_en$weightedusage =tRNA_occ_df_en$weightedusage
		sigcols <- c('abundance','weightedusage','availability','aacor_dwell_time','dwell_time')
	    sigcols <- syms(sigcols)%>%setNames(sigcols)
	    sigcol=sigcols[1]
	    sigcolx=sigcols[1]
	    sigcoly=sigcols[4]
	    fractioni='total'
	    # sigcolx <- sig
		# sigcoly <- sym('dwell_time')
		plots<-lapply(sigcols,function(sigcoly){
			lapply(sigcols,function(sigcolx){
				if(sigcolx==sigcoly){	
				    plot <- tRNA_occ_df_en%>%filter(time%in%timeset,fraction==ifraction)%>%
				    	{qplot(data=.,x=!!sigcolx,fill=time,geom='blank')+
				    	# geom_histogram(aes(y=..density..))+
				    	geom_density(aes(y=..density..),alpha=I(0.5))+
				    	scale_fill_manual(values=stagecols)}
				}else{
					str_sigx = as.character(sigcolx)
					str_sigy = as.character(sigcoly)
					dat = tRNA_occ_df_en%>%filter(time%in%timeset,fraction==ifraction)
					if(str_sigx=='aacor_dwell_time') dat = dat[dat[[str_sigx]]!=0,]
					if(str_sigy=='aacor_dwell_time') dat = dat[dat[[str_sigy]]!=0,]
					mdat = dat%>%group_by(codon)%>%summarise_at(vars(str_sigx,str_sigy),mean)
					corlabel = cor.test(mdat[[str_sigx]],mdat[[str_sigy]])%>%tidy%>%makecorlabel
				    plot <- dat%>%
				    	{qplot(data=.,x=!!sigcolx,color=time,y=!!sigcoly,geom='point',position='')+scale_color_manual(values=stagecols)+
				    	ggtitle(corlabel)} 
				}
				plot <- plot+theme_bw()
			})
		})
		plots <- ggarrange(plotlist=plots%>%purrr::flatten(.),ncol=length(sigcols),nrow=length(sigcols))
		plots <- plots+theme_bw()
		plotfile = paste0('plots/Figures/Figure3/codon_stat_grid_',ifraction,'_',paste0(timeset,collapse='-'),'.pdf')
	    pdf(plotfile,h=20,w=20)
	    print(plots)
	    dev.off()
		normalizePath(plotfile)%>%message
		}
	}


	# dat%>%ungroup%>%filter(time=='E13')%>%mutate(highdt=dwell_time>median(dwell_time,na.rm=T))%>%{split(.$abundance,.$highdt)}%>%{t.test(.[[1]],.[[2]])}
	# dat%>%ungroup%>%filter(time=='E13')%>%mutate(highdt=dwell_time>median(dwell_time,na.rm=T))%>%{split(.$availability,.$highdt)}%>%{t.test(.[[1]],.[[2]])}
	# dat%>%ungroup%>%filter(time=='E13')%>%mutate(highdt=dwell_time>median(dwell_time,na.rm=T))%>%{split(.$dwell_time,.$highdt)}%>%
	# 	ggplot(aes(x=abundance))

# # Linear modeling of tRNA - occupancy

	################################################################################
	########explorating quantative relationship - avail, AA, dwell time
	################################################################################
		


	tRNA_occ_df_mean = tRNA_occ_df_en%>%filter(fraction=='total')%>%group_by(AA,codon)%>%summarise_at(vars('availability','abundance','aacor_dwell_time','dwell_time','tAI'),mean)
	tRNA_occ_df_en%<>%group_by(time)%>%mutate(n_availability = availability - mean(availability))
	tRNA_occ_df_en%<>%group_by(time)%>%mutate(n_abundance = abundance - mean(abundance))
	tRNA_occ_df_en%<>%ungroup
	#So AA definitely explains a lot of the variaiton in dwell time - and this effect varies.
	tRNA_occ_df_mean%>%lm(data=.,dwell_time ~ AA)%>%anova
	tRNA_occ_df%>%lm(data=.,dwell_time ~ AA+AA:time)%>%anova
	#so averaging over time points, it's not really teh case that higher dwell time equates to higher
	tRNA_occ_df_mean%>%lm(data=.,dwell_time ~ abundance)%>%anova
	tRNA_occ_df_mean%>%lm(data=.,dwell_time ~ availability)%>%anova
	tRNA_occ_df_mean%>%lm(data=.,dwell_time ~ tAI)%>%anova

	cor.test(tRNA_occ_df_mean$dwell_time,tRNA_occ_df_mean$abundance)
	cor.test(tRNA_occ_df_mean$aacor_dwell_time,tRNA_occ_df_mean$tAI)
	cor.test(tRNA_occ_df_en$aacor_dwell_time,tRNA_occ_df_en$tAI)
	tRNA_occ_df_en%>%filter(fraction=='total')%>%{cor.test(.$aacor_dwell_time,.$abundance)}
	cor.test(tRNA_occ_df_mean$aacor_dwell_time,tRNA_occ_df_mean$abundance)
	cor.test(tRNA_occ_df_mean$dwell_time,tRNA_occ_df_mean$tAI)


	# tRNA_occ_df%>%lm(data=.,abundance ~ aacor_dwell_time)%>%anova

	dblcor_tbl <- tRNA_occ_df_en%>%filter(fraction=='total')%>%
		group_by(AA,time)%>%
		mutate(aacor_dwell_time = dwell_time - mean(dwell_time))%>%
		mutate(aacor_abundance = abundance - mean(abundance))%>%
		mutate(aacor_availability = availability - mean(availability))%>%
		group_by(time)%>%mutate(tcor_abundance = abundance - mean(abundance))
		
		# {cor.test(.$aacor_dwell_time,.$aacor_abundance)}

	dblcortest = dblcor_tbl%>%	{cor.test(.$aacor_dwell_time,.$aacor_availability)}

	dblcor_tbl%>%	{quicktest(.$aacor_dwell_time,.$aacor_availability)}
	dblcor_tbl%>%	{quicktest(.$aacor_dwell_time,.$aacor_abundance)}
	dblcor_tbl%>%	{quicktest(.$aacor_dwell_time,.$tcor_abundance)}

	#availibilty
	tRNA_occ_df_mean%>%lm(data=.,dwell_time ~ abundance)%>%anova
	#nor is that the case within individual timepoints
	tRNA_occ_df_en%>%filter(time=='E175',fraction=='total')%>%lm(data=.,dwell_time ~ availability)%>%anova
	tRNA_occ_df_en%>%filter(time=='P0',fraction=='total')%>%lm(data=.,dwell_time ~ abundance)%>%anova
	#and a simple correction for the AA effect doesn't help
	tRNA_occ_df_en%>%filter(time=='E175',fraction=='total')%>%lm(data=.,dwell_time ~ AA+availability)%>%anova
	tRNA_occ_df_en%>%filter(time=='E175',fraction=='total')%>%lm(data=.,dwell_time ~ AA+abundance)%>%anova
	#however if we lok at all timepoints and do a tp sepcific AA correction, now there's a link
	tRNA_occ_df_en%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ AA+availability)%>%anova
	tRNA_occ_df_en%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ AA:time+(availability))%>%anova
	tRNA_occ_df_en%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ AA:time+abundance)%>%anova
	tRNA_occ_df_en%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ AA:time+n_availability)%>%anova
	tRNA_occ_df_en%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ AA*time+n_abundance)%>%anova

	#but now we have codon identity confounding our points - does this explain the effect? Yes pretty much
	tRNA_occ_df_en%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ AA:time+codon*abundance)%>%anova

	
	#True statements:
	#There is a strong variation in mean availability among amino acides (about 1/3 of the variance in availbility)
	tRNA_occ_df_en%>%filter(time=='E175',fraction=='total')%>%lm(data=.,dwell_time ~ AA+availability)%>%anova
	#No real tendency availibilty have an average higher dwell time, once we correct for amino acid identity
	tRNA_occ_df_en%>%filter(time=='E175',fraction=='total')%>%lm(data=.,dwell_time ~ AA+availability)%>%confint
	#HOWEVER there does seem to be some effect of change in tRNA levels over time.
	tRNA_occ_df_en%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ AA+codon+codon:availability)%>%anova




	#but only a few:
	tRNA_occ_df_en%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ codon+codon:availability)%>%.$coefficients%>%sort
	tRNA_occ_df_en%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ codon+codon:availability)%>%confint
	relcods = tRNA_occ_df_en%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ codon+codon:availability)%>%confint%>%as.data.frame%>%
		rownames_to_column('term')%>%tail(-1)%>%filter(`2.5 %`>0)%>%.$term%>%str_extract('(?<=codon)...')
	# stopifnot(relcods == c('CTG','TAC','TCG'))
	#these guys are among the slowest codons	
	tRNA_occ_df%>%filter(time=='E13')%>%mutate(relcod = codon%in%relcods)%>%arrange(availability)%>%as.data.frame
	#or at E13 even more so
	tRNA_occ_df%>%filter(fraction=='total',time=='E13')%>%mutate(relcod = codon%in%relcods)%>%arrange(-dwell_time)%>%head(10)%>%as.data.frame
	#about the same in poly or 80S
	tRNA_occ_df%>%filter(fraction=='80S',time=='E13')%>%mutate(relcod = codon%in%relcods)%>%arrange(-dwell_time)%>%head(10)%>%as.data.frame

	tRNA_occ_df%>%filter(fraction=='Poly',time=='E13')%>%mutate(relcod = codon%in%relcods)%>%arrange(-availability)%>%as.data.frame
	tRNA_occ_df%>%filter(fraction=='80S',time=='E13')%>%mutate(relcod = codon%in%relcods)%>%arrange(-availability)%>%as.data.frame
	tRNA_occ_df%>%filter(fraction=='Poly',time=='P0')%>%mutate(relcod = codon%in%relcods)%>%arrange(-availability)%>%as.data.frame
	tRNA_occ_df%>%filter(fraction=='80S',time=='P0')%>%mutate(relcod = codon%in%relcods)%>%arrange(-availability)%>%as.data.frame

	tRNA_occ_df_en%>%filter(fraction=='total',)%>%lm(data=.,dwell_time ~ codon+codon:availability)%>%.$coefficients%>%sort
	

	tRNA_occ_df_en%>%filter(fraction=='total')%>%
		filter(codon%in%c('CTG','TAC','TCG'))%>%
		group_by(codon)%>%
		filter(n()==5)%>%
		group_slice(2)%>%
		mutate(codrelav=availability/mean(availability),codreldt=dwell_time/mean(dwell_time))%>%
		{quicktest(.$dwell_time,.$availability)}


	#overall relationshiop
	tRNA_occ_df_en%>%filter(fraction=='total')%>%
		#filter(codon%in%c('CTG','TAC','TCG'))%>%
		group_by(codon)%>%
		filter(n()==5)%>%
		#group_slice(1)%>%
		mutate(codrelav=availability/mean(availability),codreldt=dwell_time/mean(dwell_time))%>%
		{quicktest(.$dwell_time,.$availability)}
	#it's not just that these guys have a very high var in availibiity
	tRNA_occ_df_en%>%group_by(codon)%>%summarise(avsd = sd(availability))%>%arrange(-avsd)
	#what if I were to test cors and then multiple correct ? No way, too many corrections - only barely sig  as is
	codcor_df = tRNA_occ_df_en%>%filter(fraction=='total',codon%in%codons4occ)%>%group_by(codon)%>%filter(n()>3)%>%{split(.,.$codon)}%>%map_df(.id='codon',.%>%{tidy(cor.test(.$dwell_time,.$availability,method='pearson',alt='greater'))})%>%
		arrange(desc(codon %in% relcods))%>%
		mutate(p.value=p.adjust(p.value))

	# codcor_df%>%safe_left_join(tRNA_occ_df%>%filter(time=='E13',fraction=='total')%>%select(codon,dwell_time))%>%
		# {quicktest(.$estimate,.$dwell_time)}

	# codonoccs%>%left_join(trna_ab_df%>%mutate(anticodon=anticodon%>%str_replace('...\\-','')),by=c('time','codon'='codon'))%>%filter(fraction=='total')%>%{quicktest(.$dwell_time,.$availability)}
	tRNA_occ_df_en%>%group_by(time,AA)%>%
		filter(n()>1)%>%mutate(ishighdt = dwell_time==max(dwell_time))%>%
		mutate(ishighav = availability==max(availability))%>%
		{table(.$ishighav,.$ishighdt)}%>%
		fisher.test

	
	tRNA_occ_df_en%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ time*AA*availability)%>%anova
	tRNA_occ_df_en%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ AA+AA:availability)%>%anova
	
	################################################################################
	########Now do the above with samples as data points rather than averages
	################################################################################
	{
	trna_ab_df_samp$rep%<>%str_replace('rep','')
	trna_occ_df_samp <- trna_ab_df_samp%>%left_join(codonoccs_samp,by=c('fraction','codon','time','rep'))
	trna_occ_df_samp%<>%mutate(AA = GENETIC_CODE[codon]%>%qs('S+'))
	trna_occ_df_samp%<>%rename('dwell_time':='occupancy')
	trna_occ_df_samp_mean = trna_occ_df_samp%>%filter(fraction=='total')%>%group_by(AA,codon)%>%summarise_at(vars('abundance','availability','dwell_time'),mean)
	
	codonoccs_samp%>%group_by(codon,time,fraction)%>%group_slice(1)
	trna_ab_df_samp%>%group_by(codon,time,fraction)%>%tally
	trna_occ_df_samp%>%group_by(codon,time,fraction)%>%tally
	#So AA definitely explains a lot of the variaiton in dwell time - and this effect varies.
	trna_occ_df_samp_mean%>%lm(data=.,dwell_time ~ AA)%>%anova
	tRNA_occ_df%>%lm(data=.,dwell_time ~ AA+AA:time)%>%anova
	#so averaging over time points, it's not really teh case that higher dwell time equates to higher
	trna_occ_df_samp_mean%>%lm(data=.,dwell_time ~ abundance)%>%anova
	#availibilty
	# tRNA_occ_df_mean%>%lm(data=.,dwell_time ~ )%>%anova
	#But within individual timeopints, it is (though it's a weak realtionships)
	trna_occ_df_samp%>%filter(time=='E175',fraction=='total')%>%lm(data=.,dwell_time ~ abundance)%>%anova
	#THis remains, even after the AA effect.
	trna_occ_df_samp%>%filter(time=='E175',fraction=='total')%>%lm(data=.,dwell_time ~ AA+abundance)%>%anova
	#if we lok at all timepoints and do a tp sepcific AA correction, now there's a relatively  significant lnk
	trna_occ_df_samp%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ AA:time+abundance)%>%anova
	#but now we have codon identity confounding our points - does this explain the effect? Yes pretty much
	trna_occ_df_samp%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ AA*time+codon*abundance)%>%anova

	#True statements:
	#There is a strong variation in mean abundance among amino acides (about 1/3 of the variance in availbility)
	trna_occ_df_samp%>%filter(time=='E175',fraction=='total')%>%lm(data=.,dwell_time ~ AA+abundance)%>%anova
	#No real tendency availibilty have an average higher dwell time, once we correct for amino acid identity
	trna_occ_df_samp%>%filter(time=='E175',fraction=='total')%>%lm(data=.,dwell_time ~ AA+abundance)%>%confint
	#HOWEVER there does seem to be some effect of change in tRNA levels over time.
	trna_occ_df_samp%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ codon*abundance)%>%anova
	trna_occ_df_samp%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ codon*abundance)%>%anova
	# trna_occ_df_samp%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ AA:time+codon*abundance)%>%anova
	#but only a few:
	trna_occ_df_samp%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ codon+codon:abundance)%>%.$coefficients%>%sort
	trna_occ_df_samp%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ codon+codon:abundance)%>%confint
	relcods = trna_occ_df_samp%>%filter(fraction=='total')%>%lm(data=.,dwell_time ~ codon+codon:abundance)%>%confint%>%as.data.frame%>%
		rownames_to_column('term')%>%tail(-1)%>%filter(`2.5 %`>0)%>%.$term%>%str_extract('(?<=codon)...(?=:abund)')%>%na.omit
	stopifnot(c('CTG','TAC','TCG')%in%relcods)
	#Not really a relationship to abundanc, 
	tRNA_occ_df%>%filter(time=='E13')%>%mutate(relcod = codon%in%relcods)%>%arrange(abundance)%>%as.data.frame
	#But these guys are among the slowest codons or at E13 even more so
	tRNA_occ_df%>%filter(fraction=='total',time=='E13')%>%mutate(relcod = codon%in%relcods)%>%arrange(-dwell_time)%>%as.data.frame%>%
		select(relcod,AA,codon,dwell_time)%>%mutate(slownessrank=rank(dwell_time))%>%filter(relcod)
	tRNA_occ_df%>%filter(fraction=='total',time=='E13')%>%mutate(relcod = codon%in%relcods)%>%{split(.$abundance,.$relcod)}%>%{t.test(.[['TRUE']],.[['FALSE']])}%>%tidy
	tRNA_occ_df%>%filter(fraction=='total',time=='E175')%>%mutate(relcod = codon%in%relcods)%>%{split(.$abundance,.$relcod)}%>%{t.test(.[['TRUE']],.[['FALSE']])}%>%tidy
	tRNA_occ_df%>%filter(fraction=='total',time=='E13')%>%mutate(relcod = codon%in%relcods)%>%{split(.$dwell_time,.$relcod)}%>%{t.test(.[['TRUE']],.[['FALSE']])}%>%tidy
	tRNA_occ_df%>%filter(fraction=='total',time=='E175')%>%mutate(relcod = codon%in%relcods)%>%{split(.$dwell_time,.$relcod)}%>%{t.test(.[['TRUE']],.[['FALSE']])}%>%tidy

	# tRNA_occ_df%>%filter(fraction=='poly',time=='E13')%>%mutate(relcod = codon%in%relcods)%>%{split(.$dwell_time,.$relcod)}%>%{t.test(.[['TRUE']],.[['FALSE']])}%>%tidy

	#about the same in poly or 80S
	tRNA_occ_df%>%filter(fraction=='80s',time=='E13')%>%mutate(relcod = codon%in%relcods)%>%arrange(-dwell_time)%>%as.data.frame

	tRNA_occ_df%>%filter(fraction=='poly',time=='E13')%>%mutate(relcod = codon%in%relcods)%>%arrange(-abundance)%>%as.data.frame
	tRNA_occ_df%>%filter(fraction=='80s',time=='E13')%>%mutate(relcod = codon%in%relcods)%>%arrange(-abundance)%>%as.data.frame
	tRNA_occ_df%>%filter(fraction=='poly',time=='P0')%>%mutate(relcod = codon%in%relcods)%>%arrange(-abundance)%>%as.data.frame
	tRNA_occ_df%>%filter(fraction=='80s',time=='P0')%>%mutate(relcod = codon%in%relcods)%>%arrange(-abundance)%>%as.data.frame



	tRNA_occ_df%>%.$dwell_time%>%na.omit%>%txtdensity

	trna_occ_df_samp%>%filter(fraction=='total',)%>%lm(data=.,dwell_time ~ codon+codon:abundance)%>%.$coefficients%>%sort
	
	tRNA_occ_df%>%filter(fraction=='total')%>%group_by(time)%>%summarise(mean(abundance))
	
	slowclassdf = tRNA_occ_df%>%
		filter(fraction=='total')%>%
		group_by(codon)%>%summarise(dwell_time=mean(dwell_time))%>%
		# filter(fraction=='total',time=='E13')%>%
		# arrange(desc(dwell_time))%>%
		# group_by(dwell_time = dwell_time / mean(dwell_time))
		mutate(slowclass = case_when(
			dwell_time > quantile(dwell_time,.75,na.rm=TRUE) ~ 'Top 25% DT at E13 (slowest)',
			dwell_time < quantile(dwell_time,.25,na.rm=TRUE) ~ 'Bottom 25% DT at E13 (fastest)',
			TRUE ~ 'Middle 50%'
		))%>%select(codon,slowclass)%>%mutate(slowclass=factor(slowclass,rev(unique(slowclass))))
	codreldf = trna_occ_df_samp%>%filter(fraction=='total')%>%
		left_join(slowclassdf)%>%
		# mutate(codrelab=abundance-mean(abundance),codreldt=dwell_time/mean(dwell_time))
		group_by(time)%>%mutate(abundance = abundance - mean(na.rm=T,abundance))%>%
		# group_by(time)%>%mutate(abundance = availability - mean(na.rm=T,availability))%>%
		group_by(codon)%>%
		# mutate(codrelab=abundance-mean(abundance))%>%	
		mutate(codrelab=abundance-mean(abundance))%>%	
		mutate(codreldt=dwell_time/mean(dwell_time))
    codreltests = codreldf%>%group_by(slowclass)%>%summarise(tidy(cor.test(codreldt,codrelab)))
	codreltests = codreltests%>%mutate(acor_label=paste0('rho = ',round(estimate,3),'\n','pval = ',ifelse(p.value > 0.001,round(p.value,4),format(p.value,format='e',digits=4))))
	#now plot
	plotfile<- here(paste0('plots/','cod_rel_df','.pdf'))
	pdf(plotfile,width=5,h=9)
		codrelplot = ggplot(data=codreldf , aes(x = codrelab,y = codreldt,color=time))+geom_point()+
		facet_grid(slowclass~.,scales='free_x')+
		geom_smooth(aes(color=NA),method='lm')+
		geom_text(data=codreltests,color=I('black'),aes(label=acor_label,x=1.5,y=1.5))+
		theme_bw()
	#
	codrelplot+
		scale_color_manual(name='Time',values=stagecols)+
		scale_x_continuous(paste0('log2(Abundance / Mean Abundance [Codon])'),limits=c(-3,2))+
		scale_x_continuous(paste0('log2(Abundance / Mean Abundance [Codon])\n time corrected'),limits=c(-4,3))+
		scale_y_continuous(paste0('Dwell Time / Mean DT [Codon] (relative slowness)'))+
		ggtitle(paste0('Relationship Between tRNA abundance and Dwell time for fast/slow codons'))+
		theme_bw()
	dev.off()
	normalizePath(plotfile)


	trna_occ_df_samp%>%filter(fraction=='total')


		# {quicktest(.$dwell_time,.$abundance)}

	}

# ## Gingold Class analysis

	totstarttrandf = tRNA_occ_df%>%filter(fraction=='total',time=='E13')
	cdsseq%>%oligonucleotideFrequency(step=6,width=6)%>%colSums%>%enframe('hexan','count')%>%mutate(codon=str_extract(hexan,'...$'),p_codon=hexan%>%str_extract('^...'))%>%select(-hexan)%>%
		# filter(codon%>%is_in(relcods))%>%
		left_join(totstarttrandf,by=c('p_codon'='codon'))%>%
		group_by(codon)%>%
		summarise(pcodslowness= logSumExp(log2(count/sum(count))+dwell_time,na.rm=T))%>%
		mutate(isrelcod=codon%in%relcods)%>%arrange(pcodslowness)%>%as.data.frame
	#change in availability over time is associated with changes in dwell time (weakly though)


	# Sys.glob('pipeline/star/reports/*/*idx*')%>%str_subset('total')%>%str_subset('.bam.',neg=T)%>%setNames(.,.)%>%map_dbl(.%>%fread%>%filter(V1%>%str_detect('chr'))%>%group_by(V1=='chrM')%>%summarise(n=sum(V3))%>%{.$n[2]/sum(.$n)})%>%
	# 	txtplot

	# Sys.glob('pipeline/star/reports/*/*idx*')%>%str_subset('ribo')%>%str_subset('.bam.',neg=T)%>%setNames(.,.)%>%map_dbl(.%>%fread%>%filter(V1%>%str_detect('chr'))%>%group_by(V1=='chrM')%>%summarise(n=sum(V3))%>%{.$n[2]/sum(.$n)})%>%
	# 	txtplot



	gingold_types = 'ext_data/gingold_etal_2014_trna_types.tsv'%>%read_tsv
	gingold_types=gingold_types%>%mutate(codon = DNAStringSet(anticodon)%>%reverseComplement%>%as.character)


	if(!'type'%in%colnames(tRNA_occ_df_en))tRNA_occ_df_en%<>%left_join(gingold_types%>%select(codon,type),by='codon')

	#now plot
	plotfile<- here(paste0('plots/','gingold_type_dist','.pdf'))
	pdf(plotfile)
		tRNA_occ_df_en%>%
		filter(fraction=='total',!is.na(type))%>%
		ggplot(.,aes(fill=type,x=dwell_time))+
		geom_density(alpha=I(0.5))+
		facet_wrap(time~.)+
	    # scale_fill_manual(values = stagecols) +
		scale_x_continuous(paste0('Dwell Time'))+
		ggtitle(paste0('DT vs gingold types'))+
		theme_bw()
	dev.off()
	message(normalizePath(plotfile))

	lapply(unique(tRNA_occ_df_en$time),function(itp){
		tRNA_occ_df_en%>%
			filter(fraction=='total',!is.na(type))%>%
			filter(time==itp,type%in%c('prol','diff'))%>%
			{split(.$dwell_time,.$type)}%>%{t.test(.[[1]],.[[2]])}
	})%>%map_df(.id='time',tidy)%>%write_tsv('tables/gingold_ttest_vals.tsv')
	message(normalizePath('tables/gingold_ttest_vals.tsv'))

profvarpca <- codonprofiles%>%
	select(sample,readlen,position,occ_nonorm,occupancy,codon)%>%
	split(.,.$sample)%>%
	map_df(.id='sample',.%>%
		split(.,list(.$readlen))%>%
		map_df( .id='readlen',.%>%
			mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
			filter(position> -numreadlen+6,position < -6)%>%
			# filter(sample=='E13_ribo_1')%>%
			ungroup%>%
			select(-numreadlen,-occ_nonorm,-readlen,-sample)%>%
			group_by(codon)%>%
			spread(position,occupancy)%>%
			{set_rownames(.[,-1],.$codon)}%>%
			princomp%>%{.$loadings[,1]}%>%{./.[which.max(abs(.))]}%>%enframe('position','pca1')
		)
	)

profvarpca%<>%select(sample,readlen,position,pca1)

profvarpca%>%group_by(readlen,position)%>%summarise(pca1=mean(pca1))%>%
	mutate(pca13wind = pca1+lead(pca1)+lead(pca1,2))

library(rlang)#
plotfile<-'plots/Figures/Figure3/fppos_vs_codon_pcascore.pdf'
offsets%<>%mutate(readlen=paste0('rl',length))
pdf(plotfile,w=12,h=12)
profvarpca%>%slice_by(sample,c(1,2,3,4,5,6))%>%
ggplot(data=.,aes(y=pca1,x=as.numeric(position)))+geom_point()+
	facet_grid(readlen~sample)+
		geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=offsets,aes(xintercept= -offset-5),color=I('green'),linetype=2)
dev.off()
normalizePath(plotfile)


library(rlang)#
plotfile<-'plots/Figures/Figure3/fppos_vs_codon_pcascore_3wind.pdf'
offsets%<>%mutate(readlen=paste0('rl',length))
pdf(plotfile,w=12,h=12)
profvarpca%>%
	mutate(pca13wind = pca1+lead(pca1)+lead(pca1,2))%>%
	slice_by(sample,c(1,2,3,4,5,6))%>%
	ggplot(data=.,aes(y=pca13wind,x=as.numeric(position)))+geom_point()+
	facet_grid(readlen~sample)+
		geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=offsets,aes(xintercept= -offset-5),color=I('green'),linetype=2)
dev.off()
normalizePath(plotfile)

pcaderrivedpos = profvarpca%>%
	mutate(pca13wind = pca1+lead(pca1)+lead(pca1,2))%>%
	group_by(readlen,position)%>%summarise(pca13wind=mean(pca13wind))%>%
	group_by(readlen)%>%slice(which.max(pca13wind))
pcaderrivedpos=pcaderrivedpos%>%mutate(position = map(as.numeric(position),~ (.:(.+2))))%>%unnest

codonprofiles_pcawind = codonprofiles%>%semi_join(pcaderrivedpos%>%mutate(position=as.numeric(position)))
codonprofiles_pcawind = codonprofiles_pcawind%>%group_by(sample,codon,fraction)%>%
	summarise(occ_nonorm=sum(occ_nonorm),occupancy=sum(occupancy))%>%
	separate(sample,c('time','assay','rep'))%>%
	group_by(codon,time)%>%
	summarise(occ_nonorm=sum(occ_nonorm),occupancy=sum(occupancy))


quicktest(codonprofiles_pcawind$occ_nonorm,codonprofiles_pcawind$occupancy)

offsets%<>%mutate(readlen=paste0('rl',length))
pdf('plots/Figures/Figure3/fppos_vs_codon_variance.pdf',w=12,h=12)
#plotting variance amongst codons at each point.
codonprofiles%>%
	ungroup%>%
	group_by(sample,readlen,position)%>%
	# filter(signal==0)%>%.$codon%>%unique
	# filter(signal!=0)%>%
	# filter(position==1)%>%
	filter(!is.nan(signal))%>%
	mutate(signal=occ_nonorm)%>%
	summarise(sdsig=sd(signal,na.rm=T)/mean(signal,na.rm=T))%>%
	separate(sample,c('time','assay','rep'))%>%
	group_by(readlen,time,assay,position)%>%
	summarise(sdsig=mean(sdsig))%>%
	mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
	filter(position> -numreadlen+6,position < -6)%>%
	filter(numreadlen>=25,numreadlen<=31)%>%
	arrange(position)%>%
	{
		qplot(data=.,x=position,y=sdsig)+
		theme_bw()+
		facet_grid(readlen~time)+
		scale_y_continuous('between codon variation (meannorm)')+
		scale_x_continuous('5 read position relative to codon ')+
		geom_vline(data=offsets,aes(xintercept= -offset),color=I('blue'),linetype=2)+
		geom_vline(data=offsets,aes(xintercept= -offset-5),color=I('green'),linetype=2)+
		ggtitle("variance of 5' read occurance vs position")
	}%>%print
dev.off()
normalizePath('plots/Figures/Figure3/fppos_vs_codon_variance.pdf')



rfreqdf = tRNA_occ_df%>%distinct(AA,codon,freq)%>%group_by(AA)%>%mutate(rfreq=freq/sum(freq))
rfreqdtdf = tRNA_occ_df%>%filter(fraction=='total')%>%left_join(rfreqdf)%>%select(time,codon,AA,rfreq,dwell_time)%>%
		group_by(codon)%>%summarise(rfreq=rfreq[1],dwell_time=mean(dwell_time))
codreltests = rfreqdtdf%>%summarise(tidy(cor.test(rfreq,dwell_time)))
codreltests = codreltests%>%mutate(acor_label=paste0('rho = ',round(estimate,3),'\n','pval = ',ifelse(p.value > 0.001,round(p.value,4),format(p.value,format='e',digits=4))))
	
#now plot
plotfile<- here(paste0('plots/Figures/Figure3/','rfreq_vs_dwell_time','.pdf'))
pdf(plotfile,w=5,h=5)
rfreqdtdf%>%
	ggplot(.,aes(x=rfreq,y=dwell_time))+
	geom_point()+
	scale_x_continuous(paste0('Optimality'))+
	scale_y_continuous(paste0('Dwell Time'))+
	ggtitle(paste0('Codon optimality vs dwell time'))+
	scale_color_manual(values=stagecolsdisplay)+
	geom_text(data=codreltests,aes(label=acor_label,x=Inf,y=Inf),color=I('black'),vjust=1,hjust=1)+
	# facet_grid(.~time)+
	theme_bw()
dev.off()
message(normalizePath(plotfile))

# we can show that the dominant effect on dwell time is the amino acid coded for, in line with Riba et al. We can also show that those AA effects change significantly over time. We can show that within each amino acid coded for, there are codon-specific differences in dwell time which are consistent between time points. We can also show that there exists a relatively weak relationship between tRNA abundance and dwell time, and that if we control for the above, this becomes  controlling for these (and only if we control for them), a relationship between tRNA abundance and dwell time, with slower codons apparently protecting 
################################################################################
########GC content of codons
################################################################################
tRNA_occ_df%>%filter(fraction=='total',time=='E13')
codontable = 'tables/S3.xlsx'%>%readxl::read_xlsx(2)
codontable%>%.$codon%>%n_distinct

# codontable=tRNA_occ_df_en
codontable%<>%mutate(gcnum = str_count(codon,'G|C'))
codontable%<>%mutate(gc1 = str_count(codon,'(G|C)..'))
codontable%<>%mutate(gc2 = str_count(codon,'.(G|C).'))
codontable%<>%mutate(gc3 = str_count(codon,'..(G|C)'))

#now plot
codontable%>%filter(time==first(time))%>%lm(data=.,dwell_time ~ AA+gc1+gc2+gc3)%>%tidy%>%filter(term%>%str_detect('gc'))
codontable%>%group_by(time,AA,codon,gc1,gc2,gc3)%>%summarise(dwell_time=mean(dwell_time))%>%lm(data=.,dwell_time ~ AA+gc1+gc2+gc3)%>%tidy%>%filter(term%>%str_detect('gc'))
codontable%>%group_by(time,AA,codon,gc1,gc2,gc3)%>%summarise(dwell_time=mean(dwell_time))%>%lm(data=.,dwell_time ~ AA+gc1+gc2+gc3)%>%tidy%>%filter(term%>%str_detect('gc'))
codontable%>%group_by(time,AA,codon,gc1,gc2,gc3)%>%summarise(dwell_time=mean(dwell_time))%>%lm(data=.,dwell_time ~ AA+gc1+gc2+gc3)%>%tidy%>%filter(term%>%str_detect('gc'))


codontable%>%group_by(time,AA,codon,gcnum)%>%summarise(dwell_time=mean(dwell_time))%>%lm(data=.,dwell_time ~ gcnum)%>%tidy%>%filter(term%>%str_detect('gc'))
codontable%>%group_by(time,AA,codon,gcnum)%>%summarise(dwell_time=mean(dwell_time))%>%lm(data=.,dwell_time ~ AA+gcnum)%>%tidy%>%filter(term%>%str_detect('gc'))
codontable%>%group_by(time,AA,codon,gcnum,gc1,gc2,gc3)%>%summarise(dwell_time=mean(dwell_time))%>%lm(data=.,dwell_time ~ gc1+gc2+gc3)%>%tidy%>%filter(term%>%str_detect('gc'))
codontable%>%group_by(time,AA,codon,gcnum,gc1,gc2,gc3)%>%summarise(dwell_time=mean(dwell_time))%>%lm(data=.,dwell_time ~ AA+gc1+gc2+gc3)%>%tidy%>%filter(term%>%str_detect('gc'))


codontable%>%group_by(time)%>%nest()%>%mutate(lmres=map(data,~lm(data=.,dwell_time ~AA+ gcnum)%>%tidy%>%filter(term%>%str_detect('gc'))))%>%unnest(lmres)%>%
	select(-data)
codontable%>%group_by(time)%>%nest()%>%mutate(lmres=map(data,~lm(data=.,dwell_time ~ AA+gc1+gc2+gc3)%>%tidy%>%filter(term%>%str_detect('gc'))))%>%unnest(lmres)%>%
	select(-data)


plotfile<- here(paste0('plots/','codon_gc_dt_boxplots','.pdf'))
pdf(plotfile,w=16,h=8)
codontable%>%
	# filter(time=='E12.5')%>%
	# group_by(AA,codon,gcnum,gc1,gc2,gc3)%>%summarise(dwell_time=mean(dwell_time,na.rm=T))%>%
	ggplot(.,aes(x=as.factor(gcnum),y=dwell_time))+
	geom_boxplot()+
	geom_jitter()+
	# scale_color_discrete(name='colorname',colorvals)+
	facet_grid(.~time)+
	scale_x_discrete(paste0('G/C nucleotides'))+
	scale_y_continuous(paste0('average(dwell_time)'))+
	ggtitle(paste0('gc_content'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))

save.image('data/fig2_codon_occ_2.Rdata')
load('data/fig2_codon_occ_2.Rdata')