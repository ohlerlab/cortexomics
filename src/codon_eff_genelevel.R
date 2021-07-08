

library(rtracklayer)
library(Biostrings)
library(GenomicFeatures)
library(GenomicAlignments)
require(txtplot)
require(Rsamtools)
require(rlang)
#coverageplots

library(data.table)

base::source(here::here('src/R/Rprofile.R'))
if(!exists('cdsgrl')) base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
if(!exists('iso_tx_countdata')) load('data/1_integrate_countdata.R')
if((!exists('allcodsigmean_isomerge'))||(!'availability'%in%colnames(allcodsigmean_isomerge))){
	base::source(here('src/Figures/Figure3/3_tRNA_array_analysis.R'))
} 
stagecolsdisplay <- c(E12.5 = "#214098", E14 = "#2AA9DF", E15.5 = "#F17E22", E17 = "#D14E28", P0 = "#ED3124")
displaystageconv <- names(stagecolsdisplay)%>% setNames(c("E13", "E145", "E16", "E175", "P0"))
stagecols <- stagecolsdisplay %>% setNames(c("E13", "E145", "E16", "E175", "P0"))
stageconv <- names(stagecols) %>% setNames(c("E13", "E145", "E16", "E175", "P0"))
GENETIC_CODE<-Biostrings::GENETIC_CODE
#id dfs
trgiddf=ids_nrgname%>%select(g_id=gene_id,tr_id=transcript_id)%>%distinct
gnm_gid=ids_nrgname%>%select(gnm=gene_name,g_id=gene_id)%>%distinct



decon_genegrps = list(
	CPN=c(0,5,15,10,4),
	ScPN=c(16,7,18,11,3),
	CthPN=c(6,9,2,12,17)
)
decon_genegrps<-decon_genegrps%>%enframe('decongrp','cluster')%>%
	unnest(cluster)

tbr2_ipc_gids = read_tsv('ext_data/tbr2_ipc_basal.tsv',col_names=FALSE)[[1]]

deconglist = readxl::read_xlsx('ext_data/DECON_detected.xlsx',skip=1)%>%
	rename(gnm='gene_id...1')%>%
	select(gnm,cluster)%>%
	left_join(decon_genegrps)%>%
	filter(!is.na(decongrp))
deconglist<-deconglist%>%inner_join(gnm_gid,by='gnm')
deconglist<-deconglist%>%{split(.$g_id,.$decongrp)}

telleygdf = 'ext_data/telley_weights_comb.xlsx'%>%readxl::read_xlsx(.)
telleygnmsplit = telleygdf[,3:4]%>%
	tail(-1)%>%set_colnames(c('gnm','w'))%>%
	mutate(set=ifelse(w>0,'Neuron','AP'))%>%
	inner_join(gnm_gid)%>%
	{split(.$g_id,.$set)}

#load info on dTE genes
dte_df = readxl::read_xlsx("tables/S2.xlsx",1,col_types=c(time='text'))
dte_df%<>%filter(time=='3')
dte_df$time='E175'
teupgenes = dte_df%>%filter(adj_p_value<0.05,log2fc>0)%>%.$gene_id
tedowngenes = dte_df%>%filter(adj_p_value<0.05,log2fc<0)%>%.$gene_id
dtegenes = c(teupgenes,tedowngenes)

#telley 2016 genes
N_allage=read_tsv('ext_data/telley_2016_neuronallages.tsv',col_names='gnm')%>%inner_join(gnm_gid)%>%.$g_id
P_allage=read_tsv('ext_data/telley_2016_progenallages.tsv',col_names='gnm')%>%inner_join(gnm_gid)%>%.$g_id
allglists = c(list(tbr2=tbr2_ipc_gids),deconglist,telleygnmsplit,
	list(dTEup=teupgenes,dTEdown=tedowngenes),
	list(N_allage=N_allage,P_allage=P_allage)
)

isoformab <- iso_tx_countdata$abundance%>%
	as.data.frame%>%
	rownames_to_column('tr_id')%>%
	pivot_longer(-tr_id,names_to='sample',values_to='TPM')%>%
	# gather(contrast,TPM,-tr_id)%>%
	separate(sample,into=c('time','assay','rep'))%>%
	group_by(time,assay,tr_id)%>%
	summarise(TPM=mean(TPM))
isoformab %<>% left_join(trgiddf,by='tr_id')
isoformab %<>% filter(assay=='total')

trcodusage <- codonfreqs%>%
	as.data.frame%>%
	rownames_to_column('tr_id')%>%
	gather(codon,count,-tr_id)

if(!file.exists(here('data/gcodusagedf.rds'))){
	gcodusagedf <- isoformab%>%
		left_join(trcodusage)%>%
		group_by(g_id,time,codon)%>%
		summarise(usage=sum(count*TPM))%>%
		group_by(g_id,time)%>%
		mutate(usage = usage/sum(usage))
	saveRDS(gcodusagedf,here('data/gcodusagedf.rds'))
}else{
	gcodusagedf<-readRDS(here('data/gcodusagedf.rds'))
}

gset = telleygnmsplit[[set1nm]]
get_set_pref<-function(gcodusagedf,gset){
	codusage <- gcodusagedf%>%filter(g_id%in%gset)%>%
		group_by(time,codon)%>%summarise_at(vars(usage),mean,na.rm=T)
	codusage$AA <- GENETIC_CODE[codusage$codon]
	codpref = codusage%>%group_by(time,AA)%>%mutate(pref = usage/sum(usage))
	codpref
}


set1nm<-"dTEdown"
set2nm<-"dTEup"
set1nm<-"P_allage"
set2nm<-"N_allage"
set1 <- allglists[[set1nm]]
set2 <- allglists[[set2nm]]
set1usage <- get_set_pref(gcodusagedf,set1)
set2usage <- get_set_pref(gcodusagedf,set2)
setratiodf = left_join(set1usage,set2usage,by=c('time','codon','AA'))%>%
	mutate(prefratio=log2(pref.y/pref.x))%>%
	mutate(usageratio=log2(usage.y/usage.x))
codonoccs <- read_tsv('tables/codonoccs_orig.tsv')%>%filter(fraction=='total')%>%select(-fraction)%>%
	dplyr::rename(dwell_time=occupancy)
rat_dt_df <- setratiodf%>%
	filter(!AA=='*')%>%
	safe_left_join(codonoccs,by=c('time','codon'))%>%
		select(time,codon,AA,prefratio,usageratio,dwell_time,pref.x,pref.y)
rat_dt_df%<>%group_by(time,AA)%>%mutate(AAcorDT=ifelse(n()==1,NA,mean(dwell_time)))

#now plot
plotfile<- here(paste0('plots/',set1nm,'_',set2nm,'codon_pref_vs_dt','.pdf'))
pdf(plotfile,w=15,h=3)
col1=sym('dwell_time')
col2=sym('prefratio')
corlabel = rat_dt_df%>%filter(is.finite(!!col1),is.finite(!!col2))%>%
		group_by(time)%>%
		summarise(tidy(cor.test(!!col1, !!col2)))
corlabel = corlabel%>%
	mutate(
		pformat=format(p.value,format='e',digits=4),
		pvalstring = ifelse(p.value > 0.001,round(p.value,4),pformat),
		labl=paste0('rho = ',round(estimate,3),'\n','pval = ',pvalstring))
nlabel=rat_dt_df%>%group_by(time)%>%summarise(labl=paste0('N=',n()))
rat_dt_df%>%
	ggplot(.,aes(x=usageratio	,y=dwell_time,color=AA))+
	# ggplot(.,aes(x=prefratio	,y=AAcorDT,color=AA))+
	# ggplot(.,aes(x=prefratio,y=dwell_time,color=AA))+
	geom_point()+
	geom_smooth(method='lm',aes(color=NULL))+
	facet_grid(.~time)+
	scale_color_discrete(guide=F)+
	geom_text(show.legend=F,data=corlabel,
			hjust=1,vjust=1,x= Inf,y=Inf,aes(color=NULL,label=labl))+
	geom_text(show.legend=F,data=nlabel,
			hjust=0,vjust=1,x= -Inf,y=Inf,aes(color=NULL,label=labl))
	scale_x_continuous(paste0('Pref Ratio'))+
	scale_y_continuous(paste0('Dwell Time'))+
	ggtitle(paste0('plots/',set1nm,'_',set2nm,'codon_pref_vs_dt','.pdf'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))


#now plot
library(ggrepel)
plotfile<- here(paste0('plots/',set1nm,'_',set2nm,'codon_pref_scatter','.pdf'))
pdf(plotfile)
rat_dt_df%>%
	ungroup()%>%
	ggplot(.,aes(x=pref.x,y=pref.y,label=AA))+
	# geom_point()+
	geom_text_repel()+
	facet_grid(time~.)+
	scale_x_continuous(paste0('Usage in ',set1nm))+
	scale_y_continuous(paste0('Usage in ',set2nm))+
	ggtitle(paste0('Usage Comparisons'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))




#get codon in wide format
allgcodwide = gcodusagedf%>%
	group_by(g_id,codon)%>%summarise(usage=mean(usage,na.rm=T))%>%
	pivot_wider(names_from=codon,values_from=usage)

#pca of codon usage
codpca = gcodwide%>%ungroup%>%select(AAA:TTT)%>%
	filter(is.finite(AAA))%>%
	princomp()

#caluclate codon PREFERENCE
corallgcodlong = gcodusagedf%>%
	group_by(g_id)%>%
	group_by(g_id,codon)%>%summarise(usage=mean(usage,na.rm=T))%>%
	mutate(AA = GENETIC_CODE[codon]	)%>%
	group_by(g_id,AA)%>%mutate(usage = usage/sum(usage))%>%
	group_by(g_id)%>%mutate(usage = replace_na(usage,0))%>%
	filter(!all(usage==0))%>%
	ungroup%>%select(g_id,codon,usage)

#get this in wide form
corallgcodwide <-	corallgcodlong%>%pivot_wider(names_from=codon,values_from=usage)

#get pcas of codon prefreence
codpca = corallgcodwide%>%ungroup%>%
	{set_rownames(select(.,AAA:TTT),.$g_id)}%>%
	princomp()

#now join our glists to the codon preference pca
codonprofiles <-readRDS('data/codonprofiles.rds')
offsets <- read_tsv('ext_data/offsets_manual.tsv')

codondts = codonprofiles%>%
    filter(sample%>%str_detect('ribo'))%>%
    # filter(fraction=='total')%>%
    separate(sample,c('time','assay','rep'))%>%
    safe_left_join(offsets%>%select(readlen,offset))%>%
    filter(position== -offset-3)%>%
    group_by(time,codon)%>%
    summarise(dwell_time = mean (occupancy))%>%
    filter(time=='E13')

# openxlsx::write.xlsx( file = "tables/S3.xlsx")
codsetdf = codpca$scores[,2]%>%enframe('g_id','codpref_pca2')%>%
	inner_join(allglists%>%enframe('gset','g_id')%>%unnest(g_id))


codonstats<-read_tsv('tables/tRNA_stat_df.tsv')
codondts = codonstats%>%filter(time=='E13')%>%filter(fraction=='total')%>%select(codon,dwell_time)
dt_pref_df = corallgcodlong%>%inner_join(codondts)%>%group_by(g_id)%>%
	mutate(usage = usage/sum(usage))%>%
	summarise(dt_pref=weighted.mean(dwell_time,usage,na.rm=T))

#now plot
library(ggrepel)
library(ggExtra)
plotfile<- here(paste0('plots/','codpref_vs_dtpref_gsets','.pdf'))
pdf(plotfile)
p=codsetdf%>%left_join(dt_pref_df)%>%
	group_by(gset)%>%summarise_at(vars(dt_pref,codpref_pca1),list(mean))%>%
	ggplot(.,aes(y=codpref_pca1,x=dt_pref,color=gset,label=gset,group=gset))+
	geom_point(size=I(1))+
	geom_text_repel()+
	# scale_x_continuous(paste0('xname'))+
	# scale_y_continuous(paste0('yname'))+
	# ggtitle(paste0('title'))+
	theme_bw()
ggMarginal(p,margins='both',type='density',groupColour=TRUE)
dev.off()
message(normalizePath(plotfile))

timecodondts = codonstats%>%filter(fraction=='total')%>%select(codon,time,dwell_time)
timedt_pref_df = corallgcodlong%>%inner_join(timecodondts)%>%group_by(time,g_id)%>%
	mutate(usage = usage/sum(usage))%>%
	summarise(dt_pref=weighted.mean(dwell_time,usage,na.rm=T))
plotfile<- here(paste0('plots/','time_codpref_vs_dtpref_gsets','.pdf'))
pdf(plotfile,w=15,h=3)
p=codsetdf%>%left_join(timedt_pref_df)%>%
	group_by(time,gset)%>%summarise_at(vars(dt_pref,codpref_pca2),list(mean))%>%
	ggplot(.,aes(y=codpref_pca2,x=dt_pref,color=gset,label=gset,group=gset))+
	geom_point(size=I(1))+
	geom_text_repel()+
	facet_grid(.~time,scale='free')+
	# scale_x_continuous(paste0('xname'))+
	# scale_y_continuous(paste0('yname'))+
	# ggtitle(paste0('title'))+
	theme_bw()
# ggMarginal(p,margins='both',type='density',groupColour=TRUE)
p
dev.off()
message(normalizePath(plotfile))

#
# corallgcodlong%>%
# 	filter(g_id %in% c(goid_ribo,goid_chrom))%>%
# 	mutate(cat = ifelse(g_id %in% goid_ribo,'ribo','pattern'))%>%
# 	group_by(cat,codon)%>%
# 	summarise(usage = mean(usage))%>%
# 	spread(cat,usage)


gcodwide<-gcodwide%>%left_join(dte_df%>%select(time,g_id=gene_id,log2fc))

#load info on go terms
#definite chrom and ribo genes.
GTOGO <- 'data/GTOGO.rds'%>%readRDS%>%select(gene_name,go_id,g_id=ensembl_gene_id)
goid_ribo = 'GO:0003735'
goid_chrom = 'GO:0003682'
mphase_goid = 'GO:0000087'
pattern_goid = 'GO:0007389'
ribogenes = GTOGO%>%filter(go_id==goid_ribo)%>%.$g_id
chromgenes = GTOGO%>%filter(go_id==goid_chrom)%>%.$g_id
mphasegenes = GTOGO%>%filter(go_id==mphase_goid)%>%.$g_id
patterngenes = GTOGO%>%filter(go_id==pattern_goid)%>%.$g_id
ribchrgenes<-c(ribogenes,chromgenes)

#load info on AP/NN genes


gcodwide = gcodusagedf%>%
	filter(g_id%in%dtegenes)%>%
	filter(g_id%in%highcountgenes)%>%
	group_by(g_id,codon)%>%summarise(usage=mean(usage,na.rm=T))%>%
	pivot_wider(names_from=codon,values_from=usage)
gcodwide%<>%filter(is.finite(AAA))
codpca = princomp(gcodwide%>%ungroup%>%select(AAA:TTT))
gcodwide%<>%mutate(ischrom=g_id%in%chromgenes)
gcodwide%<>%mutate(isteup=g_id%in%teupgenes)

#gcodwide%>%ungroup%>%mutate(pca1=codpca$scores[,5])%>%mutate(ischrom=g_id%in%chromgenes)%>%
#	{split(.$pca1,.$ischrom)}%>%{t.test(.[[1]],.[[2]])}

codons4occ <- colnames(codonfreqs)
codformula = as.formula(paste0('isteup ~ ',paste0(collapse='+',codons4occ)))
glm(codformula,data=gcodwide,family='binomial')%>%summary



gingold_types = 'ext_data/gingold_etal_2014_trna_types.tsv'%>%read_tsv
gingold_types=gingold_types%>%mutate(codon = DNAStringSet(anticodon)%>%reverseComplement%>%as.character)

ging_gcodlong = gcodusagedf%>%
	# filter(g_id%in%ribchrgenes)%>%
	filter(g_id%in%dtegenes)%>%
	filter(g_id%in%highcountgenes)%>%
	left_join(gingold_types%>%select(codon,type))%>%
	mutate(type = replace_na(type,'other'))%>%
	group_by(g_id,time,type)%>%
	summarise(usage=sum(usage,na.rm=T))

ging_gcodlong%<>%mutate(ischrom=g_id%in%chromgenes)
ging_gcodlong%<>%mutate(isteup=g_id%in%teupgenes)

gnmgiddf=ids_nrgname%>%select(g_id=gene_id,gene_name=gene_name)%>%distinct


ging_gcodwide <- ging_gcodlong%>%
	pivot_wider(names_from=type,values_from=usage)


glm(isteup~diff+prol,data=ging_gcodwide,family='binomial')%>%summary

#now plot
plotfile<- here(paste0('plots/','gingc_usage_dte','.pdf'))
pdf(plotfile)
ging_gcodlong%>%
	mutate(dTE = ifelse(isteup,'up','down'))%>%
	ggplot(.,aes(x=usage,fill=dTE))+
	geom_density(alpha=I(0.5))+
	facet_grid(time~type,scale='free')+
	scale_x_continuous(paste0('Codon class usage'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))

glm(chromform,data=gcodwide,family='binomial')%>%summary


codontable = 'tables/S3.xlsx'%>%readxl::read_xlsx(2)
codontable%>%.$codon%>%n_distinct

ging_gcodlong = gcodusagedf%>%
	# filter(g_id%in%ribchrgenes)%>%
	filter(g_id%in%highcountgenes)%>%
	left_join(gingold_types%>%select(codon,type))%>%
	mutate(type = replace_na(type,'other'))%>%
	group_by(g_id,time,type)%>%
	summarise(usage=sum(usage,na.rm=T))

dpratiodf = ging_gcodlong%>%left_join(gnmgiddf)%>%group_by(g_id,gene_name,type)%>%
	summarise(usage=sum(usage))%>%
	mutate(usage=usage/sum(usage))%>%
	spread(type,usage)%>%
	mutate(dpratio=diff/prol)

#now plot
plotfile<- here(paste0('plots/','dpratiodist','.pdf'))
pdf(plotfile)
dpratiodf%>%
	ggplot(.,aes(x=log2(dpratio)))+
	geom_density()+
	scale_x_continuous(paste0('Diff_Usage/Prol_Usage'))+
	theme_bw()+
	geom_vline(xintercept=dpratiodf%>%filter(gene_name=='Satb2')%>%.$dpratio%>%log2)
dev.off()
message(normalizePath(plotfile))

gnmgiddf%>%filter(gene_name=="Satb2")
'ENSMUSG00000038331'

fname = 'plots/diff_proliff_ratio_vs_techange.pdf'
dpratiodf%>%left_join(dte_df%>%select(time,g_id=gene_id,log2fc))%>%
   make_quantcompplot(log2fc,dpratio,fname)





################################################################################
########
################################################################################
	

slowcodons = tRNA_occ_df_en%>%filter(time=='E13',fraction=='total')%>%filter(dwell_time>1.5)%>%.$codon
trcodusage$AA=NULL
discr_codons = tRNA_occ_df_en%>%distinct(codon,AA)%>%mutate(slow=codon%in%slowcodons)%>%group_by(AA)%>%filter(!all(slow),!all(!slow))
all_codons = tRNA_occ_df_en%>%distinct(codon,AA)%>%mutate(slow=codon%in%slowcodons)

slowfracdf = trcodusage%>%
	inner_join(discr_codons,by='codon')%>%
	group_by(tr_id,AA)%>%
	summarise(slowfrac=mean(slow*count)/sum(count))

slowusagemat = spread(slowfracdf,AA,slowfrac)%>%ungroup%>%select(-tr_id)%>%as.matrix
slowusagemat = slowusagemat[!apply(slowusagemat,1,function(x)any(is.na(x))),]
#looks like controlling for amino acid usage doesn't work well.
slowusagemat%>%princomp%>%.$loadings

allslowfracdf = trcodusage%>%
	group_by(tr_id)%>%
	mutate(frac=count/sum(count))
allslowfracdf%>%select(tr_id,codon,frac)%>%spread(codon,frac)
	inner_join(all_codons)%>%
	summarise(slowfrac=mean(slow*count)/sum(count))

ncodonfreqmat = codonfreqs%>%sweep(1,FUN='/',STAT=rowSums(.))

slowfrac = log2(ncodonfreqmat[,all_codons%>%filter(slow)%>%.$codon]%>%rowSums)
nonslowfrac = log2(ncodonfreqmat[,all_codons%>%filter(!slow)%>%.$codon]%>%rowSums)


longtrids = cdsgrl%>%width%>%sum%>%enframe('tr_id','length')%>%mutate(gid=trid2gid[[tr_id]])%>%group_by(gid)%>%slice(which.max(length))%>%.$tr_id

	geneslowfrac = (slowfrac[longtrids]-nonslowfrac[longtrids])%>%.[is.finite(.)]

	genemoreslow = geneslowfrac%>%{. > median(.,na.rm=T)}

	onts = c('BP','MF','CC')
	
	source('src/R/Functions/go_term_funcs.R')
stopifnot(exists("GTOGO"))
if(!'gene_id'%in%colnames(GTOGO))GTOGO%<>%mutate(gene_id=ensembl_gene_id)
	get_cluster_gos = function(clustvect){
    out=map_df(.id='ontology',onts%>%setNames(.,.),function(ont){
      lapply(unique(clustvect)%>%setNames(.,.),(function(clustval){
            gids <- names(clustvect)
            stopifnot(mean(names(clustvect)%in%GTOGO$gene_id)>.9)
            filtGTOGO <- GTOGO %>%filter(gene_id %in%names(clustvect))
             projmemoise(rungo)(
              gids[clustvect==clustval],
              filtGTOGO,
              ont,
              algo='classic'
            )
        }))%>%map('result')%>%bind_rows(.id='cluster')
    })
    out
}

	genemoreslow%>%{names(.)%<>%trid2gid[[.]];.}%>%get_cluster_gos

codpca$loadings[,1]


trcodusage%>%group_by(tr_id)%>%group_slice(1:10)%>%group_by(tr_id,AA)


countpred_df<-readRDS('data/countpred_df.rds')
gidtedf = countpred_df%>%select(g_id=gene_id,contrast,TE=logFC)%>%filter(str_detect(contrast,'TE'))%>%separate(contrast,c('time','tmp'))
gcodusagedf%>%left_join(gidtedf)




	teform = as.formula(paste0('TE ~ ',paste0(collapse='+',codons4occ)))
	gcodus_te_model = gcodusagedf%>%filter(time=='E13')%>%spread(codon,usage)%>%left_join(gidtedf)%>%lm(data=.,teform)
	temodcoefdf = gcodus_te_model%>%.$coefficients%>%.[codons4occ]%>%enframe('codon','temodel_coef')%>%mutate(time=itime)

	temodcoefdf%>%left_join(codonoccs)%>%{quicktest(log2(.$temodel_coef),.$dwell_time)}
	temodcoefdf%>%left_join(trna_ab_df%>%filter(fraction=='total'))%>%{quicktest(log2(.$temodel_coef),.$dwell_time)}
	stop()
	gcodavaildf = gcodusagedf%>%
		# filter(time=='E13')%>%
		safe_left_join(allow_missing=TRUE,tRNA_occ_df%>%filter(fraction=='total')%>%select(time,codon,availability))%>%
		filter(!is.na(availability))%>%
		group_by(g_id,time)%>%summarise(totalavail = logSumExp(log2(usage)+availability))%>%
		left_join(gidtedf)


	teavail_glm4 = gcodavaildf%>%group_by(g_id)%>%filter(!is.na(TE))%>%group_slice(1:10000)%>%filter(n()==5)%>%lm(data=.,TE ~ g_id+totalavail)
	anova(teavail_glm4)

	glmfit = glm4(count ~ 0 + gene + codon+a_codon+phase, data=sitedf,family=negative.binomial(mtheta),MXITER=400,doFit=T, sparse=T, verbose=T)
	gcodavaildf%>%lm(data=.,TE ~ g_id+totalavail)%>%
		confint


get_tidy_pval <- function(x) x %>%tidy%>%.$p.value%>%round(4)%>%{ifelse(. == 0,' p < 0.0001 ',paste0('p = ',.))}	
wtest <- function(df,contcol,factcol){
	contcol=quo_text(contcol)
	factcol=quo_text(factcol)
    stopifnot(n_distinct(df[[factcol]])==2)
    df%>%{split(.[[contcol]],.[[factcol]])}%>%{wilcox.test(.[[1]],.[[2]])}
}
get_double_pvalstring <- function(df, contcol,fct1,fct2){
	contcol = enquo(contcol)
	stopifnot(quo_text(contcol)%in%colnames(df))
	fct1 = enquo(fct1)
	fct2 = enquo(fct2)
	paste(sep='\n',
		paste0('Up: wilcox test', get_tidy_pval(wtest(df%>%filter(!(!!(fct1))), contcol, fct2 ) )),
		paste0('down: wilcox test', get_tidy_pval(wtest(df%>%filter(!(!!fct2)), contcol, fct1) ) )
	)
}
make_density_threefact_fig <- function(df,classcolname,contcol,fct1,fct2,xname,plotfile){
	fct1 = enquo(fct1)
	fct2 = enquo(fct2)
	classcolname = enquo(classcolname)
	contcol = enquo(contcol)
	fct1text = quo_text(fct1)
	fct2text = quo_text(fct2)
	stopifnot(quo_text(contcol)%in%colnames(df))
	#
	# xlabs = xbreaks[-c(1,length(xbreaks))]
	#
	lclip <- df[[quo_text(contcol)]]%>%quantile(0.01)
	rclip <- df[[quo_text(contcol)]]%>%quantile(1-0.01)
	df%<>%filter(dplyr::between(!!contcol,lclip,rclip))
	#
	xlims <- df[[quo_text(contcol)]]%>%range
	xspan <- xlims%>%dist
	tickwidth <- 10^(round(log10(xspan))-1) * 2.5
	xlimlow <- floor(xlims[1]/(tickwidth))*tickwidth
	xlimhigh <- ceiling(xlims[2]/(tickwidth))*tickwidth
	xlims <- c(xlimlow,xlimhigh)
	xbreaks <- seq(xlimlow,xlimhigh,tickwidth)
	#
	pdf(plotfile)
	{df%>%
		filter(!( (!!fct1) & (!!fct2) ))%>%
		# filter(time!='P0')%>%
		mutate(!!classcolname := case_when(
			!!fct1==1 ~ fct1text,
			!!fct2==1 ~ fct2text,
			TRUE ~ 'neither'
		))%>%
		ggplot(.,aes(x=!!contcol,color=!!classcolname))+geom_density(alpha=I(0.01))+theme_bw()+
			scale_x_continuous(name=xname,limits=xlims,labels=xbreaks,breaks=xbreaks )+
			geom_text(data=data.frame(label=get_double_pvalstring(df,!!contcol,!!fct1,!!fct2)),aes(label=label),color=I('black'),x=Inf,y=Inf,vjust=1,hjust=1)
		}%>%print
	dev.off()
	normalizePath(plotfile)
}

usedcodonfreqsnorm <- usedcodonfreqs[,]%>%{./colSums(t(.))}
timeoccscores <- map_df(.id='time',times,function(itime){
	timeoccvect <- tRNA_occ_df%>%
		filter(time==itime)%>%
		mutate(dwell_time=dwell_time-mean(na.rm=T,dwell_time))%>%
		{setNames(.$dwell_time,.$codon)[codons4occ]}
	usedcodonfreqsnorm[,names(timeoccvect)]%*%matrix(timeoccvect,ncol=1)%>%.[,1]%>%
	enframe('protein_id','elongscore')
})

fractions <- unique(tRNA_occ_df$fraction) 
timeSigscores <- map_df(.id='time',times,function(itime){
	map_df(.id='fraction',fractions%>%setNames(.,.),function(fractioni){
		timetrrnavect <- tRNA_occ_df%>%
			filter(time==itime,fraction==fractioni)%>%
			 mutate(abundance=ifelse(abundance %in% -Inf,min(abundance[is.finite(abundance)],na.rm=T)-1,abundance))%>%
			mutate(abundance=abundance-mean(na.rm=T,abundance))%>%				
			{setNames(.$abundance,.$codon)[codons4occ]}
		timetrrnavect <- timetrrnavect[!is.na(timetrrnavect)]
		
		# (t(usedcodonfreqsnorm[,names(timetrrnavect)])%*%(timetrrnavect))%>%
		# 	colMeans(na.rm=T)%>%
		# 	enframe('protein_id','abundance')

	usedcodonfreqsnorm[,names(timetrrnavect)]%*%matrix(2^timetrrnavect,ncol=1)%>%.[,1]%>%
	log2%>%
	enframe('protein_id','abundance')

	})
})

timeaAvailscores <- map_df(.id='time',times,function(itime){
	map_df(.id='fraction',fractions%>%setNames(.,.),function(fractioni){
		timetrrnavect <- tRNA_occ_df%>%
			filter(time==itime,fraction==fractioni)%>%
			 mutate(availability=ifelse(availability %in% -Inf,min(availability[is.finite(availability)],na.rm=T)-1,availability))%>%
			 mutate(availability=availability-mean(na.rm=T,availability))%>%				
			{setNames(.$availability,.$codon)[codons4occ]}
timetrrnavect <- timetrrnavect[!is.na(timetrrnavect)]

			usedcodonfreqsnorm[,names(timetrrnavect)]%*%matrix(2^timetrrnavect,ncol=1)%>%.[,1]%>%
			log2%>%
			enframe('protein_id','availability')

	})
})

#
score_techange_df<-timeoccscores%>%
	left_join(timeSigscores,by=c('time','protein_id'))%>%
	left_join(timeaAvailscores,by=c('time','fraction','protein_id'))%>%
	left_join(ms_id2protein_id%>%distinct(protein_id,gene_name))%>%
	left_join(allTEchangedf)

#Data frame showing elongation change over time for each gene
occchange_vs_te_df <- score_techange_df%>%
	group_by(up,down,protein_id)%>%
	filter(fraction=='Total')%>%
	# filter(time!='P0')%>%
	# group_slice(1:2)%>%
	nest%>%
	mutate(elongchange = map_dbl(data,~{lm(data=.,elongscore/mean(elongscore) ~ seq_along(time))$coef[2]}))%>%
	mutate(tps = map_dbl(data,nrow))%>%
	select(-data)

tRNAchange_vs_te_df <- score_techange_df%>%
filter(fraction=='Total')%>%
	group_by(fraction,up,down,protein_id)%>%
	# filter(time!='P0')%>%
	# group_slice(1:2)%>%
	nest%>%
	mutate(tRNA_score_change = map_dbl(data,~{lm(data=.,(abundance) ~ seq_along(time))$coef[2]}))%>%
	mutate(tps = map_dbl(data,nrow))%>%
	select(-data)



avail_change_vs_te_df <- score_techange_df%>%
	filter(fraction=='Total')%>%
	group_by(fraction,up,down,protein_id)%>%
	# filter(time!='P0')%>%
	# group_slice(1:2)%>%
	mutate(ntime=seq_along(time))%>%
	# select(ntime,availability)%>%
	filter(is.finite(availability))%>%
	nest%>%
	mutate(avail_change = map_dbl(data,~lm.fit(matrix(.$availability),matrix(c(rep(1,5),.$ntime),nrow=5,ncol=2))$coef[2]))%>%
	mutate(tps = map_dbl(data,nrow))%>%
	select(-data)


make_density_threefact_fig(df=occchange_vs_te_df,
	classcolname = TEchange_class,
	contcol = elongchange,
	fct1=up,fct2=down,
	plotfile='plots/figures/figure2/trna_codons/te_change_vs_occscorechange.pdf',
	xname='Predicted Elongation Rate Change - Riboseq Occ'
)

make_density_threefact_fig(df=tRNAchange_vs_te_df,
	classcolname = TEchange_class,
	contcol = tRNA_score_change,
	fct1=up,fct2=down,
	plotfile='plots/figures/figure2/trna_codons/te_change_vs_tRNAscorechange.pdf',
	xname='Predicted Change in tRNA Abundance Score'
)

make_density_threefact_fig(df=avail_change_vs_te_df,
	classcolname = TEchange_class,
	contcol = avail_change,
	fct1=up,fct2=down,
	plotfile='plots/figures/figure2/trna_codons/te_change_vs_availchange.pdf',
	xname='Predicted Change in tRNA availability Score'
)

###So is this the same genes???
library(txtplot)
# tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%{txtplot(.$tRNA_score_change,.$elongchange)}
# avail_change_vs_te_df%>%left_join(occchange_vs_te_df)%>%{cor.test(.$avail_change,.$elongchange)}

# tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%{txtplot(.$tRNA_score_change,.$elongchange)}
# avail_change_vs_te_df%>%left_join(occchange_vs_te_df)%>%{cor.test(.$tRNA_score_change,.$elongchange)}

# tRNAchange_vs_te_df%>%left_join(occchange_vs_te_df)%>%
# 	filter(elongchange>2,tRNA_score_change<4)



# testvals%>%filter(p.value<0.05)

# library(broom)


################################################################################
########Now let's do the same, but with transcriptional change rather than TE change.
################################################################################


foldchangecatdf <- readRDS(here('data/foldchangecatdf.rds'))

timeoccscores%<>%select(-matches('gene_id'))
timeoccscores%<>%safe_left_join(mcols(cds)%>%as.data.frame%>%distinct(protein_id,gene_id))

allTXNchangedf <- foldchangecatdf%>%group_by(gene_id)%>%
	mutate(up = any((transcriptional_logFC>0)&(transcriptional_adj.P.Val<0.05) ))%>%
	mutate(down = any((transcriptional_logFC<0)&(transcriptional_adj.P.Val<0.05) ))%>%
	summarise(up=any(up),down=any(down))%>%mutate_all(list(~replace_na(.,FALSE)))%>%
	filter(!(up&down))


score_txnchange_df<-timeoccscores%>%
	left_join(timeSigscores,by=c('time','protein_id'))%>%
	left_join(timeaAvailscores,by=c('time','fraction','protein_id'))%>%
	left_join(ms_id2protein_id%>%distinct(protein_id,gene_name))%>%
	left_join(allTXNchangedf)

#Data frame showing elongation change over time for each gene
occchange_vs_txn_df <- score_txnchange_df%>%
	group_by(up,down,protein_id)%>%
	filter(fraction=='Total')%>%
	# filter(time!='P0')%>%
	# group_slice(1:2)%>%
	nest%>%
	mutate(elongchange = map_dbl(data,~{lm(data=.,elongscore ~ seq_along(time))$coef[2]}))%>%
	mutate(tps = map_dbl(data,nrow))%>%
	select(-data)

tRNAchange_vs_txn_df <- score_txnchange_df%>%
	filter(fraction=='Total')%>%
	group_by(fraction,up,down,protein_id)%>%
	# filter(time!='P0')%>%
	# group_slice(1:2)%>%
	nest%>%
	mutate(tRNA_score_change = map_dbl(data,~{lm(data=.,abundance ~ seq_along(time))$coef[2]}))%>%
	mutate(tps = map_dbl(data,nrow))%>%
	select(-data)

avail_change_vs_txn_df <- score_txnchange_df%>%
	filter(fraction=='Total')%>%
	group_by(fraction,up,down,protein_id)%>%
	# filter(time!='P0')%>%
	# group_slice(1:2)%>%
	mutate(ntime=seq_along(time))%>%
	select(ntime,availability)%>%
	filter(is.finite(availability))%>%
	arrange(protein_id,fraction)%>%
	nest%>%
	# .$data%>%.[1]
	# %>%select(availability,ntime)
	mutate(avail_change = map_dbl(data,~lm.fit(matrix(.$availability),matrix(c(rep(1,length(.$ntime)),.$ntime),nrow=length(.$ntime),ncol=2))$coef[2]))%>%
	mutate(tps = map_dbl(data,nrow))%>%
	select(-data)

#

make_density_threefact_fig(df=occchange_vs_txn_df,
	classcolname = TXN_change,
	contcol = elongchange,
	fct1=up,fct2=down,
	plotfile='plots/figures/figure2/trna_codons/txn_change_vs_occscorechange.pdf',
	xname='Predicted Elongation Rate Change - Riboseq Occ'
)
make_density_threefact_fig(df=tRNAchange_vs_txn_df,
	classcolname = TXN_change,
	contcol = tRNA_score_change,
	fct1=up,fct2=down,
	plotfile='plots/figures/figure2/trna_codons/txn_change_vs_tRNAscorechange.pdf',
	xname='Predicted Change in tRNA Abundance Score'
)
make_density_threefact_fig(df=avail_change_vs_txn_df,
	classcolname = TXN_change,
	contcol = avail_change,
	fct1=up,fct2=down,
	plotfile='plots/figures/figure2/trna_codons/txn_change_vs_availchange.pdf',
	xname='Predicted Change in tRNA availability Score'
)

plotfile='plots/figures/figure2/trna_codons/txn_change_vs_tRNAscorechange.pdf'
pdf(plotfile)
tRNAchange_vs_txn_df%>%
	filter(!(up&down))%>%
	# filter(time!='P0')%>%
	mutate(Txn_change_class = case_when(
		up==1 ~ 'up',
		down==1 ~ 'down',
		TRUE ~ 'neither'
	))%>%
	ggplot(.,aes(x=tRNA_score_change,color=Txn_change_class))+geom_density(alpha=I(0.01))+theme_bw()+
		scale_x_continuous()
dev.off()
normalizePath(plotfile)
