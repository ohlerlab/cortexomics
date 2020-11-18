################################################################################
################################################################################
load('data/1_integrate_countdata.R')
base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
}
# base::source(here::here('src/R/Rprofile.R'))
# base::source(here::here('src/Figures/Figures0/1_ms_data_import.R'))

trid2gid<-load_hashmap('trid2gid.hmp')
gid2gnm<-load_hashmap('gid2gnm.hmp')
gnm2gid<-load_hashmap('gnm2gid.hmp')
trid2gnm<-load_hashmap('trid2gnm.hmp')
gnm2gid<-load_hashmap('gnm2gid.hmp')
gnm2trid<-load_hashmap('gnm2trid.hmp')

#
protgnmtrids<-readRDS('data/protgnmtrids.rds')

proteinmsdata<-readRDS('data/proteinmsdata.rds')
proteinmsdata$ms_id = proteinmsdata[[1]]
proteinmsdatamatch <- proteinmsdata
proteinmsdatamatch$cdsidx = match(proteinmsdatamatch$gene_name,fmcols(cdsgrl,gene_name))
proteinmsdatamatch <- proteinmsdatamatch%>%filter(!is.na(cdsidx))
ms_transcripts <- protgnmtrids$tr_id
#select iBAQ and input cols
sigcols <- proteinmsdatamatch%>%colnames%>%str_subset('iBAQ')%>%str_subset('input')
proteinmsdatamatch%<>%select(ms_id,gene_name,g_id,one_of(sigcols))
ms_mat <- proteinmsdatamatch[,-c(1:3)]%>%as.matrix
#
rownames(ms_mat) <- proteinmsdatamatch$ms_id
colnames(ms_mat) <- colnames(ms_mat)%>%
	str_replace('iBAQ\\.','')%>%
	str_replace('input','MS')%>%
	str_replace('rep','')%>%
	str_replace('p(\\d)','\\1')
ms_mat = proDA::median_normalization(ms_mat)

################################################################################
########Import it all to get tr length scaled counts
################################################################################
{
library(tximport)
library(tidyverse)

rnasalmonfiles = Sys.glob(here('pipeline/salmon/data/*total*/quant.sf'))
dpoutfiles = Sys.glob(here('pipeline/deepshapeprime/fakesalmonfiles/*ribo*/*'))

allquantfiles = c(rnasalmonfiles,dpoutfiles)
names(allquantfiles) <- allquantfiles%>%dirname%>%basename

dptrs = dpoutfiles[[1]]%>%fread%>%.$Name
salmontrs = allquantfiles[[1]]%>%fread%>%.$Name%>%str_extract('ENSMUST\\w+')
# inclusiontable(dptrs,salmontrs)
trs = intersect(dptrs,salmontrs)


####KEY STEP
ms_trs = intersect(trs,ms_transcripts)

ms_tx_countdata = tximport(files=allquantfiles,
	ignoreTxVersion=TRUE,
	tx2gene=mcols(cds)[,c('transcript_id','gene_id')]%>%as.data.frame%>%distinct(transcript_id,gene_id),
	type='salmon',
	countsFromAbundance='scaledTPM',
	importer=function(file){
		read_tsv(file,col_types=cols())%>%
			mutate(Name=str_extract(Name,'ENSMUST\\w+'))%>%
			filter(Name%in%ms_trs)%>%arrange(match(Name,trs))
})

randomround = function(x)floor(x)+rbinom(length(x),1,x%%1)
ms_tx_countdata$counts%<>%as.data.frame%>%
	mutate_all(randomround)%>%
	as.matrix%>%
	set_rownames(rownames(ms_tx_countdata$counts))

ms_tx_countdata$counts%>%as.data.frame%>%
	rownames_to_column('gene_id')%>%
	mutate(gene_name = gid2gnm[[gene_id]])%>%
	select(-gene_id)%>%select(gene_name,everything())%>%
	mutate_at(vars(-gene_name),list(randomround))%>%
	write_tsv('data/ms_tx_scaled_countData.tsv')
}

{
sharedgenes = intersect(ms_tx_countdata$abundance%>%rownames,tx_countdata$abundance%>%rownames)

quicktest(
	log2(ms_tx_countdata$abundance[sharedgenes,1]),
	log2(tx_countdata$abundance[sharedgenes,1])
)
ms_matgidrows <- ms_mat%>%set_rownames(proteinmsdatamatch$g_id)
mssharedgenes <- intersect(rownames(ms_matgidrows),
	intersect(
		ms_tx_countdata$abundance[,]%>%rownames,
		tx_countdata$abundance[highcountgenes,]%>%rownames)
		)


quicktest(
	log2(ms_matgidrows[mssharedgenes,1]),
	log2(tx_countdata$abundance[mssharedgenes,'E13_total_1'])
)

quicktest(
	log2(ms_matgidrows[mssharedgenes,1]),
	log2(tx_countdata$abundance[mssharedgenes,'E13_ribo_1'])
)


quicktest(
	log2(ms_matgidrows[mssharedgenes,1]),
	log2(ms_tx_countdata$abundance[mssharedgenes,'E13_total_1'])
)

quicktest(
	log2(ms_matgidrows[mssharedgenes,1]),
	log2(ms_tx_countdata$abundance[mssharedgenes,'E13_ribo_1'])
)
library(LSD)
source('Applications/LSD/R/LSD.heatscatter.R')

sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
proDAfitms<-readRDS('data/proDAfitms.rds')
sel_prodpreds<-readRDS('data/sel_prodpreds.rds')
contrdf<-readRDS('data/contrdf.rds')
sel_ms_mat<-readRDS('data/sel_ms_mat.rds')



countcontr_df <- readRDS(here('data/countcontr_df.rds'))


}

# stop()

stage_count='E13'
stage_ms='E13'


tps = sel_ms_mat%>%colnames%>%str_extract('[^_]+')%>%unique
tps %<>%setNames(.,.)
cassays = c('total','ribo')%>%setNames(.,.)

allxtail = Sys.glob('pipeline/xtail/*')%>%map_df(.id='time',fread)%>%group_by(gene_name)
allxtail$gene_id = gnm2gid[[ allxtail$gene_name]]

txnchangedf <- countcontr_df%>%filter(assay=='all')%>%group_by(gene_id)%>%
  mutate(sig = (adj.P.Val < 0.05)& (abs(logFC)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (logFC > 0))),
    down = as.numeric(any(sig & (logFC < 0)))
  )

ribochangedf <- countcontr_df%>%filter(assay=='ribo')%>%group_by(gene_id)%>%
  mutate(sig = (adj.P.Val < 0.05)& (abs(logFC)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (logFC > 0))),
    down = as.numeric(any(sig & (logFC < 0)))
  )

techangedf <- allxtail%>%group_by(gene_id,gene_name)%>%
  mutate(sig = (adj_p_value < 0.05)& (abs(log2fc)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (log2fc > 0))),
    down = as.numeric(any(sig & (log2fc < 0)))
  )

#or use ms
mschangedf = contrdf%>%
  group_by(gene_id)%>%
  mutate(sig = (adj_pval < 0.05)& (abs(diff)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (diff > 0))),
    down = as.numeric(any(sig & (diff < 0)))
  )

mschangedf%>%group_by(up,down)%>%tally
techangedf%>%group_by(up,down)%>%tally
ribochangedf%>%group_by(up,down)%>%tally
techangedf%>%group_by(up,down)%>%tally

genesets = list(
	all = techangedf$gene_id,
	te_up = txnchangedf%>%filter(up==1)%>%.$gene_id,
	te_down = txnchangedf%>%filter(down==1)%>%.$gene_id,
	te_up = techangedf%>%filter(up==1)%>%.$gene_id,
	te_down = techangedf%>%filter(down==1)%>%.$gene_id,
	ms_up = mschangedf%>%filter(up==1)%>%.$gene_id,
	ms_down = mschangedf%>%filter(down==1)%>%.$gene_id
)



avtpms = ms_tx_countdata$abundance[mssharedgenes,]%>%
	as.data.frame%>%
	rownames_to_column('gene_id')%>%
	gather(contrast,TPM,-gene_id)%>%
	separate(contrast,into=c('time','assay','rep'))%>%
	group_by(time,assay,gene_id)%>%summarise(TPM=mean(TPM))%>%
	unite(contrast,c('time','assay'))
ms_metadf <- readRDS('data/ms_metadf.rds')

itime = 'E13'

msscatterlist = lapply(cassays,function(countassay){
	lapply(genesets,function(geneset){
	lapply(tps,function(stage_count){
	lapply(tps,function(stage_ms){
	countname = paste0('log10(TPM) ',countassay,'_',stage_count)
	MSname = paste0('MS_',stage_ms)
	misscol = paste0(stage_ms,'_anymissing')
	istagepresent <- ms_metadf%>%filter(!(!!sym(misscol)))

	scatdf = avtpms%>%
		filter(!str_detect(contrast,'_TE'))%>%
		separate(contrast,c('time','assay')	)%>%
		filter(time==stage_count)%>%
		filter(assay==countassay)%>%
		filter(gene_id%in%geneset)%>%
		left_join(sel_prodpreds%>%filter(time==stage_ms)%>%select(gene_id,diff))%>%
		filter(gene_id %in% mssharedgenes)%>%
		inner_join(istagepresent)%>%
		filter(TPM>1)%>%
		mutate(TPM = log10(TPM))%>%
		select(gene_id,countvals=TPM,msvals=diff)
	#
	testtidy = tidy(cor.test(scatdf$countvals,scatdf$msvals))
	acor_label=paste0('rho = ',round(testtidy$conf.low,3),' - ',round(testtidy$conf.high,3))
	# acor_label=paste0('rho = ',round(testtidy$estimate,3),'\n','pval = ',round(testtidy$p.value,3))
	pdf('tmp.pdf')
	gplot = heatscatter((scatdf$countvals),scatdf$msvals,ggplot=TRUE)+
		scale_x_continuous(countname)+
		scale_y_continuous(MSname)+
		ggtitle(acor_label)
	dev.off()
	list(gplot,testtidy)
})
})
})
})



stop()
#now plot

plotfile<- here(paste0('plots/','ms_count_scatter','.pdf'))
pdf(plotfile)
msscatterlist[[2]][[1]][[1]][[1]]
dev.off()
normalizePath(plotfile)


totcormat = msscatterlist[[1]][['all']]%>%map_depth(2,~.[[2]]$estimate)%>%map_df(.id='stage_count',~bind_rows(.id='stage_ms',.))#%>%spread(stage_ms,cor)%>%{as.matrix(.[,-1])}%>%set_rownames(.,colnames(.))
ribocormat = msscatterlist[[2]][['all']]%>%map_depth(2,~.[[2]]$estimate)%>%map_df(.id='stage_count',~bind_rows(.id='stage_ms',.))#%>%spread(stage_ms,cor)%>%{as.matrix(.[,-1])}%>%set_rownames(.,colnames(.))


countmscors = msscatterlist%>%map_df(.id='assay',.%>%map_df(.id='geneset',.%>%map_depth(2,~.[[2]]$estimate)%>%map_df(.id='stage_count',~bind_rows(.id='stage_ms',.))))

#now plot
plotfile<- here(paste0('plots/','cortilesribo','.pdf'))
pdf(plotfile,w=10,h=15)
countmscors%>%
	ggplot(data=.,aes(fill=cor,x=stage_count,y=stage_ms))+
	geom_tile()+
	geom_text(aes(label = round(cor,2)))+
	facet_grid(geneset~assay)+
	scale_x_discrete(paste0('Stage Ribo-seq'))+
	scale_y_discrete(paste0('Stage MS'))+
	scale_fill_continuous(low='#8904B1',middle='white',high='#FF8000')+
	ggtitle(paste0('Stage-stage correlations MS vs Ribo-seq'))+
	theme_bw()
dev.off()
normalizePath(plotfile)

techangedf = read_tsv(here('tables/xtailTEchange.tsv'))



################################################################################
########Time correlations
################################################################################
cests = countpred_df%>%filter(str_detect(contrast,'total'))%>%distinct(gene_id,contrast,logFC)%>%spread(contrast,logFC)
rcests = countpred_df%>%filter(str_detect(contrast,'ribo'))%>%distinct(gene_id,contrast,logFC)%>%spread(contrast,logFC)
pests = sel_prodpreds%>%distinct(gene_id,time,diff)%>%spread(time,diff)

rnacors = pests%>%left_join(cests)%>%.[,-1]%>%apply(1,function(x)cor(x[1:5],x[6:10]))
ribocors = pests%>%left_join(rcests)%>%.[,-1]%>%apply(1,function(x)cor(x[1:5],x[6:10]))

corggdf = bind_rows(.id='assay',rnaseq=data.frame(cor=rnacors),ribo=data.frame(cor=ribocors))
#now plot
plotfile<- here(paste0('plots/','timecordensity','.pdf'))
pdf(plotfile)
corggdf%>%
	ggplot(.,aes(x=cor,color=assay))+
	geom_density()+
	# facet_grid(assay~.)+
	# scale_color_discrete(name='assay',colorvals)+
	scale_x_continuous(paste0('Temporal Correlation'))+
	ggtitle(paste0('Correlations for Riboseq vs RNAseq'))+
	theme_bw()
dev.off()
normalizePath(plotfile)


cassaynames = c('RNAseq','Riboseq')%>%setNames(c('all','ribo'))

stepcountcontrdf <- readRDS(here('data/stepcountcontrdf.rds'))
stepstepprotcontrdf <- readRDS(here('data/stepstepprotcontrdf.rds'))

for(cassay in c('all','ribo')){
steplcdf = stepcountcontrdf%>%mutate(csig = adj.P.Val<0.05)%>%select(gene_id,time,logFC,assay,csig)%>%
	filter(assay==cassay)%>%
	inner_join(stepstepprotcontrdf%>%mutate(psig = adj_pval<0.05)%>%select(gene_id,time,prot_logFC=diff,psig))%>%
	group_by(gene_id)
csigname = paste0(cassaynames[cassay],' Sig')
steplcdf%<>%mutate(sigset = case_when(
	csig & psig ~ 'Both Sig',
	psig ~ 'MS Sig',
	csig ~ csigname,
	TRUE ~ 'None'
))
lfccorlabels = steplcdf%>%group_by(time)%>%filter(is.finite(logFC),is.finite(prot_logFC))%>%summarise(tidy(cor.test(logFC,prot_logFC)))%>%
	mutate(label = paste0('rho = ',round(conf.low,3),'-',round(conf.high,3),'\np = ',format(p.value,sci=T,digits=2)%>%str_replace('.*0.0.*',' < 1.0e-300')
))

sisetcolors = c('None'='grey',csigname='blue','MS Sig'='red','Both Sig'='Purple')
names(sisetcolors)%<>%str_replace('csigname',csigname)
steplcdf%>%arrange(-match(sigset,names(sisetcolors)))%>%.$sigset%>%setdiff(names(sisetcolors))
#now plot
plotfile<- here(paste0('plots/','lfc_cors_',cassay,'.pdf'))
pdf(plotfile,w=16,h=4)
print(steplcdf%>%arrange(-match(sigset,names(sisetcolors)))%>%
	ggplot(aes(x=logFC,y=prot_logFC,color=sigset))+facet_grid(.~time)+geom_point(size=I(0.2),aes(alpha=sigset!='None'))+
	geom_text(data=lfccorlabels,aes(label=label,color=NULL),y=I(-Inf),x=I(-Inf),vjust=0,hjust=0)+
	scale_color_manual(name='colorname',values=sisetcolors)+
	scale_x_continuous(paste0('RNAseq Fold Change'),limits=c(-5,5))+
	scale_y_continuous(paste0('MS Fold Change'),limits=c(-5,5))+
	ggtitle(paste0('Fold Change Comparison - MS vs ',cassaynames[cassay]))+
	theme_bw())
dev.off()
message(normalizePath(plotfile))

}

stepcountcontrdf%>%group_by(gene_id)%>%group_slice(1)

################################################################################
########Look at binarized trajectories - concordances
################################################################################

mschangedf%>%
	left_join(techangedf,by='gene_id',suffix=c('_MS','_TE'))%>%
	left_join(txnchangedf,by='gene_id')%>%
	group_by(up,down,up_MS,down_MS,up_TE,down_TE)%>%tally%>%as.data.frame
	filter(! (up_MS&down_MS),!(up_TE & down_TE))

txntrajcatdf = mschangedf%>%
	left_join(txnchangedf,by='gene_id',suffix=c('_MS','_txn'))%>%
	mutate_at(vars(matches('up|down')),as.logical)%>%
	mutate(txngroup = case_when(
			up_txn & down_txn ~ 'txn-Both',
			up_txn ~ 'txn-Up',
			down_txn ~ 'txn-Down',
			TRUE ~ 'txn-NoChange'
		),msgroup = case_when(
			up_MS & down_MS ~ 'MS-Both',
			up_MS ~ 'MS-Up',
			down_MS ~ 'MS-Down',
			TRUE ~ 'MS-NoChange'
		))%>%
	group_by(txngroup,msgroup)
txntrajcatdf %>% tally%>%
	arrange(
		desc(txngroup%>%str_detect('Change')),desc(txngroup%>%str_detect('Both')),desc(txngroup%>%str_detect('Down')),
		desc(msgroup%>%str_detect('Change')),desc(msgroup%>%str_detect('Both')),desc(msgroup%>%str_detect('Down')))%>%
	group_by(txngroup)%>%mutate(Freq = round(n/sum(n),3))


ribontrajcatdf = mschangedf%>%
	left_join(ribochangedf,by='gene_id',suffix=c('_MS','_ribo'))%>%
	mutate_at(vars(matches('up|down')),as.logical)%>%
	mutate(ribogroup = case_when(
			up_ribo & down_ribo ~ 'Ribo-Both',
			up_ribo ~ 'Ribo-Up',
			down_ribo ~ 'Ribo-Down',
			TRUE ~ 'Ribo-NoChange'
		),msgroup = case_when(
			up_MS & down_MS ~ 'MS-Both',
			up_MS ~ 'MS-Up',
			down_MS ~ 'MS-Down',
			TRUE ~ 'MS-NoChange'
		))%>%
	group_by(ribogroup,msgroup)
ribontrajcatdf %>%	tally%>%
	arrange(
		desc(ribogroup%>%str_detect('Change')),desc(ribogroup%>%str_detect('Both')),desc(ribogroup%>%str_detect('Down')),
		desc(msgroup%>%str_detect('Change')),desc(msgroup%>%str_detect('Both')),desc(msgroup%>%str_detect('Down')))%>%
	group_by(ribogroup)%>%mutate(Freq = round(n/sum(n),3))


tetrajcatdf = mschangedf%>%
	left_join(techangedf,by='gene_id',suffix=c('_MS','_te'))%>%
	mutate_at(vars(matches('up|down')),as.logical)%>%
	mutate(tegroup = case_when(
			up_te & down_te ~ 'te-Both',
			up_te ~ 'te-Up',
			down_te ~ 'te-Down',
			TRUE ~ 'te-NoChange'
		),msgroup = case_when(
			up_MS & down_MS ~ 'MS-Both',
			up_MS ~ 'MS-Up',
			down_MS ~ 'MS-Down',
			TRUE ~ 'MS-NoChange'
		))%>%
	group_by(tegroup,msgroup)
tetrajcatdf %>%	tally%>%
	arrange(
		desc(tegroup%>%str_detect('Change')),desc(tegroup%>%str_detect('Both')),desc(tegroup%>%str_detect('Down')),
		desc(msgroup%>%str_detect('Change')),desc(msgroup%>%str_detect('Both')),desc(msgroup%>%str_detect('Down')))%>%
	group_by(tegroup)%>%mutate(Freq = round(n/sum(n),3))


#now all 3

txnte_trajdf = mschangedf%>%
	left_join(txnchangedf,by='gene_id',suffix=c('_MS','_txn'))%>%
	left_join(techangedf,by='gene_id')%>%rename('up_te':=up,'down_te':=down)%>%
	mutate_at(vars(matches('up|down')),as.logical)%>%
	mutate(txngroup = case_when(
			up_txn & down_txn ~ 'txn-Both',
			up_txn ~ 'txn-Up',
			down_txn ~ 'txn-Down',
			TRUE ~ 'txn-NoChange'
		),msgroup = case_when(
			up_MS & down_MS ~ 'MS-Both',
			up_MS ~ 'MS-Up',
			down_MS ~ 'MS-Down',
			TRUE ~ 'MS-NoChange'
		),tegroup = case_when(
			up_te & down_te ~ 'te-Both',
			up_te ~ 'te-Up',
			down_te ~ 'te-Down',
			TRUE ~ 'te-NoChange'
		))%>%
	group_by(txngroup,msgroup,tegroup)



txnte_trajdf %>% tally%>%
	arrange(
		desc(txngroup%>%str_detect('Change')),desc(txngroup%>%str_detect('Both')),desc(txngroup%>%str_detect('Down')),
		desc(msgroup%>%str_detect('Change')),desc(msgroup%>%str_detect('Both')),desc(msgroup%>%str_detect('Down')),
		desc(tegroup%>%str_detect('Change')),desc(tegroup%>%str_detect('Both')),desc(tegroup%>%str_detect('Down')))%>%
	group_by(txngroup)%>%mutate(Freq = round(n/sum(n),3))%>%
	filter(!str_detect(txngroup,'Both'),!str_detect(tegroup,'Both'),!str_detect(msgroup,'Both'))%>%
	as.data.frame








#so in general TE doesn't necessariy track txn change that much
mschangedf%>%
	left_join(techangedf,by='gene_id',suffix=c('_MS','_TE'))%>%
	left_join(txnchangedf,by='gene_id')%>%
	group_by(up,down,up_MS,down_MS,up_TE,down_TE)%>%tally%>%as.data.frame
	filter(! (up_MS&down_MS),!(up_TE & down_TE))	



if(F){

countpred_df<-readRDS('data/countpred_df.rds')
tx_countdata<-readRDS('data/tx_countdata.rds')
  
contrastlens = tx_countdata$length%>%
	as.data.frame%>%
	rownames_to_column('gene_id')%>%
	gather(contrast,length,-gene_id)%>%
	separate(contrast,into=c('time','assay','rep'))%>%
	group_by(time,assay,gene_id)%>%summarise(length=mean(length))%>%
	unite(contrast,c('time','assay'))

avtpms = tx_countdata$abundance[,]%>%
	as.data.frame%>%
	rownames_to_column('gene_id')%>%
	gather(contrast,TPM,-gene_id)%>%
	separate(contrast,into=c('time','assay','rep'))%>%
	group_by(time,assay,gene_id)%>%summarise(TPM=mean(TPM))%>%
	unite(contrast,c('time','assay'))

avcounts = tx_countdata$counts[,]%>%
	as.data.frame%>%
	rownames_to_column('gene_id')%>%
	gather(contrast,count,-gene_id)%>%
	separate(contrast,into=c('time','assay','rep'))%>%
	group_by(time,assay,gene_id)%>%summarise(count=mean(count))%>%
	unite(contrast,c('time','assay'))


avcounts%>%left_join(contrastlens)%>%{quicktest(log2(.$count),log2(.$length))}

avms = sel_ms_mat%>%
	as.data.frame%>%
	rownames_to_column('gene_id')%>%
	gather(contrast,diff,-gene_id)%>%
	separate(contrast,into=c('time','assay','rep'))%>%
	group_by(time,assay,gene_id)%>%summarise(diff=mean(diff,na.rm=TRUE))%>%ungroup
	# unite(contrast,c('time','assay'))


avms = sel_ms_mat%>%
	as.data.frame%>%
	rownames_to_column('gene_id')%>%
	gather(contrast,diff,-gene_id)%>%
	separate(contrast,into=c('time','assay','rep'))%>%
	group_by(time,assay,gene_id)%>%summarise(diff=mean(diff,na.rm=TRUE))%>%ungroup
	# unite(contrast,c('time','assay'))

countpred_df%>%filter(!str_detect(contrast,'_TE'))%>%safe_left_join(contrastlens)%>%
	group_by(contrast)%>%group_slice(5)%>%select(logFC,length)%>%{quicktest(.$logFC,log2(.$length))}

countpred_df%>%filter(!str_detect(contrast,'_TE'))%>%safe_left_join(contrastlens)%>%
	separate(contrast,c('time','assay')	)%>%
	left_join(sel_prodpreds%>%select(time,gene_id,diff))%>%
	group_by(time,assay)%>%group_slice(1)%>%select(logFC,diff)%>%{quicktest(.$logFC,.$diff)}

countpred_df%>%filter(!str_detect(contrast,'_TE'))%>%safe_left_join(contrastlens)%>%
	separate(contrast,c('time','assay')	)%>%
	left_join(sel_prodpreds%>%select(time,gene_id,diff))%>%
	filter(gene_id %in% mssharedgenes)%>%
	# mutate(logFC = logFC - log(length))%>%
	filter(assay=='ribo')%>%
	group_by(time,assay)%>%group_slice(1)%>%select(logFC,diff)%>%{quicktest(.$logFC,.$diff)}


avtpms%>%filter(!str_detect(contrast,'_TE'))%>%safe_left_join(contrastlens)%>%
	separate(contrast,c('time','assay')	)%>%
	left_join(sel_prodpreds%>%select(time,gene_id,diff))%>%
	filter(gene_id %in% mssharedgenes)%>%
	# mutate(TPM = TPM - log(length))%>%
	filter(assay=='ribo')%>%
	group_by(time,assay)%>%group_slice(1)%>%select(TPM,diff)%>%{quicktest(log2(.$TPM),.$diff)}


itime = 'E13'
misscol = paste0(itime,'_anymissing')

ms_metadf <- readRDS('data/ms_metadf.rds')
istagepresent <- ms_metadf%>%filter(!(!!sym(misscol)))

#These averages are teh ones to beast
avtpms%>%filter(!str_detect(contrast,'_TE'))%>%safe_left_join(contrastlens)%>%
	separate(contrast,c('time','assay')	)%>%
	left_join(sel_prodpreds%>%select(time,gene_id,diff))%>%
	# left_join(avms%>%select(time,gene_id,diff))%>%
	# left_join(avms%>%select(time,gene_id,diff))%>%
	filter(gene_id %in% mssharedgenes)%>%
	# mutate(TPM = TPM - log(length))%>%
	semi_join(istagepresent)%>%
	filter(assay=='total')%>%
	filter(TPM>1)%>%
	filter(time==itime)%>%
	{tgenest <<- .$gene_id;.}%>%
	group_by(time,assay)%>%group_slice(1)%>%select(TPM,diff)%>%{quicktest(log2(.$TPM),.$diff)}


avtpms%>%filter(!str_detect(contrast,'_TE'))%>%separate(contrast,c('time','assay'))%>%




countpred_df%>%
	filter(gene_id%in%tgenest)%>%
	filter(!str_detect(contrast,'_TE'))%>%safe_left_join(contrastlens)%>%
	separate(contrast,c('time','assay')	)%>%
	left_join(sel_prodpreds%>%select(time,gene_id,diff))%>%
	filter(gene_id %in% mssharedgenes)%>%
	# mutate(logFC = logFC - log(length))%>%
	filter(assay=='ribo')%>%
	group_by(time,assay)%>%group_slice(2)%>%select(logFC,diff)%>%{quicktest(.$logFC,.$diff)}


	# {quicktest(log2(.$TPM),.$diff)}



avms = sel_ms_mat%>%
	as.data.frame%>%
	rownames_to_column('gene_id')%>%
	gather(contrast,diff,-gene_id)%>%
	separate(contrast,into=c('time','assay','rep'))%>%
	group_by(time,assay,gene_id)%>%summarise(diff=mean(diff,na.rm=TRUE))%>%ungroup
	# unite(contrast,c('time','assay'))

#Hmmm so sometimes our sel_prodpreds are MUCH lower even though this should be the complete genes only...
avms%>%select(time,gene_id,diff)%>%
	left_join(sel_prodpreds%>%select(time,gene_id,diff),by=c('time','gene_id'))%>%
	inner_join(istagepresent)%>%
	filter(time==itime)%>%
	{quicktest(.$diff.x,.$diff.y)}

tgene = avms%>%select(time,gene_id,diff)%>%
	left_join(sel_prodpreds%>%select(time,gene_id,diff),by=c('time','gene_id'))%>%
	inner_join(istagepresent)%>%
	filter(time==itime)%>%
	arrange(-abs(log2(diff.x/diff.y)))%>%select(time,gene_id,diff.y,diff.x)%>%.$gene_id%>%.[1]

#gene which differs a lot
'ENSMUSG00000026701'
tmsid = ms_metadf%>%filter(gene_id==tgene)%>%.$ms_id
#This gene should agree because none of it should be missing....
proDAfitms$abundances[tmsid,]

}