library(tidyverse)
library(data.table)
library(magrittr)
library(readxl)
source('src/Figures/Figure0/0_load_annotation.R')
################################################################################
########
################################################################################
# Correlation panel
panel.cor <- function(x, y){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits=2)
    txt <- paste0("R = ", r)
    cex.cor <- 0.7/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y)
}


ixnoselongfiles = Sys.glob('../Ribotransformer/pipeline/ixnos_elong/*/*.elong.csv')
ixnoselongfiles%<>%setNames(basename(dirname(.)))
ixnoseelongdata = ixnoselongfiles%>%map_df(.id='sample',read_csv)
ixnoseelongdata$tr_id%<>%str_replace('\\.\\d+','')
ixnoseelongdata%<>%separate(sample,c('time','assay','rep'))
load('data/1_integrate_countdata.R')

ixnoseelongdata%<>%group_by(tr_id)%>%filter(between(elong,0.5,1.5))%>%ungroup
ixnoseelongspread = ixnoseelongdata%>%unite(sample,time,assay,rep)%>%spread(sample, elong)

#pretty well correlated between samples, like rho = .92
left_join(
	ixnoseelongdata%>%filter(time=='E13',rep=='1')%>%select(tr_id,elong1=elong),
	ixnoseelongdata%>%filter(time=='E16',rep=='2')%>%select(tr_id,elong2=elong),
	by='tr_id'
)%>%{quicktest(.$elong1,.$elong2)}

# Create the plots

my_fn <- function(data, mapping, ...){
      p <- ggplot(data = data, mapping = mapping) + 
        stat_density2d(aes(fill=..density..), geom="tile", contour = FALSE) +
        scale_fill_gradientn(colours=rainbow(100))
      p
}
library(GGally)
#now plot
plotfile<- here(paste0('plots/','elong_repscat','.pdf'))
pdf(plotfile)
ixnoseelongspread%>%
	head(1000)%>%
	.[,-1]%>%as.data.frame%>%ggpairs(., lower=list(continuous=my_fn),upper = list(continuous = wrap("cor", size = 3)))
dev.off()
message(normalizePath(plotfile))

ixnoseelongdata%<>%filter(!((time=='E13') & (rep == 2)))


tr_elong_df = iso_tx_countdata$abundance%>%as.data.frame%>%
	rownames_to_column('tr_id')%>%
	gather(sample, tpm, -tr_id)%>%
	separate(sample,c('time','assay','rep'))%>%
	left_join(ixnoseelongdata%>%select(-assay),by=c('tr_id','time','rep'))%>%
	left_join(ids_nrgname%>%
	select(tr_id=transcript_id, gene_id=gene_id))%>%
	filter(!is.na(elong))%>%
	filter(gene_id%in%highcountgenes)

tr_elong_df <- tr_elong_df%>%
	group_by(gene_id, tr_id, assay, time)%>%
	summarise_at(vars(tpm,elong),list(mean))

#okay so elongation rate correlates with rnaseq
tr_elong_df%>%filter(time=='E13',assay=='total')%>%filter(tpm>1)%>%{quicktest(.$elong,log2(.$tpm))}
#but riboseq more so.
tr_elong_df%>%filter(time=='E13',assay=='ribo')%>%filter(tpm>1)%>%{quicktest(.$elong,log2(.$tpm))}
#what about TE?
#now plot
elongte_df<-tr_elong_df%>%
	spread(assay,tpm)%>%
	mutate(TE = log2(ribo/total))

labltibble <-   elongte_df%>%filter(tpm>1)%>%group_by(time)%>%
	summarise(labl=tidy(cor.test(elong,TE))%>%
		{paste0('rho= ',round(.$estimate,4),'p < 10e-15')})
	
labltibble

plotfile<- here(paste0('plots/','ixnos_vs_te','.pdf'))
pdf(plotfile)
elongte_df%>%
	ggplot(.,aes(x=elong,y=TE))+
	geom_point()+geom_smooth(method='lm')+
	facet_wrap(time~.)+
	scale_x_continuous(paste0('Mean Dwell Time (Ixnos)'))+
	geom_text(show.legend=F,data=labltibble,
		hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))+
	scale_y_continuous(paste0('log2(Ribo/RNAseq)'))+
	ggtitle(paste0('Ixnos DT vs TE'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))




#actually this should usee weighting when calculating TE for genes

labltibble <-  
elongte_df%>%
	group_by(gene_id, time)%>%summarise_at(vars(elong,TE,ribo,total),mean)%>%
	group_by(time)%>%summarise(labl=tidy(cor.test(elong,TE))%>%{paste0('rho= ',round(.$estimate,4),'p < 10e-15')})


labltibble
plotfile<- here(paste0('plots/','ixnos_vs_te_gn','.pdf'))
pdf(plotfile)
elongte_df%>%
	group_by(gene_id, time)%>%summarise_at(vars(elong,TE),mean)%>%
	ggplot(.,aes(x=elong,y=TE))+
	geom_point()+geom_smooth(method='lm')+
	facet_wrap(time~.)+
	scale_x_continuous(paste0('Mean Dwell Time (Ixnos)'))+
	geom_text(show.legend=F,data=labltibble,
		hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))+
	scale_y_continuous(paste0('log2(Ribo/RNAseq)'))+
	ggtitle(paste0('Ixnos DT vs TE'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))




################################################################################
########predicts MS?
################################################################################



elongte_df_gnm <- elongte_df%>%
	group_by(gene_id, time, rep)%>%
	summarise(ribo=sum(ribo), total = sum(total),
		TE = mean(TE*(ribo/sum(ribo)),na.rm=T),
		elong = mean(elong*(ribo/sum(ribo)),na.rm=T)
	)

ms_data = readxl::read_xlsx(path = "tables/S1.xlsx",sheet='matched_ms_data')

tpconv = setNames(c("E13", "E145", "E16", "E175", "P0"), c("E12.5","E14","E15.5","E17","P0"))
summsdata = ms_data%>%gather(dset, MS, -gene_id)%>%
	separate(dset,c('time', 'assay', 'rep'), sep='_')%>%
	mutate(time=tpconv[time])%>%
	group_by(time,gene_id)%>%
	summarise_at(vars(MS),mean,na.rm=T)


summsdata%>%filter(time=='E13')%>%left_join(elongte_df_gnm)%>%
	lm(data=., MS ~ log2(ribo) + elong)%>%anova

summsdata%>%filter(time=='E13')%>%left_join(elongte_df_gnm)%>%
	lm(data=., MS ~ log2(ribo) + log2(elong))%>%anova

summsdata%>%filter(time=='E13')%>%left_join(elongte_df_gnm)%>%
	lm(data=., MS ~ log2(ribo) + log2(elong))%>%tidy

#now plot
# BiocManager::install(c('LSD'))
library(LSD)
library(ggpubr)
#
base::source('https://raw.githubusercontent.com/stineb/LSD/master/R/LSD.heatscatter.R')
#
tps<-summsdata$time%>%unique
plotfile<- here(paste0('plots/','elong_MS_comp','.pdf'))
pdf(plotfile)
ggarrange(plotlist = lapply(tps,function(tp){
	summsdata%>%filter(time==tp)%>%left_join(elongte_df_gnm)%>%{
	x=filter(., is.finite(ribo), is.finite(MS), is.finite(elong))
	# lm(data=., MS ~ log2(ribo) + log2(elong))%>%anova
	lm(data=x, MS ~ log2(ribo))$residuals%>%
	heatscatter(y=.,x=log2(x$elong),cor=TRUE,ggplot=TRUE)+
	geom_smooth(method='lm')+
	xlab('log2(MeanDT)')+ylab('Residual(log2(MS) ~ log2(ribo)')
}}))
dev.off()
message(normalizePath(plotfile))


summsdata%>%filter(time=='E13')%>%left_join(elongte_df_gnm)%>%
	lm(data=., MS ~ total + ribo + elong)%>%anova




tr_elong_df%>%filter(assay=='total')%>%
	group_by(time, tr_id)%>%
	group_slice(1)%>%
	summarise(elong = mean(elong,na.rm=T))

gnm_elong_df<-tr_elong_df%>%
	group_by(gene_id, time, rep)%>%summarise(elong=mean(elong*(tpm/sum(tpm)), na.rm=T))%>%
	select(gene_id,  time, rep,  elong)


################################################################################
########Codon based elongartion rate
################################################################################
'/fast/AG_Ohler/dharnet/cortexomics/src/Figures/Figure3/3_tRNA_array_analysis.R'

weighted_codon_usage <- lapply(times, function(itime) {
    # stopifnot(rownames(codonfreqs)==cdsexprdf$transcript_id)
    (codonfreqs[cdsexprdf$transcript_id, ] * cdsexprdf[[itime]]) %>%
    sweep(.,F='/',STAT=rowSums(.),MARGIN=1)%>%
    as.data.frame%>%rownames_to_column('tr_id')%>%gather(codon,wuse,-tr_id)
}) %>% bind_rows(.id = "time")

offsets <- read_tsv('ext_data/offsets_manual.tsv')
codonprofiles <- readRDS(here('data/codonprofiles.rds'))
codonoccs<- codonprofiles%>%
    filter(sample%>%str_detect('ribo'))%>%
    safe_left_join(offsets%>%select(readlen,offset))%>%
    filter(position== -offset-3)%>%
    separate(sample,c('time','assay','rep'))%>%
    group_by(time,codon)%>%
    summarise(dwell_time = mean (occupancy))

mymeandt = weighted_codon_usage%>%left_join(codonoccs)%>%
	group_by(time, tr_id)%>%
	summarise(meanDT = sum(dwell_time*wuse,na.rm=T))%>%
	left_join(ids_nrgname%>%select(gene_id, tr_id = transcript_id))%>%
	group_by(time, gene_id)%>%
	left_join(elongte_df_gnm)%>%
	summarise(ribo=sum(ribo), total = sum(total),
		TE = mean(TE*(ribo/sum(ribo)),na.rm=T),
		elong = mean(elong*(ribo/sum(ribo)),na.rm=T),
		meanDT = mean(meanDT*(ribo/sum(ribo)),na.rm=T)
	)

myelongdf = summsdata%>%left_join(mymeandt)
tps<-mymeandt$time%>%unique
plotfile<- here(paste0('plots/','Myelong_MS_comp','.pdf'))
pdf(plotfile)
ggarrange(plotlist = lapply(tps,function(tp){
	summsdata%>%left_join(mymeandt)%>%filter(time==tp)%>%{
	df=filter(., is.finite(ribo), is.finite(MS), is.finite(meanDT),total>1)
	res=lm(data=df, MS ~ log2(ribo))$residuals
	# lm(data=., MS ~ log2(ribo) + log2(elong))%>%anova
	labltibble <-   tidy(cor.test(res,log2(df$meanDT)))%>%
		{paste0('rho= ',round(.$estimate,4),'p < 10e-15')}%>%
		tibble(labl=.)
	heatscatter(y=res,x=log2(df$meanDT),cor=TRUE,ggplot=TRUE)+
	geom_text(show.legend=F,data=labltibble,
		hjust=0,vjust=1,x= -Inf,y=Inf,size=I(3),aes(label=labl))+
	geom_smooth(method='lm')+
	xlab('log2(MeanDT)')+ylab('Residual(log2(MS) ~ log2(ribo)')
}}))
dev.off()
message(normalizePath(plotfile))

tp=tps[1]
lapply(tps,function(tp){
	summsdata%>%left_join(mymeandt)%>%filter(time==tp)%>%
	filter(., is.finite(ribo), is.finite(MS), is.finite(meanDT))%>%
	filter(., 0<(ribo), 0<(MS), 0<(meanDT))%>%
	lm(data=., MS ~ log2(ribo) + log2(meanDT))%>%anova
})

################################################################################
########Vs change?
################################################################################
	


gnm_elong_df<-tr_elong_df%>%
	group_by(gene_id, time, rep)%>%summarise(elong=mean(elong*(tpm/sum(tpm)), na.rm=T))%>%
	select(gene_id,  time, rep,  elong)

elongchange_df = gnm_elong_df%>%
	# separate(sample, c('time', 'assay', 'rep'))%>%
	group_by(gene_id, time)%>%summarise_at(vars(elong), list(mean))%>%
	group_by(gene_id)%>%
	summarise(elongchange = log2(elong[time=='E175']/elong[time=='E13']))%>%
	select(gene_id, elongchange)

techangedf = readxl::read_xlsx(path = "tables/S2.xlsx",sheet=1)


compdf = techangedf%>%filter(time==3)%>%	
	left_join(elongchange_df,by='gene_id')%>%
	select(gene_id, log2fc, adj_p_value, elongchange)
	
compdf%>%filter(adj_p_value < 0.05)%>%{quicktest(.$log2fc,.$elongchange)}
compdf%>%filter(abs(elongchange)%>%{.>quantile(.,0.1,na.rm=T)})%>%{quicktest(.$log2fc,.$elongchange)}


