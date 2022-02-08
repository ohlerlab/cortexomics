################################################################################
################################################################################
base::source(here::here('src/Rprofile.R'))
if(!exists("cdsgrl")) {
	load('data/1_integrate_countdata.R')
	base::source("src/Preprocess/0_load_annotation.R")
}
library(scales)
library(dplyr)

# conflict_prefer("intersect", "dplyr")
# conflict_prefer("exprs", "Biobase")
intersect = dplyr::intersect
exprs = Biobase::exprs
#A4 size is 8 by 11

HIGHCOUNTTHRESH <- 32

{

allgids = trid2gidv[alltrs]%>%unique
allgnms = trid2gnmv[alltrs]%>%unique
}
#THis statement from dunn and weismann 2013:
##A4 size is 8 by 11Remarkably, we find that the range of translation efficiencies for different messages spans four orders of magnitude, a range comparable to that observed for mRNA abundance of well-counted genes (Figure 1D). Moreover, translation efficiency is uncorrelated with mRNA abundance (r2 = 8.29 × 10−5; Figure 1E) and mRNA abundance predicts only one third of the variance in the rate of protein production as measured by ribosome footprint density (Figure 1F). Translational regulation is therefore a major determinant of gene expression in the early embryo (supplementary table 1 at Dryad: Dunn et al., 2013), and ribosome profiling provides a quantitative and robust means to monitor translational regulation during development.

stagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
stageconv = names(stagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
#Get the counts at each timepoint.

# countexprdata <- readRDS(countfile)
# countfile='pipeline/exprdata/countexprset.rds'

TEs = Sys.glob('pipeline/xtail/*.txt')%>%setNames(.,basename(.))%>%map_df(.%>%fread%>%select(feature_id=gene_name,dplyr::matches('log2TE')))
TEs = bind_rows(
	TEs%>%select(feature_id,E13_log2TE)%>%group_by(feature_id)%>%summarise(val=mean(E13_log2TE))%>%mutate(time = 'E13'),
	TEs%>%select(-E13_log2TE)%>%gather(time,val,-feature_id)%>%mutate(time = str_extract(time,'[^_]+'))
)


HIGHCOUNTTHRESH =32

txi_counttbl_reps <- tx_countdata$counts%>%
	as.data.frame%>%
	rownames_to_column('gene_id')%>%
	pivot_longer(cols=-gene_id,names_to='dataset',values_to='count')%>%
	separate(dataset,c('time','assay','rep'))%>%
	pivot_wider(names_from=rep,values_from=count,names_prefix='Replicate_')%>%
	# mutate(highcount=ishighcount[gene_id])
	# mutate(highcount = allcountkeep[gene_id])%>%
	mutate(highcount = Replicate_1+Replicate_2 >= HIGHCOUNTTHRESH)

testgn = 'ENSMUSG00000064356'
txi_tetbl_reps <- tx_countdata$abundance%>%as.data.frame%>%
	rownames_to_column('gene_id')%>%
	gather(dataset,val,-gene_id)%>%
	separate(dataset,into=c('time','assay','replicate'))%>%
	spread(assay,val)%>%
	mutate(te = log2(ribo) - log2(total))%>%
	select(gene_id,time,replicate,te)%>%
	mutate(replicate = ifelse(replicate==1,'Replicate_1','Replicate_2'))%>%
	spread(replicate,te)%>%
	left_join(txi_counttbl_reps%>%filter(assay=='ribo')%>%distinct(gene_id,time,highcount))

repcors = txi_counttbl_reps%>%group_by(time,assay)%>%filter(highcount)%>%summarise(cor=cor(use='complete',Replicate_1,Replicate_2))



terepcors = txi_tetbl_reps%>%group_by(time)%>%filter(highcount)%>%filter(is.finite(Replicate_1))%>%
	filter(is.finite(Replicate_2))%>%
	filter(Replicate_1%>%between(-12,12))%>%
	filter(Replicate_2%>%between(-12,12))%>%
	summarise(cor=cor(use='complete',Replicate_1,Replicate_2))

txi_tetbl_reps%>%filter(gene_id==testgn)

plotfile <- 'plots/QC_plots/fig1b_b_tete.pdf'
cairo_pdf(plotfile,w=12,h=3)
txi_tetbl_reps%>%arrange(highcount)%>%
	filter(Replicate_1%>%between(-12,12))%>%
	filter(Replicate_2%>%between(-12,12))%>%
ggplot(.,aes(x=((2^Replicate_1)),y=((2^Replicate_2)),color=highcount))+
	scale_color_manual(values = setNames(c('black','grey'),c('TRUE','FALSE')))+
	geom_point(size=I(0.2))+
	facet_grid( ~ time)+
	scale_x_log10(name = 'TE Replicate 1',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(name = 'TE Replicate 2',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	geom_text(data=terepcors%>%mutate(text = paste0('rho = ',round(cor,3))),x=I(0),y=I(3),color=I('black'),aes(label=text))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))


#
plotfile <- 'plots/QC_plots/fig1b_rpfrpf.pdf'
cairo_pdf(plotfile,w=12,h=3)
ggplot(txi_counttbl_reps%>%filter(assay=='ribo')%>%arrange(highcount),aes(x=(1+Replicate_1),y=(1+Replicate_2),color=highcount))+
	scale_color_manual(values = setNames(c('black','grey'),c('TRUE','FALSE')))+
	geom_point(size=I(0.2))+
	facet_grid( ~ time)+
	scale_x_log10(name = 'RPF Counts Replicate 1',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(name = 'RPF Counts Replicate 2',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	geom_text(data=repcors%>%filter(assay=='ribo')%>%mutate(text = paste0('p = ',round(cor,3))),x=1,y=4,color=I('black'),aes(label=text))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))
#
plotfile <- 'plots/QC_plots/fig1b_b_rnarna.pdf'
cairo_pdf(plotfile,w=12,h=3)
ggplot(txi_counttbl_reps%>%filter(assay=='total')%>%arrange(highcount),aes(x=(1+Replicate_1),y=(1+Replicate_2),color=highcount))+
	scale_color_manual(values = setNames(c('black','grey'),c('TRUE','FALSE')))+
	# scale_color_manual(values = stagecols)+
	# scale_fill_manual(values=stagecols)+
	geom_point(size=I(0.2))+
	facet_grid( ~ time)+
	scale_x_log10(name = 'RNA Counts Replicate 1',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(name = 'RNA Counts Replicate 2',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	geom_text(data=repcors%>%filter(assay=='total')%>%mutate(text = paste0('p = ',round(cor,3))),x=1,y=4,color=I('black'),aes(label=text))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))


repavexprtbl<-
	txi_counttbl_reps%>%group_by(gene_id,time)%>%
	summarise(highcount = any(highcount[assay=='ribo']),mRNA=(Replicate_2+Replicate_1)[assay=='total']/2,ribo=(Replicate_2+Replicate_1)[assay=='ribo']/2)


mRNA_ribo_cortbl<-repavexprtbl%>%filter(highcount)%>%group_by(time)%>%summarise(cor = cor(use='complete',log2(1+mRNA),log2(1+ribo)))



'plots/QC_plots/fig1b_b_rna_rpf.pdf'%>%dirname%>%dir.create(rec=TRUE)
plotfile <- 'plots/QC_plots/fig1b_b_rna_rpf.pdf'
cairo_pdf(plotfile,w=12,h=3)
ggplot(repavexprtbl%>%arrange(highcount),aes(x=mRNA,y=ribo,color=highcount))+
	scale_color_manual(values = setNames(c('black','grey'),c('TRUE','FALSE')))+
	geom_point(size=I(0.2))+
	facet_grid( ~ time)+
	scale_x_log10(name = 'Mean RNAseq Counts',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(name = 'Mean Riboseq Counts',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	geom_text(data=	mRNA_ribo_cortbl%>%mutate(text = paste0('p = ',round(cor,3))),x=1,y=4,color=I('black'),aes(label=text))+
	theme_bw()
dev.off()
plotfile%>%normalizePath%>%message



{
mRNAvals <- tx_countdata$abundance%>%
	as.data.frame%>%rownames_to_column('gene_id')%>%
	pivot_longer(-gene_id,names_to='dataset',values_to='TPM')%>%
	separate(dataset,c('time','assay','rep'))%>%
	filter(assay=='total')%>%
	group_by(time,gene_id)%>%
	summarise(val=log2(mean(TPM)),type='mRNA')

TEs$gene_id = gnm2gidv[TEs$feature_id]
ggdf = bind_rows(mRNAvals,TEs)

#add high count col
ggdf%<>%safe_left_join(txi_counttbl_reps%>%filter(assay=='ribo')%>%select(gene_id,time,highcount)%>%distinct)
ggdf %<>% mutate(time = stageconv[time])

ggTPMvals = ggdf%>%filter(type=='mRNA',highcount)%>%group_by(time)%>%mutate(val=val-median(val))%>%arrange(highcount)
ggTEvals = ggdf%>%filter(type=='TE',highcount)%>%group_by(time)%>%mutate(val=val-median(na.omit(val)))
library(scales)

l10scales =   scale_x_log10(name = 'TE / TPM',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))
plotfile <- 'plots/QC_plots/fig1b_c_TE_rna_dists.pdf'
cairo_pdf(plotfile,w=12,h=3)
ggplot(
	data = ggTPMvals,
	aes(x=2^val))+
	geom_density(aes(fill=time),alpha=I(0.5))+
	geom_density(
		data=ggTEvals,alpha=I(1),
	aes(color=time))+
	facet_grid(~time)+
	scale_color_manual(values = stagecols)+
	scale_fill_manual(values = stagecols)+
	# scale_x_log10()+
	l10scales+
	theme_bw()
dev.off()
plotfile%>%normalizePath%>%message

}

{

mRNA_TE_tbl <- allvoom$E%>%as.data.frame%>%
	rownames_to_column('gene_id')%>%
	gather(dataset,val,-gene_id)%>%
	separate(dataset,into=c('time','assay','replicate'))%>%
	spread(assay,val)%>%
	group_by(gene_id,time)%>%
	summarise(TE = mean(ribo) -  mean(total) )
ribo_TE_tbl <- tx_countdata$abundance%>%as.data.frame%>%
	rownames_to_column('gene_id')%>%
	gather(dataset,val,-gene_id)%>%
	separate(dataset,into=c('time','assay','replicate'))%>%
	spread(assay,val)%>%
	group_by(gene_id,time)%>%
	summarise(TE = log2(mean(ribo)) -  log2(mean(total) ))
ribotpms <- tx_countdata$abundance%>%
	as.data.frame%>%rownames_to_column('gene_id')%>%
	pivot_longer(-gene_id,names_to='dataset',values_to='TPM')%>%
	separate(dataset,c('time','assay','rep'))%>%
	filter(assay=='ribo')%>%
	group_by(time,gene_id)%>%
	summarise(TPM=log2(mean(TPM)))
# mRNA_TE_tbl%<>%safe_left_join(x=.,y=txi_counttbl_reps%>%filter(assay=='total')%>%select(gene_id,time,highcount)%>%distinct)
mRNA_TE_tbl%<>%safe_left_join(x=.,y=txi_counttbl_reps%>%filter(assay=='ribo')%>%distinct(gene_id,time,highcountribo = highcount))
mRNA_TE_tbl%<>%safe_left_join(x=.,y=txi_counttbl_reps%>%filter(assay=='total')%>%distinct(gene_id,time,highcountrna = highcount))
mRNA_TE_tbl$highcount = mRNA_TE_tbl$highcountribo & mRNA_TE_tbl$highcountrna
mRNA_TE_tbl%<>%safe_left_join(mRNAvals%>%select(gene_id,time,TPM=val))
mRNA_TE_tbl_cor<-mRNA_TE_tbl%>%filter(highcount)%>%group_by(time)%>%filter(is.finite(TPM+TE))%>%summarise(cor = cor(TPM,TE,use='complete'))
ribo_TE_tbl%<>%safe_left_join(x=.,y=txi_counttbl_reps%>%filter(assay=='ribo')%>%distinct(gene_id,time,highcountribo = highcount))
ribo_TE_tbl%<>%safe_left_join(x=.,y=txi_counttbl_reps%>%filter(assay=='total')%>%distinct(gene_id,time,highcountrna = highcount))
ribo_TE_tbl$highcount = ribo_TE_tbl$highcountribo & ribo_TE_tbl$highcountrna
ribo_TE_tbl%<>%safe_left_join(ribotpms)
ribo_TE_tbl_cor<-ribo_TE_tbl%>%filter(highcount)%>%group_by(time)%>%filter(is.finite(TPM+TE))%>%summarise(cor = cor(TPM,TE,use='complete'))

TErange = mRNA_TE_tbl%>%filter(highcount)%>%.$TE%>%2^.%>%range
TPMrange = mRNA_TE_tbl%>%filter(highcount)%>%.$TPM%>%2^.%>%range

'plots/QC_plots/fig1b_d_te_rna.pdf'%>%dirname%>%dir.create(rec=TRUE)
plotfile <- 'plots/QC_plots/fig1b_TE_rna.pdf'
cairo_pdf(plotfile,w=12,h=3)
ggplot(mRNA_TE_tbl%>%arrange(highcount),aes(x=2^TPM,y=2^TE,color=highcount))+
	scale_color_manual(values = setNames(c('black','grey'),c('TRUE','FALSE')))+
	geom_point(size=I(0.2))+
	facet_grid( ~ time)+
	scale_x_log10(name = 'Mean RNAseq TPM',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(name = 'Mean TE',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	geom_text(data=	mRNA_TE_tbl_cor%>%mutate(text = paste0('p = ',round(cor,3))),x=1,y=4,color=I('black'),aes(label=text))+
	coord_cartesian(xlim=TPMrange,ylim=TErange)+
	theme_bw()
dev.off()
plotfile%>%normalizePath%>%message

}

{




#
ribo_TE_tbl_cor<-ribo_TE_tbl%>%filter(highcount)%>%group_by(time)%>%summarise(cor = cor(TPM,TE))

TErange = ribo_TE_tbl%>%filter(highcount)%>%.$TE%>%2^.%>%range
TPMrange = ribo_TE_tbl%>%filter(highcount)%>%.$TPM%>%2^.%>%range

'plots/QC_plots/fig1b_d_te_ribo.pdf'%>%dirname%>%dir.create(rec=TRUE)
plotfile <- 'plots/QC_plots/fig1b_d_TE_ribo.pdf'
cairo_pdf(plotfile,w=12,h=3)
ggplot(ribo_TE_tbl%>%arrange(highcount),aes(x=2^TPM,y=2^TE,color=highcount))+
	scale_color_manual(values = setNames(c('black','grey'),c('TRUE','FALSE')))+
	geom_point(size=I(0.2))+
	facet_grid( ~ time)+
	scale_x_log10(name = 'Mean Riboseq TPM',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(name = 'Mean TE',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	geom_text(data=	ribo_TE_tbl_cor%>%mutate(text = paste0('p = ',round(cor,3))),x=1,y=4,color=I('black'),aes(label=text))+
	coord_cartesian(xlim=TPMrange,ylim=TErange)+
	theme_bw()
dev.off()
plotfile%>%normalizePath%>%message

}
# countexprdata %>% saveRDS(here('data/fig1countexprdata_w_high.rds'))

# tx_countdata$abundance%>%as.data.frame%>%
# 	rownames_to_column('gene_id')%>%
# 	gather(dataset,val,-gene_id)%>%
# 	separate(dataset,into=c('time','assay','replicate'))%>%
# 	spread(assay,val)%>%
# 	group_by(gene_id,time)%>%
# 	filter(highcount)%>%
# 	filter(gene_id=='ENSMUSG00000000003')


ribo_TE_tbl%>%filter(TE==Inf)

'plots/QC_plots/TE_outliers.pdf'%>%dirname%>%dir.create(rec=TRUE)
plotfile <- 'plots/QC_plots/TE_outliers.pdf'
cairo_pdf(plotfile,w=12,h=3)
ribo_TE_tbl = ribo_TE_tbl%>%arrange(highcount)%>%mutate(isTEoutlier = (2^TE)>50)
ribo_TE_tbl_cor<-ribo_TE_tbl%>%filter(highcount,!isTEoutlier)%>%group_by(time)%>%summarise(cor = cor(TPM,TE))
ggplot(ribo_TE_tbl,aes(x=2^TPM,y=2^TE,color=!isTEoutlier))+
	scale_color_manual(values = setNames(c('black','grey'),c('TRUE','FALSE')))+
	geom_point(size=I(0.2))+
	facet_grid( ~ time)+
	scale_x_log10(name = 'Mean Riboseq TPM',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(name = 'Mean TE',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	geom_text(data=	ribo_TE_tbl_cor%>%mutate(text = paste0('p = ',round(cor,3))),x=1,y=4,color=I('black'),aes(label=text))+
	theme_bw()
dev.off()
plotfile%>%normalizePath%>%message


source("src/R/Functions/go_term_funcs.R")
outliersumtbl = ribo_TE_tbl%>%group_by(gene_id)%>%summarise(isTEoutlier=any(isTEoutlier))
outliergidvect = outliersumtbl%>%.$isTEoutlier%>%as.numeric%>%as.factor%>%setNames(outliersumtbl$gene_id)
library(topGO)
minnodesize=20
my_stat = 'Fisher.elim'
for(myontology in c('BP','MF','CC')){
    go_data <- new("topGOdata",
                   ontology = myontology,
                   allGenes = outliergidvect,
                   nodeSize = minnodesize,
                   # annotationFun = annFUN.db
                   annotationFun = annFUN.org,
                   mapping = "org.Mm.eg",
                   ID = "ensembl"
                  )
    #
    results <- runTest(go_data, algorithm = "elim", statistic = "fisher")
    #
    results.tab <- GenTable(object = go_data, elimFisher = results,topNodes = 100)
    results.tab%<>%mutate(Enrichment = Significant / Expected )
    results.tab%<>%mutate(elimFisher = as.numeric(elimFisher) )
    results.tab%<>%dplyr::mutate(gene_ratio = Significant/Annotated)
    #
    results.tab%>%write_tsv(str_interp('tables/TEoutlier_${myontology}.'))

    goplotfile=str_interp('plots/QC_plots/TEoutlier_goplot_minsize${minnodesize}_Ontology_${myontology}_stat_${my_stat}.pdf')
    cairo_pdf(goplotfile,w=6)
    print(plot_go_enrich(results.tab,sort_var = 'elimFisher',goplotfile))
    dev.off()
    message(normalizePath(goplotfile))
    
}

results.tab<-read_tsv(str_interp('tables/TEoutlier_MF.pdf'))
results.tab%>%as.data.frame%>%arrange(elimFisher)%>%head
