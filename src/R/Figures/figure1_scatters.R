library(scales)
library(dplyr)




library(conflicted)
conflict_prefer("intersect", "dplyr")
conflict_prefer("exprs", "Biobase")
#A4 size is 8 by 11

HIGHCOUNTTHRESH <- 32

#THis statement from dunn and weismann 2013:
##A4 size is 8 by 11Remarkably, we find that the range of translation efficiencies for different messages spans four orders of magnitude, a range comparable to that observed for mRNA abundance of well-counted genes (Figure 1D). Moreover, translation efficiency is uncorrelated with mRNA abundance (r2 = 8.29 × 10−5; Figure 1E) and mRNA abundance predicts only one third of the variance in the rate of protein production as measured by ribosome footprint density (Figure 1F). Translational regulation is therefore a major determinant of gene expression in the early embryo (supplementary table 1 at Dryad: Dunn et al., 2013), and ribosome profiling provides a quantitative and robust means to monitor translational regulation during development.

stagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
stageconv = names(stagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
#Get the counts at each timepoint.

# countexprdata <- readRDS(countfile)
# countfile='pipeline/exprdata/countexprset.rds'

TEs = Sys.glob('pipeline/xtail/*.txt')%>%setNames(.,basename(.))%>%map_df(.%>%fread%>%select(feature_id,dplyr::matches('log2TE')))

TEs = bind_rows(
	TEs%>%select(feature_id,E13_log2TE)%>%group_by(feature_id)%>%summarise(val=mean(E13_log2TE))%>%mutate(time = 'E13'),
	TEs%>%select(-E13_log2TE)%>%gather(time,val,-feature_id)%>%mutate(time = str_extract(time,'[^_]+'))
)

TEs%<>%inner_join(.,
	as.data.frame(fData(countexprdata))%>%filter(is_gid_highest)%>%select(feature_id=gene_id,gene_name,protein_id)%>%distinct
)


TEs%<>%select(-feature_id)
TEs%<>%select(feature_id = gene_name,everything())

# TEs <- itimecountebayes$coef%>%as.data.frame%>%select(dplyr::matches('riboTRUE'))%>%
# 	gather(time,val)%>%
# 	mutate(time = str_replace(time,'riboTRUE:?',''))%>%
# 	mutate(time = ifelse(time=='','timeE13',time))%>%
# 	mutate(time = str_replace(time,'time',''))


TEs %<>%mutate(type='TE')

# countexprdata%>%saveRDS(here('data/fig1countexprdata.rds'))
countexprdata <- readRDS(here('data/fig1countexprdata.rds'))

fData(countexprdata)$gene_id%>%n_distinct

mscountrows <- (fData(countexprdata)$protein_id %in% ms_id2protein_id$protein_id) & fData(countexprdata)$is_gid_highest

conflict_prefer('rowMedians','Biobase')

exprs(countexprdata)[mscountrows,'E13_ribo_1',drop=F]%>%rowMedians%>%add(1)%>%log2%>%txtdensity


mRNAvals = exprs(countexprdata)[fData(countexprdata)$is_gid_highest,]
mRNAvals <- countvoom$E[fData(countexprdata)$protein_id,][fData(countexprdata)$protein_id[fData(countexprdata)$is_gid_highest],]%>%as.data.frame%>%
	rownames_to_column('protein_id')%>%
	select(protein_id,dplyr::matches('total'))%>%
	gather(dataset,val,-protein_id)%>%
	mutate(time = str_extract(dataset,'[^_]+'),type='mRNA')

counttbl_reps<-exprs(countexprdata)[fData(countexprdata)$is_gid_highest,]%>%
	as.data.frame%>%
	rownames_to_column('feature_id')%>%
	gather(dataset,count,-feature_id)%>%
	separate(dataset,c('time','assay','replicate'))%>%
	mutate(replicate = paste0('Replicate_',letters[as.numeric(replicate)]))%>%
	spread(replicate,count)%>%
	group_by(feature_id)

counttbl_reps$highcount = HIGHCOUNTTHRESH <= counttbl_reps$Replicate_a+counttbl_reps$Replicate_b

tetbl_reps <- countvoom$E[fData(countexprdata)$protein_id[fData(countexprdata)$is_gid_highest],]%>%as.data.frame%>%
	rownames_to_column('protein_id')%>%
	gather(dataset,val,-protein_id)%>%
	mutate(time = str_extract(dataset,'[^_]+'),type='mRNA')%>%
	separate(dataset,into=c('time','assay','replicate'))%>%
	spread(assay,val)%>%
	mutate(te = ribo - total)%>%
	select(protein_id,time,replicate,te)%>%
	mutate(replicate = ifelse(replicate==1,'Replicate_a','Replicate_b'))%>%
	spread(replicate,te)

tetbl_reps%<>%safe_left_join(x=.,y=counttbl_reps%>%filter(assay=='total')%>%distinct(protein_id=feature_id,time,highcount))


# tetbl_reps<-data.frame(
# 	counttbl_reps[counttbl_reps$assay=='ribo',]%>%select(-dplyr::matches('Replic')),
# 	Replicate_a = counttbl_reps[counttbl_reps$assay=='ribo',]$Replicate_a / counttbl_reps[counttbl_reps$assay=='total',]$Replicate_a,
# 	Replicate_b =counttbl_reps[counttbl_reps$assay=='ribo',]$Replicate_b / counttbl_reps[counttbl_reps$assay=='total',]$Replicate_b
# )

repcors = counttbl_reps%>%group_by(time,assay)%>%filter(highcount)%>%summarise(cor=cor(use='complete',Replicate_a,Replicate_b))
terepcors = tetbl_reps%>%group_by(time)%>%filter(highcount)%>%filter(is.finite(Replicate_a))%>%
	filter(is.finite(Replicate_b))%>%summarise(cor=cor(use='complete',Replicate_a,Replicate_b))

plotfile <- 'plots/figures/figure1/fig1b_rpfrpf.pdf'
cairo_pdf(plotfile,w=12,h=3)
ggplot(counttbl_reps%>%filter(assay=='ribo')%>%arrange(highcount),aes(x=(1+Replicate_a),y=(1+Replicate_b),color=highcount))+
	scale_color_manual(values = setNames(c('black','grey'),c('TRUE','FALSE')))+
	# scale_color_manual(values = stagecols)+
	# scale_fill_manual(values=stagecols)+
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

plotfile <- 'plots/figures/figure1/fig1b_b_rnarna.pdf'
cairo_pdf(plotfile,w=12,h=3)
ggplot(counttbl_reps%>%filter(assay=='total')%>%arrange(highcount),aes(x=(1+Replicate_a),y=(1+Replicate_b),color=highcount))+
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


repavexprtbl<-counttbl_reps%>%group_by(feature_id,time)%>%summarise(highcount = any(highcount[assay=='ribo']),mRNA=(Replicate_b+Replicate_a)[assay=='total']/2,ribo=(Replicate_b+Replicate_a)[assay=='ribo']/2)


mRNA_ribo_cortbl<-repavexprtbl%>%filter(highcount)%>%group_by(time)%>%summarise(cor = cor(use='complete',log2(1+mRNA),log2(1+ribo)))


'plots/figures/figure1/fig1b_b_rna_rpf.pdf'%>%dirname%>%dir.create(rec=TRUE)
plotfile <- 'plots/figures/figure1/fig1b_b_rna_rpf.pdf'
pdf(plotfile,w=12,h=3)
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





plotfile <- 'plots/figures/figure1/fig1b_b_tete.pdf'
pdf(plotfile,w=12,h=3)
ggplot(tetbl_reps%>%arrange(highcount)%>%arrange(highcount),aes(x=(2^Replicate_a),y=(2^Replicate_b),color=highcount))+
	scale_color_manual(values = setNames(c('black','grey'),c('TRUE','FALSE')))+
	# scale_color_manual(values = stagecols)+
	# scale_fill_manual(values=stagecols)+
	geom_point(size=I(0.2))+
	facet_grid( ~ time)+
	scale_x_log10(name = 'TE Replicate 1',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(name = 'TE Replicate 2',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	geom_text(data=terepcors%>%mutate(text = paste0('p = ',round(cor,3))),x=1,y=4,color=I('black'),aes(label=text))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))


ggdf = bind_rows(mRNAvals,TEs)


#add high count col
ggdf%<>%safe_left_join(counttbl_reps%>%filter(assay=='total')%>%select(protein_id=feature_id,time,highcount)%>%distinct)


ggdf %<>% mutate(time = stageconv[time])

library(scales)
l10scales =   scale_x_log10(name = 'TE / CPM',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))
ggdf%>%filter(type=='TE',highcount)%>%.$time%>%table

plotfile <- 'plots/figures/figure1/fig1b_c_TE_rna_dists.pdf'
pdf(plotfile,w=12,h=3)
ggplot(ggdf%>%filter(type=='mRNA',highcount)%>%group_by(time)%>%mutate(val=val-median(val))%>%arrange(highcount),aes(x=2^val))+
	geom_density(aes(fill=time),alpha=I(0.5))+
	geom_density(data=ggdf%>%filter(type=='TE',highcount)%>%group_by(time)%>%mutate(val=val-median(na.omit(val))),alpha=I(1),aes(color=time))+facet_grid(~time)+
	scale_color_manual(values = stagecols)+
	scale_fill_manual(values = stagecols)+
	# scale_x_log10()+
	l10scales+
	theme_bw()
dev.off()
message(normalizePath('tmp.pdf'))
# 'plots/figures/fig1.pdf'%>%dirname%>%dir.create
file.copy('tmp.pdf',plotfile,overwrite=T)
plotfile%>%normalizePath%>%message


mRNA_TE_tbl <- countvoom$E[fData(countexprdata)$protein_id[fData(countexprdata)$is_gid_highest],]%>%as.data.frame%>%
	rownames_to_column('protein_id')%>%
	gather(dataset,val,-protein_id)%>%
	mutate(time = str_extract(dataset,'[^_]+'),type='mRNA')%>%
	separate(dataset,into=c('time','assay','replicate'))%>%
	spread(assay,val)%>%
	group_by(protein_id,time)%>%
	summarise(mRNA = mean(total), TE = mean(ribo) -  mean(total) )

mRNA_TE_tbl%<>%safe_left_join(x=.,y=counttbl_reps%>%filter(assay=='total')%>%select(protein_id=feature_id,time,highcount)%>%distinct)
mRNA_TE_tbl%<>%safe_left_join(x=.,y=counttbl_reps%>%filter(assay=='ribo')%>%distinct(protein_id=feature_id,time,highcountribo = highcount))
# mRNA_TE_tbl$highcount = mRNA_TE_tbl$highcount & mRNA_TE_tbl$highcountribo
mRNA_TE_tbl$highcount =  mRNA_TE_tbl$highcountribo
mRNA_TE_tbl%<>%safe_left_join(cds%>%split(.$protein_id)%>%width%>%sum%>%enframe('protein_id','length'),by='protein_id')
mRNA_TE_tbl%<>%mutate(RPKM = mRNA - log2(length/1e3))


mRNA_TE_tbl

mRNA_TE_tbl_cor<-mRNA_TE_tbl%>%filter(highcount)%>%group_by(time)%>%summarise(cor = cor(RPKM,TE))

'plots/figures/figure1/fig1b_d_te_rna.pdf'%>%dirname%>%dir.create(rec=TRUE)
plotfile <- 'plots/figures/figure1/fig1b_TE_rna.pdf'
pdf(plotfile,w=12,h=3)
ggplot(mRNA_TE_tbl%>%arrange(highcount),aes(x=2^RPKM,y=2^TE,color=highcount))+
	scale_color_manual(values = setNames(c('black','grey'),c('TRUE','FALSE')))+
	geom_point(size=I(0.2))+
	facet_grid( ~ time)+
	scale_x_log10(name = 'Mean RNAseq RPKM',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(name = 'Mean TE',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	geom_text(data=	mRNA_TE_tbl_cor%>%mutate(text = paste0('p = ',round(cor,3))),x=1,y=4,color=I('black'),aes(label=text))+
	theme_bw()
dev.off()
plotfile%>%normalizePath%>%message

ribo_TE_tbl <- countvoom$E[fData(countexprdata)$protein_id[fData(countexprdata)$is_gid_highest],]%>%as.data.frame%>%
	rownames_to_column('protein_id')%>%
	gather(dataset,val,-protein_id)%>%
	mutate(time = str_extract(dataset,'[^_]+'),type='Ribo')%>%
	separate(dataset,into=c('time','assay','replicate'))%>%
	spread(assay,val)%>%
	group_by(protein_id,time)%>%
	summarise(Ribo = mean(ribo), TE = mean(ribo) -  mean(total) )

ribo_TE_tbl%<>%safe_left_join(x=.,y=counttbl_reps%>%filter(assay=='total')%>%select(protein_id=feature_id,time,highcount)%>%distinct)
ribo_TE_tbl%<>%safe_left_join(x=.,y=counttbl_reps%>%filter(assay=='ribo')%>%distinct(protein_id=feature_id,time,highcountribo = highcount))
ribo_TE_tbl$highcount = ribo_TE_tbl$highcountribo
ribo_TE_tbl%<>%safe_left_join(cds%>%split(.$protein_id)%>%width%>%sum%>%enframe('protein_id','length'),by='protein_id')
ribo_TE_tbl%<>%mutate(RPKM = Ribo - log2(length/1e3))

ribo_TE_tbl_cor<-ribo_TE_tbl%>%filter(highcount)%>%group_by(time)%>%summarise(cor = cor(RPKM,TE))

'plots/figures/figure1/fig1b_d_te_ribo.pdf'%>%dirname%>%dir.create(rec=TRUE)
plotfile <- 'plots/figures/figure1/fig1b_d_TE_ribo.pdf'
pdf(plotfile,w=12,h=3)
ggplot(ribo_TE_tbl%>%arrange(highcount),aes(x=2^RPKM,y=2^TE,color=highcount))+
	scale_color_manual(values = setNames(c('black','grey'),c('TRUE','FALSE')))+
	geom_point(size=I(0.2))+
	facet_grid( ~ time)+
	scale_x_log10(name = 'Mean Riboseq RPKM',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	scale_y_log10(name = 'Mean TE',breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
	geom_text(data=	ribo_TE_tbl_cor%>%mutate(text = paste0('p = ',round(cor,3))),x=1,y=4,color=I('black'),aes(label=text))+
	theme_bw()
dev.off()
plotfile%>%normalizePath%>%message

fData(countexprdata)$highcount <- ribo_TE_tbl$highcount[fData(countexprdata)%>%.$protein_id%>%match(ribo_TE_tbl$protein_id)]

# countexprdata %>% saveRDS(here('data/fig1countexprdata_w_high.rds'))


# save.image('/fast/groups/ag_ohler/work/dharnet_m/cortexomics/data/figure1_scatters.Rdata')
# load('/fast/groups/ag_ohler/work/dharnet_m/cortexomics/data/figure1_scatters.Rdata')