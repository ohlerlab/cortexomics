

library(tidyverse)
library(data.table)
library(here)
library(magrittr)

pdf <- cairo_pdf
chekhits <- fread('ext_data/chekulaeva_etal_ebp1hits.tsv')
neurosomtbl <- fread('ext_data/neurites_zappulo_etal_2017.csv')

#lots of data in this table, fold changes an pvalues and some clusteirng I can't remember
#the point of. 
neurites <- here('ext_data/neurites_zappulo_etal_2017.csv')
neurites%<>%fread(skip=2)
ids = fread('pipeline/ids.txt')

#match by name or gene id

neuriteidmatch <- match(neurites$gene_id, ids$gene_id) %>% 
	ifelse(!is.na(.),.,match(neurites$gene_name,ids$gene_name))
stopifnot(mean(is.na(neuriteidmatch))<0.06)
neurites$gene_id <- ids$gene_id[neuriteidmatch]
neurites%<>%filter(!is.na(gene_id))

neurite_ribo_gids<-neurites%>%filter(RiboSeq_padj_Neurite_Soma<0.05)%>%.$gene_id
neurite_ribo_gids<-neurites%>%filter(RiboSeq_padj_Neurite_Soma<0.05)%>%.$gene_id
neurite_te_gids<-neurites%>%filter(`Prot log2FC >1, RNA log2FC <0`==1|`Prot log2FC >1, 1> RNA log2FC >0`==1)%>%.$gene_id
neurite_te_gids%>%colnames

neurites%>%colnames


rnaseqres <- neurites%>%
	mutate(neurite_soma_diff = rzRNA_padj_Neurite_Soma < 0.05)%>%
	select(neurite_soma_diff,gene_id,rzRNA_log2FC_Neurite_Soma,rzRNA_padj_Neurite_Soma)%>%
	left_join(chekhits%>%select(gene_id)%>%mutate(is_ebp1_hit=TRUE))



rnaseqresshort <- rnaseqres%>%	group_by(neurite_soma_diff)%>%
	summarise(is_ebp1_hit = sum(is_ebp1_hit %in% TRUE),total=n())%>%
	head(2)%T>%
	fisher.test(c(.$is_ebp1_hit,.$total))


library(ggrepel)

plotfile<-'plots/ebp1/checkscatter.pdf'

neurites$is_ebp1_hit <- neurites$gene_id %in% chekhits$gene_id
p <- ggplot(neurites,aes(y=RiboSeq_log2FC_Neurite_Soma,x=(RiboSeq_mean_Neurite+RiboSeq_mean_Soma)/2,label=gene_name,color=is_ebp1_hit))+
	geom_point(data=neurites%>%filter(!is_ebp1_hit))+
	# geom_point(data=neurites)+
	geom_point(data=neurites%>%filter(is_ebp1_hit))+
	scale_color_manual(guide=F,values=c('TRUE'='red','FALSE'='grey'))+
	ggrepel::geom_text_repel(data=neurites%>%filter(is_ebp1_hit))+
	geom_hline(yintercept=0,linetype='dashed')+
	theme_bw()

normalizePath(plotfile)%>%message

p$scales$scales[[1]]$range%>%str

plotfile2<-'plots/ebp1/checkscatter_rna.pdf'

neurites$is_ebp1_hit <- neurites$gene_id %in% chekhits$gene_id
p2<-ggplot(neurites,aes(y=rzRNA_log2FC_Neurite_Soma,x=(rzRNA_log2_mean_RPKM_Neurite+rzRNA_log2_mean_RPKM_Soma)/2,color=is_ebp1_hit,label=gene_name))+
	geom_point(data=neurites%>%filter(!is_ebp1_hit))+
	# geom_point(data=neurites)+
	geom_point(data=neurites%>%filter(is_ebp1_hit))+
	ggrepel::geom_text_repel(guide=F,data=neurites%>%filter(is_ebp1_hit))+
	scale_color_manual(guide=F,values=c('TRUE'='red','FALSE'='grey'))+
	geom_hline(yintercept=0,linetype='dashed')+
	theme_bw()
# p2<-change_plot_ylimits(p,ylims)


plotfile3<-'plots/ebp1/checkscatter_prot.pdf'

neurites$is_ebp1_hit <- neurites$gene_id %in% chekhits$gene_id
p3<-ggplot(neurites,aes(y=Proteomics_log2FC_Neurite_Soma,x=(Proteomics_log2_mean_LFQ_Neurite+Proteomics_log2_mean_LFQ_Soma)/2,color=is_ebp1_hit,label=gene_name))+
	geom_point(data=neurites%>%filter(!is_ebp1_hit))+
	# geom_point(data=neurites)+
	geom_point(data=neurites%>%filter(is_ebp1_hit))+
	ggrepel::geom_text_repel(guide=F,data=neurites%>%filter(is_ebp1_hit))+
	scale_color_manual(guide=F,values=c('TRUE'='red','FALSE'='grey'))+
	geom_hline(yintercept=0,linetype='dashed')+
	theme_bw()


breakint <- 0.5
commonyspanplots  <- function(trajectoryplots,breakint=0.5){
	trajectoryplots_ranges <- map(trajectoryplots,~ggplot_build(.)$layout$panel_scales_y[[1]]$range$range)
	trajectoryplots_range_centers <- trajectoryplots_ranges%>%map_dbl(mean)
	centeredranges <- map2(trajectoryplots_ranges,trajectoryplots_range_centers,`-`)
	centbreakrangemin<-centeredranges%>%map_dbl(1)%>%min%>%divide_by(breakint)%>%floor%>%multiply_by(breakint)
	centbreakrangemax<-centeredranges%>%map_dbl(2)%>%max%>%divide_by(breakint)%>%ceiling%>%multiply_by(breakint)
	trajrangeminsnap <- trajectoryplots_ranges%>%map_dbl(1)%>%map(divide_by,breakint)%>%map(floor)%>%map(multiply_by,breakint)
	trajrangemaxsnap <- trajrangeminsnap%>%map_dbl(add,centbreakrangemax-centbreakrangemin)
	for(i in seq_along(trajectoryplots)){
		trajectoryplots[[i]] = 	trajectoryplots[[i]] + coord_cartesian(ylim=c(trajrangeminsnap[[i]],trajrangemaxsnap[[i]]))
	}
	trajectoryplots
}


# normalizePath(plotfile2)%>%message

# c(p,p2) %<-% commonyspanplots(list(p,p2))

pdf(plotfile)
print(p+scale_y_continuous(limits = c(-6.5,6.5))+scale_x_continuous(limits=c(-5,25)))
dev.off()
pdf(plotfile2)
print(p2+scale_y_continuous(limits = c(-6.5,6.5))+scale_x_continuous(limits=c(-5,25)))
dev.off()
pdf(plotfile3)
print(p3)
dev.off()

normalizePath(plotfile3)


# neurites%>%
# 	mutate(sighit = RiboSeq_padj_Neurite_Soma < 0.05 & (abs(RiboSeq_log2FC_Neurite_Soma) > 0))%>%
# 	select(sighit,gene_id)%>%
# 	left_join(chekhits%>%select(gene_id)%>%mutate(is_ebp1_hit=TRUE))%>%
# 	group_by(sighit)%>%
# 	summarise(is_ebp1_hit = sum(is_ebp1_hit %in% TRUE),total=n())%>%
# 	head(2)%T>%print%>%
# 	fisher.test(c(.$is_ebp1_hit,.$total))

# neurites%>%
# 	mutate(sighit = (abs(QuaNCAT_log2FC_Neurite_Soma) > 0.25))%>%
# 	select(sighit,gene_id)%>%
# 	left_join(chekhits%>%select(gene_id)%>%mutate(is_ebp1_hit=TRUE))%>%
# 	group_by(sighit)%>%
# 	summarise(is_ebp1_hit = sum(is_ebp1_hit %in% TRUE),total=n())%>%
# 	head(2)%T>%print%>%
# 	fisher.test(c(.$is_ebp1_hit,.$total))


# # neurites%>%
# # 	mutate(sighit = (abs(pSILAC_log2FC_Neurite_Soma) > 0))%>%
# # 	select(sighit,gene_id)%>%
# # 	left_join(chekhits%>%select(gene_id)%>%mutate(is_ebp1_hit=TRUE))%>%
# # 	group_by(sighit)%>%
# # 	summarise(is_ebp1_hit = sum(is_ebp1_hit %in% TRUE),total=n())%T>%print%>%head(2)%>%
# # 	fisher.test(c(.$is_ebp1_hit,.$total))


neurites%>%
	mutate(sighit = Proteomics_padj_Neurite_Soma < 0.05 & ((Proteomics_log2FC_Neurite_Soma) > 0.5))%>%
	select(sighit,gene_id)%>%
	left_join(chekhits%>%select(gene_id)%>%mutate(is_ebp1_hit=TRUE))%>%
	group_by(sighit)%>%
	summarise(is_ebp1_hit = sum(is_ebp1_hit %in% TRUE),total=n())%>%
	head(2)%T>%print%>%
	fisher.test(c(.$is_ebp1_hit,.$total))

