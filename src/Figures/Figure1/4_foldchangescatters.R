timepoints <- c("E13", "E145", "E16", "E175", "P0")
################################################################################
########We want to show changes in translation, transcription and color them appropriately
################################################################################
gnm2gid <- ids_nrgname%>%distinct(gene_name,gene_id)%>%{setNames(.$gene_id,.$gene_name)}
multi_spread <- function(df, key, value) {
    # quote key
    keyq <- rlang::enquo(key)
    # break value vector into quotes
    valueq <- rlang::enquo(value)
    s <- rlang::quos(!!valueq)
    df %>% gather(variable, value, !!!s) %>%
        unite(temp, !!keyq, variable) %>%
        spread(temp, value)
}
fr_stepcountcontrdf <- readRDS(here('data/fr_stepcountcontrdf.rds'))
#
alllimmares$contrast%>%table
fr_stepcountcontrdf$time%>%table
stepcountcontrdf <- stepcountcontrdf
#
foldchangetbl<-bind_rows(
	stepcountcontrdf%>%
		filter(assay=='all')%>%
		mutate(changetype='transcriptional'),
	stepcountcontrdf%>%
		filter(assay=='TE')%>%
		mutate(changetype='translational')
)
foldchangetbl%<>%select(-assay)
foldchangetbl$gene_name <- gid2gnm[foldchangetbl$gene_id]
foldchangetbl$time%>%unique	
limmaobcols<-c("time", "logFC", "CI.L", "CI.R", "AveExpr",
"t", "P.Value", "adj.P.Val", "B", "gene_name", "gene_id",
"changetype")
stopifnot(all(foldchangetbl%>%colnames%>%setequal(limmaobcols)))
xtailcols <- c('time',"gene_id", "mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1",
"pvalue_v1", "E13_log2TE", "E145_log2TE", "log2FC_TE_v2", "pvalue_v2",
"logFC", "p_value", "adj.P.Val", "base_mean", "gene_name",
"E16_log2TE", "E175_log2TE", "P0_log2TE")
stopifnot(all(foldchangetbl%>%colnames%>%setequal(c("time", "logFC", "CI.L", "CI.R", "AveExpr",
"t", "P.Value", "adj.P.Val", "B", "gene_name", "gene_id",
"changetype"))))
xtailfoldchange<-Sys.glob('pipeline/xtail/xtail_*')%>%setNames(timepoints[-1])%>%
	map_df(.id='time',fread)%>%
	mutate(gene_id = gnm2gid[gene_name])%>%
	mutate(changetype='translational_xtail')
foldchangetblall<-bind_rows(
	foldchangetbl,
	xtailfoldchange%>%select(changetype,time,gene_id,logFC=log2fc,adj.P.Val=adj_p_value)
)

foldchangetblall_spread<-foldchangetblall%>%
	select(changetype,time,gene_id,logFC,adj.P.Val)%>%
	multi_spread(changetype,c(logFC,adj.P.Val))

foldchangetblall_spread%>%mutate(xtailsig=translational_xtail_adj.P.Val<0.05,lTEsig = translational_adj.P.Val<0.05)%>%
	group_by(xtailsig,lTEsig)%>%tally

inclusiontable(foldchangetblall_spread$gene_id,fData(countexprdata)$gene_id)
inclusiontable(foldchangetblall_spread$gene_id,highcountgenes)
inclusiontable(foldchangetblall_spread$gene_id,highcountgenes)

'plots/figures/figure1/foldchangecomp_te_xtail_limma.pdf'%>%dirname%>%dir.create(rec=TRUE)

plotfile <- 'plots/figures/figure1/foldchangecomp.pdf'
cairo_pdf(plotfile,w=24,h=6)
foldchangetblall_spread%>%
	mutate(lTEsig = translational_adj.P.Val<0.05)%>%
	arrange(lTEsig)%>%
	ggplot(aes(x=translational_logFC,y=translational_xtail_logFC,color=lTEsig))+geom_point(size=1)+theme_bw()+
	ggplot(aes(x=translational_logFC,y=translational_xtail_logFC,color=lTEsig))+geom_point(size=1)+theme_bw()+
	facet_grid( ~ time,scale='free')+
	scale_y_continuous(limits=c(-5,5))+
	scale_x_continuous(limits=c(-5,5))	
foldchangetblall%>%filter(changetype=='translational_xtail')%>%head
dev.off()
plotfile%>%normalizePath%>%message
foldchangetblall_spread%>%colnames

FCchangecols <- c('No Sig Change' = '#F2F2F2',
	'Transcriptional Only'= '#6E6E6E',
'Translational Only' = '#04B404',
'Concurrent Change' = '#B4045F',
'Compensating Change' = '#FFF800')

foldchangecatdf <- foldchangetblall_spread%>%
	mutate(lTEsig = translational_adj.P.Val<0.05)%>%
	mutate(lTRsig = transcriptional_adj.P.Val<0.05)%>%
	mutate(xtailsig = translational_adj.P.Val<0.05)%>%
	mutate(sigstatus = case_when(
		(lTRsig)  & (lTEsig) & (sign(translational_logFC)==sign(transcriptional_logFC)) ~ 'Concurrent Change',
		(lTRsig) & (lTEsig) & (sign(translational_logFC)!=sign(transcriptional_logFC)) ~ 'Compensating Change',
		(lTRsig) & (!lTEsig) ~ 'Transcriptional Only',
		(!lTRsig) & (lTEsig) ~ 'Translational Only',
		TRUE ~ 'No Sig Change'
	 ))%>%
	{assert_that(all(.$sigstatus%in%names(FCchangecols)));.}%>%
	# )%>%filter(lTEsig)%>%filter(lTRsig)%>%group_by(sign(translational_logFC),sign(translational_logFC))
	arrange(match(sigstatus,names(FCchangecols)))%>%
	identity

foldchangecatdf%>%group_by(gene_id)%>%summarise(sig=any(xtailsig))%>%group_by(sig)%>%tally

plotfile <- 'plots/figures/figure1/foldchangecomp_limma.pdf'
plotfile%>%dirname%>%dir.create(rec=TRUE)
cairo_pdf(plotfile,w=24*.6,h=6*.6)
foldchangecatdf%>%	ggplot(aes(x=transcriptional_logFC,y=translational_logFC,color=sigstatus))+
	geom_point(size=0.5)+theme_bw()+
	facet_grid( ~ time,scale='free')+
	scale_y_continuous(limits=c(-5,5))+
	scale_x_continuous(limits=c(-5,5))+
	scale_color_manual(values = FCchangecols)+
	guides(colour = guide_legend(override.aes = list(size=10)))
dev.off()
plotfile%>%normalizePath%>%message

stop()


catlevels=c("Translational Only", "Transcriptional Only",
"Concurrent Change", "Compensating Change","No Sig Change")

plotfile <- 'plots/figures/figure1/foldchange_numbers.pdf'
plotfile%>%dirname%>%dir.create(rec=TRUE)
cairo_pdf(plotfile,w=24*.6,h=6*.6)
foldchangecatdf%>%
	mutate(sigstatus = sigstatus%>%as_factor%>%fct_relevel(catlevels))%>%
	group_by(time,sigstatus)%>%
	filter(!sigstatus%in%c('No Sig Change'))%>%tally%>%ggplot(aes(x=time,y=n,fill=sigstatus))+
	# filter(!sigstatus%in%c('No Sig Change','Transcriptional Only'))%>%tally%>%ggplot(aes(x=sigstatus,y=n,fill=sigstatus))+
	stat_identity(geom='bar')+
	theme_bw()+
	facet_grid( ~ sigstatus)+
	scale_fill_manual(values = FCchangecols)+
	theme(axis.text.x = element_text(angle=45))+
	guides(fill = guide_legend(override.aes = list(size=10)))
dev.off()
plotfile%>%normalizePath%>%message

foldchangecatdf%>%saveRDS(here('data/foldchangecatdf.rds'))






