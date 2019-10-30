

################################################################################
########We want to show changes in translation, transcription and color them appropriately
################################################################################

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

bestmscountebayes

techangetbl  <- lapply(tps[-1],function(testtp){
	eBayes(contrasts.fit(lmFit(mscountvoom[best_uprotein_ids,]),contrasts = timeTEeffect[,2:5]))%>%topTable(coef=which(tps==testtp)-1,number=1e9,confint=0.95)%>%rownames_to_column('uprotein_id')%>%
	safe_left_join(ms_id2protein_id%>%distinct(ms_id,gene_name,gene_id,uprotein_id))
})

allchangetbl <- lapply(tps[-1],function(testtp){
	eBayes(contrasts.fit(lmFit(mscountvoom[best_uprotein_ids,]),contrasts = alltimeeff[,2:5]))%>%topTable(coef=which(tps==testtp)-1,number=1e9,confint=0.95)%>%rownames_to_column('uprotein_id')%>%
	safe_left_join(ms_id2protein_id%>%distinct(ms_id,gene_name,gene_id,uprotein_id))
})

foldchangetbl<-bind_rows(
	allchangetbl%>%setNames(timepoints[-1])%>%bind_rows(.id='time')%>%mutate(changetype='transcriptional'),
	techangetbl%>%setNames(timepoints[-1])%>%bind_rows(.id='time')%>%mutate(changetype='translational')
)

stopifnot(all(foldchangetbl%>%colnames%>%`==`(c("time", "uprotein_id", "logFC", "CI.L", "CI.R", "AveExpr",
"t", "P.Value", "adj.P.Val", "B", "ms_id", "gene_name", "gene_id",
"changetype"))))

xtailcols <- c('time',"gene_id", "mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1",
"pvalue_v1", "E13_log2TE", "E145_log2TE", "log2FC_TE_v2", "pvalue_v2",
"logFC", "p_value", "adj.P.Val", "base_mean", "gene_name",
"E16_log2TE", "E175_log2TE", "P0_log2TE")

xtailfoldchange<-Sys.glob('pipeline/xtail/xtail_*')%>%setNames(timepoints[-1])%>%map_df(.id='time',fread)%>%set_colnames(xtailcols)%>%
	select(one_of(colnames(foldchangetbl)))%>%mutate(changetype='translational_xtail')
foldchangetblall<-bind_rows(foldchangetbl,xtailfoldchange)

foldchangetblall_spread<-foldchangetblall%>%select(changetype,time,gene_id,logFC,adj.P.Val)%>%
	multi_spread(changetype,c(logFC,adj.P.Val))

foldchangetblall_spread%>%mutate(xtailsig=translational_xtail_adj.P.Val<0.05,lTEsig = translational_adj.P.Val<0.05)%>%
	group_by(xtailsig,lTEsig)%>%tally

'plots/figures/figure1/foldchangecomp_te_xtail_limma.pdf'%>%dirname%>%dir.create(rec=TRUE)
plotfile <- 'plots/figures/figure1/foldchangecomp.pdf'
pdf(plotfile,w=24,h=6)
foldchangetblall_spread%>%
	mutate(lTEsig = translational_adj.P.Val<0.05)%>%
	arrange(lTEsig)%>%
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
	mutate(xtailsig = translational_xtail_adj.P.Val<0.05)%>%
	mutate(sigstatus = case_when(
		(lTRsig)  & (xtailsig) & (sign(translational_xtail_logFC)==sign(transcriptional_logFC)) ~ 'Concurrent Change',
		(lTRsig) & (xtailsig) & (sign(translational_xtail_logFC)!=sign(transcriptional_logFC)) ~ 'Compensating Change',
		(lTRsig) & (!xtailsig) ~ 'Transcriptional Only',
		(!lTRsig) & (xtailsig) ~ 'Translational Only',
		TRUE ~ 'No Sig Change'
	 ))%>%
	{assert_that(all(.$sigstatus%in%names(FCchangecols)));.}%>%
	# )%>%filter(xtailsig)%>%filter(lTRsig)%>%group_by(sign(translational_logFC),sign(translational_logFC))
	arrange(match(sigstatus,names(FCchangecols)))%>%
	identity

foldchangecatdf%>%group_by(gene_id)%>%summarise(sig=any(xtailsig))%>%group_by(sig)%>%tally

plotfile <- 'plots/figures/figure1/foldchangecomp_limma.pdf'
plotfile%>%dirname%>%dir.create(rec=TRUE)
pdf(plotfile,w=24*.6,h=6*.6)
foldchangecatdf%>%	ggplot(aes(x=transcriptional_logFC,y=translational_xtail_logFC,color=sigstatus))+
	geom_point(size=0.5)+theme_bw()+
	facet_grid( ~ time,scale='free')+
	scale_y_continuous(limits=c(-5,5))+
	scale_x_continuous(limits=c(-5,5))+
	scale_color_manual(values = FCchangecols)+
	guides(colour = guide_legend(override.aes = list(size=10)))
dev.off()
plotfile%>%normalizePath%>%message


plotfile <- 'plots/figures/figure1/foldchange_numbers.pdf'
plotfile%>%dirname%>%dir.create(rec=TRUE)
pdf(plotfile,w=24*.6,h=6*.6)
foldchangecatdf%>%group_by(time,sigstatus)%>%
	# filter(!sigstatus%in%c('No Sig Change','Transcriptional Only'))%>%tally%>%ggplot(aes(x=sigstatus,y=n,fill=sigstatus))+
	filter(!sigstatus%in%c('No Sig Change'))%>%tally%>%ggplot(aes(x=time,y=n,fill=sigstatus))+
	stat_identity(geom='bar')+
	theme_bw()+
	facet_grid( ~ sigstatus)+
	scale_fill_manual(values = FCchangecols)+
	theme(axis.text.x = element_text(angle=45))+
	guides(fill = guide_legend(override.aes = list(size=10)))
	# coord_flip()
dev.off()
plotfile%>%normalizePath%>%message



foldchangetblall_spread%>%filter(translational_xtail_adj.P.Val<0.05)

foldchangecatdf%>%head









