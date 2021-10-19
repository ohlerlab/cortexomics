################################################################################
########
################################################################################
	
source('src/Figures/load_annotation.R')


ctobs <- get(load('data/1_integrate_countdata.R'))


#first get the gc content for every transcripts 3' UTR
fputrseq = fputrs%>%extractTranscriptSeqs(x=fafileob)
fputrgc = (fputrseq%>%str_count('G|C')) / (fputrseq%>%nchar)
fputrgc %<>% setNames(names(fputrseq))

#Then weight that by abundance to get mean gc content per gene
abtxs = iso_tx_countdata$abundance%>%rownames

fputrgcdf = iso_tx_countdata$abundance%>%as.data.frame%>%
	rownames_to_column('tr_id')%>%
	gather(sample, tpm, -tr_id)%>%
	left_join(enframe(fputrgc, 'tr_id', 'gc'))%>%
	left_join(ids_nrgname%>%select(tr_id=transcript_id, gene_id=gene_id))%>%
	filter(!is.na(gc))%>%
	filter(gene_id%in%highcountgenes)%>%
	group_by(gene_id, sample)%>%summarise(gc=mean(gc*(tpm/sum(tpm)), na.rm=T))%>%
	select(gene_id,  sample,  gc)

fpgcchange = fputrgcdf%>%separate(sample, c('time', 'assay', 'rep'))%>%
	filter(assay=='total')%>%
	group_by(gene_id, time)%>%summarise_at(vars(gc), list(mean))%>%
	group_by(gene_id)%>%
	summarise(gcchange = log2(gc[time=='E175']/gc[time=='E13']))%>%
	select(gene_id, gcchange)


#Now correlate this with TE changes.
techange = readxl::read_xlsx(  "tables/S2.xlsx", 1, col_types=c(time='text'))
techange$log2fc%<>%as.numeric

gc_te_changedf <- fpgcchange%>%
	left_join(techange%>%filter(time==3)%>%select(gene_id, log2fc, adj_p_value))%>%
	select(gene_id, gcchange, log2fc, adj_p_value)

gc_te_changedf$gcchange%>%abs%>%.[is.finite(.)]%>%.[.!=0]%>%cut_number(5)%>%table




gc_te_changedf%>%
	filter(abs(gcchange)!=0)%>%
	{quicktest(.$gcchange, .$log2fc)}

gc_te_changedf%>%
	filter(abs(gcchange)!=0)%>%
	filter(adj_p_value<0.05)%>%
	{cor.test(.$gcchange, .$log2fc, method='spearman')}

gc_te_changedf%>%
	filter(abs(gcchange)>0.1)%>%
	filter(adj_p_value<0.05)%>%
	{quicktest(.$gcchange, .$log2fc)}


#now plot
plotfile<- here(paste0('plots/','te_fputrgc_changescat','.pdf'))
pdf(plotfile)
gc_te_changedf%>%
	filter(abs(gcchange)>0.1)%>%
	ggplot(.,aes(x=gcchange,y=log2fc))+
	geom_point()+
	scale_x_continuous(paste0('xname'))+
	scale_y_continuous(paste0('yname'))+
	ggtitle(paste0('title'))+
	theme_bw()
dev.off()
message(normalizePath(plotfile))

##############################################################################
########5' TOPs
##############################################################################
{
fptopgr = rtracklayer::import('pipeline/fpTOP_scan/fpTOP_scan/fimo.gff')
fptopgr%<>%subset(score>43)
fptopgr%<>%subset(str_detect(seqnames,'\\(\\+\\)'))
fptopgr$tr_id = fptopgr@seqnames%>%str_extract('\\w+')
gid2trid = ids_nrgname%>%distinct(gene_id,transcript_id)%>%
	{setNames(.[[1]],.[[2]])}
fptopgr$g_id = gid2trid[fptopgr$tr_id]%>%setNames(NULL)
fptopgr$tr_id%>%setdiff(ids_nrgname$transcript_id)
fptopgr%<>%subset(!is.na(g_id))
fptopgr%<>%as.data.frame%>%
	select(tr_id, gene_id=g_id)
#
te_top_df = techange%>%filter(time==3)%>%left_join(fptopgr)%>%
	mutate(has_top = gene_id %in% fptopgr$gene_id)%>%
	select(has_top, gene_id, log2fc, adj_p_value)

topggdf = te_top_df%>%	
	# filter(adj_p_value<0.05)%>%
	filter(gene_id %in% highcountgenes)

topggdf%>%{split(.$log2fc,.$has_top)}%>%{t.test(.[[1]],.[[2]])}

testlbl = topggdf%>%{split(.$log2fc,.$has_top)}%>%{wilcox.test(.[[1]],.[[2]])}%>%
	tidy%>%{str_interp('wilox test on all logfc, p = ${round(.$p.value, 5)}')}

topggdf%>%filter(adj_p_value<0.05)%>%mutate(negchange=log2fc<0)%>%
	{fisher.test(.$has_top,.$negchange)}



#now plot
plotfile<- here(paste0('plots/','tophist','.pdf'))
pdf(plotfile)
print(topggdf%>%
	# filter(adj_p_value<0.05)%>%
	mutate(log2fc = pmax(log2fc,-3))%>%
	mutate(log2fc = pmin(log2fc,3))%>%
	{
		# ggplot(filter(.,adj_p_value<0.05),aes(fill=has_top,x=log2fc))+
		ggplot(filter(.),aes(fill=has_top,x=log2fc))+
	geom_histogram()+
	# geom_histogram(data=.,alpha=I(0.5))+
	facet_grid(has_top~.,scale='free')+
	scale_x_continuous(paste0('delta-TE'))+
	geom_text(show.legend=F,data=tibble(has_top=TRUE,labl=testlbl),
		hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))+
	# coord_cartesian(ylim=c(0,1000))+
	# scale_y_continuous(paste0('yname'))+
	ggtitle(paste0('TE Fold Change Distribution(E175/E13)\n5\' TOP genes vs others'))+
	theme_bw()})
dev.off()
message(normalizePath(plotfile))



#now plot
plotfile<- here(paste0('plots/','topdens','.pdf'))
pdf(plotfile)
print(topggdf%>%
	mutate(log2fc = pmax(log2fc,-3))%>%
	mutate(log2fc = pmin(log2fc,3))%>%
	ggplot(.,aes(color=has_top,x=log2fc))+
	geom_density()+
	# facet_grid(has_top~.,scale='free')+
	# scale_x_continuous(paste0('xname'))+
	# scale_y_continuous(paste0('yname'))+
	ggtitle(paste0('title'))+
	theme_bw())
dev.off()
message(normalizePath(plotfile))

}
