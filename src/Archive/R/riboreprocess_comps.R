################################################################################
########
################################################################################
library(GenomicAlignments)
library(GenomicFeatures)

get_cds_vects<-function(bamfile,cdsgrl,offsets,expn=30){
    #read psites from a bam file using our offsets
    psites = get_genomic_psites(bamfile,cdsgrl%>%unlist,offsets,mapqthresh=20,comps=c('chrM'='chrM'))
    #now convert to a cocverage vector
    psitecovrles = mapToTranscripts(psites,cdsgrl)%>%coverage
}

cdsgrl <- cds%>%split(.,.$protein_id)
cdsexonsgrl <- exons%>%split(.$transcript_id)%>%.[(fmcols(cdsgrl,transcript_id))]
names(cdsexonsgrl)=names(cdsgrl)
expcdsgrl <- get_exp_cds(cdsgrl,cdsexonsgrl,0)
expcdsgrl <- expcdsgrl%>%split(.,names(.))
expcdsgrl <- expcdsgrl%>%dropSeqlevels('chrM',pruning='coarse')

if(!file.exists('data/psitecovrles.rds')) {	
	psitecovrles = Sys.glob('pipeline/star/data/*/*ribo*.bam')%>%
		setNames(.,basename(dirname(.)))%>%
		lapply(FUN=get_cds_vects,cdsgrl=expcdsgrl[names(bestcds)],offsets=offsets)

	psitecovrles%>%saveRDS('data/psitecovrles.rds')

}else{
	psitecovrles<-readRDS('data/psitecovrles.rds')	
}
if(!file.exists(here('data/spec_vals.rds'))){
	spec_vals = psitecovrles%>%mclapply(.%>%lapply(.%>%as.vector%>%ftestvect)%>%map(1))
	saveRDS(spec_vals,here('data/spec_vals.rds'))
}else{
	spec_vals<-readRDS(here('data/spec_vals.rds'))
}
cdslens = expcdsgrl%>%width%>%sum%>%.[names(bestcds)]

#now use our psite coverage vector to get occupancies for each codon
spec_mags = spec_vals%>%lapply(.%>%flatten_dbl%>%multiply_by(cdslens))

spec_mag_df = spec_mags%>%map_df(.id='sample',enframe,'protein_id','spec_mag')

df2matrix <- function(df){
	set_rownames(as.matrix(df[,-1]),df[[1]])
}

ribodesign=spec_mag_df%>%distinct(sample)%>%separate(sample,c('time','assay','rep'),remove=F)

spec_voom <- spec_mag_df%>%
	mutate(spec_mag = replace_na(spec_mag,0))%>%
	spread(sample,spec_mag)%>%
	df2matrix%>%
 	{limma::voom(.,design= model.matrix(~time,ribodesign))}

count_voom <- psitecovrles%>%lapply(sum)%>%map_df(.id='sample',enframe,'protein_id','count')%>%
	mutate(count = replace_na(count,0))%>%
	spread(sample,count)%>%
	df2matrix%>%
 	{limma::voom(.,
 		design= model.matrix(~time,
 		data.frame(sample=colnames(.))%>%separate(sample,c('time','assay','rep'),remove=F)
 	))}

spec_fcs = eBayes(contrasts.fit(lmFit(count_voom),coef='timeP0'))%>%topTable(n=Inf)%>%as.data.frame%>%rownames_to_column("protein_id")
count_fcs = eBayes(contrasts.fit(lmFit(spec_voom),coef='timeP0'))%>%topTable(n=Inf)%>%as.data.frame%>%rownames_to_column("protein_id")

fcompdf = prodalfcs%>%mutate(protein_id=msid2pid[[name]])%>%inner_join(spec_fcs)%>%inner_join(count_fcs,by='protein_id',suffix=c('_spec','_count'))

quicktest(fcompdf$diff,fcompdf$logFC_count)
quicktest(fcompdf$diff,fcompdf$logFC_spec)
quicktest(fcompdf$logFC_count,fcompdf$logFC_spec)

quicktest
quicktest(
	topTable(spec_fcs,n=Inf)%>%as.data.frame%>%rownames_to_column("protein_id")%>%arrange(protein_id)%>%.$logFC,
	topTable(count_fcs,n=Inf)%>%as.data.frame%>%rownames_to_column("protein_id")%>%arrange(protein_id)%>%.$logFC
)

#now plot
scattertitle = fcompdf%>%{cor.test(use='complete',.[['logFC_count']],.[['logFC_spec']])}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n ${.$conf.low} - ${.$conf.high}')}
plotfile<- here(paste0('plots/','fold_change_comparisons_spec_vs_count','.pdf'))
pdf(plotfile)
	ggplot(data=fcompdf,aes(x=logFC_count,y=logFC_spec))+
	geom_point()+
	scale_x_continuous(paste0('Limma log2FC - count'))+
	scale_y_continuous(paste0('Limma log2FC - Spectral coeffiecient'))+
	ggtitle(paste0('Fold Change in RPF count vs Spectral coefficient - E13 vs P0'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)

#now plot
scattertitle = fcompdf%>%{cor.test(use='complete',.[['logFC_count']],.[['diff']])}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n ${.$conf.low} - ${.$conf.high}')}
plotfile<- here(paste0('plots/','fold_change_comparisons_count_v_ms','.pdf'))
pdf(plotfile)
	ggplot(data=fcompdf,aes(x=logFC_count,y=diff))+
	geom_point()+
	scale_x_continuous(paste0('Limma log2FC - count'))+
	scale_y_continuous(paste0('Limma log2FC - Mass Spec (proDA)'))+
	ggtitle(paste0('Fold Change in RPF count vs Spectral coefficient - E13 vs P0'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)

#now plot
scattertitle = fcompdf%>%{cor.test(use='complete',.[['logFC_spec']],.[['diff']])}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n ${.$conf.low} - ${.$conf.high}')}
plotfile<- here(paste0('plots/','fold_change_comparisons_spec_v_ms','.pdf'))
pdf(plotfile)
	ggplot(data=fcompdf,aes(x=logFC_spec,y=diff))+
	geom_point()+
	scale_x_continuous(paste0('Limma log2FC - Spectral coeffiecient'))+
	scale_y_continuous(paste0('Limma log2FC - Mass Spec (proDA)'))+
	ggtitle(paste0('Fold Change in RPF count vs Spectral coefficient - E13 vs P0'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)

#now use this to implement the model.

