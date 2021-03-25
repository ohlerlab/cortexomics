library(tidyverse)
library(Gviz)
library(magrittr)


# devtools::load_all('~/work/Applications/SaTann')
annofile <-'data/test_arabidopsis.gtf.gz_Rannot'
anno <- get(load(annofile))
# plastidres <- riboseqcoutput%>%map(~subset(.,as.character(seqnames) %in% c('Chr1'))%>%subset(start > 1e6)%>%subset(start < 2e6))
# save(plastidres,file='data/root_plastid_for_SaTAnn')

#run_SaTAnn('data/root_for_SaTAnn',annotation_file='data/test_arabidopsis.gtf.gz_Rannot',genome_seq = "data/test_arabidopsis.fa",n_cores=20)
# run_SaTAnn('data/root_plastid_for_SaTAnn',annotation_file='data/test_arabidopsis.gtf.gz_Rannot',genome_seq = "data/test_arabidopsis.fa",n_cores=20)

# plastid_satann <- get(load('data/root_plastid_for_SaTAnn_final_SaTAnn_results'))

# selgene <- genetxnum%>%sort%>%tail(1)%>%names

# genetrnums <- anno$txs_gene%>%elementNROWS
# geneseltxnum <- plastid_satann$selected_txs%>%str_replace('\\.\\d+$','')%>%table
# txdiscardednum <- genetrnums[names(geneseltxnum)] - geneseltxnum

# geneorfnums <- anno$cds_txs_coords%>%as.data.frame%>%group_by(gene_id)%>%tally%>%arrange(desc(n))%>%{setNames(.$n,.$gene_id)}
# names(plastid_satann$ORFs_tx)<-NULL

# ORFs_picked <- mcols(plastid_satann$ORFs_tx)%>%.[,c('ORF_id_tr','gene_id')]%>%as.data.frame%>%group_by(gene_id)%>%tally%>%arrange(desc(n))%>%{setNames(.$n,.$gene_id)}

# numorfs_discarded <- geneorfnums[names(ORFs_picked)] - ORFs_picked
# orfs_discarded <- which(numorfs_discarded > 1)
# # plastid_satann%>%names
# # plastid_satann$ORFs_tx%>%mcols%>%.$pval%>%`<`(0.05)

# txs_discarded <- which(txdiscardednum >2 & (geneseltxnum > 2))





riboseqcoutput<-get(load('data/root_for_SaTAnn'))

#So get a gene for which both transcripts and ORFs were discarded
# selgene <- names(txs_discarded) %>% intersect(names(orfs_discarded)) %>% head(1)


genes_withmulttx_junc_sup_df<-riboseqcoutput$junctions%>%mcols%>%as.data.frame%>%dplyr::select(tx_name,gene_id,reads)%>%mutate(tx_name=map(tx_name,~ as.list(.)))%>%mutate(gene_id=as.character(gene_id))%>%unnest%>%identity%>%
	filter(reads>0)%>%
	group_by(gene_id)%>%
	summarise(n_sup_tx = n_distinct(tx_name))	

altsplicejuncs<-riboseqcoutput$junctions%>%{d=mcols(.);d$start=start(.);d$end=end(.);.}%>%as.data.frame%>%dplyr::select(tx_name,gene_id,reads,start,end)%>%
	mutate(tx_name=map(tx_name,~ as.list(.)))%>%
	mutate(gene_id=as.character(gene_id))%>%
	unnest%>%identity%>%
	distinct(gene_id,start,end,reads)%>%
	filter(reads>0)%>%
	group_by(gene_id)%>%
	filter(any(table(start)>1),any(table(end)>1))

# selgene <- altsplicejuncs$gene_id%>%unique%>%.[3]
selgene<-genes_withmulttx_junc_sup_df%>%arrange(desc(n_sup_tx))%>%.$gene_id%>%head(1)%>%tail(1)




{
selgenerange <- anno$genes%>%subset(gene_id==selgene)
subsetres <- riboseqcoutput%>%map(~subset(selgenerange))
save(subsetres,file='data/root_subset_for_SaTAnn')
#run_SaTAnn('data/root_for_SaTAnn',annotation_file='data/test_arabidopsis.gtf.gz_Rannot',genome_seq = "data/test_arabidopsis.fa",n_cores=20)
run_SaTAnn('data/root_for_SaTAnn',annotation_file='data/test_arabidopsis.gtf.gz_Rannot',genome_seq = "data/test_arabidopsis.fa",n_cores=20,gene_id=selgene)
subset_satann <- get(load('data/root_for_SaTAnn_final_SaTAnn_results'))
}

{

#orfs_quantified_track <-
orfs_quantified_gen <-  subset_satann$ORFs_gen%>%subset(.,str_detect(names(.),selgene))

mcols(orfs_quantified_gen) <- mcols(subset_satann$ORFs_tx)[match(names(orfs_quantified_gen),subset_satann$ORFs_tx$ORF_id_tr),]
seltxs <- orfs_quantified_gen$transcript_id%>%unique

orfs_quantified_gen$feature <- 'CDS'
orfs_quantified_gen$transcript=orfs_quantified_gen$transcript_id
orfs_quantified_tr <- anno$exons_tx%>%.[unique(orfs_quantified_gen$transcript_id)]

#get the orfs, then get the negative coverage
#add transcript info to the ORFs_tx object
seqinf <- Seqinfo(names(anno$exons_tx),anno$exons_tx%>%width%>%sum)
#Now get the negatives for each ORF
utrs <- subset_satann$ORFs_tx%>%subset(gene_id==selgene)%>%keepSeqlevels(seltxs)%>%{seqinfo(.)<-seqinf[seltxs];.}%>%
	coverage%>%
	as('GRanges')%>%subset(score==0)%>%mapFromTranscripts(anno$exons_tx)%>%
	{.$transcript <- names(anno$exons_tx)[.$transcriptsHits];.}%>%
	{.$feature='utr';.}

orfquantgr <- c(
	orfs_quantified_gen[,c('feature','transcript')]
	,utrs[,c('feature','transcript')]
)
orfquantgr$feature[orfquantgr$feature=='CDS'] <- names(orfquantgr)[orfquantgr$feature=='CDS']

orfscores<- orfquantgr$feature%>%unique%>%setNames(match(.,orfs_quantified_gen$ORF_id_tr)%>%orfs_quantified_gen$TrP_pNpM[.],.)
orfcols <- orfscores%>%{./max(na.omit(.))}%>%
	map_chr(~possibly(rgb,'white')(0,.,0))%>%setNames(orfquantgr$feature%>%unique)

###Define non selected
disctxs<-anno$txs_gene[selgene]%>%unlist%>%.$tx_name%>%unique%>%setdiff(seltxs)

disc_orfquantgr <- anno$cds_txs_coords%>%keepSeqlevels(disctxs,'coarse')%>%{seqinfo(.)<-seqinf[disctxs];.}%>%	coverage%>%
	as('GRanges')%>%
	subset(.$score==0)%>%
	{	
		txgr = .
		out = mapFromTranscripts(txgr,anno$exons_tx)
		out$score = txgr$score[out$xHits]
		out
	}%>%
	{.$transcript <- names(anno$exons_tx)[.$transcriptsHits];.}%>%
	{.$feature=ifelse(.$score==0,'utr','CDS');.}
disc_orfquantgr %<>% c(.,anno$cds_txs[disctxs]%>%unlist%>%{.$feature=rep('CDS',length(.));.$transcript=names(.);.})


fakejreads <- riboseqcoutput$junctions%>%subset(any(gene_id==selgene))%>%resize(width(.)+2,'center')%>%
	{.$cigar <- paste0('1M',width(.)-2,'N','1M');.}

fakejreads <- fakejreads[map2(seq_along(fakejreads[]),fakejreads$reads,rep)%>%unlist]



# anno$cds_txs_coords%>%

# library(GenomicFeatures)

# tx_mapped_cds <- anno$cds_txs%>%unlist%>%mapToTranscripts(.,anno$exons_tx)%>%reduce
# seltx_mapped_cds<-tx_mapped_cds%>%subset(str_detect(seqnames,selgene))%>%resize(width(.)-3)
# isquanted <- seltx_mapped_cds %in% orfs_quantified
# stopifnot(any(isquanted))
# discarded_orfs <- seltx_mapped_cds[!isquanted]
# discarded_orfs_gen <- discarded_orfs %>% mapFromTranscripts(anno$exons_tx)

}





{
orfquantgrscores = orfs_quantified_gen$TrP_pNpM[match(names(orfquantgr),orfs_quantified_gen$ORF_id_tr)]

orfscores_nona <- orfscores[!is.na(orfscores)]

make_fakelegend_ranges <-function(selgenerange,orfscores){
	seq(start(selgenerange),end(selgenerange),length.out = (length(na.omit(orfscores))*2)+2 )%>%
	head(-1)%>%tail(-1)%>%
	floor%>%
	matrix(byrow=T,nrow=length(na.omit(orfscores)))%>%
	as.data.frame%>%
	transmute(seqnames=as.character(seqnames(selgenerange)[1]),start=V1,end=V2,feature=names(sort(na.omit(orfscores))))%>%
	GRanges(.)
}

plotfile <- 'summary_plots_satann.pdf'
pdf(plotfile,width=14,h=7)
SaTAnn::plot_SaTAnn_results(
	for_SaTAnn_file='data/root_for_SaTAnn',
	SaTAnn_output_file='data/root_for_SaTAnn_final_SaTAnn_results',
	annotation_file = 'data/test_arabidopsis.gtf.gz_Rannot')
dev.off()
normalizePath(plotfile)

library(Gviz)
plotfile <- 'locusplot.pdf'
par(mfrow=c(1,2))
pdf(plotfile,width=14,h=7)
plottitle <- paste0('ORFquant: ',selgene)
p<-plotTracks(main=plottitle,cex.main=2,legend=TRUE,
	from=start(selgenerange),to=end(selgenerange),#zoomed in on the orf in question
	sizes=c(1,1,1,1,1,1,1,1),
	c(
		GenomeAxisTrack(range=selgenerange),
		# rnaseqtrack, # plot the riboseq signal
		# txs_discarded_Track,
		# txs_selected_track,
		# DataTrack(riboseqcoutput$P_sites_all%>%subsetByOverlaps(selgenerange),type='hist'),
		GeneRegionTrack(name='discarded\ntranscripts',anno$exons_tx[disctxs]%>%unlist%>%{.$transcript=names(.);.$feature=rep('exon',length(.));.},fill='#F7CAC9'),
		GeneRegionTrack(exon='green',name='selected\ntranscripts',anno$exons_tx[seltxs]%>%unlist%>%{.$transcript=names(.);.$feature=rep('exon',length(.));.},fill='#F7CAC9'),
		DataTrack(legend=TRUE,name='P Sites',col.histogram='green',riboseqcoutput$P_sites_all%>%subsetByOverlaps(selgenerange),type='hist'),
		AlignmentsTrack(name='\n\n\n\n\n',col.sashimi='green',fakejreads[,'cigar'],type='sashimi',sashimiNumbers=TRUE),
		# GeneRegionTrack(discarded_orfs_gen),
		GeneRegionTrack(name='Discarded\nORFs',disc_orfquantgr,collapse=FALSE,thinBoxFeature='utr',CDS='blue',utr='white'),
		GeneRegionTrack(box.legend=TRUE,name='Selected\nORFs',orfquantgr,collapse=FALSE,
			thinBoxFeature='utr',CDS='red',utr='white',
			showFeatureId=TRUE,
			id=orfquantgrscores)%>%
		{displayPars(.)[names(orfcols)]<-orfcols;.},
		DataTrack(type='heatmap',showColorBar=TRUE,
			gradient = c(orfcols[which.min(orfscores)],orfcols[which.max(orfscores)]),
			name='legend',data=matrix(na.omit(orfscores),nrow=3,ncol=1),
			selgenerange%>%resize(width(.),'end')%>%resize(1)
		)
	),
	col.labels='black',
	chr=seqnames(selgenerange)
)
dev.off()
normalizePath(plotfile)

plotfile <- 'colscale.pdf'
pdf(plotfile,width=4,height=4)
data.frame(orf = unique(names(orfquantgr)),score=unique(orfquantgrscores))%>%filter(!is.na(score))%>%
	mutate(col=orfcols[as.character(orf)])%>%
	ggplot(aes(x=orf,fill=I(col),color=I('black'),y=score,ymin=0))+stat_identity(geom='bar')+scale_y_continuous(limits= c(0,max(na.omit(orfquantgrscores))))+
	theme(axis.text.x=element_text(angle=45,vjust=0.5))
dev.off()
normalizePath(plotfile)


orfquantgrscores = orfs_quantified_gen$TrP_pNpM[match(,orfs_quantified_gen$ORF_id_tr)]
}