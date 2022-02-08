library(GenomicFeatures)
MAPQTHRESH <- 200

base::source('src/Rprofile.R')
ribo_TE_tbl <- allvoom$E%>%as.data.frame%>%
	rownames_to_column('gene_id')%>%
	gather(dataset,val,-gene_id)%>%
	separate(dataset,into=c('time','assay','replicate'))%>%
	spread(assay,val)%>%
	group_by(gene_id,time)%>%
	summarise(TE = mean(ribo) -  mean(total) ,Ribo=mean(ribo))
#
# ribo_TE_tbl%<>%safe_left_join(x=.,y=counttbl_reps%>%filter(assay=='total')%>%select(protein_id=feature_id,time,highcount)%>%distinct)
ribo_TE_tbl%<>%safe_left_join(x=.,y=txi_counttbl_reps%>%filter(assay=='ribo')%>%distinct(gene_id,time,highcountribo = highcount))
ribo_TE_tbl$highcount = ribo_TE_tbl$highcountribo
ribo_TE_tbl%<>%safe_left_join(ribotpms)
ribo_TE_tbl_cor<-ribo_TE_tbl%>%filter(highcount)%>%group_by(time)%>%summarise(cor = cor(TPM,TE))

################################################################################
########This will compare out TE differential genes to others using a few different measures
################################################################################

allTEchangedf<-read_tsv('tables/manuscript/go_all_highcount_updown.tsv')

techanges <- Sys.glob('pipeline/xtail/*.txt')%>%setNames(.,basename(.))%>%map_df(.%>%fread%>%select(log2fc,adj_p_value,gene_name),.id='time')
techanges$time%<>%str_replace('xtail_','')%>%str_replace('.txt','')
techanges%<>%filter(!is.na(adj_p_value))

techanges%>%group_by(time)%>%group_slice(1)%>%filter(adj_p_value<0.05)

cairo_pdf('plots/QC_plots/TEchange_horizontaldist.pdf',w=12,h=6)
techanges%>%group_by(time)%>%arrange(log2fc)%>%
	# filter(adj_p_value<0.05)%>%
	mutate(xpos = 1:n())%>%
	ggplot(data=.,aes(alpha = adj_p_value < 0.05,color=adj_p_value < 0.05, ymax=pmax(log2fc,0),ymin=pmin(log2fc,0),x=xpos))+
	geom_linerange()+
	scale_color_manual(values = c('grey','black'))+
	scale_y_continuous(limits=c(-3,3))+
	facet_wrap(time~.,scales='free_x',ncol=1)+
	theme_bw()
dev.off()
normalizePath('plots/QC_plots/TEchange_horizontaldist.pdf')%>%message


fcdirnums<-techanges%>%filter(gene_name%in%highcountgnms)%>%group_by(time,sign(log2fc))%>%
filter(adj_p_value<0.05)%>%
filter(abs(log2fc)>1)%>%
tally
fcdirnums%>%tally
fcdirnums %>% summarise(ratio=n[1]/n[2])

#define genes changing their TE
alltechangeids <- techanges %>% filter(adj_p_value<0.05)

library(rtracklayer)
#Define various continous predictors.
tr_gid_df<-ids_nrgname%>%distinct(transcript_id,gene_name)
transcript_fp_utrlen <- 
	'pipeline/fputrs.gtf'%>%import%>%split(.,.$transcript_id)%>%width%>%sum%>%enframe('transcript_id','fputr_length')
transcript_tp_utrlen <- 
	'pipeline/tputrs.gtf'%>%import%>%split(.,.$transcript_id)%>%width%>%sum%>%enframe('transcript_id','tputr_length')



stops = cds%>%
	split(.,.$protein_id)%>%sort %>%
	{.[as(ifelse( any(strand(.)=='-') , 1, elementNROWS(.)),'IntegerList')]} %>%
	resize(3,'end')

starts = cds%>%
	split(.,.$protein_id)%>%sort %>%
	{.[as(ifelse( any(strand(.)=='+') , 1, elementNROWS(.)),'IntegerList')]} %>%
	resize(3,'start')



################################################################################
########Not currently used - gave negative results
################################################################################
	
if(F){
bams <- here('pipeline/star/data/*/*.bam')%>%Sys.glob%>%
	str_subset(neg=TRUE,'transcript')%>%
	str_subset('_ribo_|total')

startstoppsites <- mclapply(bams,F=get_genomic_psites,c(stops,starts)%>%unlist)
startstoppsites %<>% setNames(bams%>%dirname%>%basename)

#get dfs with the counts over codons
startcountsdf <- startstoppsites%>%map_df(.id='sample',.%>%countOverlaps((unlist(starts)),.)%>%data.frame(startcount=.,protein_id=(unlist(starts)$protein_id)))
stopcountsdf <- startstoppsites%>%map_df(.id='sample',.%>%countOverlaps((unlist(stops)),.)%>%data.frame(stopcount=.,protein_id=(unlist(stops)$protein_id)))

startstopcountdf<-allsegcounts_nz%>%
	select(sample,protein_id,total,centercount,length)%>%
	left_join(startcountsdf)%>%
	left_join(stopcountsdf)%>%
	safe_left_join(protid_gid_df)

cor(use='complete',tmp$count,tmp$startcount)
txtplot(tmp$count,tmp$startcount)

startstopcountdf%<>%mutate(startcountfrac = startcount / (total/length))

startstopcountdf%>%filter(sample=='P0_ribo_1')%>%filter(total>32)%>%{txtplot(log2((1+.$total)/.$length),log2(1+.$startcount))}

pdf('tmp.pdf')
samp2plot = 'P0_ribo_2'
startstopcountdf%>%filter(sample==samp2plot)%>%filter(total>32)%>%filter(startcount<total)%>%
	group_by(gene_id)%>%slice(which.max(startcount))%>%
	# mutate(get_density(log10((1+total)/length),log10(startcount+1),15)) + 
	ggplot(data=.,aes(x=log10((1+total)/length),y=log10(startcount+1))) + geom_point(size=I(0.5),alpha=I(0.5))+
	theme_bw()+ggtitle(samp2plot)
dev.off()
normalizePath('tmp.pdf')


pdf('tmp.pdf')
samp2plot = 'E13_ribo_1'
startstopcountdf%>%filter(sample==samp2plot)%>%filter(total>32)%>%filter(stopcount<total)%>%
	group_by(gene_id)%>%slice(which.max(stopcount))%>%
	# mutate(get_density(log10((1+total)/length),log10(stopcount+1),15)) + 
	ggplot(data=.,aes(x=log10((1+total)/length),y=log10(stopcount+1))) + geom_point(size=I(0.5),alpha=I(0.5))+
	theme_bw()+ggtitle(samp2plot)
dev.off()
normalizePath('tmp.pdf')


################################################################################
########DEXseq analysis of the stop codons
################################################################################
	

stopdexdf<-startstopcountdf%>%
	group_by(gene_id)%>%
	filter(gene_id %in% highcountgenes)%>%
	filter(protein_id == protein_id[which.max(stopcount)])%>%
	filter(any(stopcount>10))%>%
	select(sample,gene_id,centercount,stopcount)%>%
	gather(region,count,-sample,-gene_id)

dexseqwidedf<-stopdexdf%>%
	# filter(sample%>%str_detect(c('E13','E16')))%>%
	filter(sample %>%str_detect('ribo'))%>%spread(sample,count)
dexseqcoldata<-stopdexdf%>%
	# filter(sample%>%str_detect(c('E13','E16')))%>%
	filter(sample %>%str_detect('ribo'))%>%ungroup%>%distinct(sample)%>%separate(sample,into=c('time','assay','rep'),remove=F)%>%as.data.frame

dxd  <- DEXSeqDataSet(
	dexseqwidedf%>%ungroup%>%select(-gene_id,-region)%>%mutate_at(vars(everything()),list(~replace_na(.,0)))%>%as.matrix,
	featureID = paste0(dexseqwidedf$gene_id,dexseqwidedf$region	),
	groupID = paste0(dexseqwidedf$gene_id),
	design = ~ sample + exon + time:exon,
	sampleData =dexseqcoldata
)

sizeFactors(dxd) <- dxd%>%counts%>%{.[rownames(.)%>%str_detect('center'),]}%>%estimateSizeFactorsForMatrix
dxd%<>%estimateDispersions
dxd = testForDEU(dxd)
dxd = estimateExonFoldChanges( dxd, fitExpToVar="time")

mcols(dxd)%>%head(2)
dxr1 = DEXSeqResults( dxd )

stopchangedf<-data.frame(
	gene_id=rownames(dxr1)%>%str_extract('[^:]+'),
	stop_change = dxr1$log2fold_P0_E13[rownames(dxr1)%>%str_detect('stop')]
)

# (dxr1%>%.$log2fold_E16_E13%>%`>`(2)) & ()




cairo_pdf(h=5,w=5,'plots/QC_plots/cumdist_TEchange_stopchange.pdf')
plot1=allTEchangedf%>%
	left_join(stopchangedf)%>%
	ungroup%>%
	arrange(stop_change)%>%
	filter(!is.na(stop_change))%>%
	mutate(TE_change = case_when(up==1 ~ 'up',down==1~ 'down',TRUE ~ 'None'))%>%
	group_by(TE_change)%>%
	mutate(`Cumulative Fraction`= (1:n() / n()) )%>%
	ggplot(data=.,aes(y=`Cumulative Fraction`,color=TE_change,x=stop_change))+
	geom_line()+
	theme_bw()
# ggarrange(plot1,)
plot1
dev.off()
normalizePath('plots/QC_plots/cumdist_TEchange_stopchange.pdf')

cairo_pdf(h=5,w=5,'plots/QC_plots/stop_madist.pdf')
plotMA(dxr1%>%subset(featureID%>%str_detect('stop')))
dev.off()
normalizePath('plots/QC_plots/stop_madist.pdf')

changesgid<-allTEchangedf%>%filter(up|down)%>%.$gene_id

library(scales)
stop_ma_dist_df<-dxr1%>%subset(!featureID%>%str_detect('stop'))
stop_ma_dist_df$TEcol <- ifelse(stop_ma_dist_df$groupID%in%changesgid,'red','grey')%>%alpha(.,0.5)
library(ggExtra)
# 
cairo_pdf(h=5,w=5,'plots/QC_plots/stop_madist.pdf')
p=ggplot(data=as.data.frame(stop_ma_dist_df)%>%arrange(TEcol),aes(y= -log2fold_P0_E13,x=log2(exonBaseMean),color=TEcol))+geom_point(size=I(0.1),alpha=I(0.5))+
	theme_bw()+coord_cartesian(ylim=c(-2,2))+
	scale_y_continuous(name='Estimated Pre Stop Codon Fold Change, P0 to E13')+
	scale_x_continuous(name='CDS Base Mean')+
	scale_color_manual(name='TE Change',values = stop_ma_dist_df$TEcol%>%unique%>%rev,labels=c('No TE Change','TE Change'))
ggMarginal(p,margins='y',type='hist',groupFill=TRUE,groupColour=TRUE)
dev.off()
normalizePath('plots/QC_plots/stop_madist.pdf')

abstedf <- bestmscountebayes$coef[,'TE']%>%enframe('uprotein_id','TE')%>%left_join(ms_id2protein_id%>%distinct(uprotein_id,gene_id))

cairo_pdf('plots/QC_plots/abs_TE_vs_start_occ.pdf')
stop_ma_dist_df%>%as.data.frame%>%left_join(abstedf,by=c('groupID'='gene_id'))%>%qplot(data=.,log2fold_P0_E13,TE)%>%print
dev.off()
normalizePath('plots/QC_plots/abs_TE_vs_start_occ.pdf')


as.data.frame(stop_ma_dist_df)%>%{split(.$log2fold_P0_E13,.$TEcol)}%>%{wilcox.test(.[[1]],.[[2]])}

# cairo_pdf(h=5,w=5,'plots/QC_plots/stop_madist_hist.pdf')
# %>%{.$TEchange<-.$groupID%in%changesgid;.}%>%as.data.frame%>%filter(padj<0.05)%>%qplot(data=.,x=log2fold_E16_E13,fill=TEchange,geom='histogram')
# dev.off()
# normalizePath('plots/QC_plots/stop_madist_hist.pdf')

dxr1%>%as.data.frame%>%filter(padj<0.05)%>%filter(log2fold_E16_E13< -5)%>%.$featureID




################################################################################
########DEXseq analysis of the start codons
################################################################################
		

startdexdf<-startstopcountdf%>%
	group_by(gene_id)%>%
	filter(gene_id %in% highcountgenes)%>%
	filter(protein_id == protein_id[which.max(startcount)])%>%
	filter(any(startcount>10))%>%
	select(sample,gene_id,centercount,startcount)%>%
	gather(region,count,-sample,-gene_id)

dexseqwidedf<-startdexdf%>%
	# filter(sample%>%str_detect(c('E13','E16')))%>%
	filter(sample %>%str_detect('ribo'))%>%spread(sample,count)
dexseqcoldata<-startdexdf%>%
	# filter(sample%>%str_detect(c('E13','E16')))%>%
	filter(sample %>%str_detect('ribo'))%>%ungroup%>%distinct(sample)%>%separate(sample,into=c('time','assay','rep'),remove=F)%>%as.data.frame

startdxd  <- DEXSeqDataSet(
	dexseqwidedf%>%ungroup%>%select(-gene_id,-region)%>%mutate_at(vars(everything()),list(~replace_na(.,0)))%>%as.matrix,
	featureID = paste0(dexseqwidedf$gene_id,dexseqwidedf$region	),
	groupID = paste0(dexseqwidedf$gene_id),
	design = ~ sample + exon + time:exon,
	sampleData =dexseqcoldata
)

startdxd%<>%estimateSizeFactors
startdxd%<>%estimateDispersions
startdxd = testForDEU(startdxd)
startdxd = estimateExonFoldChanges( startdxd, fitExpToVar="time")

mcols(startdxd)%>%head(2)
startdxr1 = DEXSeqResults( startdxd )

startchangedf<-data.frame(
	gene_id=rownames(startdxr1)%>%str_extract('[^:]+'),
	start_change = startdxr1$log2fold_P0_E13[rownames(startdxr1)%>%str_detect('start')]
)




cairo_pdf(h=5,w=5,'plots/QC_plots/cumdist_TEchange_startchange.pdf')
plot1=allTEchangedf%>%
	left_join(startchangedf)%>%
	ungroup%>%
	arrange(start_change)%>%
	filter(!is.na(start_change))%>%
	mutate(TE_change = case_when(up==1 ~ 'up',down==1~ 'down',TRUE ~ 'None'))%>%
	group_by(TE_change)%>%
	mutate(`Cumulative Fraction`= (1:n() / n()) )%>%
	ggplot(data=.,aes(y=`Cumulative Fraction`,color=TE_change,x=start_change))+
	geom_line()+
	theme_bw()
# ggarrange(plot1,)
plot1
dev.off()
normalizePath('plots/QC_plots/cumdist_TEchange_startchange.pdf')

library(scales)
start_ma_dist_df<-startdxr1%>%subset(!featureID%>%str_detect('start'))
start_ma_dist_df$TEcol <- ifelse(start_ma_dist_df$groupID%in%changesgid,'red','grey')%>%alpha(.,0.5)
library(ggExtra)
# 
cairo_pdf(h=5,w=5,'plots/QC_plots/start_madist.pdf')
p=ggplot(data=as.data.frame(start_ma_dist_df)%>%arrange(TEcol),aes(y= -log2fold_P0_E13,x=log2(exonBaseMean),color=TEcol))+geom_point(size=I(0.1),alpha=I(0.5))+
	theme_bw()+coord_cartesian(ylim=c(-2,2))+
	scale_y_continuous(name='Estimated Start Codon Fold Change, P0 to E13')+
	scale_x_continuous(name='CDS Base Mean')+
	scale_color_manual(name='TE Change',values = start_ma_dist_df$TEcol%>%unique%>%rev,labels=c('No TE Change','TE Change'))
ggMarginal(p,margins='y',type='hist',groupFill=TRUE,groupColour=TRUE)
dev.off()
normalizePath('plots/QC_plots/start_madist.pdf')



# get_density <- function(x, y, ...) {
#   dens <- MASS::kde2d(x, y, ...)
#   ix <- findInterval(x, dens$x)
#   iy <- findInterval(y, dens$y)
#   ii <- cbind(ix, iy)
#   return(dens$z[ii])
# }
}

