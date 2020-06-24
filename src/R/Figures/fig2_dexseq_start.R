library(GenomicFeatures)
source(here::here('src/R/Rprofile.R'))

#DEXseq position specific changes

#Get counts at the start


library(GenomicFeatures)
fmcols <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)[start(grl@partitioning)]
}
fmcols_List <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)%>%split(grl@partitioning)
}
get_genomic_psites <- function(bam,windows,offsets,mapqthresh=200,comps=c('chrM'='chrM')) {
  require(GenomicAlignments)
  riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=windows)
  reads <- readGAlignments(bam,param=riboparam)
  mcols(reads)$length <- qwidth(reads)
  stopifnot('comp' %in% colnames(offsets))
  #
  if(is.null(offsets)){
    mcols(reads)$offset <- floor(qwidth(reads)/2)
  }else{
    useqnms <- as.character(unique(seqnames(reads)))%>%setdiff(names(comps))

    compmap = safe_hashmap(c(useqnms,names(comps)),c(rep('nucl',length(useqnms)),comps))
    stopifnot(all(compmap$values()%in%offsets$comp))
    
    mcols(reads)$offset <-
      data.frame(length=mcols(reads)$length,
        comp=compmap[[as.character(seqnames(reads))]])%>%
      safe_left_join(offsets%>%select(offset,length,comp),allow_missing=TRUE,by=c('length','comp'))%>%.$offset
  }
  #
  reads <- reads%>%subset(!is.na(mcols(reads)$offset))
  #
  # mcols(reads)$length <- width(reads)
  reads%<>%subset(!is.na(offset))
  psites <- apply_psite_offset(reads,c('offset'))%>%as("GRanges")
  mcols(psites)$length <- mcols(reads)$length
  psites
}

anno <- projmemoise(function(...){rtracklayer::import(...)})(here('../Ebp1_ribo/pipeline/gencode.vM12.annotation.gtf'))
#We want a function that 
cds <- anno%>%subset(type=='CDS')%>%split(.,.$protein_id)%>%sort_grl_st
exons <- anno%>%subset(type=='exon')%>%split(.,.$transcript_id)%>%sort_grl_st 
# sorted_is_dj <- function(grl){
#   starts <- start(grl)%>%{.[ IntegerList(as.list(rep(-1,length((.)))) ) ]}
#   ends = start(grl)%>%{.[ IntegerList(as.list(-elementNROWS(.)) ) ]}
#   all(all(starts < ends))
# }
# sorted_is_dj(cds)


pid2gidmap <- safe_hashmap(fmcols(cds,protein_id),fmcols(cds,gene_id)%>%str_replace('\\.\\d+$',''))
#NOTE - this was the end of trim_cds before, might have dependencies.

#extend exons so we can always make our windows
{
cutoffs <- 'ext_data/offsets_manual.tsv'%>%read_tsv

bakcds <- cds
cds <- bakcds
bams <- Sys.glob(here('pipeline/star/data/*/*ribo*.bam'))%>%setNames(.,basename(dirname(.)))
bams%<>%str_subset(neg=TRUE,'transcript')%>%str_subset(neg=TRUE,'E145|E16|P0')
bamseqnames <- seqinfo(Rsamtools::BamFile(bams[1]))@seqnames
sample = names(bams)[1]
chr = 'chrM'

spl_mapFromTranscripts<-function(trspacegr,exons_grl){
  exons_tr<-exons_grl%>%unlist%>%mapToTranscripts(exons_grl)%>%.[names(.)==seqnames(.)]
  ov <- findOverlaps(trspacegr,exons_tr)
  trspacegr_spl <- suppressWarnings({trspacegr[queryHits(ov)]%>%pintersect(exons_tr[subjectHits(ov)])})
  genomic_trspacegr <- mapFromTranscripts(
    trspacegr_spl,
    exons_grl
  )
  genomic_trspacegr$xHits <- queryHits(ov)[genomic_trspacegr$xHits]
  genomic_trspacegr
}

countdfs <- Sys.glob('pipeline/cdscounts/*/*trimcounts*')%>%
  setNames(.,basename(dirname(.)))%>%
  map_df(.id='sample',fread)

if(!exists('metainfo')){
  highcountgenes <- countdfs%>%group_by(gene_id)%>%filter(any(totalcount>32))%>%.$gene_id
bestcds <- countdfs%>%filter(gene_id%in%highcountgenes)%>%group_by(gene_id)%>%slice(which.max(trimcount))%>%.$protein_id
bestcds <- bestcds%>%setdiff(cds%>%unlist%>%.[seqnames(.)%in%'chrM']%>%.$protein_id)
bestcds <- cds[bestcds]
}else{
  names(cds) %<>% str_replace('\\.\\d+$','')
  bestcds = cds[metainfo%>%filter(isbest)%>%.$protein_id%>%unique]
}
#cds<-cdsbak
#cdsbak <- cds
windback <- 50
windforw <- 150
exonsexp <- exons[fmcols(bestcds,transcript_id)]%>%resize_grl(sum(width(.))+windback,'end',check=FALSE)%>%
	resize_grl(sum(width(.))+windforw,'start',check=FALSE)%>%
	.[!any(is_out_of_bounds(.))]


#let's try, 25 to 50, -15 to 5, end-30 to end, and the rest

exonsexp <- exons[fmcols(bestcds,transcript_id)]%>%resize_grl(sum(width(.))+windback,'end',check=FALSE)%>%
	resize_grl(sum(width(.))+windforw,'start',check=FALSE)%>%
	.[!any(is_out_of_bounds(.))]


rangestart=25
rangeend=50

get_cdsrange_genomic <- function(bestcds,exons,rangestart,rangeend,origin='start'){
	if(origin=='start'){
		origins <- pmapToTranscripts (bestcds%>%resize_grl(1),exons[fmcols(bestcds,transcript_id)])%>%unlist

	}else{
		origins <- pmapToTranscripts (bestcds%>%resize_grl(1,'end'),exons[fmcols(bestcds,transcript_id)])%>%unlist	
	}
	origins = origins %>% shift(rangestart) %>% resize(rangeend-rangestart+1,ignore.strand=TRUE)


	origins %>% {spl_mapFromTranscripts(.,exons_grl = exons[seqnames(.)])} %>%split(.,names(.))%>%.[names(bestcds)]
}

ribodexdf <- mclapply(bams,mc.cores=6,function(bam){
	message(bam)
	psites <- get_genomic_psites(bam,exonsexp%>%unlist%>%GenomicRanges::reduce(),cutoffs)

	startcount = countOverlaps(subject=psites,get_cdsrange_genomic(bestcds,exonsexp,-18,6))
	midcount = countOverlaps(subject=psites,get_cdsrange_genomic(bestcds,exonsexp,24,54))
	endcount = countOverlaps(subject=psites,get_cdsrange_genomic(bestcds,exonsexp,-18,6,'end'))
	restcount = countOverlaps(subject=psites,bestcds%>%
				resize_grl(sum(width(.))-6,'end')%>%
				resize_grl(sum(width(.))-18,'start'))

	ribodexdf <- data.frame(s_m18_6=startcount,s_24_54=midcount,e_m18_6=endcount,rest=restcount)%>%rownames_to_column('protein_id')%>%
			gather(region,count,-protein_id)
	message('done')
	ribodexdf
})



dexseqwidedf <- ribodexdf%>%bind_rows(.id='sample')%>%spread(sample,count)
dexseqwidedf%>%write_tsv(here('data/dexseqwidedf.tsv'))
bams%<>%setNames(.,basename(dirname(.)))
dexseqcoldata <- names(bams)%>%str_extract('E13|E175')%>%{data.frame(condition=.,sample=names(bams))}

library(DEXSeq)

dexseqwidedf%>%select(-protein_id,-region)%>%mutate_at(vars(everything()),list(~replace_na(.,0)))%>%as.matrix


ribodexdf%>%bind_rows(.id='sample')%>%group_by(region,sample)%>%summarise(mean(count==0))%>%as.data.frame

dxd  <- DEXSeqDataSet(
	dexseqwidedf%>%select(-protein_id,-region)%>%mutate_at(vars(everything()),list(~replace_na(.,0)))%>%as.matrix,
	featureID = paste0(dexseqwidedf$protein_id,dexseqwidedf$region	),
	groupID = paste0(dexseqwidedf$protein_id),
	design = ~ sample+condition*exon,
	sampleData =dexseqcoldata
)

sizeFactors(dxd) <- dxd%>%counts%>%{.[rownames(.)%>%str_detect('rest'),]}%>%estimateSizeFactorsForMatrix
dxd%<>%estimateDispersions
dxd = testForDEU(dxd)
dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")

dxd%>%rownames%>%length

mcols(dxd)%>%head(2)
dxr1 = DEXSeqResults( dxd )

library(txtplot)
dxr1$padj%>%na.omit%>%txtdensity
}
dxr1[rownames(dxr1)%>%str_detect('s_m18'),]$padj%>%na.omit%>%txtdensity

dxr1[rownames(dxr1)%>%str_detect('e_m18'),]$padj%>%na.omit%>%txtdensity
dxr1[rownames(dxr1)%>%str_detect('s_24'),]$padj%>%na.omit%>%txtdensity
dxr1[rownames(dxr1)%>%str_detect('rest'),]$padj%>%na.omit%>%txtdensity


dxr1[rownames(dxr1)%>%str_detect('s_m18'),]$padj%>%na.omit%>%`<`(0.05)%>%table


dxr1[rownames(dxr1)%>%str_detect('s_m18'),]%>%as.data.frame%>%filter(padj<0.05)%>%{pid2gidmap[[.$groupID]]}

highchanges <- dxr1[rownames(dxr1)%>%str_detect('s_m18'),]%>%subset(padj<0.05)%>%subset(rank(-log2fold_E175_E13) < 10)%>%.$groupID

dexseqwidedf%>%filter(protein_id%in%highchanges[1])%>%arrange(region)

dxr1%>%subset(groupID%in%highchanges[1])


pdf<-grDevices::pdf

pid2gnm <- safe_hashmap(trim_ids(fmcols(cds,protein_id)),fmcols(cds,gene_name))
pid2gidmap <- safe_hashmap(trim_ids(fmcols(cds,protein_id)),trim_ids(fmcols(cds,gene_id)))


dir.create(here('tables'))
startchangetbl<-tibble(
  protein_id = dxrstart$groupID,
  gene_id = pid2gidmap[[dxrstart$groupID]],
  gene_name = pid2gnm[[dxrstart$groupID]],
  log2fold_E175_E13 = dxrstart$log2fold_E175_E13,
  padj = dxrstart$padj
)%>%mutate(
  start_up = (padj<0.05) & (log2fold_E175_E13>0)
)%T>%write_tsv(.,here('tables/start_change_tbl.tsv')%T>%message)

startchangetbl$start_up%>%table


kostartchangetable <- read_tsv('../Ebp1_ribo/tables/start_change_tbl.tsv')%>%mutate(protein_id=trim_ids(protein_id))

startchangetbl%>%inner_join(kostartchangetable,by=c('protein_id','gene_id','gene_name'))%>%
  {fisher.test(table(.$padj.x<0.05,.$padj.y<0.05))}

startchangetbl%>%inner_join(kostartchangetable,by=c('protein_id','gene_id','gene_name'))%>%
  {fisher.test(table(.$padj.x<0.05,.$padj.y<0.05))}

startchangetbl%>%inner_join(kostartchangetable,by=c('protein_id','gene_id','gene_name'))%>%
  filter(padj.x<0.05,padj.y<0.05)%>%
  {cor.test(.$log2fold_si_mock,.$log2fold_E175_E13)}


startchangetbl%>%inner_join(kostartchangetable,by=c('protein_id','gene_id','gene_name'))%>%
  filter(padj.x<0.05,padj.y<0.05)%>%
  filter(log2fold_si_mock<5,log2fold_E175_E13<5)%>%
  {cor.test(.$log2fold_si_mock,.$log2fold_E175_E13)}


#now plot
plotfile<- here(paste0('plots/','startupdown_fc_comparison','.pdf'))
pdf(plotfile)
startchangetbl%>%inner_join(kostartchangetable,by=c('protein_id','gene_id','gene_name'))%>%
  filter(padj.x<0.05,padj.y<0.05)%>%
  filter(log2fold_si_mock<5,log2fold_E175_E13<5)%>%
  ggplot(aes(x=log2fold_E175_E13,y=log2fold_si_mock))+
  geom_point()+
  # scale_color_discrete(name='colorname',colorvals)+
  # scale_x_continuous(paste0('xname'))+
  # scale_y_continuous(paste0('yname'))+
  ggtitle(paste0('Start Codon Occupancy change over time in Cortex vs\n In response to Ebp1 KO'))+
  theme_bw()
dev.off()
normalizePath(plotfile)

################################################################################
########Plot of the fold change distribution.
################################################################################
  
#now plot
plotfile<- here(paste0('plots/','start_codon_usage','.pdf'))
pdf(plotfile)
dxr1[rownames(dxr1)%>%str_detect('s_m18'),]%>%
	# subset(padj<0.05)%>%
	# subset(log2fold_E175_E13 > 0)%>%
	as.data.frame%>%
	ggplot(.,aes(x=log2fold_E175_E13,fill=(padj<0.05)%in%TRUE))+
	# geom_density(alpha=I(0.8))+
	geom_histogram(alpha=I(0.8))+
	scale_fill_manual(name='Significant Change (p<0.05)',values=c('TRUE'='red','FALSE'='grey'))+
	scale_x_continuous(paste0('log2fc(Start Codon Relative Density)'),limits=c(-5,5))+
	# scale_y_continuous(paste0('yname'))+
	ggtitle(paste0('Change in start codon usage (-18 : +6 bp) \n as per DEXseq'))+
	theme_bw()
dev.off()
normalizePath(plotfile)

sigchangers <-dxr1[rownames(dxr1)%>%str_detect('s_m18'),]%>%
	subset(padj<0.05)%>%
	subset(log2fold_E175_E13 > 0)%>%
	.$groupID

endname='start'

################################################################################
########Metaplots
################################################################################
  
catnum = sigchangers%>%intersect(names(bestcds))%>%n_distinct
posdf <- meatplotdfs[[endname]]%>%
    filter(protein_id%in%sigchangers)%>%
    safe_left_join(normdf)
    # filter(totalcount>32)%>%
plotdf<-posdf   %>%mutate(count = count / density)%>% group_by(sample,pos)%>% summarise(count=sum(count))%>%
    mutate(condition = condition[sample])%>%
    group_by(sample)%>%mutate(count = count/sum(count))%>%
    # group_by(condition,pos)%>%summarise(count=mean(count))
    group_by(condition,pos)%>%
    group_modify(~mean_cl_boot(.$count))%>%
    rename('count':=y)
myymax=plotdf$count%>%max
catnm<-'Start Usage Up'
  #
plotfile<- here(paste0('plots/','metaplot_',ovnm,'_',catnm,'_',endname,'_','.pdf'))
pdf(plotfile,h=5,w=7)
#
p <- plotdf%>%    
  rename('upper':=ymax,'lower':=ymin)%>%
  ggplot(.,aes(x=pos,y=count,color=condition))+
  geom_line()+
  geom_ribbon(alpha=I(0.5),aes(y=count,ymin=lower,ymax=upper,fill=condition),color=NA)+
  scale_x_continuous(condition,breaks = profbrklist[[endname]])+
  scale_y_continuous(paste0('Psite Count Density (w/ 95% CI)'),limit=c(0,myymax),breaks = seq(0,myymax,0.03))+
  ggtitle(str_interp(paste0('Psite distribution Metaplot\nCategory ${catnm} \n(${catnum}) genes) ')))+
  # coord_cartesian(xlim = c(0,profbrklist[[1]][4]))+
  facet_grid(. ~ condition )+
  theme_bw()
print(p)
dev.off()
normalizePath(plotfile)%>%message
#now plot

notsigchangers <-dxr1[rownames(dxr1)%>%str_detect('s_m18'),]%>%
	# subset(padj<0.05)%>%
	subset(abs(log2fold_E175_E13) < 0.2)%>%
	.$groupID

catnum = notsigchangers%>%intersect(names(bestcds))%>%n_distinct
posdf <- meatplotdfs[[endname]]%>%
    filter(protein_id%in%notsigchangers)%>%
    safe_left_join(normdf)
    # filter(totalcount>32)%>%
plotdf<-posdf   %>%mutate(count = count / density)%>% group_by(sample,pos)%>% summarise(count=sum(count))%>%
    mutate(condition = condition[sample])%>%
    group_by(sample)%>%mutate(count = count/sum(count))%>%
    # group_by(condition,pos)%>%summarise(count=mean(count))
    group_by(condition,pos)%>%
    group_modify(~mean_cl_boot(.$count))%>%
    rename('count':=y)
myymax=plotdf$count%>%max
catnm<-'Start Usage Not Up'
  #
plotfile<- here(paste0('plots/','metaplot_',ovnm,'_',catnm,'_',endname,'_','.pdf'))
pdf(plotfile,h=5,w=7)
#
p <- plotdf%>%    
  rename('upper':=ymax,'lower':=ymin)%>%
  ggplot(.,aes(x=pos,y=count,color=condition))+
  geom_line()+
  geom_ribbon(alpha=I(0.5),aes(y=count,ymin=lower,ymax=upper,fill=condition),color=NA)+
  scale_x_continuous(condition,breaks = profbrklist[[endname]])+
  scale_y_continuous(paste0('Psite Count Density (w/ 95% CI)'),limit=c(0,myymax),breaks = seq(0,myymax,0.03))+
  ggtitle(str_interp(paste0('Psite distribution Metaplot\nCategory ${catnm} \n(${catnum}) genes) ')))+
  # coord_cartesian(xlim = c(0,profbrklist[[1]][4]))+
  facet_grid(. ~ condition )+
  theme_bw()
print(p)
dev.off()
normalizePath(plotfile)%>%message
#now plot

catnm=names(catlist[4])
bestcdsgns <- pid2gidmap[[names(bestcds)]]
sigchangersgns <- pid2gidmap[[sigchangers]]
nm<-function(x)setNames(x,x)
catnms <- nm(names(catlist)%>%setdiff('all'))
catnm = 'mitochondria'

fatr <- c('FALSE','TRUE')
catovtests <- lapply(catnms,function(catnm){
	fisher.test(table(bestcdsgns%in%sigchangersgns,bestcdsgns%in%catlist[[catnm]])[c(fatr),(fatr)])
})
catovtests%>%map_df(.id='category',tidy)%>%arrange(p.value)%>%select(-method,-alternative)
#exluding mitochondria
mitobestcdsgns <- pid2gidmap[[names(bestcds)]]%>%setdiff(catlist$mitochondria)
sigchangersgns <- pid2gidmap[[sigchangers]]

mitocatovtests <- lapply(catnms,function(catnm){
	if(! any(mitobestcdsgns%in%catlist[[catnm]])) return(NULL)
	fisher.test(table(as.factor(mitobestcdsgns%in%sigchangersgns),as.factor(mitobestcdsgns%in%catlist[[catnm]])))
})
mitocatovtests%>%map_df(.id='category',tidy)%>%arrange(p.value)%>%select(-method,-alternative)

catovtests%>%map_df(.id='category',tidy)%>%arrange(p.value)%>%select(-method,-alternative)


dxrstart <- dxr1[rownames(dxr1)%>%str_detect('s_m18'),]

pulldowndata <- fread('ext_data/no_DSP.csv')

pulldowndata$Gene_name %in% fmcols(cds,gene_name)

startchangetbl%>%inner_join(pulldowndata,by=c(gene_name='Gene_name'))%>%
  nrow
library(ggplot2)
#now plot

startchangetbl%>%inner_join(pulldowndata,by=c(gene_name='Gene_name'))%>%
  filter(CDS_reads_IP>32,CDS_reads_total>32)%>%
  filter(gene_id%in%str_replace(highcountgenes,'\\.\\d+',''))%>%
  filter(between(log2fold_E175_E13,-5,5))%>%
  filter(padj<0.05)%>%
  {cor.test(.$log2fold_E175_E13,.$RPKM_ratio)}


startchangetbl%>%inner_join(pulldowndata,by=c(gene_name='Gene_name'))%>%
  filter(CDS_reads_IP>32,CDS_reads_total>32)%>%
  filter(gene_id%in%str_replace(highcountgenes,'\\.\\d+',''))%>%
  filter(between(log2fold_E175_E13,-5,5))%>%
  {fisher.test(table(.$padj<0.05,.$RPKM_ratio>1.1))}

nm <-.%>% setNames(.,.)
source('../Ebp1_ribo/src/gofuncs.R')
GTOGO%<>%rename('gene_id':=ensembl_gene_id)

gotbls <-map_df(.id='direction',c('up','notup')%>%nm,function(dir){
  map_df(.id='ontology',c('MF','BP','CC')%>%nm,function(ont){
  startchangetbl%>%filter(gene_id%in%trim_ids(highcountgenes))%>%{setNames((.$padj<0.05&(.$log2fold_E175_E13>0)),.$gene_id)}%>%{if(dir=='notup') (! .) else . }%>%rungo(GTOGO%>%filter(gene_id%in%trim_ids(highcountgenes)),ont)
})})

gotbls%>%group_by(ontology)%>%arrange(elimFisher)%>%slice(1:5)
gotbls%>%group_by(ontology)%>%arrange(elimFisher)%>%slice(1:5)%>%as.data.frame

topterms <- gotbls%>%group_by(direction,ontology)%>%arrange(elimFisher)%>%slice(1:5)%>%mutate(name=Term%>%str_replace(' ','_'))
topterms%>%write_tsv('tables/topterms.tsv')


####dotplot of start vs pulldown

plotfile<- here(paste0('plots/','starteffect_vs_pulldown','.pdf'))
pdf(plotfile)
startchangetbl%>%inner_join(pulldowndata,by=c(gene_name='Gene_name'))%>%
  filter(CDS_reads_IP>32,CDS_reads_total>32)%>%
  filter(gene_id%in%str_replace(highcountgenes,'\\.\\d+',''))%>%
  filter(log2(RPKM_ratio)>0.4)%>%
  ggplot(.,aes(color = padj < 0.05,group= padj<0.05,x=(log2fold_E175_E13),y=(RPKM_ratio)))+
  geom_point()+
  geom_smooth(method='lm')+
  scale_color_discrete(name='colorname',c('TRUE'='red','FALSE'='grey'))+
  scale_x_continuous(paste0('Log2 Fold Change siRNA/Mock'))+
  scale_y_continuous(paste0('RPKM ratio - seRP'))+
  ggtitle(paste0('siRNA knockdown effect vs EBP1 pulldown enrichment'))+
  theme_bw()
dev.off()
normalizePath(plotfile)


plotfile<- here(paste0('plots/','starteffect_vs_pulldown','.pdf'))
pdf(plotfile)
startchangetbl%>%inner_join(pulldowndata,by=c(gene_name='Gene_name'))%>%
  filter(CDS_reads_IP>32,CDS_reads_total>32)%>%
  filter(gene_id%in%str_replace(highcountgenes,'\\.\\d+',''))%>%
  filter(log2(RPKM_ratio)>0.4)%>%
  # ggplot(.,aes(color = padj < 0.05,group= padj<0.05,x=(log2fold_E175_E13),y=(RPKM_ratio)))+
  ggplot(.,aes(color = padj < 0.05,group= padj<0.05,x=(log2fold_E175_E13)))+
  geom_density()+
  scale_color_discrete(name='colorname',c('TRUE'='red','FALSE'='grey'))+
  scale_x_continuous(paste0('Log2 Fold Change siRNA/Mock'))+
  scale_y_continuous(paste0('RPKM ratio - seRP'))+
  ggtitle(paste0('siRNA knockdown effect vs EBP1 pulldown enrichment'))+
  theme_bw()
dev.off()
normalizePath(plotfile)


trim_ids <- .%>%str_replace('\\.\\d+$','')



plotfile<- here(paste0('plots/','starteffect_vs_pulldown_boxplot','.pdf'))
pdf(plotfile)
startchangetbl%>%inner_join(pulldowndata,by=c(gene_name='Gene_name'))%>%
  filter(CDS_reads_IP>32,CDS_reads_total>32)%>%
  filter(gene_id%in%str_replace(highcountgenes,'\\.\\d+',''))%>%
  # filter(between(log2fold_E175_E13,-5,5))%>%
  ggplot(.,aes(color = padj<0.05,x= padj<0.05,y=RPKM_ratio))+
  geom_boxplot()+
  scale_color_discrete(name='colorname',c('TRUE'='red','FALSE'='grey'))+
  scale_x_discrete(paste0('Log2 Fold Change siRNA/Mock'))+
  scale_y_continuous(paste0('RPKM ratio - seRP'))+
  ggtitle(paste0('siRNA knockdown effect vs EBP1 pulldown enrichment'))+
  theme_bw()
dev.off()
normalizePath(plotfile)

dir.create(here('tables'))
startchangetbl<-tibble(
	protein_id = dxrstart$groupID,
	gene_id = pid2gidmap[[dxrstart$groupID]],
	gene_name = pid2gnm[[dxrstart$groupID]],
	log2fold_E175_E13 = dxrstart$log2fold_E175_E13,
	padj = dxrstart$padj
)%>%mutate(
	start_up = (padj<0.05) & (log2fold_si_mock>0)
)%T>%write_tsv(.,here('tables/start_change_tbl.tsv')%T>%message)

startchangetbl$start_up%>%table

################################################################################
########Now compare our start effect to the transcripitonal and translational effects
################################################################################
resultslist <- readRDS('data/resultslist.rds')

#now plot
plotfile<- here(paste0('plots/','fc_comp_TE_vs_Start','.pdf'))
pdf(plotfile)
resultslist$xtailTE_fraction%>%
	mutate(gene_id=str_replace(gene_id,'\\.\\d+',''))%>%
	inner_join(startchangetbl,by='gene_id')%>%
	rename('start_log2fc':=log2fold_si_mock)%>%
	rename('TE_log2fc'=log2FoldChange)%>%
	filter(padj.y<0.05)%>%
	ggplot(aes(x=TE_log2fc,y=start_log2fc))+
	geom_point()+
	# scale_color_discrete(name='colorname',colorvals)+
	# scale_x_continuous(paste0('xname'))+
	scale_y_continuous(paste0('Relative Start Usage'),limits=c(-5,5))+
	ggtitle(paste0('TE change vs Start Usage'))+
	theme_bw()
dev.off()
normalizePath(plotfile)

#now plot
plotfile<- here(paste0('plots/','fc_comp_mRNA_vs_Start','.pdf'))
pdf(plotfile)
resultslist$fraction_si_vs_mock%>%as.data.frame%>%
	mutate(gene_id=str_replace(gene_id,'\\.\\d+',''))%>%
	inner_join(startchangetbl,by='gene_id')%>%
	rename('start_log2fc':=log2fold_si_mock)%>%
	rename('mRNA_log2fc'=log2FoldChange)%>%
	filter(padj.y<0.05)%>%
	ggplot(aes(x=mRNA_log2fc,y=start_log2fc))+
	geom_point()+
	# scale_color_discrete(name='colorname',colorvals)+
	# scale_x_continuous(paste0('xname'))+
	scale_y_continuous(paste0('Relative Start Usage'),limits=c(-5,5))+
	ggtitle(paste0('mRNA change vs Start Usage'))+
	theme_bw()
dev.off()
normalizePath(plotfile)

save.image()
