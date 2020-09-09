N_SPLINE = 4
#Now clip the ends
ibrary(GenomicFeatures)

spl_mapFromTranscripts<-function(tmptest,tmpexons){
  exons_tr<-tmpexons%>%unlist%>%setNames(paste0('exon_',seq_along(.)))%>%mapToTranscripts(tmpexons)
  ov <- findOverlaps(tmptest,exons_tr)
  tmptest_spl <- suppressWarnings({tmptest[queryHits(ov)]%>%pintersect(exons_tr[subjectHits(ov)])})
  genomic_tmptest <- mapFromTranscripts(
    tmptest_spl,
    # exons_tr[subjectHits(ov)]%>%split(.,seqnames(.))
    tmpexons
  )
  genomic_tmptest$xHits <- queryHits(ov)[genomic_tmptest$xHits]
  genomic_tmptest
}


#get offset values
offsets<-read_tsv(here('ext_data/offsets_manual.tsv'))

#get the bams for each library
bams <- Sys.glob(here('pipeline/star/data/*/*.bam'))%>%str_subset(negate=TRUE,'transcript')

dexseqbams<-bams%>%str_subset('_(ribo)')
dexseqbamsrna<-bams%>%str_subset('_(total)')


STARTCLIP = 15
ENDCLIP = 5


metainfo<-read_tsv(here('pipeline/exprdata/limma_genemetadata.tsv'))

#select
library(GenomicAlignments)
library(GenomicFeatures)
metainfo$dTE%>%table

cds_counttrs <- cds%>%subset(protein_id %in% (metainfo%>%filter(highcount,is_gid_highest)%>%.$protein_id))
cdssplit <- cds_counttrs%>%split(.,.$transcript_id)
trlens<-cdssplit%>%width%>%sum%>%enframe('transcript_id','length')%>%filter(length>(STARTCLIP*3+ENDCLIP*3))
trregions<-GRanges(trlens$transcript_id,IRanges(STARTCLIP*3+1,trlens$length-(ENDCLIP*3)))
trregions_genome<-trregions%>%spl_mapFromTranscripts(.,cdssplit)
trregions_genome$transcript_id <- seqnames(trregions)[trregions_genome$xHits]
trregions_genome$gene_id <- Rle(cds_counttrs$gene_id)[match(trregions_genome$transcript_id,cds_counttrs$transcript_id)]




#
dexseqexons<-trregions_genome%>%subset(gene_id%in%(unique(dexseqgenes)))
#
#get exon level counts with our offsets
exonribopsites <- mclapply(dexseqbams,function(bam) get_genomic_psites(bam,dexseqexons%>%IRanges::reduce(.),offsets=offsets)%>%
	countOverlaps(dexseqexons,.))
exonribopsites%>%length
exonribopsitesna <- mclapply(mc.cores=1,dexseqbamsrna,function(bam) get_genomic_psites(bam,invertStrand(dexseqexons)%>%IRanges::reduce(.),offsets=NULL)%>%
	countOverlaps(invertStrand(dexseqexons),.))
exonribopsitesna%>%length
#
dexseqsamples <- c(dexseqbams,dexseqbamsrna)%>%dirname%>%basename
save.image('data/fig2_dexseq.Rdata')
#
dexseqsampledf<-data.frame(tmp=c(dexseqsamples))%>%separate(tmp,c('time','assay','rep'))
#
library(DEXSeq)
dexseqcountdf <- c(exonribopsites,exonribopsitesna)%>%setNames(dexseqsamples)%>%as.data.frame 
exonids <- data.frame(gid=dexseqexons$gene_id)%>%group_by(gid)%>%mutate(n=1:n())%>%{paste0(.$gid,'_',.$n)}
rownames(dexseqcountdf)<-exonids

#sampledf
dexseqsampledf%<>%
	mutate(time = factor(time,sort(unique(time))))%>%
	mutate(assay = factor(assay,c('total','ribo')))%>%
  mutate(ntime = as.numeric(as.factor(time)))
#
N_SPLINE=4
dexseqsampledf %<>% bind_cols(ns(dexseqsampledf$ntime,N_SPLINE)%>%set_colnames(paste0('spl',1:N_SPLINE))%>%as.data.frame)
#
#create DEXseq object with 
te_dxd <- DEXSeqDataSet(
	dexseqcountdf,
	featureID =  exonids,
	groupID =  paste0(dexseqexons$gene_id),
	 ~ exon*(spl1+spl2+spl3+spl4)*assay,
   # ~ exon*ns(ntime,N_SPLINE)*assay,
	 sampleData = dexseqsampledf
)
#
te_dxd%>%counts%>%colMedians
#
#create 
sizeFactors(te_dxd) <- te_dxd%>%counts%>%estimateSizeFactorsForMatrix
te_dxd%<>%estimateDispersions
#
te_dxd = testForDEU(te_dxd,
	fullModel = 	 ~ exon*(spl1+spl2+spl3+spl4)*assay,
	reducedModel = 	 ~ exon+(spl1+spl2+spl3+spl4)+assay+exon:(spl1+spl2+spl3+spl4)+exon:assay+(spl1+spl2+spl3+spl4):assay
)


mcols(te_dxd)%>%head(2)
te_dxr1 = DEXSeqResults( te_dxd )


te_dxd = testForDEU(te_dxd,
	fullModel = 	 ~ exon*time*assay,
	reducedModel = 	 ~ exon+time+assay+exon:time+exon:assay+time:assay
)


mcols(te_dxd)%>%head(2)
te_dxr1 = DEXSeqResults( te_dxd )

te_dxr1$padj%>%vals%>%{table(.<0.05)}

te_dxr1$groupID%>%n_distinct

dtegenes <- metainfo%>%filter(dTE)%>%.$gene_id%>%unique
dtegenes%>%length
exon_dtegenes <- te_dxr1%>%subset(padj<0.05)%>%.$groupID
allgids <- metainfo$gene_id%>%unique
table(TE_change = allgids%in%dtegenes,exon_TE_change = allgids%in%exon_dtegenes)%T>%print%>%
fisher.test()




metainfo%>%subset(gene_id%in%c(exon_dtegenes))%>%.$gene_name%>%unique%>%sort


mcols(te_dxd)%>%colnames
mcols(te_dxd)%>%colnames
colData(te_dxd)$timeassay <- colData(te_dxd)%>%as.data.frame%>%paste0(.$time,'_',.$assay)
te_dxd = estimateExonFoldChanges( te_dxd, fitExpToVar=c('timeassay'))

#to get the dexseq for TE.

#Select genes with TE exons

#properotion of genes with TE exons

#GO enrichment of genes with TE exons

# #Proprotion of TE diff which is TE exons
# get_cds_offsets<-	function(reads,offsets,compartments){
# 	#
#     compartmentcol <- if('compartment'%in%colnames(offsets)) 'compartment' else NULL
#     #
#     if(!'compartment' %in% colnames(mcols(reads))){
#       mcols(reads)$compartment <- as(compartments,"List")[seqnames(reads)]%>%unlist(use.names=F)
#     }
#     if(!'length' %in% colnames(mcols(reads))){
#       mcols(reads)$length <- qwidth(reads)
#     }
#     #
#     phasecol <- if('phase'%in%colnames(offsets)) 'phase' else NULL
#     joincols <- c('length',compartmentcol,phasecol)
#     #
# 	assert_that(all(has_name(mcols(reads),joincols)))
# 	#
#   reads%>%
#     mcols%>%
#     as_tibble%>%
#     safe_left_join(offsets%>%
#       ungroup%>%
#       # select(-phase,-score),by=c('length','compartment'),
#       select(offset,one_of(joincols)),by=joincols,
#       allow_missing=TRUE
#     )%>%
#     .$offset
# }

