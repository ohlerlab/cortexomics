base::source(here::here('src/R/Rprofile.R'))

#get exons

# presalmontrdata = 'grep -e ">" /fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.pc_transcripts.fa'%>%
  # pipe%>%readLines
# salmontrdata = 'grep -e ">" /fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.pc_transcripts_filter.fa'%>%
  # pipe%>%readLines
# salmontrs = salmontrdata%>%str_extract('ENSMUST\\w+')
# presalmontrs = presalmontrdata%>%str_extract('ENSMUST\\w+')
salmontrs = here('pipeline/salmon/data/E13_total_1/quant.sf')%>%fread%>%.[[1]]%>%str_extract('ENSMUST\\w+')
dptrs = Sys.glob(here('pipeline/deepshapeprime/*ribo*/run*'))%>%head(1)%>%fread%>%.[[1]]%>%str_extract('ENSMUST\\w+')


stopifnot(dptrs%>%setdiff(salmontrs)%>%length%>%`==`(273))
stopifnot(salmontrs%>%setdiff(dptrs)%>%length%>%`==`(0L))
alltrs = salmontrs

gtf = here('pipeline/my_gencode.vM12.annotation.gtf')
fafile = Rsamtools::FaFile('pipeline/my_GRCm38.p5.genome.chr_scaff.fa')

# rm(gtf_gr,exons,cds)
if(!exists('gtf_gr')) gtf_gr<-rtracklayer::import(con=gtf,format='gtf')%>%subset(transcript_id%in%alltrs)
if(!exists('exons',where=baseenv())) exons <- gtf_gr%>%subset(type=='exon')%>%subset(transcript_id%in%alltrs)
if(!exists('cds',where=baseenv())) cds <- gtf_gr%>%subset(type=='CDS')%>%subset(transcript_id%in%alltrs)

cdsgrl = cds%>%split(.,.$transcript_id)%>%.[alltrs]
exonsgrl = exons%>%split(.,.$transcript_id)%>%.[alltrs]




# rm(gtftxdb,tputrs,fputrs)
library(GenomicFeatures)
if(!exists('gtftxdb')) gtftxdb = makeTxDbFromGRanges(gtf_gr)
if(!exists('tputrs')) tputrs = threeUTRsByTranscript(gtftxdb,use.names=TRUE)
if(!exists('fputrs')) fputrs = fiveUTRsByTranscript(gtftxdb,use.names=TRUE)

inclusiontable(gtf_gr$transcript_id,alltrs)
stopifnot(setdiff(names(tputrs),alltrs)%>%length%>%`==`(0))

tputrtrs = names(tputrs)
fputrtrs = names(fputrs)
utrtrs = intersect(fputrtrs ,tputrtrs)
is3bpmult = cdsgrl%>%width%>%`%%`(3)%>%sum%>%`==`(0)

inclusiontable(names(fputrs),alltrs)
inclusiontable(names(tputrs),alltrs)
inclusiontable(utrtrs,alltrs)

 #make gene names non redundant - put a _2 after the couple of repeated ones
if(!exists('ids_nrgname')){
	ids_nrgname <- mcols(cds)%>%as.data.frame%>%distinct(gene_name,gene_id)%>%group_by(gene_name)%>%
	mutate(new_gene_name=paste0(gene_name[1],c('',paste0('_',2:n())))[1:n()])%>%
	ungroup%>%
	select(gene_name=new_gene_name,gene_id)%>%
	left_join(cds%>%mcols%>%as.data.frame%>%distinct(transcript_id,gene_id,protein_id))
	stopifnot(ids_nrgname%>%group_by(gene_id)%>%filter(n_distinct(gene_name)>1)%>%nrow%>%identical(0L))
}

################################################################################
########Make ID conversion hashmaps
################################################################################
	
if(!file.exists('data/gnm2trid.hmp')){

trid2gid = cds%>%mcols%>%as.data.frame%>%select(transcript_id,gene_id)%>%{safe_hashmap(.[[1]],.[[2]])}
gid2gnm = ids_nrgname%>%distinct(gene_id,gene_name)%>%{safe_hashmap(.[[2]],.[[1]])}
gnm2gid = ids_nrgname%>%distinct(gene_id,gene_name)%>%{safe_hashmap(.[[1]],.[[2]])}
trid2gnm = ids_nrgname%>%select(transcript_id,gene_name)%>%{hashmap(.[[1]],.[[2]])}
gnm2gid = ids_nrgname%>%distinct(gene_name,gene_id)%>%{safe_hashmap(.[[1]],.[[2]])}
gnm2trid = ids_nrgname%>%distinct(gene_name,transcript_id)%>%{safe_hashmap(.[[1]],.[[2]])}

allgids = trid2gid[[alltrs]]%>%unique
allgnms = trid2gnm[[alltrs]]%>%unique


trid2gid%>% save_hashmap(here('data/trid2gid.hmp'))
gid2gnm%>%save_hashmap(here('data/gid2gnm.hmp'))
gnm2gid%>%save_hashmap(here('data/gnm2gid.hmp'))
trid2gnm%>%save_hashmap(here('data/trid2gnm.hmp'))
gnm2gid%>%save_hashmap(here('data/gnm2gid.hmp'))
gnm2trid%>%save_hashmap(here('data/gnm2trid.hmp'))


}else{

	trid2gid<-load_hashmap(here('data/trid2gid.hmp'))
	gid2gnm<-load_hashmap(here('data/gid2gnm.hmp'))
	gnm2gid<-load_hashmap(here('data/gnm2gid.hmp'))
	trid2gnm<-load_hashmap(here('data/trid2gnm.hmp'))
	gnm2gid<-load_hashmap(here('data/gnm2gid.hmp'))
	gnm2trid<-load_hashmap(here('data/gnm2trid.hmp'))

}

