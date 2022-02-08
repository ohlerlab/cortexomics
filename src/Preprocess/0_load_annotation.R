#set input files, verify they exist
base::source(here::here('src/Rprofile.R'))
gtf <- here('pipeline/my_gencode.vM12.annotation.gtf')
fafile <- here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa')




stopifnot(file.exists(fafile))
salmonfiles <- Sys.glob(here('pipeline/salmon/data/*total*/quant.sf'))
salmonfiles%T>%{stopifnot(file.exists(.))}
salmonfile <- salmonfiles[1]
stopifnot(length(salmonfile)>0)
salmontrs <- here(salmonfile)%>%fread%>%.[[1]]%>%str_extract('ENSMUST\\w+')
alltrs = salmontrs
#parse out the gtf
fafileob = Rsamtools::FaFile(fafile)


if(!file.exists(paste0(fafile,'.fai')))Rsamtools::indexFa(fafile)
if(!exists('gtf_gr')) gtf_gr<-rtracklayer::import(con=gtf,format='gtf')%>%
	subset(transcript_id%in%alltrs)
if(!exists('exons',where=baseenv())) exons <- gtf_gr%>%subset(type=='exon')%>%
	subset(transcript_id%in%alltrs)
if(!exists('cds',where=baseenv())) cds <- gtf_gr%>%subset(type=='CDS')%>%
	subset(transcript_id%in%alltrs)
#
cdsgrl <- cds%>%split(.,.$transcript_id)%>%.[alltrs]
exonsgrl <- exons%>%split(.,.$transcript_id)%>%.[alltrs]

#get utrs
library(GenomicFeatures)
if(!exists('gtftxdb')) gtftxdb <- makeTxDbFromGRanges(gtf_gr)
if(!exists('tputrs')) tputrs <- threeUTRsByTranscript(gtftxdb, use.names = TRUE)
if(!exists('fputrs')) fputrs <- fiveUTRsByTranscript(gtftxdb, use.names = TRUE)

inclusiontable(gtf_gr$transcript_id,alltrs)
stopifnot(setdiff(names(tputrs),alltrs)%>%length%>%`==`(0))

tputrtrs <- names(tputrs)
fputrtrs <- names(fputrs)
utrtrs <- intersect(fputrtrs ,tputrtrs)
is3bpmult <- cdsgrl%>%width%>%`%%`(3)%>%sum%>%`==`(0)

 #make gene names non redundant - put a _2 after the couple of repeated ones
if(!exists('ids_nrgname')){
	ids_nrgname <- mcols(cds)%>%as.data.frame%>%distinct(gene_name,gene_id)%>%
	group_by(gene_name)%>%
	mutate(new_gene_name=paste0(gene_name[1],c('',paste0('_',2:n())))[1:n()])%>%
	ungroup%>%
	select(gene_name=new_gene_name,gene_id)%>%
	left_join(cds%>%mcols%>%as.data.frame%>%
		distinct(transcript_id,gene_id,protein_id))
	mod_gnms_unique = ids_nrgname%>%group_by(gene_id)%>%
		filter(n_distinct(gene_name)>1)%>%nrow%>%identical(0L)
	stopifnot(mod_gnms_unique)
}

################################################################################
########Make ID conversion vectors
################################################################################
trid2gidv = cds%>%mcols%>%as.data.frame%>%select(transcript_id,gene_id)%>%
	{setNames(.$gene_id,.$transcript_id)}
gid2gnmv = ids_nrgname%>%distinct(gene_id,gene_name)%>%
	{setNames(.$gene_name,.$gene_id)}
gnm2gidv = ids_nrgname%>%distinct(gene_id,gene_name)%>%
	{setNames(.$gene_id,.$gene_name)}
trid2gnmv = ids_nrgname%>%select(transcript_id,gene_name)%>%
	{setNames(.$gene_name,.$transcript_id)}
gnm2gidv = ids_nrgname%>%distinct(gene_name,gene_id)%>%
	{setNames(.$gene_id,.$gene_name)}
genesofinterest=c('Nes','Tle4','Flna','Satb2','Bcl11b')


dir.create(showWarn=F,rec=TRUE,'plots/QC_plots')
dir.create(showWarn=F,rec=TRUE,'plots/Positional_Analysis')
dir.create(showWarn=F,rec=TRUE,'plots/Motif_Analysis')
dir.create(showWarn=F,rec=TRUE,'plots/Kinetic_Modeling')
dir.create(showWarn=F,rec=TRUE,'plots/Hierarch_clust')
dir.create(showWarn=F,rec=TRUE,'plots/GO_plots')
dir.create(showWarn=F,rec=TRUE,'plots/Variance_Decomp')
dir.create(showWarn=F,rec=TRUE,'plots/Trajectory_plots')
