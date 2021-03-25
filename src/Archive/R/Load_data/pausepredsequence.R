
pauspredfiles = Sys.glob('pipeline/pauspred/*ribo*.pause')


#get the pausfiles as dfs
pausepredsites <- pauspredfiles%>%setNames(.,.)%>%map_df(.id='sample',.%>%fread()%>%
	rename('start':=coordinate_position)%>%
	rename('seqname':=gene_name)%>%
	mutate(end=start)%>%
	filter(start!='NA'))%>%
	{GRanges(.)}
pausepredsites$sample%<>%basename
pausepredsites$sample%<>%str_replace('.pause','')

gtf = here('pipeline/my_gencode.vM12.annotation.gff3')
#get exons
if(!exists('gtf_gr')) gtf_gr<-rtracklayer::import(con=gtf,format='gtf')
exons <- gtf_gr%>%subset(type=='exon')
cds <- 'pipeline/my_gencode.vM12.annotation.cds.gtf'%>%rtracklayer::import(.)

fastafile='pipeline/my_GRCm38.p5.genome.chr_scaff.fa'

library(rtracklayer)

starts = cds%>%split(.,.$protein_id)%>%sort_grl_st%>%resize_grl(3)
stops = cds%>%split(.,.$protein_id)%>%sort_grl_st%>%resize_grl(3,'end')


pausepredsites%<>% subset(Pause_score > 20 )

pausepredsites$startoverlap = pausepredsites%>%resize(width=100,'center')%>%
	countOverlaps(starts)%>%`!=`(0)

pausepredsites$stopoverlap = pausepredsites%>%resize(width=100,'center')%>%
	countOverlaps(stops)%>%`!=`(0)

pausepredsitessplit = pausepredsites%>%split(.,.$sample)
pausepredsitessplit%<>%.[sort(names(.))]

pausepredsitessplit%<>%.[names(.)%in%c(samplereps)&names(.)%in%c(samples)]


samples = pausepredsites$sample%>%unique
samplereps= samples%>%str_replace('1$','2x')%>%str_replace('2$','1x')%>%str_replace('x$','')
samplereps%<>%setNames(samples)

for(sample in names(pausepredsitessplit)){
	pausepredsitessplit[[sample]]$replicates <- 
		pausepredsitessplit[[sample]] %>% 
		resize(.,13,'center')%>%
		countOverlaps(pausepredsitessplit[samplereps[sample]])%>%
		`!=`(0)
}

filteredpauses <- pausepredsitessplit%>%
	lapply(function(gr){
		gr <- gr[order(gr$Pause_score)]%>%subset(replicates)%>%subset(!startoverlap)%>%subset(!stopoverlap);head(gr,1000)
	})%>%
	GRangesList%>%
	unlist%>%
	identity


sampi = samples[1]

for(sampi in unique(filteredpauses$sample)){
	sampesites <- filteredpauses%>%subset(sample==sampi)
	seqs = sampesites%>%{paste0(mcols(.)[['50_upstream_seq']],mcols(.)[['50_downstream_seq(including_pause_position)']])}
	seqs%<>%DNAStringSet
	names(seqs) = paste0(seqnames(sampesites),'_',start(sampesites))
	seqfile = here('pipeline','pauspred',paste0(sampi,'.seq.fa'))
	writeXStringSet(seqs,seqfile)
}
list.files(here('pipeline','pauspred'))


offsets='ext_data/offsets_manual.tsv'%>%read_tsv
bams <- Sys.glob(here('pipeline/star/data/*/*ribo*.bam'))%>%setNames(.,basename(dirname(.)))
bams%<>%str_subset(neg=TRUE,'transcript')%>%str_subset('ribo')

testwinds = filteredpauses%>%
	subset(sample==sampi)%>%
	subset(Pause_score > 20)
testwinds%<>%setNames(as.character(.))
library(GenomicFeatures)

psites <- testwinds%>%resize(51,'center')%>%get_genomic_psites(bams[6],.,offsets)

psites%>%mapToTranscripts(testwinds%>%resize(51,'center'))%>%
        coverage%>%
        {set_rownames(as.matrix(.),names(.))}%>%
        colSums%>%txtplot

indexFa(FaFile('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'))
testwinds[1]%>%resize(7,'start')%>%getSeq(x=FaFile('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'),.)

#now as GRanges objects

#calculate overlap between successive positions

#calculate overlap with start codons