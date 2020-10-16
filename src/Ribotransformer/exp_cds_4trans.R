
base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
library(GenomicFeatures)
windback <- 60
windforw <- 60
outprefix = 'ribotrans_'

#get expanded exons - so that we never hit the edge of the exons when expanding
#our cds by the necessary amount
is3bpmult = cdsgrl%>%width%>%sum%>%`%%`(3)%>%`==`(0)
filtcdsgrl = cdsgrl[is3bpmult]
filtcdsgrl = filtcdsgrl%>%sort_grl_st
cdsseq = extractTranscriptSeqs(x=fafile,filtcdsgrl)
validORF = cdsseq%>%translate%>%str_detect('^M[^*]+\\*$')
filtcdsgrl = filtcdsgrl[validORF]

#get exons for our cds
cdsexonsgrl <- exonsgrl[fmcols(filtcdsgrl,transcript_id)]
cdsexonsgrl%<>%sort_grl_st
#get an object representing the CDS In transript space
cdstrspace = pmapToTranscripts(filtcdsgrl,cdsexonsgrl[fmcols(filtcdsgrl,transcript_id)])
#ensure all cds map cleanly to the exons
stopifnot(cdstrspace%>%elementNROWS%>%`==`(1))

cdsstartpos = setNames(start(cdstrspace@unlistData),names(cdstrspace))
endpos = sum(width(cdsexonsgrl))-end(cdstrspace@unlistData)

#expand our first exon when needed
startposexpansion = pmax(0,windback - cdsstartpos + 1)
#expand/trim the 5' end of the exons
cdsexonsgrl@unlistData[start(cdsexonsgrl@partitioning)]%<>%resize(width(.)+startposexpansion,'end')
#expand or trim the last exon when needed
endposexpansion = pmax(0,windforw - endpos)
cdsexonsgrl@unlistData[cdsexonsgrl@partitioning@end]%<>%resize(width(.)+endposexpansion,'start')

#now map our cds to that
cds_exptrspc = pmapToTranscripts(filtcdsgrl,cdsexonsgrl)
stopifnot(cds_exptrspc%>%elementNROWS%>%`==`(1))

expcds_exptrspc=cds_exptrspc
#expand our cds exons
expcds_exptrspc%<>%resize(width(.)+windback,'end',ignore.strand=TRUE)
#and expand the 3' ends
expcds_exptrspc%<>%resize(width(.)+windforw,'start',ignore.strand=TRUE)
#now back to genome space
expcdsgenspace = spl_mapFromTranscripts(expcds_exptrspc%>%unlist,cdsexonsgrl)
expcdsgenspace = split(expcdsgenspace,names(expcdsgenspace))
#check they are intact, the above worked.
splitcdstr = names(expcds_exptrspc)[head(which(expcds_exptrspc%>%elementNROWS%>%`>`(1)),1)]
stopifnot(length(splitcdstr)==0)
#get the sequences
expcdsgenspaceseq <- 
	expcdsgenspace%>%
	sort_grl_st%>%
	extractTranscriptSeqs(x=fafile)


fastanames <- paste(sep='|',
	fmcols(cdsexonsgrl,transcript_id),
	fmcols(cdsexonsgrl,gene_id),
	fmcols(cdsexonsgrl,havana_gene),
	fmcols(cdsexonsgrl,havana_transcript),
	fmcols(cdsexonsgrl,transcript_name),
	fmcols(cdsexonsgrl,gene_name),
	sum(width(cdsexonsgrl)),
	paste0('UTR5:1-',start(cds_exptrspc)-1),
	paste0('CDS:',start(cds_exptrspc),'-',end(cds_exptrspc)),
	paste0('UTR3:',1+end(cds_exptrspc),'-',sum(width(cdsexonsgrl))),
	'|')

names(expcdsgenspaceseq) <- fastanames

library(rtracklayer)

#write the expanded exons sequences to disk
writeXStringSet(expcdsgenspaceseq,paste0(outprefix,'_gencodeM12.fa'))

#also write our cds coordinates to disk in the new trspace
new_trspc_anno <- c(
	GRanges(as.character(seqnames(cds_exptrspc)),IRanges(1,sum(width(cdsexonsgrl))))%>%{.$type='exon';.},
	cds_exptrspc%>%unlist%>%{.$type='CDS';.}
)
new_trspc_anno%>%export(paste0(outprefix,'_trspaceanno.gtf'))
as.data.frame(cds_exptrspc)%>%select(seqnames,start,end)%>%export(paste0(outprefix,'_trcds.txt'))

#as.data.frame(cds_exptrspc)%>%select(seqnames,start,end)%>%mutate(phase=((end-start)%%3))%>%.$phase%>%table

if(FALSE){
# cds4salmonseq%>%writeXStringSet('pipeline/')
expcdsgenspaceseq%>%writeXStringSet(here('pipeline/gencode.vM12.trimtrs.fa'))

expcdata = getSeq(x=Rsamtools::FaFile('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/ext_data/gencode.vM12.pc_transcripts_filter.fa'))

exheader = expcdata%>%names%>%.[1000]
exttr = exheader%>%str_extract('ENSMUST\\d+')
cdsgrl[exttr]%>%pmapToTranscripts(exonsgrl[exttr])
exonsgrl[exttr]%>%width%>%sum
}
# current strategy - quantify with salmon and dp
# but this means that we transcripts that hardly have any UTR
# This makes the counts reflect differences in UTR lenght rather than just TE
# I should quantify TE using the 




 