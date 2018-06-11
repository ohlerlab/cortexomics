library(Gviz)
library(magrittr)
library(Rsamtools) 
library(tidyverse)
library(rtracklayer)
library(stringr)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm10)
select=dplyr::select

options(ucscChromosomeNames=FALSE)
        
mapToTranscriptsGR <- function(gr,transcriptgr,trcol = 'transcript_id'){
  
  names(ranges(transcriptgr)) = mcols(transcriptgr)[[trcol]]
  transcriptgr %<>% split(.,names(.))
  out = GenomicFeatures::mapToTranscripts(  
     gr,
    transcripts = transcriptgr
  )
  mcols(out) = mcols(gr)[mcols(out)[['xHits']],]
  out
}
###input files

########
annofile <- '/fast/projects/cubit/0.12.0/static_data/annotation/GENCODE/M12/GRCm38/gencode.vM12.annotation.gtf'
assert_that(file.exists(annofile))
annotation <- annofile%>% {rtracklayer::import(.)}
annotation = annofile %>% import
annotation$transcript_id %<>% str_replace('\\.[0-9]+$','')
  
#get our main transcript 
vtr =
  annotation%>%
  subset(gene_name%>%str_detect('Pa2g4'))%>%subset(type=='transcript')%>%
  subset(transcript_name%>%str_detect('Pa2g4-001'))
vtr$transcript_id<-vtr$transcript_id%>%str_replace('\\.[0-9]+$','')

#get the annotation for our gene and transfer it onto the transcript coordinates
startcod<-annotation%>%
  subset(gene_id ==vtr$gene_id)%>%
  subset(type=='start_codon')%>%
  unique

inframewith<-function(gr1,gr2){
  fp2 = ifelse(strand(gr2)=='-',end(gr2),start(gr2))
  fp1 = ifelse(strand(gr1)=='-',end(gr1),start(gr1))
  (fp2 - fp1) %% 3 ==0
}

gchr = seqnames(vtr)%>%as.character
gstart = start(vtr)
gend = end(vtr)
vtrname = vtr$transcript_id

#get exons of transcript for liftover
vtrexons = annotation%>%
  subset(transcript_id%>%str_detect(vtrname))%>%subset(type%in%c('exon'))

startcodons<-matchPattern(reverseComplement(DNAString('ATG')),Mmusculus[[gchr%>%as.character]][gstart:gend])%>%
  as('IRanges')%>%
  shift(gstart-1)%>%
  GRanges(gchr,.,strand='-')%>%
  mapToTranscriptsGR(vtrexons)%>%
  keep(.,inframewith(.,mapToTranscriptsGR(startcod,vtrexons)))%>%
  {strand(.)='+';.}%>%
  sort

codtrack =     Gviz::AnnotationTrack(startcodons,feature='start_codon',chr=vtr$transcript_id,shape='box')

exontrack =
  annotation%>%
  subset(transcript_id%>%str_detect(vtrname))%>%
  subset(type%in%c('UTR','CDS'))%>%
  .[,'type']%>%
  {.$feature<-as.character(.$type);.}%>%
  mapToTranscriptsGR(vtrexons)%>%
  Gviz::GeneRegionTrack(.,thinBoxFeature=c("UTR"))
  
  # {Gviz::AnnotationTrack(.,thinBoxFeature=c("utr"),chromosome=vtr$transcript_id,shape='box',col='black',fill='yellow')}

#create the datatrack for gviz
options(ucscChromosomeNames=FALSE)
selection =   GRanges(vtr$transcript_id,IRanges(1,sum(width(vtrexons))))


#get the data track names
import.bw.neg <- function(file,selection){
  data <- import.bw(file,sel=selection)
  data$score <- data$score * -1
  data
}

make_track<-function(file,trackname,isneg=FALSE,...){
  if(str_detect(trackname,'total')) isneg=TRUE
  if(isneg){ 
    importfunction = import.bw.neg
    ylims = c(-30,0)
  }else{
    importfunction= function(file,selection) import.bw(file,sel=selection)
    ylims = c(0,10)
  }
  DataTrack(file,name = trackname,chromosome=vtr$transcript_id,stream=TRUE,importFunction = importfunction,ylim = ylims)
}
datafiles <- Sys.glob(file.path(root,'data/mergedbigwigs/*/*/*transcript*'))[c(1:4,21:24)]
dTrack <- datafiles%>%
  data_frame(file=.)%>%
  mutate(trackname = file%>%basename%>%str_replace('(.transcript)?.bw','') )%>%
  filter(trackname%>%str_detect('ribo.*pos|(total.*neg)'))%>%
  {map2(.$file,.$trackname,.f = make_track)}

# testTrack<-DataTrack(file='tmp.foo',genome='mm10',name='test',chromosome=vtr$transcript_id,stream=TRUE,importFunction = testimportfun)
# plotTracks(c(testTrack,codtrack,exontrack),chromosome=vtr$transcript_id,from=1,to=300, type="hist",ylim = c(-10,10))


stop()

pdf(file.path(root,'exploration','mRNA_covplot.pdf'),h=12,w=7)
plotTracks(c(dTrack[1:4],c(codtrack,exontrack)),chromosome = vtr$transcript_id,from=1,to=107+88, type="hist")
plotTracks(c(dTrack[-c(1:4)],c(codtrack,exontrack)),chromosome = vtr$transcript_id,from=1,to=107+88, type="hist")
plotTracks(c(codtrack,exontrack),chromosome = vtr$transcript_id,from=1,to=107+88, type="hist")
dev.off()


plotTracks(c(dTrack[3],codtrack,exontrack),chr=vtr$transcript_id,from=1,to=107+88, type="hist")

pdf('~/projects/cortexomics/mRNA_covplot_exon1.pdf',h=18,w=7)
plotTracks(c(dTrack[3],codtrack,exontrack),chr=vtr$transcript_id,from=1,to=107+88, type="hist")
dev.off()

plotTracks(c(dTrack,codtrack,exontrack),chr=vtr$transcript_id,from=1,to=sum(width(vtrexons)), type="hist")

plotTracks(c(dTrack,codtrack,exontrack),chr=vtr$transcript_id,from=1200,to=2500, type="hist")

####################################################
#create annotation tracks for exons
#read counts on pa2g4 sections

# startcodons<-
secondcod = matchPattern(reverseComplement(DNAString('ATG')),Mmusculus[[gchr%>%as.character]][(gstart):gend])%>%
  as('IRanges')%>%
  shift(gstart-1)%>%
  GRanges(gchr,.,strand='-')%>%
  subsetByOverlaps(startcod%>%resize(100))%>%
  .[1]


zone1 <- GRanges('chr10',IRanges(start(startcod),end(startcod+26)))
# zone1 <- GRanges('chr10',IRanges(start(startcod),gend))
zone2 <- GRanges('chr10',IRanges(end(secondcod)+1,start(zone1)-1))
zone3 <- GRanges('chr10',IRanges(gstart,end(secondcod)))
zones <- c(zone1,zone2,zone3)
zones$Gene_ID = c("5'UTR-AUG1","AUG1-AUG2","AUG2-STOP")

anndf = data.frame(
  GeneID= zones$Gene_ID,
  Chr = seqnames(zones)%>%as.character,
  Start = zones%>%start,
  End = zones%>%end,
  Strand = zones%>%strand
)

# zones%>%export('data/pa2g4_zones.gff3')

# get bam files
bam_files <- 
        #get bam files from star folder
        list.files(file.path(root,"data",'star/data/'),
                       pattern=".*\\.bam$", full.names = T, recursive = T)%>%
        #exclude transcript bams
        grep(val=TRUE,inv=TRUE,patt='star_transcript')%>%
        #full paths
        paste('realpath',.)%>%
        map_chr(.%>%system(.,intern=TRUE))
stopifnot(bam_files%>%file.exists()%>%all)
#also get library sizes
samples <- bam_files%>%dirname%>%basename%>%setNames(.,.)
#library size
libsizes <- map_dbl( samples, .%>%
    sprintf('data/star/reports/%s/%s.bam.bamstats.txt',.,.)%>%
    readLines%>%
    str_extract(regex('(?<=reads mapped:\\s{0,80})\\d+$'))%>%
    keep(Negate(is.na))%>%
    as.numeric
)
is_ribo <- grepl(x = samples , patt = 'ribo')
#now count reads on our areas
getcounts = partial(Rsubread::featureCounts,
                                             annot.ext = anndf,
                                             useMetaFeatures = F,
                                             isPairedEnd = 0,
                                             nthreads = 4,
                                             allowMultiOverlap=TRUE,
                                             tmpDir='/tmp'
)
ribocountres = quietly(getcounts)(files = bam_files[is_ribo],strandSpecific=2)
totalcountres = quietly(getcounts)(files = bam_files[!is_ribo],strandSpecific=1)

ribocounts = ribocountres$result$counts
totalcounts = totalcountres$result$counts
ribocounts %<>% set_colnames(samples%>%str_subset('ribo'))
totalcounts %<>% set_colnames(samples%>%str_subset('total'))
#bind them and then get counts per million
counts <- cbind(ribocounts,totalcounts)[,samples]
counts %<>% sweep(MARGIN = 2,STATS = libsizes/1e6,FUN='/')

counts<-
  bind_cols(zone=rownames(counts),as_tibble(counts))%>%
    gather(dataset,CPM,-zone)%>%
    separate(dataset,c('time','assay','rep'))

scale_y_log2_linlabels<-scale_y_continuous(
  breaks = function(lims) {
    nb = 7
    step = ceiling((lims[2] - lims[1])/7)
    lims[1] = nb*(floor(lims[1]/nb))
    lims[2] = nb*(ceiling(lims[2]/nb))
    seq(from=lims[1],to=lims[2],by=step);
  },
  labels = function(x) format((2^x),digits = 2)
)

for(zone_i in unique(counts$zone)){
  
  counts %>% 
    filter(zone==zone_i)%>%
    mutate(CPM = CPM)%>%
    group_by(zone,assay)%>%
    mutate(CPM = log2(CPM/mean(CPM[time=='E13'])))%>%
    { 
      qplot(data=.,y=CPM,x = time, geom = 'point' , ylab =  'T/T0')+
        facet_grid( scale = 'free',assay ~ . )+
        expand_limits(y=c(floor(min(.$CPM)),ceiling(max(.$CPM))))+
        expand_limits(y=c(-2,1))+
        scale_y_log2_linlabels+
        ggtitle(str_interp("CPM for ${zone_i}"))+
        theme_bw()+
        theme(text = element_text(size = 16))
      
    }%>%print
  
}

markerTE <-  
  counts%>%
  # filter(name==name_i)%>%
  dplyr::select(assay,zone,time,rep,CPM)%>%
  spread(assay,CPM)%>%
  mutate(TE = ribo/total)%>%
  dplyr::select(-ribo,-total)

markerTE<-
  markerTE%>%
  mutate(refTE=mean(TE[time==time[1]]))%>%
  group_by(zone)%>%mutate(TE = log2(TE/refTE))%>%
  dplyr::select(-refTE)%>%
  mutate(TE =  replace(TE,is_in(TE,c(Inf)),NA))


for(zone_i in (markerTE$zone%>%unique)){
  
  
  markerTE %>% 
    
    filter(zone==zone_i)%>%
  { 
    ymax = ungroup(.)%>%select(.,`TE`)%>%.[[1]]%>%keep(is.finite)%>%keep(Negate(is.na))%>%keep(Negate(is.nan))%>%max%>%multiply_by(1.2)
    
 
      qplot(data=.,y=`TE`,x = time, geom = 'point' ,ylab = 'TE - T/T0')+
        expand_limits(y=c(floor(min(.$TE)),ceiling(max(.$TE))))+
        scale_y_log2_linlabels+    
        ggtitle(str_interp("Translational Efficiency for ${zone_i}"))+
        theme_bw()+
        theme(text = element_text(size = 16))
    
  }%>%print

  }


# 
# # datafiles <- Sys.glob('data/star/data/P0_ribo_2/P0_ribo_2.bam'
# # annofile = '~/projects/cortexomics/data/transcripts.gff3'
# 
# # #biomart import gtrack code
# # mart = biomaRt::useMart('ENSEMBL_MART_MOUSE')
# # listDatasets(mart)
# # mart = biomaRt::useMart('ENSEMBL_MART_MOUSE',dataset = 'mc57bl6nj_gene_ensembl')
# # biomTrack <- BiomartGeneRegionTrack(mart,
# #                                     chromosome = as.character(gchr), start = gstart, end = gend,
# #                                     name = "ENSEMBL")
# # 
# 
# ########
# annofile <- '~/projects/cortexomics/data/static_local/gencode.vM12.annotation.gtf'
# transcripts <- annofile%>% {rtracklayer::import(.)}
# transcripts = annofile %>% import
# vgene = transcripts%>%subset(gene_name%>%str_detect('Pa2g4') & (type=='gene'))
# gchr = vgene%>%seqnames
# gstart = (vgene%>%start)
# gend = (vgene%>%end)
# #create the gene track for gviz
# gtrack = transcripts%>%{
#   transcripts = .
#   transsubset=transcripts%>%subsetByOverlaps(vgene)
#   transsubset$feature %<>% as.character
#   transsubset%>%subset(feature=='UTR')
#   transsubset<-transsubset[ transsubset$feature %in% c("CDS","UTR")]
#   transsubset$exon = transsubset$exon_id
#   transsubset$gene = transsubset$gene_name
#   transsubset$gene = transsubset$transcript_name
#   transsubset$symbol = transsubset$transcript_name
#   transsubset$feature %<>% as.character
#   transsubset$transcript= transsubset$transcript_id
#   # transsubset$group = transsubset$transcript_id 
#   gtrack = GeneRegionTrack(transsubset,thinBoxFeature=c("UTR"),showId=TRUE,geneSymbol=TRUE)
# }
# 
# #create the datatrack for gviz
# datafiles = Sys.glob("~/projects/cortexomics/data/star/data/E13_ribo_1/E13_ribo_1.bam"
# dTrack <- DataTrack(range=datafile, genome="mm10", name="Coverage", 
#                     window=-1, chromosome=gchr, importFunction=strandedBamImport, 
#                     stream=TRUE) 
# 
# plotTracks(list(dTrack,gtrack),chr = as.character(gchr), from = gend - 500, to = gend, col=c("red", "blue"), 
#            groups=c("+", "-"), type="hist", col.histogram=NA)
# 
