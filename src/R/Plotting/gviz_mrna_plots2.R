library(Gviz)
library(assertthat)
library(magrittr)
library(Rsamtools) 
library(tidyverse)
library(rtracklayer)
library(stringr)
library(biomaRt)
library(BSgenome.Mmusculus.UCSC.mm10)
select=dplyr::select

#todo - should probalby just get the data and then put that in a track
#slow otherwise.

source = function(path,...){
  require(stringr)
  require(magrittr)
  path%<>%str_replace('/home/dharnett','~')
  message(path)
  base::source(path,...)
}
root='~/projects/cortexomics'

datafiles <- Sys.glob(file.path(root,'data/bigwigs/*/*/*transcript*bw'))
assert_that(length(datafiles)>4)
# datafiles %<>% str_subset('E13|P0')

annofile <- '/fast/projects/cubit/0.12.0/static_data/annotation/GENCODE/M12/GRCm38/gencode.vM12.annotation.gtf'
assert_that(file.exists(annofile))

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

inframewith<-function(gr1,gr2){
  fp2 = ifelse(strand(gr2)=='-',end(gr2),start(gr2))
  fp1 = ifelse(strand(gr1)=='-',end(gr1),start(gr1))
  (fp2 - fp1) %% 3 ==0
}
#get the data track names
import.bw.neg <- function(file,selection){
  data <- import.bw(file,sel=selection)
  data$score <- data$score * -1
  data
}

make_track<-function(file, trackname, sizefactor, vtr,isneg=FALSE,...){

  selection =   GRanges(vtr$transcript_id,IRanges(1,sum(width(vtrexons))))

  if(str_detect(trackname,'total')) isneg=TRUE
  if(isneg){ 
    importfunction = import.bw.neg
    # ylims = c(-30,0)
  }else{
    importfunction= function(file,selection) import.bw(file,sel=selection)
    # ylims = c(0,10)
  }
  DataTrack(file,name = trackname,chromosome=vtr$transcript_id,stream=TRUE,importFunction = importfunction)

}


###input files

# ########
if(!exists('annotation')) annotation <- annofile%>% {rtracklayer::import(.)}
annotation = annofile %>% import
annotation$transcript_id %<>% str_replace('\\.[0-9]+$','')
  
genes_to_check = list('Pa2g4','Satb2','Nes','Flna')
#get our main transcript 
annotation%>%subset(gene_name %in% genes_to_check)%>%.$transcript_id%>%table

vtrs =
  annotation%>%
  subset(gene_name %in% genes_to_check)%>%
  subset(type=='transcript')%>%
  {assert_that(length(.)>1);.}

vtr = vtrs%>%.[1]


#get exons of transcript for liftover
vtrexonlist <- annotation%>%
    subset(type%in%c('exon'))%>%
    split(.,.$transcript_id)


gene_featurecounts<- Sys.glob(file.path(root,'data/feature_counts/data/*/*genefeature_counts'))
genecounts <-
  gene_featurecounts%>%
  setNames(.,basename(.))%>%
  map(data.table::fread,skip=1,select=c(1,7))%>%
  Reduce(f=partial(left_join,by=c('Geneid')))
sample_size_factors<-DESeq2::estimateSizeFactorsForMatrix(as.matrix(genecounts[,-1]))
sample_size_factors%<>%setNames(sample_size_factors%>%names%>%dirname%>%basename)


make_coverage_plot<-function(datafiles, vtr, annotation){
  options(ucscChromosomeNames=FALSE)


  gchr = seqnames(vtr)%>%as.character
  gstart = start(vtr)
  gend = end(vtr)
  vtrname = vtr$transcript_id

  vtrexons = vtrexonlist[[vtr$transcript_id]]
  totlength <- width(vtrexons)%>%sum

  #get the annotation for our gene and transfer it onto the transcript coordinates
  startcod<-annotation%>%
    subset(transcript_id==vtrname)%>%
    subset(type=='start_codon')%>%
    unique

  startcod<- mapToTranscriptsGR(startcod,vtrexons)
  if(as.logical(strand(vtr)=='-')){
    revcompfun = reverseComplement
  }else{
    revcompfun = identity
  }
  startcodons<-
    matchPattern(revcompfun(DNAString('ATG')),Mmusculus[[gchr%>%as.character]][gstart:gend])%>%
    as('IRanges')%>%
    shift(gstart-1)%>%
    GRanges(gchr,.,strand='-')%>%
    mapToTranscriptsGR(vtrexons)%>%
    keep(.,inframewith(.,startcod))%>%
    {strand(.)='+';.}%>%
    sort

  codtrack = Gviz::AnnotationTrack(startcodons,feature='start_codon',chr=vtr$transcript_id,shape='box')

  exontrack =
    annotation%>%
    subset(transcript_id%>%str_detect(vtrname))%>%
    subset(type%in%c('UTR','CDS'))%>%
    .[,'type']%>%
    {.$feature<-as.character(.$type);.}%>%
    mapToTranscriptsGR(vtrexons)%>%
    Gviz::GeneRegionTrack(.,thinBoxFeature=c("UTR"))

  fileregex <- ifelse(strand(vtr)=='-','ribo.*pos|(total.*neg)','ribo.*neg|(total.*pos)')

  dfiles <- datafiles%>%
    data_frame(file=.)%>%
    mutate(trackname = file%>%basename%>%str_replace('(.transcript)?.bw','') )%>%
    filter(trackname%>%str_detect(fileregex))%>%
    mutate(sample = trackname%>%str_replace('\\.[a-z]+',''))%>%
    mutate(sizefactor = sample_size_factors[sample])%>%
    mutate(isneg = str_detect(file,'total'))%>%
    separate(sample,c('time','assay','rep'), remove=FALSE)%>%
    arrange(time,rep,assay)
  

  get_sig_grs<-function(file, trackname, sizefactor, vtr,isneg=FALSE,...){
    if(isneg){ 
      importfunction = import.bw.neg
    }else{
      importfunction= function(file,selection) import.bw(file,sel=selection)
    }
    importfunction
  }

  dTrack<- dfiles%>%
    as.list%>%
    {pmap(.l=., .f = get_sig_grs, vtr = vtr)}

  #get the region we're selecting
  vtrgr <-GRanges(vtr$transcript_id,IRanges(1,totlength))
  #get signal over this region in gr objects
  siggrs <- dfiles$file%>%setNames(.,.)%>%mclapply(import.bw,selection=vtrgr)
  #reverse some of them (total rna)
  siggrs[dfiles$isneg]  %<>%map(function(x) {x$score %<>% multiply_by(-1);x})
  #scale by size factors
  siggrs%<>%map2(dfiles$sizefactor, function(x,y){x$score%<>%divide_by(y);x})
  #now average replicates
  dfiles$sampe_norep <- dfiles$sample%>%str_replace('_\\d+','')
  siggrs <- siggrs%>%split(.,  dfiles$sampe_norep)
  siggrs <- mclapply(siggrs,function(x){  
    as((coverage(x[[1]],weight='score')[vtr$transcript_id] + 
    coverage(x[[2]],weight='score')[vtr$transcript_id])/2,'GRanges')
  })
  dfiles = dfiles%>%group_by(sampe_norep)%>%dplyr::slice(1)
  #now get ymax for all of these
  allylims <- siggrs%>%map(~.$score)%>%range%>%{c(floor(min(0,.[1])),ceiling(max(.[2],0)))}

  trackylims <- list(c(0,allylims[2]), c(allylims[1], 0))[as.numeric(dfiles$isneg)+1]
  trackylims <- rep(trackylims[which.max(abs(trackylims))],2)
  # trackylims <- list( c(allylims[1], 0), c(0,allylims[2]))[as.numeric(dfiles$isneg)+1]

  dTrack <-  lapply(seq_along(siggrs),function(i) {
    tname <- dfiles$trackname[i]%>%str_replace('\\.\\w+','')
    DataTrack(siggrs[[i]], name = tname, ylim = trackylims[[i]], cex.title=.5)
  })

  #length of transript to plot

  plotTracks(c(dTrack,c(codtrack,exontrack)),chromosome = vtr$transcript_id,from=1,to=totlength, type="hist", 
    main = paste0(vtr$gene_name,'_',vtr$transcript_name,'_',vtr$transcript_id)
  )


  plotTracks(c(dTrack,c(codtrack,exontrack)),chromosome = vtr$transcript_id,from=start(startcod)-60,to=start(startcod)+60, type="hist", 
    main = paste0(vtr$gene_name,'_',vtr$transcript_name,'_',vtr$transcript_id)
  )

}


file.path(root,'plots','mRNA_covplots')%>%dir.create


for(i in seq_along(vtrs)){
  tname <- vtrs$transcript_name[i]
  gname <- vtrs$gene_name[i]
  covplot <- file.path(root,'plots','mRNA_covplots',paste0(gname,'_',tname,'.pdf'))
  pdf(covplot)
  make_coverage_plot(datafiles, vtrs[i], annotation)
  dev.off()
  message(covplot)
}



