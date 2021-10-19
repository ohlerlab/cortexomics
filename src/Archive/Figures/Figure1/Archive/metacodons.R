################################################################################
########This script produces metacodon plots
################################################################################
library(rlang)
library(tidyverse)
library(GenomicRanges)
library(magrittr)
library(ggplot2)
library(riboWaltz)

`%>%` <- magrittr::`%>%`
`%<>%` <- magrittr::`%<>%`
GRanges <- GenomicRanges::GRanges

#TODOS
#remove the requirement for expression
library(R.utils)
args <- R.utils::commandArgs(trailingOnly = TRUE, asValues = TRUE,
    defaults = list(
        bamparentfolder = "../Liuetal_pipeline/pipeline/star/ORFext/data/", #folder with subfolders containing bam files
        fafile = '/fast/work/groups/ag_ohler/dharnet_m/Splicing_Lausanne/ext_data/annotation/gencode.v24lift37.annotation.orfext.fa',
        gtf = '../Liuetal_pipeline/pipeline/gencode.v24lift37.annotation.gtf',
        riboexprfolder = '../Liuetal_pipeline/pipeline/ribotrans_process/',
        FLANKCODS = 14, 
        minreadlen=25,
        maxreadlen=35,
        offsetfile='../Liuetal_pipeline/pipeline/riboseqc/data/ribo_0h/_P_sites_calcs',
        outputfolder='../Liuetal_pipeline/pipeline/metacodon_offset_analysis'
        foo = 'bar'

))
# Turn arguments into R variables
keys <- attachLocally(args)
cat("Command-line arguments attached to global environment:\n");
print(keys);
str(mget(keys,  envir = globalenv()))

base::source(here::here(paste0('../cortexomics/src/Rprofile.R')))
dir.create(here('data'),showWarn=FALSE)

dir.create(outputfolder,recursive=TRUE,showWarnings=FALSE)

riboexpr = riboexprfolder%>%
    {Sys.glob(paste0(., '*/ribotrans_expr.tr_expr.tsv'))}%>%
    setNames(., basename(dirname(.)))%>%
    map_df(.id = 'sample', fread)
riboexpr%<>%mutate(msample = sample%>%str_replace('ribo', 'r'))

highcountcovtrs = riboexpr%>%
    filter(cds_len > 100)%>%
    group_by(tr_id)%>%summarise(fdens = mean(RPF_dens))%>%
    arrange(desc(fdens))%>%
    slice(1:5000)%>%
    .$tr_id

# allbamtbls = Sys.glob(paste0('../Liuetal_pipeline/pipeline/deepshapebamdata/*.bam.reformat'))
# stopifnot(length(allbamtbls)>0)
# names(allbamtbls) = allbamtbls%>%str_extract('[^/]*?(? = .bam.ref)')
# mainbamtbls <- allbamtbls
# bamtbl = mainbamtbls[1]
# mainsamps = names(mainbamtbls)

library(riboWaltz)
annotation_dt <- create_annotation(gtf)
annotation_dt$transcript %<>%str_replace('\\.\\d+$', '')
stopifnot(gtf%>%str_detect('Liu.*'))#note this is for


seqs <- readDNAStringSet(fafile)
seqheaders <- names(seqs)
names(seqs) <- str_extract(names(seqs), '^[^|]+')
cdslens = seqheaders%>%str_match('CDS:(\\d+)\\-(\\d+)')%>%.[, -1]%>%apply(2, as.numeric)%>%{.[, 2]-.[, 1]+1}
seqlens = nchar(seqs)%>%setNames(names(seqs))
names(cdslens) = names(seqlens)

names(seqlens)%<>%str_extract('^[^|]+')

annotation_dt <- data.table(transcript = names(seqlens), l_utr5 = 60,
    l_cds = cdslens, l_utr3 = 57)%>%
    mutate(l_tr = l_utr5 + l_cds + l_utr3)
hctr_lens = annotation_dt$l_tr[match(highcountcovtrs, annotation_dt$transcript)]

# testtrs = annotation_dt

bfolders = Sys.glob(paste0(bamparentfolder, "/*/"))
names(bfolders) = basename(bfolders)
stopifnot(!anyDuplicated(names(bfolders)))

#get allcodlist granges object descxribing codon positions in the transcripts
message('reading coverage from bam files')
if (!file.exists(here("data/fpcovlist.rds"))) {
    fpcovlist  <-  lapply(bfolders, function(bfolder) {
        reads_list <- bamtolist(
                bamfolder = bfolder,
                annotation = annotation_dt)
        stopifnot(length(reads_list) == 1)
        reads_list = reads_list[[1]]
        GRanges(reads_list$transcript, IRanges(reads_list$end5, w = 1),
            readlen = reads_list$length) %>%
            subset(seqnames %in% highcountcovtrs) %>%
            subset(between(readlen, 25,35)) %>%
            {seqlevels(.) <- highcountcovtrs ;.} %>%
            {seqlengths(.) <- hctr_lens ;.} %>%
            {split(.,.$readlen)} %>%
            lapply(coverage)
    })
    fpcovlist%<>%setNames(names(bfolder))
    saveRDS(fpcovlist, here("data/fpcovlist.rds"))
}else{
    fpcovlist <- readRDS(here("data/fpcovlist.rds"))
    stopifnot(allcodlist@seqinfo@seqnames %>% setequal(highcountcovtrs))
}

################################################################################
########Verify offsets with metacodon plots
################################################################################


i <- 1
allcodons <- names(Biostrings::GENETIC_CODE) %>% setNames(., .)
# seqlengths <- GenomicRanges::seqlengths
.='foo'

trspacecds = annotation_dt%>%
    {GRanges(.$transcript,IRanges(
        .$l_utr5 + 1,
        .$l_utr5 + .$l_cds
    ))}%>%setNames(.,seqnames(.))

get_codon_gr <- function(codon, seqs, trspacecds,
    startbuff=60, endbuff=60, flankcods=FLANKCODS) {
    message(codon)
    #exclude the start ccodon
    codmatches <- Biostrings::vmatchPattern(pattern = codon, seqs)
    #
    matchgr <- codmatches %>% unlist %>% GRanges(names(.), .)
    cdsstarts <- start(trspacecds[as.vector(GenomicRanges::seqnames(matchgr))])
    matchgr$cdspos <- start(matchgr) - cdsstarts
    matchgr %<>% subset(cdspos %% 3 == 0)
    seqlengths(matchgr) <- setNames(nchar(seqs),names(seqs)) %>%
        .[seqlevels(matchgr)]
    innercds <- trspacecds %>%
        subset(width > (3 + startbuff + endbuff)) %>%
        resize(width(.) - startbuff, "end") %>%
        resize(width(.) - endbuff, "start")
    matchgr <- matchgr %>% subsetByOverlaps(innercds)
    codmatchwindows <- matchgr %>%
        resize(width(.) + (2 * (3 * flankcods)), "center")
    codmatchwindows <- codmatchwindows[!is_out_of_bounds(codmatchwindows)]
    codmatchwindows %<>% subsetByOverlaps(innercds)
    codmatchwindows
}

#get allcodlist granges object descxribing codon positions in the transcripts
if (!file.exists(here("data/allcodlist.rds"))) {
    allcodlist <- lapply(allcodons, F = get_codon_gr,
        seqs = seqs[highcountcovtrs],
        trspacecds = trspacecds)
    allcodlist <- allcodlist %>% GRangesList %>% unlist
    saveRDS(allcodlist, here("data/allcodlist.rds"))
}else{
    allcodlist <- readRDS(here("data/allcodlist.rds"))
    stopifnot(allcodlist@seqinfo@seqnames %>% setequal(highcountcovtrs))
}


if (!file.exists(here("data/fprustprofilelist.rds"))) {
    message('collecting RUST metacodon profiles')
    fprustprofilelist <- imap(fpcovlist ,function(sampfpcov, sampname){
        #sum over counts for that transcript
        trsums <- sampfpcov %>% map(sum) %>% purrr::reduce(., `+`)
        sampfpcov %>% lapply(function(rlfpcov) {
            rlfpcov <- rlfpcov > mean(rlfpcov)
            nztrnames <- names(trsums)[trsums!=0]
            allcodlistnz <- allcodlist %>% subset(seqnames %in% nztrnames)
            cods <- names(allcodlistnz) %>% str_split("\\.") %>% map_chr(1)
            message('.')
            out <- rlfpcov[allcodlistnz] %>%
                split(cods) %>%
                lapply(as.matrix) %>%
                map(colMeans)
            out
        })
    })
    saveRDS(fprustprofilelist,here('data/fprustprofilelist.rds'))
}else{
    fprustprofilelist<-readRDS(here('data/fprustprofilelist.rds'))
}

rustprofiledat = fprustprofilelist%>%
        map_depth(3,.%>%enframe('position','count'))%>%
        map_df(.id='sample',.%>%
            map_df(.id='length',.%>%
                bind_rows(.id='codon')
            )
        )
rustprofiledat%<>%mutate(position = position - 1 - (FLANKCODS*3))
rustprofiledat%<>%group_by(sample,length,codon)%>%mutate(count= count / median(count))
rustprofiledat%<>%filter(!codon %in% c('TAG','TAA','TGA'))
rustprofiledat$length%<>%as.numeric


################################################################################
########testing the pca based a-site calls
################################################################################

get_metacodon_var_offsets<-function(rustprofiledat){
    profdat = rustprofiledat%>%select(position,length,sample,count)
    profdat$length= profdat$length%>%as.numeric
    profdat$phase = (-profdat$position)%%3
    vardf = profdat%>%
        ungroup%>%
        group_by(sample,length,position,phase)%>%
        summarise(sdsig=sd(count,na.rm=T)/median(count,na.rm=T))
    vardf <- 
        vardf%>%
        group_by(sample,length,phase)%>%
        arrange(position)%>%
        mutate(sdsigpair = sdsig+lag(sdsig))
    vardf <- vardf%>%
        mutate(ismode=(sdsig>lag(sdsig)) & ((sdsig)<lead(sdsig)))%>%
        filter(position> -length+6,position <  -6)
    bestmode = vardf%>%group_by(sample,phase)%>%slice(which.max(sdsigpair))
    vardf%>%
        group_by(sample,length,phase)%>%
        slice(which.max(sdsigpair))%>%
        mutate(offset=-position)%>%
        select(length,phase,offset)
}
varoffsets <- get_metacodon_var_offsets(rustprofiledat)

#now plot the variation in occupancy at each position
#for each readlength
{
plotfile <- paste0(outfolder,'/fppos_vs_codon_variance.pdf')
grDevices::pdf(plotfile,w=12,h=12)
#plotting variance amongst codons at each point.
rustprofiledat%>%
    ungroup%>%
    group_by(sample,length,position)%>%
    filter(!is.nan(count))%>%
     mutate(phase = -position %%3)%>%
    summarise(sdsig=sd(count,na.rm=T)/median(count,na.rm=T))%>%
    # separate(sample,c('time','assay','rep'))%>%
    group_by(length,sample,position)%>%
    summarise(sdsig=mean(sdsig))%>%
    # filter(position> -numreadlen+6,position < -6)%>%
    filter(position> -length+(-6),position < - (-6))%>%
    filter(length>=minreadlen,length<=maxreadlen)%>%
    arrange(position)%>%
    {
        qplot(data=.,x=position,y=sdsig)+
        theme_bw()+
        facet_grid(length~sample,scale='free_y')+
        scale_y_continuous('between codon variation (meannorm)')+
        scale_x_continuous('5 read position relative to codon ')+
        geom_vline(data=varoffsets,aes(xintercept = -offset),
            color=I('blue'),linetype=2)+
        geom_vline(data=varoffsets%>%mutate(offset = offset+3),
            aes(xintercept= -offset),
            color=I('red'),linetype=2)+
        ggtitle("variance of 5' read occurance vs position")
    }%>%print
dev.off()
normalizePath(plotfile)
}

varoffsets%>%write_tsv(outfolder,'/variance_offsets.txt')

################################################################################
########Can also plot the pca to maybe seperate P and A site?
################################################################################

newcodonoccs = rustprofiledat%>%
    left_join(offsets)%>%
    filter(position == -offset-3)%>%
    identity
newcodonoccs%<>%separate(sample,c('time','assay','rep'))
newcodonoccs%<>%group_by(codon,time)%>%summarise(count=mean(count))

getpca<-function(pcadat, num=1){
    pcadat%>%
    princomp%>%
    {.$loadings[,num] * .$sdev[num]}%>%
    enframe('position','score')
}

#
profvarpca = rustprofiledat%>%
    split(.,.$sample)%>%
    map_df(.id="sample",.%>%
        split(.,list(.$readlen))%>%
        map_df(.id="length", function(profdata) {
            pcadat = profdata%>%
                mutate(numreadlen=str_extract(readlen,'\\d+')%>%as.numeric)%>%
                filter(position> -numreadlen+6,position < -6)%>%
                ungroup%>%
                select(-numreadlen,-readlen,-sample)%>%
                spread(position,count)%>%
                {set_rownames(.[,-1],.$codon)}
            bind_rows(.id='pca',
                pca1=pcadat%>%getpca(1),
                pca2=pcadat%>%getpca(2),
                pca3=pcadat%>%getpca(3)
            )
        })
    )
#
profvarpca%<>%select(sample,length,position,pca,score)
#Now plot the pcas vs the offsets
plotfile<-paste0(outfolder,'fppos_vs_codon_pcascore.pdf')
grDevices::pdf(plotfile,w=12,h=12);print(
profvarpca%>%
    ggplot(data=.,aes(y=score,x=as.numeric(position),color=pca))+
    geom_point()+
    geom_line()+
    facet_grid(length~sample)
    geom_vline(data=varoffsets,aes(xintercept = -offset),
        color=I('blue'),linetype=2)+
    geom_vline(data=varoffsets%>%mutate(offset = offset+3),
        aes(xintercept= -offset),
        color=I('red'),linetype=2)+
);dev.off()
normalizePath(plotfile)
