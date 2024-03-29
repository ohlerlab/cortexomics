#options(repos = BiocManager::repositories());packrat::init(options = list(ignored.packages = c('sleuth','xtail','SaTAnn','RiboseQC','ORFquant','rseq','rseqdata','proDD','riboWaltz','colorout')))

try(silent=T,{library(colorout)})
library(Biostrings)
library(checkmate)
library(memoise)
library(assertthat)
library(stringr)
library(tidyverse)
library(magrittr)
library(checkmate)
# library(conflicted)
message('loading libraries')
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))
suppressMessages(library(data.table))
suppressMessages(library(DESeq2))
suppressMessages(library(assertthat))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))
suppressMessages(library(here))
suppressMessages(library(biomaRt))
suppressMessages(library(testthat))
library(zeallot)
library(splines)
library(GenomicRanges)
library(limma)
library(broom)
library(txtplot)
library(multitaper)

# #
# conflict_prefer('setdiff','BiocGenerics')
# conflict_prefer('rowMedians','Biobase')
# conflict_prefer('setequal','S4Vectors')
# conflict_prefer("between", "dplyr")
# conflict_prefer("intersect", "BiocGenerics")
# conflict_prefer("lag", "dplyr")
matches <- dplyr::matches
filter<-dplyr::filter
select<-dplyr::select
slice<-dplyr::slice
qs<-checkmate::qassert


# ###memoise
project_cache=here::here('R_cache')
if(!exists('project_cache'))project_cache=tempdir()
message(system(str_interp('du -h ${project_cache}'),intern=T))
mycache=memoise::cache_filesystem(project_cache)

myclearcache=function() system(str_interp('rm -rf ${project_cache}'))


mymemoise <- function(f){
  if(!is.memoised(f)){
    memoise(f,cache=mycache)
    } else{ 
      f
  }
}
projmemoise<-mymemoise
addfileinf <- function(file){
  attr(file,'fileinfo')<-file.info(file)
  file
}

if(!interactive()) mymemoise=identity
  gigsused <- function(x)system(paste0("cat /proc/",Sys.getpid(),"/status | grep VmSize"),intern=TRUE)%>%str_extract('\\d+')%>%as.numeric%>%divide_by(1e6)
  # message('memory in use ',gigsused())


# rm(foomat)
# gc(reset=TRUE,full=TRUE)
# message('memory in use ',gigsused())


# #
# foo<-function(x){message('foooooobar called')}
# myf <- mymemoise(foo)
# myf(2)

# Q
#  2 %>% mymemoise(function(x,.foo=foo){.foo();message('foooo'); x +1})(.)
# # e <- as.call(c(as.name("{"),quote(message('foo_ins')),body(foo)[-1]))
# body(foo) <- e

# foo

# mymemoise(function(x,.foo=foo)

safe_filter <- function(...){
  filtered = filter(...)
  assert_that(nrow(filtered)>0)
  filtered
}

#' safe_left_join
#' @description left join that fails if a row in x is either duplicated or
#'   unmatched.
#' @param x table to join
#' @param y table to join
#' @param by a character vector of column names to join by.
#' @param verbose Default is TRUE.
#' @export


safe_left_join = function (x, y, by = NULL, verbose = TRUE,allow_missing=FALSE,allow_dups=FALSE) {
  rows_start = nrow(x)

  if (is.null(by)) {
    by = intersect(names(x), names(y))
  } else {
    by = as.character(by)
  }

  y[["..1.."]] = 1
  x = left_join(x, y, by)

  if(!allow_dups){
    if (nrow(x) > rows_start) {
      stop("Rows have been duplicated in 'safe' left join")
    }
  }
  if(!allow_missing){
    if (any(ind <- is.na(x[["..1.."]]))) {
      sample = sample(which(ind), min(10, sum(ind)))
      examples = distinct(x[sample, by, drop = FALSE])
      if (verbose) print(examples)
      stop(sprintf("Failed to match %d rows in x.", sum(ind)))
    }
  }

  x[["..1.."]] = NULL

  x

}

is_offchr<-function(gr,si){
  if(is(gr,'GenomicRangesList')){
   (end(gr) > split(seqlengths(gr)[as.character(unlist(seqnames(gr)))],gr@partitioning) ) %in% TRUE
  }else{
    seqinfo(gr)<-si
    end(gr) > seqlengths(gr)[as.character(seqnames(gr))]
  }
}
is_out_of_bounds <- function(gr,si = seqinfo(gr)){
  start(gr)<1 | is_offchr(gr,si) 
}


get_all_obsizes <- function(){.GlobalEnv%>%names%>%discard(~is.function(get(.)))%>%setNames(.,.)%>%map(~object.size(get(.)))}

# allobjsizes<-get_all_obsizes()
# allobjsizes%<>%enframe
# allobjsizes$value%<>%unlist
# allobjsizes$value%<>%divide_by(1e6)
# allobjsizes%>%arrange(desc(value))
# obsizes <- get_all_obsizes()

# obsizes%>%enframe('object','size')%>%mutate(size=as.numeric(size)/1e6)%>%arrange(desc(size))


purely <- function(fun,throw_error=TRUE,allow_functions=FALSE){
  
  funname <- rlang::quo_text(enquo(fun))

  globalobs <- .GlobalEnv%>%names

  globalobs <-   discard(globalobs,~is.function(get(.,envir=.GlobalEnv)) && (!identical(environment(get(.,envir=.GlobalEnv)),.GlobalEnv )))

  if(!allow_functions){
    globalobs <-   discard(globalobs,~is.function(get(.,envir=.GlobalEnv)) )
  }
 

  keep(globalobs,~is.function(get(.,envir=.GlobalEnv)) && (!identical(environment(get(.,envir=.GlobalEnv)),.GlobalEnv )))

  assert_that(! '%>%' %in% globalobs)

  funargs <- names(formals(fun))

  environment(fun) <- new.env()

  if(throw_error) messagefun = stop else messagefun = warning

  for(obnm in globalobs) {

    warnenv<-new.env()
    warnenv[['obnm']]<-obnm
    warnenv[['val']]<- force(get(obnm,envir=.GlobalEnv))
    delayedAssign(obnm,
      {
        messagefun(force(paste0('Object: ',obnm,' is being evaluated, but is not a formal argument of ',funname)));
        val
      },
      eval.env=warnenv,
      assign.env=environment(fun))
  }
  
  return(fun)
}

# myfun<-function(a){ a+1+b}

# a<-1
# b<-10
# myfun(a)
# pf<-purely(myfun)

# purely(myfun)(1)

# environment(pf)$get_frame_entropy



read_compressed_gfile <- function(annofile,annotype,fformat='gtf'){
  f=tempfile();
  stopifnot(file.exists(annofile))
  catbin = ifelse(tools::file_ext(annofile)=='gz','zcat','cat')
  system(str_interp('${catbin} ${annofile} | grep -e "\t${annotype}\t" > ${f}'));
  out = import(f,format=fformat) 
  file.remove(f)
  out
}
load_objs <- function(f){
    env <- new.env()
    nms <- load(f, env)
    map(nms,message)
    as.list(env)
}


extract_oneof <- function(strings,ids){
  
  matchlist <- map(strings,~str_extract(pattern = ids,string = .))
  
  matchnum <- matchlist%>%map(~sum(!is.na(.)))
  
  stopifnot(all(matchnum < 2 ))
  stopifnot(all(matchnum > 0))

  matches <- matchlist%>%map_chr(keep,Negate(is.na))

  matches
}

extract_id <- function(strings,ids){
	
	matchlist <- map(strings,~str_extract(pattern = sampleids,string = .))
	
	matchnum <- matchlist%>%map(~sum(!is.na(.)))
	stopifnot(all(matchnum < 2 ))
	stopifnot(all(matchnum > 0))

	matches <- matchlist%>%map_chr(keep,Negate(is.na))

	matches
}



read_compressed_gfile <- function(annofile,annotype,fformat='gtf'){
  f=tempfile();
  stopifnot(file.exists(annofile))
  catbin = ifelse(tools::file_ext(annofile)=='gz','zcat','cat')
  system(str_interp('${catbin} ${annofile} | grep -e "\t${annotype}\t" > ${f}'));
  out = import(f,format=fformat) 
  file.remove(f)
  out
}
load_objs <- function(f){
    env <- new.env()
    nms <- load(f, env)
    map(nms,message)
    as.list(env)
}


extract_oneof <- function(strings,ids){
  
  matchlist <- map(strings,~str_extract(pattern = ids,string = .))
  
  matchnum <- matchlist%>%map(~sum(!is.na(.)))
  
  stopifnot(all(matchnum < 2 ))
  stopifnot(all(matchnum > 0))

  matches <- matchlist%>%map_chr(keep,Negate(is.na))

  matches
}

DT2GR = function(dt,seqinf=si,checksi=TRUE){

  if(is(dt,'GenomicRanges')) {
    warning('already a GRanges Object')
    return(dt)
  }


  stopifnot(c('seqnames','start')%in%colnames(dt))
  stopifnot(c('width')%in%colnames(dt)|c('end')%in%colnames(dt))
  if(checksi){stopifnot(dt[['seqnames']] %in% seqlevels(seqinf))
  }else{seqinf=NULL}
  
  hasend=FALSE
  haswidth=FALSE

  if('end' %in% colnames(dt) ){
    stopifnot (dt[['end']] %>% `>`(0) %>%all)
    hasend=TRUE
  }
  if('width' %in% colnames(dt) ){
    stopifnot (dt[['width']] %>% `>`(0) %>%all)
    haswidth=TRUE
  }
  
  stopifnot(dt[['start']] %>% is.numeric)
  stopifnot(hasend|haswidth )
  
  if(haswidth & ! hasend ){
    dt[['end']]  = dt[['start']]+dt[['width']]-1 
  } 
  if(hasend ){

  } 

  #strand
  if(! 'strand' %in% colnames(dt)){
    dt[['strand']] ='*'
  }

  stopifnot(dt[['strand']] %in% c('+','-','*'))
  



  mdatcols = colnames(dt) %>% setdiff(c('seqnames','start','width','strand','end')) 
  #create basic granges object
  if(checksi){
    gr=GRanges(dt[['seqnames']],IRanges(dt[['start']],dt[['end']]),strand=dt[['strand']],seqinfo=seqinf)
  }else{    gr=GRanges(dt[['seqnames']],IRanges(dt[['start']],dt[['end']]),strand=dt[['strand']])}

  #add in the metadata if we need to
  if(length(mdatcols)){
    if(is.data.table(dt)){ mcols(gr) = dt[,mdatcols,with=FALSE]%>%as("DataFrame")
    }else{ mcols(gr) = dt[,mdatcols]%>%as("DataFrame")}
  }

    stopifnot(all(colnames(dt) %in% c(colnames(mcols(gr)),'seqnames','start','end','width' ,'strand')))

  gr
}


mergeGR2DT = function(mergegr){
  grcols = mergegr%>%vapply(is,TRUE,'GRanges')%>%which
  mergedt=cbind(
    mergegr[[grcols[[1]]]]%>%GR2DT,
    mergegr[[grcols[[2]]]]%>%GR2DT
  )
  cnames = colnames(mergedt)
  cnames[duplicated(cnames)]%<>%paste0('.1')
  setnames(mergedt,cnames)
  mergedt
}

GR2DT = function(gr){
  if(is.data.table(gr)) {
    warning('already a data table')
    return(gr)
  }
  #all columns must be vectors
  for(i in seq_len(ncol(mcols(gr)))){
    mcols(gr)[[i]] = as.vector(mcols(gr)[[i]])
  }

  dt = as.data.frame(gr,row.names=NULL,optional=FALSE)%>%as.data.table
  dt$strand= as.character(dt$strand)
  # setkey(dt,seqnames,strand,start,end)

  stopifnot( all(colnames( mcols(gr)) %in% colnames(dt) )  )
  dt
}



load_objs <- function(f){
    env <- new.env()
    nms <- load(f, env)
    map(nms,message)
    as.list(env)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


#very similiar to above but reutrns the spec value as welll as the ftest
take_Fvals_spect<-function(x,n_tapers,time_bw,slepians_values){
     if(length(x)<25){
          remain<-50-length(x)
          x<-c(rep(0,as.integer(remain/2)),x,rep(0,remain%%2+as.integer(remain/2)))
     }
     if(length(x)<1024/2){padding<-1024}
     if(length(x)>=1024/2){padding<-"default"}

     #this calculates a frequency object from a numeric vector using the 
     resSpec1 <- spec.mtm(as.ts(x), k=n_tapers, nw=time_bw, nFFT = padding, centreWithSlepians = TRUE, Ftest = TRUE, maxAdaptiveIterations = 100,returnZeroFreq=F,plot=F,dpssIN=slepians_values)

     #this selects only the frequency range surrounding 0.33 - the one we are concerned with - but isn't used
     resSpec2<-dropFreqs(resSpec1,0.29,0.39)
     
     freq_max_3nt<-resSpec1$freq[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]
     
     Fmax_3nt<-resSpec1$mtm$Ftest[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]

     spect_3nt<-resSpec1$spec[which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))]

     return(c(Fmax_3nt,spect_3nt))
     
}



read_compressed_gfile <- function(annofile,annotype,fformat='gtf'){
  f=tempfile();
  stopifnot(file.exists(annofile))
  catbin = ifelse(tools::file_ext(annofile)=='gz','zcat','cat')
  system(str_interp('${catbin} ${annofile} | grep -e "\t${annotype}\t" > ${f}'));
  out = import(f,format=fformat) 
  file.remove(f)
  out
}


get_riboprofdata<- function(exons_tr,bigwigpair){
  
  message( names(bigwigpair))
  stopifnot(c('+','-') %in% names(bigwigpair))

  for(i in ls()) assign(i,get(i),envir=.GlobalEnv)



  profilegrange <- 
    suppressWarnings(lapply(bigwigpair,import,which = unlist(exons_tr)))%>%
    {for(strandi in names(.)) strand(.[[strandi]]) <- strandi;.}%>%
    Reduce(f=c,.)%>%
    subset(score>0)

  seqlevs <-list(profilegrange,exons_tr)%>%unique%>%as.character

  shared_seqinfo <- suppressWarnings(intersect(seqinfo(BigWigFile(bigwigpair[[1]])),seqinfo(exons_tr)))
  
  trseqinfo <- Seqinfo(seqnames=names(exons_tr),seqlengths=as.vector(sum(width(exons_tr))))

  #now map our prifle data to the exons space
  rle <- suppressWarnings(mapToTranscripts(profilegrange,exons_tr))
  rle$score<-profilegrange$score[rle$xHits];
  seqinfo(rle)<-trseqinfo[seqinfo(rle)@seqnames]
  rle <- coverage(rle, weight='score')

  rle %>%
    #selecting the transcirpt we want
    lapply(FUN=.%>%
    #turn that into a dataframe
      { 
        pos = which(.!=0)
        data_frame(pos=pos, score = as.vector(.[pos]))
      }
    )%>%
    bind_rows(.id='tid')%>%
    mutate(frame=as.factor(pos %% 3))
}



# reads_tr<-oldenv$cdsread_trmap
fp<-function(gr) ifelse(strand(gr)=='-',end(gr),start(gr))
tp<-function(gr) ifelse(strand(gr)=='-',start(gr),end(gr))
fpend<-function(x)resize(x,1,'start')
tpend<-function(x)resize(x,1,'end')

# gr1<-GRanges(c('chr1:5-6:+'))
# gr2<-GRanges(c('chr1:50-51:+','chr1:40-51:+'))

downstream_dist_till<-function(gr1,gr2){
  (fp(gr2)[precede(gr1,gr2)] - tp(gr1)) * ( ifelse(strand(gr1)=='+',1,-1))
}

upstream_dist_till<-function(gr1,gr2){
  (tp(gr2)[follow(gr1,gr2)] - fp(gr1)) * ( ifelse(strand(gr1)=='+',-1,1))
}

istpmost<-function(cds,groupvar='transcript_id'){ 
  ids <- seq_along(cds)
  tpmostids<-data_frame(id=ids,end=end(cds),strand=as.vector(strand(cds)),groupvar=mcols(cds)[[groupvar]])%>%group_by(groupvar)%>%slice(which.max(end*ifelse(strand=='-',-1,1)))%>%.$id
  ids %in% tpmostids
}
isfpmost<-function(cds,groupvar='transcript_id'){
  ids <- seq_along(cds)
  fpmostids<-data_frame(id=ids,start=start(cds),strand=as.vector(strand(cds)),groupvar=mcols(cds)[[groupvar]])%>%group_by(groupvar)%>%slice(which.max(start*ifelse(strand=='-',1,-1)))%>%.$id
  ids %in% fpmostids
}



clip_start <- function(x,n) resize(x,width(x)-n,fix='end')
clip_end <- function(x,n) resize(x,width(x)-n,fix='start')
setstrand<-function(x) {strand(x)<-Rle('+') 
  x}


# BiocManager::install('GenomicRanges')
testthat::expect_equal(downstream_dist_till(GRanges(c('chr1:10:-')),GRanges(c('chr1:5:-','chr1:40-51:+'))),5)
testthat::expect_equal(downstream_dist_till(GRanges(c('chr1:10:+')),GRanges(c('chr1:14:+','chr1:40-51:+'))),4)

testthat::expect_equal(upstream_dist_till(GRanges(c('chr1:10:-')),GRanges(c('chr1:11:+','chr1:12:-','chr1:40-51:+'))),2)

testthat::expect_equal(upstream_dist_till(GRanges(c('chr1:10:+')),GRanges(c('chr1:11:+','chr1:6:+','chr1:40-51:+'))),4)
testthat::expect_equal(upstream_dist_till(GRanges(c('chr1:10:-')),GRanges(c('chr1:11:+','chr1:14:-','chr1:40-51:+'))),4)




apply_cds_offsets<- function(reads_tr,bestscores){

  #if you pass null to this, we can just take the middle of each read
  if(is.null(bestscores)) return(floor(width(reads_tr)/2))

  compartmentcol <- if('compartment'%in%colnames(bestscores)) 'compartment' else NULL
  phasecol <- if('phase'%in%colnames(bestscores)) 'phase' else NULL
  joincols <- c('length',compartmentcol,phasecol)
  assert_that(all(has_name(mcols(reads_tr),joincols)))

  reads_tr%>%
    mcols%>%
    as_tibble%>%
    safe_left_join(bestscores%>%
      ungroup%>%
      # select(-phase,-score),by=c('length','compartment'),
      select(offset,one_of(joincols)),by=joincols,
      allow_missing=TRUE
    )%>%
    .$offset
}



take_Fvals_spect<-function(x,n_tapers,time_bw,slepians_values){
  if(length(x)<25){
      remain<-50-length(x)
      x<-c(rep(0,as.integer(remain/2)),x,rep(0,remain%%2+as.integer(remain/2)))
  }
  #
  if(length(x)<1024/2){padding<-1024}
  if(length(x)>=1024/2){padding<-"default"}
  #
  resSpec1 <- spec.mtm(as.ts(x), k=n_tapers, nw=time_bw, nFFT = padding, centreWithSlepians = TRUE, Ftest = TRUE, maxAdaptiveIterations = 100,returnZeroFreq=F,plot=F,dpssIN=slepians_values)
  #
  resSpec2<-dropFreqs(resSpec1,0.29,0.39)
  #
  closestfreqind <- which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))
  # 
  freq_max_3nt<-resSpec1$freq[closestfreqind]
  Fmax_3nt<-resSpec1$mtm$Ftest[closestfreqind]
  spect_3nt<-resSpec1$spec[closestfreqind]
  return(c(Fmax_3nt,spect_3nt))
}

ftestvect<-function(psit,k=24,bw=12){
  sl<-dpss(n=length(psit)%>%ifelse(.<25,50,.),k=k,nw=bw)
  vals<-take_Fvals_spect(x = psit,n_tapers = k,time_bw = bw,slepians_values = sl)
  pval <- pf(q=vals[1],df1=2,df2=(2*24)-2,lower.tail=F)
  return(c(vals[2],pval))
}


#could include the lengths and variables it needs
#could include a method for getting the data
#and a method for shiftng a granges object
get_predictedshifts<-function(seqshiftmodel,data4readshift,probcutoff=0.5){
  require(ranger)
  outpred <- predict(seqshiftmodel,data=data4readshift)$prediction

  if(is.matrix(outpred)){
    maxcol <- max.col(outpred)
    is_highprob <- outpred[matrix(c(1:nrow(outpred),maxcol),ncol=2)] > probcutoff
    outpred=ifelse(is_highprob,colnames(outpred)[maxcol],0)
  }

  outpred <- outpred%>%as.character%>%as.numeric

  assert_that(all((outpred %% 1)==0))

  outpred

}

get_cds_offsets = function(reads_tr,offsets,compartments){

    compartmentcol <- if('compartment'%in%colnames(offsets)) 'compartment' else NULL

    if(!'compartment' %in% colnames(mcols(reads_tr))){
      mcols(reads_tr)$compartment <- as(compartments,"List")[seqnames(reads_tr)]%>%unlist(use.names=F)
    }
    if(!'length' %in% colnames(mcols(reads_tr))){
      mcols(reads_tr)$length <- qwidth(reads_tr)
    }


    phasecol <- if('phase'%in%colnames(offsets)) 'phase' else NULL
    joincols <- c('length',compartmentcol,phasecol)


  assert_that(all(has_name(mcols(reads_tr),joincols)))
  

  reads_tr%>%
    mcols%>%
    as_tibble%>%
    safe_left_join(offsets%>%
      ungroup%>%
      # select(-phase,-score),by=c('length','compartment'),
      select(offset,one_of(joincols)),by=joincols,
      allow_missing=TRUE
    )%>%
    .$offset
}

setGeneric('apply_psite_offset',function(offsetreads,offset) strandshift(offsetreads,offset))
setMethod('apply_psite_offset','GAlignments',function(offsetreads,offset){
  if(is.character(offset)){
    offset = rowSums(as.matrix(mcols(offsetreads)[,offset]))
  } 
  if(length(offset)==1) offset = rep(offset,length(offsetreads))
  isneg <-  as.logical(strand(offsetreads)=='-')
  offsetreads[!isneg] <- qnarrow(offsetreads[!isneg],start=offset[!isneg]+1,end = offset[!isneg]+1)
  ends <- qwidth(offsetreads[isneg])-offset[isneg]
  offsetreads[isneg] <- qnarrow(offsetreads[isneg],start=ends,end = ends  ) 
  offsetreads
})

fp <-function(gr)ifelse(strand(gr)=='-',end(gr),start(gr))
tp <-function(gr)ifelse(strand(gr)=='-',start(gr),end(gr))
strandshift<-function(gr,shift) shift(gr , ifelse( strand(gr)=='-',- shift,shift))



get_seqforrest_data <- function(trainreads,seq,nbp=2,trim=TRUE){

  stopifnot('length' %in% colnames(mcols(trainreads)))
  stopifnot(all(seqnames(trainreads)%in%seqinfo(seq)@seqnames))

  fp<-trainreads%>%as("GRanges")%>%resize(nbp,'start')%>%resize(nbp*2,'end')

  tp <- trainreads%>%as("GRanges")%>%resize(nbp,'end')%>%resize(nbp*2,'start')

  
  startseq <- getSeq(seq,fp)
  endseq <- getSeq(seq,tp)

  seqmat <- cbind(dnaseq2onehot(startseq,'fp.'),dnaseq2onehot(endseq,'tp.'))

  seqmat%<>%cbind(length=mcols(trainreads)$length)%>%as.data.frame

  if(has_name(mcols(trainreads),'dist')){
    seqmat %<>% cbind(dist=as.numeric(as.character(trainreads$dist)))
    #note that we need to negate the phase to get it from the 'dist', so e.g. -1 means the 2nd base of the stop,
    #while 1 means the 3rd base of the one before AUG
    seqmat$phase = as.factor(abs(-(as.numeric(as.character((seqmat)$dist))%%3)))
  }
  
  seqmat

}



#' group_slice
#' 
#' get N groups.
#' #' 
#' @param dt - a grouped tbl
#' @param v - a vector for slicing
#' @examples e.g. - an example
#' this 
#' @return returns 
#' @export 
group_slice<-function(dt,v){
  stopifnot('data.frame'%in%class(dt))
  stopifnot(v>0,(v%%1)==0)
  if(n_groups(dt)==1) warning('Only 1 group to slice')
  #get the grouping variables
  groups=
    unique(dt%>%select())%>%
    group_by()%>%
    slice(v)
  out=inner_join(dt,groups)
  out
}



#' slice by - common operation to get the first group according ot a set of variables
#' 
#' get N groups.
#' #' 
#' @param dt - a grouped tbl
#' @param v - a vector for slicing
#' @examples e.g. - an example
#' this 
#' @return returns 
#' @export 

slice_by<-function(dt,...,n=1){
  dots = enquos(...)
  n= eval_tidy(dots[[length(dots)]])
  dots = dots[-length(dots)]
  ogrps = group_vars(dt)
  ogrpsq<-syms(ogrps)
  dt%>%ungroup%>%inner_join(ungroup(dt)%>%distinct(!!!dots)%>%slice(n),by=map_chr(dots,quo_text))%>%group_by(!!!ogrpsq)
}





make_lfc_trajplot<-function(lfc_tbl_plot){
  
  assert_that(all(  lfc_tbl_plot %has_name% c('log2fc','lmax','lmin','gene_name','assay','time')))
  assert_that(nrow(lfc_tbl_plot)!=0)

  lfc_tbl_plot %<>% filter(time %>% {.==unique(.)[1]}) %>%mutate(log2fc=0,time=factor('E13',levels=levels(time)))%>%rbind(lfc_tbl_plot)

  opac <- 1

    ggplot(data=lfc_tbl_plot,aes(y=log2fc,x = as.numeric(time),group=gene_name,color=class))+
      geom_line(data=lfc_tbl_plot%>%filter(class%in%c('other','No sig TE change')),position='identity',alpha = I(1)) +
      geom_line(data=lfc_tbl_plot%>%filter(!class%in%c('No sig TE change')),position='identity',alpha=I(opac)) +
      # geom_line(position='identity',alpha='0.5') +
      # expand_limits(y=c(floor(min(lfc_tbl_plot$log2fc)),ceiling(max(lfc_tbl_plot$log2fc))))+
      # coord_cartesian(y=c(-3,3))+
      ggtitle(str_interp("TE trajectory - all genes"))+
      theme_bw()+
      scale_color_manual(values=class_color_dict)+
      guides(colour = guide_legend(override.aes = list(alpha = 1)))+
      scale_x_continuous(labels=levels(time),name='Stage')+
      theme(text = element_text(size = 16))

      # geom_ribbon(data=lfc_tbl_plot,aes(ymax=lmax,ymin=lmin,x=ntime),alpha=0.5,fill=I('darkgreen'))+
        # geom_line(data=lfc_tbl_plot,aes(y=log2fc,x=ntime))+
        
}

inclusiontable<-function(a,b){
  library(rlang)
  aname = rlang::quo_text(enquo(a))
  bname = rlang::quo_text(enquo(b))
  all = BiocGenerics::union(a,b)
  outtab = table(all %in% a,all %in% b)
  dimnames(outtab)%<>%setNames(c(aname,bname))
  outtab
}



#define methods so I can 'apply_psite_offset' either to GAlignments objets or just GRanges objects
setGeneric('apply_psite_offset',function(offsetreads,offset) strandshift(offsetreads,offset))
setMethod('apply_psite_offset','GAlignments',function(offsetreads,offset){
  if(is.character(offset)){
    offset = rowSums(as.matrix(mcols(offsetreads)[,offset]))
  } 
  if(length(offset)==1) offset = rep(offset,length(offsetreads))
  isneg <-  as.logical(strand(offsetreads)=='-')
  offsetreads[!isneg] <- qnarrow(offsetreads[!isneg],start=offset[!isneg]+1,end = offset[!isneg]+1)
  ends <- qwidth(offsetreads[isneg])-offset[isneg]
  offsetreads[isneg] <- qnarrow(offsetreads[isneg],start=ends,end = ends  ) 
  offsetreads
})
#function to get GRanges objects of the psites
get_genomic_psites <- function(bam,windows,offsets,mapqthresh=200,comps=c('chrM'='chrM')) {
  require(GenomicAlignments)
  riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=windows)
  reads <- readGAlignments(bam,param=riboparam)
  mcols(reads)$length <- qwidth(reads)
  #
  #if no 'offsets' data then just take read midpoints (e.g. for RNAseq)
  if(is.null(offsets)){
    mcols(reads)$offset <- floor(qwidth(reads)/2)
  }else{
    #else check we have offsets
  stopifnot('length' %in% colnames(offsets))
  stopifnot('offset' %in% colnames(offsets))
  stopifnot('comp' %in% colnames(offsets))
    #
    #define our compartment for all of the seqnames in the object
    comps = comps[names(comps)%in%seqnames(reads)]
    useqnms <- as.character(unique(seqnames(reads)))%>%setdiff(names(comps))
    compmap = setNames(c(rep('nucl',length(useqnms)),comps),c(useqnms,names(comps)))
    stopifnot(all(compmap$values()%in%offsets$comp))
    #now fetch an offset for all of our 
    mcols(reads)$offset <-
      data.frame(length=mcols(reads)$length,
        comp=compmap[[as.character(seqnames(reads))]])%>%
      safe_left_join(offsets%>%select(offset,length,comp),allow_missing=TRUE,by=c('length','comp'))%>%.$offset
  }
  #get rid of the reads that have no offset
  reads <- reads%>%subset(!is.na(mcols(reads)$offset))
  psites <- apply_psite_offset(reads,c('offset'))%>%as("GRanges")
  #annotate psites with length
  mcols(psites)$length <- mcols(reads)$length
  psites
}

# mapqthresh=200
# riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=testregion%>%unlist)
# reads <- readGAlignments(bams[1],param=riboparam)



GR <- GRanges
ext_grl<-function(grl,nbp,fixend='end'){
  stopifnot(identical(grl,sort(grl)))
  if(is.null(names(grl))) names(grl)<-seq_along(grl)
  #if we fix the start we want to select the end
  if(fixend=='start'){
    strandtorev <- '-' 
    endinds <-  grl%>%{as(ifelse(any(strand(.)%in%strandtorev),1,elementNROWS(.)),'IntegerList')}
  }else if (fixend=='end'){
    #and if the end then we want the start
    strandtorev <- '+'
    endinds <-  grl%>%{as(ifelse(any(strand(.)%in%strandtorev),1,elementNROWS(.)),'IntegerList')}
  }else{stop('fixend needs to be start or end')}
  cuminds <- cumsum(elementNROWS(grl))

  ulgrl <- unlist(grl,use.names=TRUE)
  endinds <- unlist(lag(cuminds,1,0)+endinds)
  #but this isn't safe since the gene might end up out of bounds
  #add the difference between it and the limit, rounded down to the nearest 3, up to 15
  stopifnot(all(start(ulgrl[endinds]) > nbp))
  ulgrl[endinds] %<>% resize(width(.)+nbp,fixend)
  grl <- split(setNames(ulgrl,NULL),names(ulgrl))
  stopifnot(is(grl,'GRangesList'))
  grl
}

stanpars_to_list <- function(stanparvect){
  #get the dimensions of each of the parameters
  iparnames <- stanparvect%>%names
  parnames <- iparnames%>%str_remove('(\\[|\\.)[0-9\\,]+(\\]|\\.)$')
  parinds <- iparnames%>%str_extract('(\\[|\\.)[0-9\\,]+(\\]|\\.)$')%>%str_extract_all(regex('\\d+'))%>%map(as.numeric)
  parinds <- parinds%>%split(parnames)%>%map(.%>%simplify2array%>%t)
  stanlist <- parinds%>%map(.%>%apply(2,max)%>%replace_na(1)%>%array(NA,dim=.))
  #
  for(parname in unique(parnames)){
    stanlist[[parname]][parinds[[parname]]] <- stanparvect[parnames==parname]
  }
  stanlist
}

vals <- function(x){
  x[is.finite(x)]

}

GRanges(1:2,1:2)%>%split(1:2)

spl_mapFromTranscripts <- function(trspacegr,exons_grl){
  #
  exons_tr<-exons_grl%>%unlist%>%
    mapToTranscripts(exons_grl)%>%.[names(.)==seqnames(.)]
  ov <- findOverlaps(trspacegr,exons_tr)
  #
  trspacegr_spl <- suppressWarnings({
    trspacegr[queryHits(ov)]%>%
    pintersect(exons_tr[subjectHits(ov)])
  })
  genomic_trspacegr <- mapFromTranscripts(
  trspacegr_spl,
  # exons_tr[subjectHits(ov)]%>%split(.,seqnames(.))
  exons_grl
  )
  genomic_trspacegr$xHits <- queryHits(ov)[genomic_trspacegr$xHits]
  # mcols(genomic_trspacegr) <- mcols(trspacegr)[genomic_trspacegr$xHits]
  genomic_trspacegr
}

resize_grl_startfix<-function(grl,width){
  #what follows is some slightly black magic using S4 vectors
  #Integerlist which showings how much we'd need to trim that exon to get to to the desired transcript length
  trim =  cumsum(width(grl)) - width 
  #Where trim is greater than the exon width, we drop it
  drop = trim >=  width(grl)
  grl = grl[!drop]
  #vector showing location of the new 3' end of each transcript
  newends = cumsum(elementNROWS(grl))
  #vector with the amount we need to trim each new 3' end by
  endtrims=trim[IntegerList(as.list(elementNROWS(grl)))]@unlistData
  #finally, use these to trim
  grl@unlistData[newends] <- resize(grl@unlistData[newends], width(grl@unlistData[newends]) - endtrims  )
  grl
}
str_order_grl<-function(grl){order( start(grl)*(((strand(grl)!='-')+1)*2 -3) )}
sort_grl_st <- function(grl)grl[str_order_grl(grl),]
resize_grl_endfix <- function(grl,width){
  grl = invertStrand(grl)%>%sort_grl_st
  # 
  grl = resize_grl_startfix(grl,width)
  invertStrand(grl)%>%sort_grl_st
}
resize_grl <- function(grl,width,fix='start',check=TRUE){
  stopifnot(all(width>0))
  assert_that(all(all(diff(str_order_grl(grl))==1) ),msg = "grl needs to be 5'->3' sorted")
  if(fix=='start'){
    grl = resize_grl_startfix(grl,width)
  }else if(fix=='end'){
    grl = resize_grl_endfix(grl,width)
  }else if(fix=='center'){
    grlwidths = sum(width(grl)) 
    diffs = (width - grlwidths)
    # 
    grl = resize_grl_startfix(grl,grlwidths + ceiling(diffs/2))
    grl = resize_grl_endfix(grl,grlwidths + diffs)
  }
  if(check){
    startstoolow <- any(start(grl)<=0)
    if(any(startstoolow)){
      stop(str_interp("${sum(startstoolow)} ranges extended below 1 .. e.g. ${head(which(startstoolow,1))}"))
    }
    grlseqs <- as.vector(unlist(use.names=F,seqnames(grl)[IntegerList(as.list(rep(1,length(grl))))]))
    endstoohigh <- any((end(grl)>seqlengths(grl)[grlseqs])%in%TRUE)
    if(any(endstoohigh)){
      stop(str_interp("${sum(endstoohigh)} ranges extended below above seqlength .. e.g. ${head(which(endstoohigh,1))}"))
    }
  }
  grl
}
trim_grl <- function(grl,bp,end='tp'){
  if(end=='tp'){
    resize_grl(grl,sum(width(grl)) - bp,fix='start')
  }else if(end=='fp'){
    resize_grl(grl,sum(width(grl)) - bp,fix='end')
  }else {
    stop("end should be fp or tp")
  }
}

setGeneric('apply_psite_offset',function(offsetreads,offset) strandshift(offsetreads,offset))
setMethod('apply_psite_offset','GAlignments',function(offsetreads,offset){
  if(is.character(offset)){
    offset = rowSums(as.matrix(mcols(offsetreads)[,offset]))
  } 
  if(length(offset)==1) offset = rep(offset,length(offsetreads))
  isneg <-  as.logical(strand(offsetreads)=='-')
  offsetreads[!isneg] <- qnarrow(offsetreads[!isneg],start=offset[!isneg]+1,end = offset[!isneg]+1)
  ends <- qwidth(offsetreads[isneg])-offset[isneg]
  offsetreads[isneg] <- qnarrow(offsetreads[isneg],start=ends,end = ends  ) 
  offsetreads
})

fp <-function(gr)ifelse(strand(gr)=='-',end(gr),start(gr))
tp <-function(gr)ifelse(strand(gr)=='-',start(gr),end(gr))
strandshift<-function(gr,shift) shift(gr , ifelse( strand(gr)=='-',- shift,shift))


#define methods so I can 'apply_psite_offset' either to GAlignments objets or just GRanges objects
setGeneric('apply_psite_offset',function(offsetreads,offset) strandshift(offsetreads,offset))
setMethod('apply_psite_offset','GAlignments',function(offsetreads,offset){
  if(is.character(offset)){
    offset = rowSums(as.matrix(mcols(offsetreads)[,offset]))
  } 
  if(length(offset)==1) offset = rep(offset,length(offsetreads))
  isneg <-  as.logical(strand(offsetreads)=='-')
  offsetreads[!isneg] <- qnarrow(offsetreads[!isneg],start=offset[!isneg]+1,end = offset[!isneg]+1)
  ends <- qwidth(offsetreads[isneg])-offset[isneg]
  offsetreads[isneg] <- qnarrow(offsetreads[isneg],start=ends,end = ends  ) 
  offsetreads
})
#function to get GRanges objects of the psites
get_genomic_psites <- function(bam,windows,offsets,mapqthresh=200,comps=c('chrM'='chrM')) {
  require(GenomicAlignments)
  riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=windows)
  reads <- readGAlignments(bam,param=riboparam)
  # browser()
  mcols(reads)$length <- qwidth(reads)

  #if no 'offsets' data then just take read midpoints (e.g. for RNAseq)
  if(is.null(offsets)){
    mcols(reads)$offset <- floor(qwidth(reads)/2)
  }else{
    #else check we have offsets
  stopifnot('length' %in% colnames(offsets))
  stopifnot('offset' %in% colnames(offsets))
  stopifnot('comp' %in% colnames(offsets))

    #define our compartment for all of the seqnames in the object
    useqnms <- as.character(unique(seqnames(reads)))%>%setdiff(names(comps))
    compmap = setNames(c(rep('nucl',length(useqnms)),comps),c(useqnms,names(comps)),)
    stopifnot(all(compmap$values()%in%offsets$comp))
    #now fetch an offset for all of our 
    mcols(reads)$offset <-
      data.frame(length=mcols(reads)$length,
        comp=compmap[as.character(seqnames(reads))])%>%
      safe_left_join(offsets%>%select(offset,length,comp),allow_missing=TRUE,by=c('length','comp'))%>%.$offset
  }
  #get rid of the reads that have no offset
  reads <- reads%>%subset(!is.na(mcols(reads)$offset))
  psites <- apply_psite_offset(reads,c('offset'))%>%as("GRanges")
  #annotate psites with length
  mcols(psites)$length <- mcols(reads)$length
  psites
}

#devtools::install_git(c('https://github.com/cran/ifultools'))
#devtools::install_git(c('https://github.com/cran/wmtsa'))
DWPT         =       function(signal){
  library(wmtsa);
  TR_lth          =       length(signal);
  if(TR_lth       <       64){
          signal  =       c(signal,rep(0,64-length(signal)));
  }
  W1              =       wavMODWPT(signal, wavelet="s4",n.levels=6);
  W2              =       wavShift(W1);
  bands           =       32:51;# use the 0.2~0.5 Hz components only.
  mx              =       matrix(0,nrow=length(bands),ncol=length(signal));
  for(i in 1:length(bands)){
          tmp     =       paste("w6.",bands[i],sep="");
          mx[i,]  =       W2$data[[tmp]];
  }
  ID_signal       =       which(signal>0);        #the positions with signal;
  mx[,-ID_signal] =       0;                      #remove noise;
  minus3nt        =       mx[-11,];
  only3nt         =       mx[11,];
  BKgrnd          =       apply(minus3nt,2,max);
  ID1             =       which(only3nt > BKgrnd);#the positions with 3nt energy higher than other frequency;
  higher3nt       =       signal;
  higher3nt[-ID1] =       0;
  out             =       higher3nt[1:TR_lth];
}


library(R6)

l=1000
s=10

# txtplot(1:10,sapply(1:10, function(s) c(3,2,1)%>%multiply_by(s)%>%rep(l)%>%{.[666]=100;.}%>%ftestvect%>%.[1]%>%round(3)%>%divide_by(l)%>%sqrt))

# txtplot(1:10,sapply(1:10, function(s) c(3,2,1)%>%multiply_by(s)%>%rep(l)%>%{.[667]=100;.}%>%DWPT%>%sum%>%divide_by(l)))



fmcols <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)[start(grl@partitioning)]
}
fmcols_List <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)%>%split(grl@partitioning)
}
quicktest<-function(x,y){
  require(txtplot)
  complete = is.finite(x) & is.finite(y)
  message(str_interp('{sum(!complete)} (round(mean(!complete)*100,2) %) missing values'))
  txtplot(x[complete],y[complete])
  cor.test(x[complete],y[complete])
}
trimids = function(str) str_replace(str,'\\.\\d+','')
factnum = function(fact) as.numeric(as.character(fact))

str_split_fast = function(x,sep=';') x %>% {str_split(.,sep,n = str_count(.,';')%>%max%>%add(1))}
sep_element_in<-function(colonlist,ridssplit,sep=';'){
  assert_that(is.character(colonlist))
  assert_that(is.character(ridssplit))
  values<-colonlist%>%str_split_fast

  inds <- rep(seq_along(colonlist),lengths(values))

  values<-flatten_chr(values)

  data_frame(inds,values)%>%
    mutate(match = values %in% ridssplit)%>%
    group_by(inds)%>%
    summarize(match=any(match))%>%
    .$match

}

make_quantcompplot <- function(compdf, col1, col2, fname){
  require(LSD)
  base::source(here('Applications/LSD/R/LSD.heatscatter.R'))
  require(broom)
  col1<-enquo(col1)
  col2<-enquo(col2)
  corlabel = compdf%>%filter(is.finite(!!col1),is.finite(!!col2))%>%
    summarise(tidy(cor.test(!!col1, !!col2)))
  corlabel = corlabel%>%
    mutate(
      pformat=format(p.value,format='e',digits=4),
      pvalstring = ifelse(p.value > 0.001,round(p.value,4),pformat),
      labl=paste0('rho = ',round(estimate,3),'\n','pval = ',pvalstring))
  #
  nlabel=tibble(labl=paste0('N=',nrow(compdf)))
  pdf(fname)
  gplot = heatscatter(ggplot=TRUE,
      compdf[[quo_name(col1)]],compdf[[quo_name(col2)]])+
    scale_x_continuous(quo_name(col1))+
    scale_y_continuous(quo_name(col2))+
    ggtitle(basename(fname))+
    geom_text(show.legend=F,data=corlabel,
      hjust=1,vjust=1,x= Inf,y=Inf,aes(label=labl))+
    geom_text(show.legend=F,data=nlabel,
      hjust=0,vjust=1,x= -Inf,y=Inf,aes(label=labl))
  dev.off()
  pdf(fname)
  gplot
  print(gplot)
  dev.off()
  message(normalizePath(fname))
}

get_cds_bin_counts <- function(
    cov,
    trspacecds,
    STOPWINDSIZE,
    STARTWINDSIZE
  ){
    ltrspacecds = trspacecds
    longcdstrs = names(ltrspacecds)[ltrspacecds%>%width%>%`>`(MINCDSSIZE)]
    ltrspacecds = ltrspacecds[longcdstrs]

    eltrspacecds <- ltrspacecds[longcdstrs]
    seqlengths(eltrspacecds)[longcdstrs] %<>% add(fputrext+tputrext)
    eltrspacecds %<>% shift(fputrext)
    eltrspacecds%<>%resize(width(.)+FPEXT,'end')%>%resize(width(.)+TPUTREXT,'start')

    midwind = eltrspacecds%>%
      resize(width(.)-(STOPWINDSIZE),'end')%>%
      resize(width(.)-(STARTWINDSIZE),'start')
    midmat = cov[midwind]%>%
      sum%>%#compress these variable length Rles down to 1 value per gene
      matrix#make a 1 column matrix
    stwind = eltrspacecds%>%resize(STARTWINDSIZE,ignore.strand=TRUE)
    startmat = cov[stwind]%>%as.matrix
    endwind = eltrspacecds%>%resize(STOPWINDSIZE,fix='end',ignore.strand=TRUE)
    endmat = cov[endwind]%>%as.matrix
    out = cbind(startmat,midmat,endmat)
  }