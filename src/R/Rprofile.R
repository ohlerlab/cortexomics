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
library(hashmap)

safe_hashmap<-setRefClass("Safe_Rcpp_Hashmap",
      contains="Rcpp_Hashmap",
      inheritPackage=TRUE
)
setMethod('[[','Safe_Rcpp_Hashmap',function (x, i, j, default,...){
  hashmapname = ''
    .local <- function (x, i, j, ..., exact = TRUE)
    {
        x$`[[`(i)
    }
    out <- .local(x, i, j, ...)
    if(missing(default)){
      if(any(is.na(out))){
        keymissingtxt = as.character(i[is.na(out)])%>%head%>%{ifelse(nchar(.)>15,paste0(substr(.,0,13),'...'),.)}%>%paste0(collapse='...')
        stop(paste0('Keys missing from safe hashmap: ',hashmapname,':',keymissingtxt))        
      }
    } else if(length(default)==1){
      out[is.na(out)] <- default
    }else{
      out[is.na(out)] <- default[is.na(out)]
    }
    return(out)
})

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
  message('memory in use ',gigsused())


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
  seqinfo(gr)<-si
  end(gr) > seqlengths(gr)[as.character(seqnames(gr))]
}
is_out_of_bounds <- function(gr,si = seqinfo(gr)){
  start(gr)<1 | is_offchr(gr,si) 
}

get_all_obsizes <- function(){.GlobalEnv%>%names%>%discard(~is.function(get(.)))%>%setNames(.,.)%>%map(~object.size(get(.)))}

allobjsizes<-get_all_obsizes()
allobjsizes%<>%enframe
allobjsizes$value%<>%unlist
allobjsizes$value%<>%divide_by(1e6)
allobjsizes%>%arrange(desc(value))
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

     if(length(x)<1024/2){padding<-1024}
     if(length(x)>=1024/2){padding<-"default"}

     resSpec1 <- spec.mtm(as.ts(x), k=n_tapers, nw=time_bw, nFFT = padding, centreWithSlepians = TRUE, Ftest = TRUE, maxAdaptiveIterations = 100,returnZeroFreq=F,plot=F,dpssIN=slepians_values)
     
     resSpec2<-dropFreqs(resSpec1,0.29,0.39)
     
     closestfreqind <- which(abs((resSpec1$freq-(1/3)))==min(abs((resSpec1$freq-(1/3)))))
     
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


library(R6)

#could include the lengths and variables it needs
#could include a method for getting the data
#and a method for shiftng a granges object
get_predictedshifts<-function(seqshiftmodel,data4readshift,probcutoff=0.5){

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


Psite_model<-R6Class("Psite_model",
  public = list(
    offsets=NULL,
    seqshiftmodel=NULL,
    compartments=NULL,
    referencefasta=NULL,
    rl_sep_seqmodel=NULL,
    initialize=function(offsets=NULL,seqshiftmodel=NULL,compartments=NULL,referencefasta=NULL){

      if(is.character(referencefasta)) referencefasta <- FaFile(referencefasta)
      #checks - compartments has nrow
      #seqshiftmodel is a ranger object,
      #compartments is a vecotr with 'nucl' in it
      if(!is.null(offsets)){
        is.data.frame(offsets)
        stopifnot(c('offset','length') %in% colnames(offsets)) 
      }
      if(!is.null(seqshiftmodel)) {
        stopifnot(!is.null(referencefasta))
        chrs <- seqinfo(referencefasta)@seqnames%>%unique
        if(is(seqshiftmodel,'ranger')){
          self$rl_sep_seqmodel <- FALSE
        }else{
          assert_that(is.list(seqshiftmodel))
          assert_that(all(map_lgl(seqshiftmodel,is,'ranger')))
          assert_that(all(names(seqshiftmodel) %in% (1:100)))
          self$rl_sep_seqmodel <- TRUE
        }

      }
      self$offsets <- offsets
      self$seqshiftmodel <- seqshiftmodel
      self$compartments <- as(compartments,'List')
      self$referencefasta <- referencefasta
  },


  get_cds_offsets = function(reads_tr){

    if('compartment'%in%colnames(self$offsets)) {
        compartmentcol <- 'compartment'
        if(!'compartment' %in% colnames(mcols(reads_tr))){      
          mcols(reads_tr)$compartment <- self$compartments[seqnames(reads_tr)]%>%unlist(use.names=F)
  
      } 
    }else{
        compartmentcol <- NULL
    }
    if(! 'length' %in% colnames(mcols(reads_tr))) {
      mcols(reads_tr)$length <- possibly(qwidth,otherwise=width(reads_tr))(reads_tr)
    }
    
    phasecol <- if('phase'%in%colnames(self$offsets)) 'phase' else NULL
    joincols <- c('length',compartmentcol,phasecol)

    assert_that(all(has_name(mcols(reads_tr),joincols)))
    

    reads_tr%>%
      mcols%>%
      as_tibble%>%
      safe_left_join(self$offsets%>%
        ungroup%>%
        # select(-phase,-score),by=c('length','compartment'),
        select(offset,one_of(joincols)),by=joincols,
        allow_missing=TRUE
      )%>%
      .$offset
    },  


  get_seq_offsets = function(seqoffreads,nbp=2){
      if(is.null(self$seqshiftmodel)) return(rep(0,length(seqoffreads)))
      
      if(! 'length' %in% colnames(mcols(seqoffreads))) {
        mcols(seqoffreads)$length <- possibly(qwidth,otherwise=width(seqoffreads))(seqoffreads)
      }
      
      data4readshift <- get_seqforrest_data(seqoffreads,self$referencefasta,nbp)
      browser()
      if(self$rl_sep_seqmodel){
        seqshift <- rep(NA,nrow(data4readshift))
        dlens <-unique(data4readshift[,'length'])
        stopifnot(all(dlens %in% names(self$seqshiftmodel)))
        for(ilen in unique(dlens)){
          imodel <- self$seqshiftmodel[[as.numeric(names(self$seqshiftmodel)) == ilen]]
          stopifnot(!is.null(imodel))
          stopifnot(!is(imodel,'ranger'))
          
          irows <- dlens==ilen

          seqshift[irows] <- self$get_predictedshifts(imodel,data4readshift[irows,],prob=0.5)
        }
        stopifnot(!any(is.na(seqshift)))

      }

      seqshift <- get_predictedshifts(self$seqshiftmodel,data4readshift,prob=0.5)
      return(seqshift)

  },

  get_offsets = function(seqoffreads,nbp=2){
      self$get_cds_offsets(seqoffreads) + self$get_seq_offsets(seqoffreads)
  },

  get_psites = function(offsetreads,filterreads=FALSE){


    offset <- self$get_offsets(offsetreads)
    mcols(offsetreads)$psite <- !is.na(offset)
    #we can also just use a pre-made shift that's in the reads' mcols
    if(is.character(offset)){
      offset = rowSums(as.matrix(mcols(offsetreads)[,offset]))
    } 
    #or a single number
    if(length(offset)==1) offset = rep(offset,length(offsetreads))
    isneg <-  as.logical(strand(offsetreads)=='-')
    offsetreads[!isneg] <- qnarrow(offsetreads[!isneg],start=offset[!isneg]+1,end = offset[!isneg]+1)
    ends <- qwidth(offsetreads[isneg])-offset[isneg]
    offsetreads[isneg] <- qnarrow(offsetreads[isneg],start=ends,end = ends  ) 

    if(filterreads){ 
      offsetreads <- offsetreads[mcols(offsetreads)$psite]
    }
    offsetreads
  })
)

# if(exists('psite_model')){
#   psite_model <- Psite_model$new(
#     offsets=psite_model$offsets,
#     seqshiftmodel=psite_model$seqshiftmodel,
#     compartments=psite_model$compartments,
#     referencefasta=psite_model$referencefasta
# )
# }


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

  # if(is(trainreads,"GAlignments")){
    # trainreads<-trainreads[njunc(trainreads)==0]
    # trainreads <- GRanges(trainreads)
  # }

  fp<-trainreads%>%resize(nbp,'start')%>%as("GRanges")%>%resize(nbp*2,'end')
  # fp_inbounds <- !is_out_of_bounds(fp,seqinfo(seq))

  tp <- trainreads%>%resize(nbp,'end')%>%as("GRanges")%>%resize(nbp*2,'start')
  # tp_inbounds <- !is_out_of_bounds(tp,seqinfo(seq))
  
  # inbounds <- fp_inbounds & tp_inbounds
  # inbounds%>%table
  
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



# conflict_prefer('has_name','assertthat')

get_genomic_psites <- function(bam,windows,mapqthresh=200,psite_model) {
  require(GenomicAlignments)
  riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=windows)
  reads <- readGAlignments(bam,param=riboparam)

  mcols(reads)$cdsshift <- get_cds_offsets(reads,psite_model$offsets,psite_model$compartments)

  reads <- reads%>%subset(width %in% psite_model$offsets$length)
  
  mcols(reads)$length <- width(reads)
  reads%<>%subset(!is.na(cdsshift))
  psites <- apply_psite_offset(reads,c('cdsshift'))%>%as("GRanges")
  mcols(psites)$length <- mcols(reads)$length   
  psites
}





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
    
    grl = resize_grl_startfix(grl,grlwidths - floor(diffs/2))
    grl = resize_grl_endfix(grl,grlwidths - ceiling(diffs/2))
    
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

get_genomic_psites <- function(bam,windows,offsets,mapqthresh=200) {
  require(GenomicAlignments)
  riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=windows)
  reads <- readGAlignments(bam,param=riboparam)
  #
  if(is.null(offsets)){
  mcols(reads)$offset <- floor(qwidth(reads)/2)
  }else{
    mcols(reads)$offset <- 
      data.frame(length=qwidth(reads),compartment='nucl')%>%
      safe_left_join(offsets,allow_missing=TRUE)%>%.$offset
  }
  #
  reads <- reads%>%subset(!is.na(mcols(reads)$offset))
  # 
  mcols(reads)$length <- width(reads)
  reads%<>%subset(!is.na(offset))
  psites <- apply_psite_offset(reads,c('offset'))%>%as("GRanges")
  mcols(psites)$length <- mcols(reads)$length   
  psites
}
#metaplots
