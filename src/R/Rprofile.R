
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

