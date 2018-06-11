#short script to make 3' and 5' utr files

require(quietly=TRUE,warn=FALSE,tidyverse)
require(quietly=TRUE,warn=FALSE,stringr)
require(quietly=TRUE,warn=FALSE,rtracklayer)
require(quietly=TRUE,warn=FALSE,data.table)


GR2DT = function(gr){
  require(data.table)
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

DT2GR = function(dt){
  require(GenomicRanges)
  if(is(dt,'GenomicRanges')) {
    warning('already a GRanges Object')
    return(dt)
  }


  stopifnot(c('seqnames','start')%in%colnames(dt))
  stopifnot(c('width')%in%colnames(dt)|c('end')%in%colnames(dt))
  
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
  gr=GRanges(dt[['seqnames']],IRanges(dt[['start']],dt[['end']]),strand=dt[['strand']])

  #add in the metadata if we need to
  if(length(mdatcols)){
    if(is.data.table(dt)){ mcols(gr) = dt[,mdatcols,with=FALSE]%>%as("DataFrame")
    }else{ mcols(gr) = dt[,mdatcols]%>%as("DataFrame")}
  }

    stopifnot(all(colnames(dt) %in% c(colnames(mcols(gr)),'seqnames','start','end','width' ,'strand')))

  gr
}






annot<-'static_local/gencode.vM12.annotation.gtf'%>%import
annot<-snakemake@input[[1]]%>%import

annot <- annot%>%GR2DT

# annottest = annot %>%group_by(transcript_id)%>%filter(any(type=='UTR'))%>%ungroup%>% filter(transcript_id %in% (transcript_id%>%na.omit%>%unique%>%.[1:4])) 



fp_utrs <- annot%>%
	group_by(transcript_id)%>%
	filter( 

		(strand=='-') & (start > end[type=='start_codon']) & (type=='UTR') |
		(strand=='+') & (end < start[type=='start_codon']) & (type=='UTR')

		)

fp_utrs%>%DT2GR%>%export(snakemake@output[[1]])

3p_utrs <- 
	annot%>%
	group_by(transcript_id)%>%
	filter( 

		(strand=='-') & (end <= end[type=='stop_codon']) & (type=='UTR') |
		(strand=='+') & (start >= start[type=='stop_codon']) & (type=='UTR')
		)

fp_utrs%>%DT2GR%>%export(snakemake@output[[2]])


