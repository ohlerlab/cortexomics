tps <- c('E13','E145','E16','E175','P0')

################################################################################
########functions for working with limma fit objects
################################################################################
  
get_limmafit_stdevs <- function(limmafit,designmatrix){
  (((limmafit$stdev.unscaled)^2) %*% t(limmafit$design))%>%
    set_colnames(designmatrix$dataset)%>%
    as.data.frame%>%
    rownames_to_column('gene_name')%>%
    gather(dataset,var_signal,-gene_name)%>%
    mutate(sd_signal = sqrt(var_signal))%>%
    left_join(designmatrix)%>%
    distinct(gene_name,time,assay,.keep_all = TRUE)
}

getlimmapredictions <- function(modfit,modname,designmatrix){
  predictedvals <- get_limmafit_predvals(modfit,designmatrix)
  predictedstdevs <- get_limmafit_stdevs(modfit,designmatrix)
  predictedvals%<>%left_join(predictedstdevs)
  predictedvals%<>%select(-matches('var_'))
  predictedvals%<>%filter(rep==1)
  predictedvals%<>%select(-rep)
  colnames(predictedvals) %<>%str_replace('predicted_signal',paste0('predicted_signal_',modname))
  colnames(predictedvals) %<>%str_replace('sd_signal',paste0('sd_signal_',modname))
  predictedvals
}


get_limmafit_predvals <- function(limmafit,designmatrix){
  (limmafit$coefficients %*% t(limmafit$design))%>%
    set_colnames(designmatrix$dataset)%>%
    as.data.frame%>%
    rownames_to_column('gene_name')%>%
    gather(dataset,predicted_signal,-gene_name)%>%
    left_join(designmatrix)%>%
    distinct(gene_name,time,assay,.keep_all = TRUE)
}
################################################################################
########`Functions for working with stan objects
################################################################################


#stan function to extract certain parameters from a stanfit object.


get_prot_samples<-function(fit) fit %>%as.data.frame%>%select(matches('prot'))%>%mutate(sample=1:nrow(.))%>%gather(par,value,-sample)%>%mutate(ppars=parse_stan_pars(par))%>%unnest%>%
  filter(parameter=='prot')%>%select(time,value,sample,gene)

parse_stan_pars<-function(stanpars,indnames=c('time','gene')){
  if(any(str_detect(stanpars,'\\['))){
    stopifnot(max(str_count(stanpars,','))==1)
    parsedpars<-stanpars%>%str_match('([^\\[]+)\\[?(\\d*),?(\\d*)\\]?')%>%as.data.frame%>%
    .[,-1]
  }else{
    # browser()
    stopifnot(any(str_detect(stanpars,'\\.')))
    stopifnot(max(str_count(stanpars,'\\.'))==2)
    parsedpars<-stanpars%>%str_match('([^\\.]+)\\.?(\\d*)\\.?(\\d*)\\.?')%>%as.data.frame%>%
    .[,-1]
  }



  n_inds <- length(colnames(parsedpars))-1

    parsedpars[,-1]%>%.[,]%>%apply(1,function(x){k = keep(x,~ . !='');c(rep(NA,n_inds-length(k)),k)})%>%t%>%
            set_colnames(indnames)%>%
            as.data.frame%>%
            map_df(.,as.integer)%>%
          mutate(parameter=parsedpars[,1])%>%
      select(parameter,!!!indnames)%>%
      split(.,seq_len(nrow(.)))

}



vparse_stan_pars<-function(stanpars,indnames=c()){

  if(any(str_detect(stanpars,'\\['))){
    stopifnot(max(str_count(stanpars,','))==1)
    parsedpars<-stanpars%>%str_match('([^\\[]+)\\[?(\\d*),?(\\d*)\\]?')%>%as.data.frame(stringsAsFactors=F)%>%
    .[,-1]
  }else{
    # browser()
    stopifnot(any(str_detect(stanpars,'\\.')))
    stopifnot(max(str_count(stanpars,'\\.'))==2)
    parsedpars<-stanpars%>%str_match('([^\\.]+)\\.?(\\d*)\\.?(\\d*)\\.?')%>%as.data.frame(stringsAsFactors=F)%>%
    .[,-1]
  }

  n_inds <- length(colnames(parsedpars))-1

  parsedpars$hasgene<- parsedpars[,2]!=''
  parsedpars$hastime<- parsedpars[,3]!=''

  parsedpars$gene = ifelse(!parsedpars$hasgene,NA,
    ifelse(!parsedpars$hastime,parsedpars$V3,parsedpars$V4)
  )
  parsedpars$time = ifelse(!parsedpars$hastime,NA,parsedpars$V3)
  colnames(parsedpars)[1]<-'parameter'
  
  parsedpars%>%mutate(par=stanpars)%>%select(par,parameter,geneind=gene,timeind=time)

}

#get maximim likelihood (kinda) and mean values from our two models, with these functions
get_ml_stanfit <- function(fit){fit%>%as.data.frame%>%slice(which.max(lp__))%>%t%>%as.data.frame%>%rownames_to_column('par')%>%mutate(ppars=parse_stan_pars(par))%>%unnest%>%select(parameter,val=V1,time,gene)}
get_parsed_summary<-function(fit) fit %>%summary%>%.$summary%>%as.data.frame%>%rownames_to_column('par')%>%mutate(ppars=parse_stan_pars(par))%>%unnest

