

#get combination initvals with that
stanpars_to_list<-function(stanparvect){
  #get the dimensions of each of the parameters
  iparnames <- stanparvect%>%names
  parnames <- iparnames%>%str_remove('(\\[|\\.)[0-9\\,]+(\\]|\\.)$')
  parinds <- iparnames%>%str_extract('(\\[|\\.)[0-9\\,]+(\\]|\\.)$')%>%str_extract_all(regex('\\d+'))%>%map(as.numeric)
  #
  parinds <- parinds%>%split(parnames)%>%map(.%>%simplify2array%>%t)
  stopifnot(all(unique(parnames)%in%c(c("sigma2", "l_st", "lcM", "l_pihalf", "prot0", "cv", "prot",
		"mRNA", "Kd", "lKd", "lKs", "zetastar", "cM", "lp__"))))
  if('zetastar' %in% names(parinds) )parinds[['zetastar']]%<>%t
  #
  stanlist <- parinds%>%map(.%>%apply(2,max)%>%replace_na(1)%>%array(NA,dim=.))
  #
  parname='zetastar'
  for(parname in unique(parnames)){
    stanlist[[parname]][parinds[[parname]]] <- stanparvect[parnames==parname]
  }
  stanlist
}


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
  #
  n_inds <- length(colnames(parsedpars))-1
  #
    parsedpars%>%.[,-1]%>%.[,]%>%
            set_colnames(indnames)%>%
            as.data.frame%>%
          mutate(parameter=parsedpars[,1])%>%
      select(parameter,!!!indnames)%>%
      mutate_at(vars(!!!indnames),list(~as.numeric(as.character(.))))%>%
      split(.,seq_len(nrow(.)))
}


tidy_opt_prots<-function(pars,gnames){
	protdf<-pars$prot%>%set_rownames(gnames)%>%set_colnames(tps)%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%
		gather(time,estimate,-uprotein_id)%>%
		filter(!str_detect(uprotein_id,'\\.\\d$'))%>%
		safe_left_join(metainfo%>%distinct(uprotein_id,gene_name))
	protdf%<>%mutate(
		CI.L = estimate - 0.0001,
		CI.R = estimate - 0.0001,
		ntime = match(time,tps),
		assay='MS'
	)
	protdf
}

parse_sampling_pars<-function(sgenesampling,uids){
	out = sgenesampling%>%rstan::summary(.)%>%.[[1]]%>%as.data.frame(stringsAsFactors=FALSE)%>%select(estimate = `50%`, CI.L = `2.5%`,CI.R = `97.5%`)%>%
	rownames_to_column('par')%>%mutate(params=par%>%parse_stan_pars(c('gene','time')))%>%unnest%>%
	select(parameter,ntime=time,gene,estimate,CI.L,CI.R)
		singlegene = out%>%filter(parameter=='prot')%>%.$gene%>%is_in(1)%>%all
	if(singlegene){	out$gene%<>%replace_na(1)}
	out%<>%mutate(uprotein_id=uids[gene])%>%
	safe_left_join(allow_missing=TRUE,metainfo%>%filter(!is.na(uprotein_id))%>%distinct(uprotein_id,gene_name))%>%
	# filter(parameter=='prot')%>%
	mutate(assay='MS')
	out%>%as.data.frame
	out
}


################################################################################
########Functions to manipulate Rstan models
################################################################################
parse_stanlines <- function(stanmodel){
  stanlines <- stanmodel@model_code%>%str_extract_all(regex('(data|parameters|transformed parameters|model|generated quantities) ?\\{\n.*?\n\\}\n',dotall=TRUE))%>%.[[1]]
  names(stanlines) <- stanlines%>%str_match('(.*?)\ *(?=\\{)')%>%.[,2]
  stopifnot(c('data','parameters','model')%in% names(stanlines))
  stanlines%<>%map(str_split,'\n')%>%map(1)
  stanlines%<>%map_df(.id='block',function(l) l %>%str_match('\\s*(real|int|vector|matrix|array)\\s*(<.*>)?\\s*(\\[.*\\])?\\s*(\\w+)\\s*.*;.*')%>%
    set_colnames(c('line','class','limits','dimensions','varname'))%>%as_tibble%>%mutate(line=l))  
  # stanlines%>%filter(line%>%str_detect('sigma2'))%>%as.data.frame
  # assert_that(!anyDuplicated(na.omit(stanlines$varname)))
  # stanlines%>%group_by(varname)%>%filter(!is.na(varname))%>%filter(n()>1)
  datavars <- stanlines%>%filter(block=='data')%>%.$varname
  paramvars <- stanlines%>%filter(block%in%c('parameters|transformed parameters'))%>%.$varname
  assert_that(!any(datavars)%in%paramvars)
  assert_that(!any(paramvars)%in%datavars)
  #
  stanlines%<>%mutate(isdef = !is.na(varname))
  # stanlines%<>%select(-matches('(lr)limit'))%>%filter(!is.na(name))%>%.$limits%>%str_match('(?<=\\<)\\s*(lower=)?(.*),(upper=)?(.*)>')%>%.[,c(3,5)]%>%as.data.frame%>%set_colnames(c('llimit','ulimit'))
  stanlines%<>%select(-matches('(l|u)limit'))%>%
    # filter(!is.na(name))%>%
    cbind(.,str_match(.$limits,'(?<=\\<)\\s*(lower=)?([^,]*)?,?(upper=)?(.*)>')%>%.[,c(3,5)]%>%as_tibble%>%set_colnames(c('llimit','ulimit')))
  # stanlines%>%filter(!is.na(name))%>%cbind(.,str_match(.$dimensions,'(?<=\\[)\\s*(.*),(.*)]')%>%.[,c(2,3)]%>%as.data.frame%>%set_colnames(c('llimit','ulimit')))
  # stanlines$limits%>%str_subset('lower')%>%str_match()
  rowvarnames <- stanlines$varname
  isdup = rowvarnames%in%rowvarnames[duplicated(rowvarnames)]
  rowvarnames[isdup]<-NA
  assert_that(all(! rowvarnames%>%na.omit%>%duplicated))
  #
  rownames(stanlines)[!is.na(rowvarnames)] <- rowvarnames[!is.na(rowvarnames)]
  stanlines
}


pars_atlimit <- function(parlist,stanmodel){
  stanlines <- parse_stanlines(stanmodel)
  parlist <- combinitvals
  assert_that(all(names(parlist)%in%stanlines$varname))
  par = 'l_st'
  for(par in names(parlist)){
    stopifnot(par %in% stanlines$varname)
    parulim = stanlines%>%filter(varname==par)%>%.$ulimit
    stopifnot(length(parulim)==1)
    parllim = stanlines%>%filter(varname==par)%>%.$llimit
    stopifnot(length(parllim)==1)
    # if(length())
    islowlim = parlist[[par]] ==parllim 
    isuplim = parlist[[par]] ==parulim 
    parlist[[par]][] <- ifelse((islowlim%in%TRUE),-1,ifelse((isuplim%in%TRUE),1,0))
  }
  parlist
}
# pars_atlimit(combinitvals,dp_model)%>%unlist%>%table

# fixlist <- list(pi_half = -3)

# fix_parameters()


# #this works
# stan_model(model_code = stanlines$line%>%paste0(collapse='\n'))

# fix_parameters(dp_model,list(pi_half=1))

# combinitvals




extract_paramlims <- function(stanmodel,pars=TRUE){
  blocks <- extract_blocks_as_df(stanmodel,c('parameters','transformed parameters'))
  parlines <- extract_pars(parlines,blocks)
  lims <- extract_lims(parlines)
  lims[pars]
}

# dp_model%>%extract_paramlims('l_pihalf')


fix_param <- function(stanmodel,parlist){
  modeltext <- modeltext(stanmodel)
  paramblock <- extract_block(modeltext,c('parameters'))
  tparamblock <- extract_block(modeltext,c('transformed parameters'))  
  datablock <-  extract_block(modeltext,c('data'))
  extparams <- extract_pars(names(parlist),list(paramblock,tparamblock))
  extparams <- make_datalines(extparams,parlist)
  datablock <- add_decs(datablock,extparams)
  outtext <- insert_blocks(modeltext,
    blocklist=list(data=datablock,parameters=paramblock,'transformed parameters'=tparamblock)
  )
  outtext
}


  

