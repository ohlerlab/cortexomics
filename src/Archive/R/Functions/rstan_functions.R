

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
  stanlines <- stanmodel@model_code%>%str_extract_all(regex('(data|parameters|transformed parameters|model|generated quantities) ?\\{\n.*?\n\\}(\n|$)',dotall=TRUE))%>%.[[1]]
  names(stanlines) <- stanlines%>%str_match('(.*?)\ *(?=\\{)')%>%.[,2]
  stopifnot(c('data','parameters','model')%in% names(stanlines))
  stanlines%<>%map(str_split,'\n')%>%map(1)
  # stanlines%<>%map_df(.id='block',function(l) l %>%str_match('\\s*(real|int|vector|matrix|array)\\s*(<.*>)?\\s*(\\w+)\\s*(\\[.*\\])?\\s*.*;.*')%>%
    # set_colnames(c('line','class','limits','dimensions','varname'))%>%as_tibble%>%mutate(line=l))
  # stanlines%<>%map_df(.id='block',function(l) l %>%str_match('\\s*(real|int|vector|matrix|array)\\s*(<.*>)?\\s*(\\w+)\\s*(\\[.*\\])?\\s*.*;.*')%>%
  #initally we just parse out the class, limits, and the varname + dimensions
  stanlines%<>%map_df(.id='block',function(l) l %>%str_match('^\\s*(real|int|vector|matrix|array)\\s*(<.*>)?\\s*([^\\s].*);.*')%>%
                        set_colnames(c('line','class','limits','dimensionsvarname'))%>%as_tibble%>%mutate(line=l))
  #the dimesnions can be before or after the varname
  stanlines$dimensions <- stanlines$dimensionsvarname%>%str_extract('\\[.*\\]')
  stanlines$varname <- stanlines$dimensionsvarname%>%str_replace('\\s*\\[.*\\]\\s*','')
  stanlines%>%filter(str_detect(line,'zeta'))%>%head(10)
  # stanlines%>%filter(line%>%str_detect('sigma2'))%>%as.data.frame
  # assert_that(!anyDuplicated(na.omit(stanlines$varname)))
  # stanlines%>%group_by(varname)%>%filter(!is.na(varname))%>%filter(n()>1)
  datavars <- stanlines%>%filter(block=='data')%>%.$varname%>%.[!is.na(.)]
  paramvars <- stanlines%>%filter(block%in%c('parameters','transformed parameters'))%>%.$varname%>%.[!is.na(.)]
  assert_that(!any(datavars%in%paramvars))
  assert_that(!any(paramvars%in%datavars))
  #
  stanlines%<>%mutate(isdef = !is.na(varname))
  # stanlines%<>%select(-matches('(lr)limit'))%>%filter(!is.na(name))%>%.$limits%>%str_match('(?<=\\<)\\s*(lower=)?(.*),(upper=)?(.*)>')%>%.[,c(3,5)]%>%as.data.frame%>%set_colnames(c('llimit','ulimit'))
  stanlines%<>%select(-matches('(l|u)limit'))%>%
    # filter(!is.na(name))%>%
    cbind(.,str_match(.$limits,'(?<=\\<)\\s*(lower=)?([^,]*)?,?(upper=)?(.*)>')%>%.[,c(3,5)]%>%as_tibble%>%set_colnames(c('llimit','ulimit')))
  # stanlines%>%filter(!is.na(name))%>%cbind(.,str_match(.$dimensions,'(?<=\\[)\\s*(.*),(.*)]')%>%.[,c(2,3)]%>%as.data.frame%>%set_colnames(c('llimit','ulimit')))
  # stanlines$limits%>%str_subset('lower')%>%str_match()
  stanlines$ulimit[stanlines$ulimit==""]<-NA
  stanlines$llimit[stanlines$llimit==""]<-NA
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

# stan_model(model_code='
# data{
# vector X = {1,2,3};
# }
# parameters{
# real Y
# }
# ')

extract_paramlims <- function(stanmodel,pars=TRUE){
  blocks <- extract_blocks_as_df(stanmodel,c('parameters','transformed parameters'))
  parlines <- extract_pars(parlines,blocks)
  lims <- extract_lims(parlines)
  lims[pars]
}

# dp_model%>%extract_paramlims('l_pihalf')
#
#stanmodel=dp_model

fix_param <- function(stanmodel,vars2fix){
  #parse our lines as a df
  stanlines <- parse_stanlines(stanmodel)
  stanlines%>%filter(!is.na(dimensions))
  #remember names
  linenames <- rownames(stanlines)
  #remember what order the blocks should be in
  blockorder <- stanlines$block%>%factor(unique(.))%>%levels
  stanlines$linenum = 1:nrow(stanlines)
  #indicate the first line of each block for sorting later
  stanlines%<>%group_by(block)%>%mutate(isdec = linenum==min(linenum) ) 
  #now reassing the line names
  stanlines <- stanlines%>%as.data.frame%>%set_rownames(linenames)
    #reassign key variables for a data line 
  for(varname in vars2fix){
    message(str_interp('fixing ${varname}'))
    stopifnot(sum(rownames(stanlines) %in% varname  )==1)
    stanlines[varname,]$block = 'data'
    stanlines[varname,]$ulimit = NA
    stanlines[varname,]$llimit = NA
    stanlines[varname,]$linenum = NA
    # if(stanlines[varname,]$class%>%is_in('vector')){
      # stanlines[varname,]$class%<>% str_replace('vector','real')#vectors switched to 'reals' for data
      # stanlines[varname,]$dimensionsvarname%<>%str_replace('(\\[.*?\\])(.*)','\\2\\1')#put dim after varname for reals
    # }
  }
  #now for each block, find the declarations at the start,
  #and show the new vars in there
  iblock='data'
  for(iblock in unique(stanlines$block)){
    defblockend = max(stanlines$linenum[stanlines$block==iblock & (stanlines$isdef) & (!stanlines$varname%in%vars2fix)],na.rm=T)
    stanlines$linenum[stanlines$block==iblock & is.na(stanlines$linenum)] <- defblockend + 0.5
  }
  #deal withpresence or abscence of limits
  stanlines%<>%mutate(limits=case_when(
      (!is.na(llimit))&(!is.na(ulimit)) ~   paste0('<','lower=',llimit,',','upper=',ulimit,'>'),
      (!is.na(ulimit)) ~   paste0('<','upper=',ulimit,'>'),
      (!is.na(llimit)) ~   paste0('<','lower=',llimit,'>'),
    TRUE ~ ''
  ))
  
  stanlines%<>%mutate_at(vars(class,limits,dimensionsvarname),list(~ replace(.,is.na(.),'')))
  stanlines%<>%mutate(line = ifelse(isdef, paste0(class,' ',limits,' ',dimensionsvarname,';'),line))
      #inspect 
  # stanlines%>%filter(isdef)%>%select(block,line,limits,ulimit,llimit,linenum)%>%mutate(line = str_replace(line,'//.*',''))
  #Now we want to make all of our fixed parameters come after the declaration of the codeblocks
  modeltext = stanlines %>%
      arrange(match(block,blockorder),linenum)%>%
      .$line%>%  
      paste0(collapse='\n')
  modeltext
}

get_stanpars <- .%>%parse_stanlines%>%filter(block=='parameters')%>%.$varname%>%.[!is.na(.)]
get_standata <- .%>%parse_stanlines%>%filter(block=='data')%>%.$varname%>%.[!is.na(.)]
get_stantrandata <- .%>%parse_stanlines%>%filter(block=='transformed data')%>%.$varname%>%.[!is.na(.)]

# parse_stanlines(dp_model)%>%filter(block=='parameters')%>%.$varname%>%unique

# fixparammodel <- stan_model(model_code = fix_param(dp_model,c('l_pihalf','lcM','l_st','prot0')))

# optimize(iter=1,)
#now print the whole thing stuck together 
# stanlines$line%>%paste0(collapse='\n')
