
read_covariates<-function(covariates_file){
  covariates <- read_tsv(covariates_file)

  assertthat::assert_that(all(c("sample_id","LABEL","CONDITION","REPL") %in% colnames(covariates)))
  assertthat::assert_that(!any(duplicated(covariates$sample_id)))
  assertthat::assert_that(!any(duplicated(covariates$LABEL)))

  covariates

}


create_metadata_file<-function(covariates,metadata_file){
  require(assertthat)
  message('creating metadata file')

  # covariates<-covariates[rev(1:nrow(covariates)),]

  covariate_names <-
    colnames(covariates)%>%
    setdiff(c("sample_id","LABEL","CONDITION","REPL"))
  fact_ordering <- rep(list(),length(covariate_names))

  covar<-covariate_names[[1]]

  covariate_names <- setNames(covariate_names,covariate_names)

  #return factor levels
  covar_levels <- lapply( covariate_names,function(covar){
    vals <- covariates[[covar]]
    #if the column is numeric order by number
    if(is.numeric(vals)) return(NULL)
    #if the column contains 'wt','WT','wild type' then make that first
    iswild <- grepl(vals,pattern='wt$|(WT$)|wild')
    if(sum(iswild==1)) return(unique(c(vals[iswild],vals[!iswild])))
    #else order as they appear in the file
    unique(vals)
  })
  covar_levels <- covar_levels[!map_lgl(covar_levels,is.null)]
  cat(yaml::as.yaml(covar_levels),file=(metadata_file))

  create_metadata_file

}

#iterate over our subsets and return lists of sampls
get_sample_vect<-function(samples_i,covariates){
  #if it's a single string interpret as a regex on the labels
  if(is.character(samples_i) & length(samples_i)==1){
    return(grep(covariates[['sample_id']],value=TRUE,pattern=samples_i))
  }else if(is.list(samples_i)){
    #else interpret as a list of columns to be filtered on and permitted values
    assert_that(is.list(samples_i))
    filtervars <- names(samples_i)
    assert_that(all(filtervars %in% colnames(covariates)))
    filters <- lapply(filtervars,function(filtervar){
      covariates[[filtervar]] %in% samples_i[[filtervar]]
    })
    filter <- Reduce('*',filters)#aggregate column filters
    return(covariates[['sample_id']][filter])
  }else {
    assert_that(all(samples_i %in% covariates[['sample_id']]))
    samples_i
  }
}


make_test_dds<-function(covariates,samples,modelstring){
  tryCatch({
  #use deseq
  require(DESeq2)
  #make a very small, fake DESeq object with our design and sample annot
  subcovars = covariates[covariates[['sample_id']] %in% samples, ]
  ex_dds <- makeExampleDESeqDataSet(m=nrow(subcovars))
  dds=DESeq2::DESeqDataSetFromMatrix(
    counts(ex_dds)[1:20,],
    colData = subcovars,
    design = as.formula(modelstring),
  )
  #fit the model
  dds <- {purrr::quietly(DESeq2::DESeq)(dds,betaPrior = FALSE )$result }
  dds},
  error = function(err){
    warning(paste0('Cannot create Test DEseq object:\n',err))
  })
}

is_subanalysis<-function(){
  'subanalyses' == (getwd()%>%dirname%>%dirname%>%dirname%>%dirname)

}


##This function traverses and checks the subfolder structure
traverse_folders <- function(subsets, subanalysisfolder, covariates, render_reports = FALSE, child=FALSE) {
  
  for(subsetname in names(subsets)){
    
    subset <- subsets[[subsetname]]
    message('Checking subset ',subsetname)
    #set up subset folder
    subsetfolder <- file.path(subanalysisfolder,subsetname) 
    dir.create(subsetfolder,showWarnings=FALSE,rec=TRUE)
    #create the sample file required
    write_tsv(data.frame(sample_id = subset$samples),file.path(subsetfolder,'samples.tsv'))
    
    for(modelname in names(subset$models)){
      
      model = subset$models[[modelname]]
      modfolder <- file.path(subanalysisfolder,subsetname,modelname)
      dir.create(modfolder,showWarnings=FALSE)
      cat(model$modelstring,file = file.path(modfolder,'modelstring.txt'))
      
      #check by creating a small DESeq object
      message('\tChecking model ',modelname,'\t',model$modelstring)
      dds <-make_test_dds(covariates,subset$samples,model$modelstring)
      resnames <- resultsNames(dds)
      
      message('\tModel ResultsNames',paste0(paste0('\n\t\t',resnames),collapse=''))
      message('\tModel Contrasts')
      
      for(contrastname in names(model$contrasts)){
        if(is.null(contrastname)) {
          message('no contrasts found for',subsetname,' - ',modelname)
          next
        }
        #use the dds object to check if the contrast specification works
        contrast = model$contrasts[[contrastname]]
        contrast_spec = contrast$contrast_spec
        
        contrastfolder <- file.path(subanalysisfolder,subsetname,modelname,contrastname)
        dir.create(contrastfolder,showWarnings=FALSE,rec=TRUE)
        #export contrast as 
        yaml::as.yaml(contrast$contrast_spec)%>%cat(file = file.path(contrastfolder,'contrastspec.txt'))
        # yaml::yaml.load_file(paste0(contrastfolder,'contrastspec.txt'))
        
        message('\t\tChecking contrast ',contrastname,':\t',contrast_spec)	
        if((length(contrast_spec)==1) &&(is.character(contrast_spec))) contrast_spec = list(contrast_spec)
        invisible(results(dds,contrast=contrast_spec))
        
      }
      file.copy(file.path(rmdfold,'model.Rmd'),modfolder,overwrite=TRUE) 
      # if the render flag is on, execute our Rmd
      if(render_reports){
        rmarkdown::render(file.path(rmdfold,'model.Rmd'),
                          knit_root_dir = modfolder,
                          output_file = 		file.path(modfolder,'model.html')
        )
      }
      if(child){
        knitr::knit_child(file.path(rmdfold,'model.Rmd'))
      }
    }
  }
  
}



