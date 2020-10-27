## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = TRUE)


## ----preparation---------------------------------------------------------
message('Generating qc report in')
message(getwd())


library(assertthat,quietly=TRUE)
library(rseqdata,quietly=TRUE)
library(stringr)
library(tidyverse)

# # source(file.path(rmdfold,'build_project.R'))
# # config <- read_config(file.path(root,'config')))
# covariates_file <- file.path(root,'covariates.txt')
# design_file<- file.path(root,'design.yaml')
# covariates <- read_covariates(covariates_file)
# #read in and use the meta data to configure factor levels
# metadata <- yaml::yaml.load_file(metadata_file)




## ----overview------------------------------------------------------------

covariates %>%
  dplyr::arrange(group) %>%
  dplyr::select(sample_id,group,sample_name)%>%
  knitr::kable()




## ----qc,results="asis",eval=TRUE-----------------------------------------
if(exists('aln_stats') && (nrow(aln_stats)>0)){
    cat(knitr::knit_child(file.path(rmdfold,"qc.Rmd"), quiet = TRUE))
    cat("No align_stats found the align stats are currently only implemented for star")
}



## ----makedds for all,results="asis",eval=TRUE,message=FALSE,warning=FALSE----
#load contrasts
   dds <-
      create_dds(counts_list_prot_coding,
                 num_cores = 8,
                 design = ~1,
                 useBetaPrior = useBetaPrior)
    
 
    rld <- DESeq2::vst(dds, blind = TRUE)



## ----counts_diagnostics,results="asis",eval=TRUE-------------------------
  modeldir = root
  cat(knitr::knit_child(file.path(rmdfold,"counts_diagnostics.Rmd"), quiet = TRUE))

