library(purrr)
library(ggExtra)
library(ggpubr)
library(tidyverse)
library(Rcpp)
library(doMC)
library(rstan)
#BiocManager::install('rstan')
library(assertthat)
library(magrittr)
library(data.table)
library(stringr)
library(magrittr)
library(splines)
library(parallel)
#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(assertthat))
suppressMessages(library(limma))
message('...done')



################################################################################
########This script loads up the expression data and then runs a simple spline 
########for prediction and imputation
################################################################################
  
root <- '/fast/work/groups/ag_ohler/dharnet_m/cortexomics/'


# message('temp commented out')

transformdexprfile=file.path('exprdata/transformed_data.txt')
designmatrixfile=file.path('exprdata/designmatrix.txt')



# save.image();stop('imagesaved')

#and export
dir.create('exprdata',showWarnings = FALSE)
exprtbl <- read_tsv(transformdexprfile) 
exprtbl %<>% select(gene_name, everything())
assert_that(map_chr(exprtbl,class)[1] == 'character')
assert_that(all(map_chr(exprtbl,class)[-1] == 'numeric'))

exprmatrix <- exprtbl  %>% { set_rownames(as.matrix(.[,-1]),.[[1]]) }


designmatrix <- read_tsv(designmatrixfile)

levels(designmatrix$assay) <- c('total','ribo','MS')


########Now, let's compare models of different complexity levels
designmatrix$ribo <- designmatrix$assay %in% c('ribo','MS')
designmatrix$MS <- designmatrix$assay %in% c('MS')

design = model.matrix( ~ time + assay + time:assay , designmatrix, xlev = list(assay = c('total','ribo','MS')) )

exprdata <- exprmatrix%>%
  set_colnames(designmatrix$dataset)%>%
  as.data.frame%>%
  rownames_to_column('gene_name')%>%
  gather(dataset,signal,-gene_name)%>%
  left_join(designmatrix)

increasinggenes <- exprdata%>%
  group_by(gene_name)%>%filter(assay=='MS')%>%
  mutate(increasingMS = (mean(signal[time=='P0']) - mean(signal[time=='E13'])) > 2 )%>%
  filter(increasingMS)%>%.$gene_name

decreasinggenes <- exprdata%>%
  group_by(gene_name)%>%filter(assay=='MS')%>%
  mutate(decreasingMS = (mean(signal[time=='P0']) - mean(signal[time=='E13'])) < -2 )%>%
  filter(decreasingMS)%>%.$gene_name


#Check we aren't using hte wrongly scaled data.
stopifnot(exprdata%>%filter(gene_name=='Satb2',assay=='ribo',time=='E13')%>%.$signal%>%`<`(10))

limmafits <- list()
designmatrix$time %<>% as_factor

#fit the full model
limmafits[['full']] = limma::lmFit(exprmatrix,
                                   design=model.matrix( ~ time*(ribo+MS), designmatrix)
)

#with 4 splines
limmafits[['spline_4']] = limma::lmFit(exprmatrix,
                                        design=model.matrix( ~ ns(as.numeric(time),4)*(ribo+MS), designmatrix)
)

#with fewer splines
limmafits[['spline_3']] = limma::lmFit(exprmatrix,
                                       design=model.matrix( ~ ns(as.numeric(time),3)*(ribo+MS), designmatrix)
)

preddf <- exprdata %>%select(gene_name,dataset,signal)%>%
  left_join(getlimmapredictions(limmafits[['full']],'full',designmatrix)) %>%
  left_join(getlimmapredictions(limmafits[['spline_4']],'spline_4',designmatrix))%>%
  left_join(getlimmapredictions(limmafits[['spline_3']],'spline_3',designmatrix))


####Quick and dirty imputation with splines
exprdata%<>%left_join(preddf%>%select(gene_name,time,assay,predicted_signal_spline_3))%>%
  mutate(signal = ifelse(is.na(signal),predicted_signal_spline_3,signal))%>%
  select(-predicted_signal_spline_3)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



#sketch
#we need a 
library(zoo)
# 
#
gnamei <- 'Satb2'
#
#increasing genes
gnamesi <- c('Satb2')
#decreasing genes
'Isoc1'
# gnamesi <- c('Satb2')
gnamesi <- c('Isoc1')
#
gnamesi <- exprdata%>%.$gene_name%>%sample(7)
#
gnamesi <- 'Satb2'
g2fit <- gnamesi



#named vector containg gene legnths
genelengths <- fread(here('pipeline/feature_counts/data/E13_ribo_1/feature_counts'))%>%
  select(Geneid,Length)%>%
  left_join(fread(here('pipeline/ids.txt')),by=c(Geneid='gene_id'))%>%
  arrange(sample(1:n()))%>%
  group_by(gene_name)%>%
  slice(1)%>%
  {setNames(.$Length,.$gene_name)}


max_cores <- 32
g2fit<-exprdata$gene_name%>%unique%>%head(3)
lengthnorm=TRUE



g2fit <- c("Acadvl", "Ak1", "Asna1", "Cald1", "Cst3", "Ctnnd2", "Cul2",
"Dclk1", "Dpysl3", "Epb41l1", "Flna", "Hspa12a", "Igsf21", "Nos1",
"Orc3", "Phactr4", "Rap1b", "Rasa3", "Rps6ka5", "Satb2", "Srcin1",
"Stx1b", "Tbc1d24", "Trmt1l", "Zc2hc1a", "Zmym3")
pars=NA
g2fit<-exprdata$gene_name%>%unique
# g2fit<-'Satb2'
