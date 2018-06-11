library(here)
library(tidyverse)
library(stringr)

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
    unique(dt%>%dplyr::select())%>%
    group_by()%>%
    slice(v)
  out=inner_join(dt,groups)
  out
}


ms_data=
  read_tsv(here('./data/MS_Ribo+Total_Proteome_raw data.txt'))%>%
  {colnames(.) = str_replace_all(colnames(.),' ','_');.}%>%
  select(gene_name=Gene_names,everything())
           

tall_ms_data <- 
  ms_data%>%
  select(gene_name,matches('iBAQ_'))%>%
  filter(!is.na(gene_name))%>%
  gather(dataset,signal,-gene_name)
  
sum_tall_ms_data <-
  tall_ms_data%>%
  group_by(dataset)%>%
  group_slice(1)%>%
  group_by(gene_name)%>%
  summarise(
    dmean = mean(signal),
    dsd = sd(signal),
    ldmean = log(dmean),
    ldsd = log(dsd),
    CV = dsd / dmean,
    nondrate = sum(signal==0),
    has_zeros = nondrate>0
)
  
qplot(data=sum_tall_ms_data,x=dmean,y=dsd)

qplot(data=sum_tall_ms_data,x=ldmean,y=CV)

qplot(data=sum_tall_ms_data,x=ldmean,y=ldsd)
qplot(data=filter(sum_tall_ms_data,has_zeros),x=ldmean,y=ldsd)


tall_ms_id <- 
  ms_data%>%
  select(id,matches('iBAQ_'))%>%
  gather(dataset,signal,-id)

sum_tall_ms_id <-
  tall_ms_id%>%
  group_by(dataset)%>%
  # group_slice(1)%>%
  group_by(id)%>%
  summarise(
    dmean = mean(signal),
    dsd = sd(signal),
    ldmean = log(dmean),
    ldsd = log(dsd),
    CV = dsd / dmean,
    nondrate = sum(signal==0),
    has_zeros = nondrate>0
  )

sig = function(x)  exp(x)/(1+exp(x))

qplot(data=filter(sum_tall_ms_id,ldmean>10),y=sig(nondrate),x=ldmean,geom=c('point','smooth'))


rnaseq_data<-
  feature_counts_data$counts%>%
    dplyr::select(id=feature_id,matches('total'))%>%
    gather(dataset,signal,-id)%>%
    group_by(id)%>%
    summarise(
      dmean = mean(signal),
      dsd = sd(signal),
      ldmean = log(dmean),
      ldsd = log(dsd),
      CV = dsd / dmean,
      nondrate = sum(signal==0),
      has_zeros = nondrate>0
    )



qplot(data=filter(rnaseq_data),y=sig(nondrate),x=ldmean,geom=c('point','smooth'))


#single prot signal is norm distributed, get mean


data%>%arrange(id)%>%.$signal
dnorm(data$signal,log = TRUE)


#wrap optim yourself
#our function would take data
#and some parameters, which are different vectors
#.. but actually you probably want optim with hardcoded 
#parameter numbers to begin with

lfunc = function(pars,data){
  
  sigpar_a = pars[1]
  sigpar_b = pars[2]
  cvpar = pars[3]
  
  means = pars[-c(1,2,3)]
  means = means[data$id]
  
  sigfit = (sigpar_a * means) + sigpar_b
  lambdadrop = 1 / ( 1 + ( sigfit  ^ (-1) ))
  lambdadrop = lambdadrop[data$id]
  pdropout = (lambdadrop) * (data$signal==0)
  
  sds = cvpar * means
  pcont = (1-lambdadrop) * dnorm(data$signal,means,sds)  
  
  L = pdroput + pcont

}

group_slice<-function(dt,v){
  require(dplyr)
  stopifnot('data.frame'%in%class(dt))
  stopifnot(v>0,(v%%1)==0)
  if(n_groups(dt)==1) warning('Only 1 group to slice')
  #get the grouping variables
  groups=
    dt%>%dplyr::select()%>%
    unique()%>%
    group_by()%>%
    .[v,]
  out=inner_join(dt,groups)
  out
}

data=tall_ms_id%>%group_by(id,dataset)%>%group_slice(1:4)

wrongmeans = data%>%summarise(m=(mean(signal*rnorm(n()))+1))%>%.$m
stopifnot(length(wrongmeans)==n_distinct(data$id))

pars = c(1,1,1,wrongmeans)

lfunc(pars,data)

optim


