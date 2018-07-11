#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(limma))
# suppressMessages(library(DESeq2))
message('...done')


defaultargs <- c(
  transformdexprfile=file.path('exprdata/transformed_data.txt'),
  designmatrixfile=file.path('exprdata/designmatrix.txt'),
  foldchangesfile='exprdata/limma_fold_changes.txt'
)

if(!exists('cargs')) cargs<- commandArgs(trailingOnly=TRUE)

cargs <- coalesce(
  cargs[1:length(defaultargs)]%>%setNames(names(defaultargs)),
  defaultargs
)
for(i in names(cargs)) assign(i, cargs[i])


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


design = model.matrix( ~ time + assay + time:assay , designmatrix, xlev = list(assay = c('total','ribo','MS')) )
design%>%colnames

limmafit = limma::lmFit(exprmatrix,design=design)
bayesfit = limma::eBayes(limmafit,trend=TRUE, robust=TRUE)
coefs<-bayesfit$p.value%>%colnames
allcoefs<-coefs%>%setNames(.,.)%>%
  map(~topTable(bayesfit,number=nrow(exprmatrix),coef=.,confint=0.95)%>%
  {cbind(gene=rownames(.),.)})%>%
  bind_rows(.id='coefficient')%>%
  as_data_frame

allcoefs %>%  write_tsv(paste0(foldchangesfile,'full.txt'))


coeffstoexport <- allcoefs %>%
  # filter(! coefficient=='(Intercept)')%>%
  filter(coefficient%>%str_detect('time'))
  identity

stopifnot(n_distinct(coeffstoexport$coefficient)>11)

coeffstoexport%>%
  select(gene_name=gene,logFC,coefficient)%>%
  spread(coefficient,logFC)%>%
  write_tsv(paste0(foldchangesfile)%T>%message)


# #plot the stochiometry heatmap
# catcolors = data_frame(color=c("Red","Black","Green"),pcat = c("Ribosomal","Translation Associated","Ebp1"))
# #colors for fold changes
# colors = c(seq(-15,-log2(1.25),length=100),seq(-log2(1.25),log2(1.25),length=100),seq(log2(1.25),15,length=100))
# my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
# #plot our heatmaps




# #this will need to change if/when we bring in additional model coefficients
# predictdf <- allcoefs%>%
#   # filter(gene=='Rps3')%>%
#   group_by(gene)%>%
#   mutate(intcept     = coefficient=='(Intercept)')%>%
#   mutate(coefficient = ifelse(intcept,'timeE13',coefficient))%>%
#   mutate(signal = ifelse(intcept,logFC,logFC+logFC[which(intcept)]))%>%
#   mutate(upper = ifelse(intcept,CI.R,CI.R+logFC[which(intcept)]))%>%
#   mutate(lower = ifelse(intcept,CI.L,CI.L+logFC[which(intcept)]))%>%
#   transmute(time=coefficient,gene_name_simp=gene,signal,upper,lower)%>%
#   mutate(time = as.factor(str_replace(time,'time','')))

# datadf<-
#   allexprdatadf %>% 
#   filter(fraction=='total')%>%
#   select(gene_name_simp,time,signal)%>%
#   mutate(signal=log2(signal))%>%
#   mutate(time=as.factor(time))

# samplerbps<-predictdf$gene_name_simp%>%sample(20)
# plotdatadf<-datadf%>%ungroup%>%filter(gene_name_simp%in%samplerbps)
# predictpredictdf<-predictdf%>%ungroup%>%filter(gene_name_simp%in%samplerbps)
# plotdatadf$gene_name_simp


# #creating trajectory plots from fold changes
# traj.plots<-
#   ggplot(data=plotdatadf,aes(x=as.numeric(time),y=signal)) +
#   geom_point(size=I(3)) +
#   geom_ribbon(data=predictpredictdf%>%plotdatafilt,aes(x=as.numeric(time),y=signal,ymax=upper,ymin=lower),
#     alpha=I(0.5)) +
#   scale_y_continuous(name='log2 iBAQ') +
#   ggtitle('cyto proteomic trajectory - 95% CI')+
#   facet_wrap(~gene_name_simp,ncol=4,scale='free')+
#   theme_bw()

# traj.plots%>%ggsave(file='~/projects/cortexomics/figures/all_ms_trajplots.pdf')



