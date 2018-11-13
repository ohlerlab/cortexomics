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


argv <- c(
  transformedexprfile=file.path('exprdata/transformed_data.txt'),
  designmatrixfile=file.path('exprdata/designmatrix.txt'),
  foldchangesfile='exprdata/limma_fold_changes.txt'
)

argv[]<- commandArgs(trailingOnly=TRUE)

for(i in names(argv)) assign(i, argv[i])


# save.image();stop('imagesaved')

#and export
dir.create('exprdata',showWarnings = FALSE)
exprtbl <- read_tsv(transformedexprfile) 
exprtbl %<>% select(gene_name, everything())
assert_that(map_chr(exprtbl,class)[1] == 'character')
assert_that(all(map_chr(exprtbl,class)[-1] == 'numeric'))

designmatrix <- read_tsv(designmatrixfile)
designmatrix %<>% filter(dataset %>%str_detect('total|ribo|MS'))
designmatrix %<>% arrange(desc(assay),time,rep)
exprmatrix <- exprtbl  %>% { set_rownames(as.matrix(.[,-1]),.[[1]]) } 

exprmatrix %<>% .[,designmatrix$dataset]
exprmatrix%>%head
levels(designmatrix$assay) <- c('total','ribo','MS')



exprmatrix[,totalscols]%>%rowMeans%>%cut_number(10)%>%table%>%names%>%str_split(pattern='[^0-9\\-\\.]')%>%map_chr(3)%>%as.numeric
exprmatrix[,riboscols]%>%rowMeans%>%cut_number(10)%>%table%>%names%>%str_split(pattern='[^0-9\\-\\.]')%>%map_chr(3)%>%as.numeric



#in terms of variance decomp what does this amoutn to....
#seperate terms for assay difference, time, and changes between them....
design = model.matrix( ~ time + assay + time:assay , designmatrix, xlev = list(assay = c('total','ribo','MS')) )

totalscols<-colnames(exprmatrix)%>%str_subset('total')
riboscols<-colnames(exprmatrix)%>%str_subset('ribo')
mscols<-colnames(exprmatrix)%>%str_subset('MS')

exprmatrix[1,totalscols] rep(c(15,15,18,18,18)i,each=2)

limmafit = limma::lmFit(exprmatrix,design=design)
bayesfit = limma::eBayes(limmafit,trend=TRUE, robust=TRUE)


#how much of total protein variation over time is accounted for by mRNA
# mRNA var / protein var 
# = time var / time var + ribo time var + MS time var 








lcoefficients <- limmafit$coefficients%>%colnames
basetimecoefs <- lcoefficients %>% str_subset('time[^:]+$')
ribotimecoeffs <- lcoefficients %>% str_subset('time.*assayribo')
mstimecoeffs <- lcoefficients %>% str_subset('time.*assayMS')

apply%>%args

timevarbase <- apply(limmafit$coefficients[TRUE,basetimecoefs],1,var)

timevarall<-apply(
    limmafit$coefficients[TRUE,basetimecoefs]+
    limmafit$coefficients[TRUE,ribotimecoeffs]+ 
    limmafit$coefficients[TRUE,mstimecoeffs] ,1,var)

pdf('../plots/variance_explained_histogram.pdf'%>%normalizePath%T>%message)
hist(timevarbase/timevarall,n=50)
dev.off()


timevarbase <- 
  limmafit$coefficients[TRUE,basetimecoefs]

protfcdf <- (limmafit$coefficients[TRUE,basetimecoefs]+
# limmafit$coefficients[TRUE,ribotimecoeffs]+ 
limmafit$coefficients[TRUE,mstimecoeffs]) %>%
  as.data.frame%>%
  rownames_to_column('gene')%>%
  gather(coeff,value,-gene)%>%
  group_by(gene)%>%arrange(gene)%>%
  rename(prot = value)

mrnafcdf <- (limmafit$coefficients[TRUE,basetimecoefs])%>%
  as.data.frame%>%
  rownames_to_column('gene')%>%
  gather(coeff,value,-gene)%>%
  group_by(gene)%>%arrange(gene)%>%
  rename(mRNA =value)

pdf('../plots/mrna_prot_cor_dist.pdf'%>%normalizePath%T>%message)
left_join(mrnafcdf,protfcdf)%>%group_by(gene)%>%summarise(mRNA_prot_cor = cor(prot,mRNA) )%>% {hist(.$mRNA_prot_cor,50)}
dev.off()


rnameans<-exprmatrix[,totalscols]%>%rowMeans%>%enframe%>%set_colnames(c('gene','meanRNA'))
ribomeans<-exprmatrix[,riboscols]%>%rowMeans%>%enframe%>%set_colnames(c('gene','meanRNA'))

rnameans$meanRNA%>%hist(20)
ribomeans$meanRNA%>%hist(20)

pdf('../plots/mrna_prot_cor_mean_scatter.pdf'%>%normalizePath%T>%message)
left_join(mrnafcdf,protfcdf)%>%group_by(gene)%>%summarise(mRNA_prot_cor = cor(prot,mRNA))%>% left_join(rnameans)%>%{qplot(data=.,y=mRNA_prot_cor,x=meanRNA,geom='point')}
dev.off()

datafilename="http://personality-project.org/r/datasets/R.appendix1.data"
data.ex1=read.table(datafilename,header=T)   #read the data into a table

aov.ex1 = aov(Alertness~Dosage,data=data.ex1)  #do the analysis of variance
summary(aov.ex1)                                    #show the summary table
print(model.tables(aov.ex1,"means"),digits=3)
  #report the means and the number of subjects/cell
boxplot(Alertness~Dosage,data=data.ex1)
  #graphical summary appears in graphics window



genesample <- sample(rownames(exprmatrix),300)

aov <- exprmatrix %>% .[genesample,]%>% as.data.frame %>% rownames_to_column('gene')%>%    
  gather(dataset,signal,-gene)%>%
  separate(dataset,c('time','assay','rep'))%>%
  filter(assay=='MS' | assay =='ribo')%>%
  aov(data=.,signal ~ gene + gene:assay + time:gene + time:gene:assay )

summary(aov)[[1]]$`Sum Sq`%>%
  setNames(rownames(summary(aov)[[1]])%>%str_replace(' ',''))%>%
  {./sum(.)}

summary(aov)[[1]]$`Sum Sq`%>%
  setNames(rownames(summary(aov)[[1]])%>%str_replace(' ',''))%>%
  {./sum(.)}%>%.[names(.)%>%str_detect('time')]%>% {. / sum(.) }

mydataframe <- rbind


exprmatrix %>% .[genesample,riboscols]%>%rowMeans%>%cut_number(10)%>%table
exprmatrix %>% .[genesample,riboscols]%>%rowMeans%>%cut_number(10)%>%table


limmafit$coefficients[1:3,c(basetimecoefs,ribotimecoeffs)]%>%
  as.data.frame%>%
  rownames_to_column('gene')%>%
  gather(coeff,value,-gene)%>%
  group_by(gene)

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



