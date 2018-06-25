
#source('../exploration/modeling/get_limma_fold_changes.R')
#!/usr/bin/env Rscript
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(limma))
suppressMessages(library(DESeq2))

defaultargs <- c(
	countfile='feature_counts/all_feature_counts',
	msfile=file.path('ms_tables/ms_LFQ_total_ms_tall.tsv'),
	outfolder=file.path('/fast/groups/cubi/projects/2017-10-10_cortexomics','pipeline','limma_fold_changes')
)


args <- coalesce(
	commandArgs(trailingOnly=TRUE)[1:length(defaultargs)]%>%setNames(names(defaultargs)),
	defaultargs
)
for(i in names(args)) assign(i,args[i])


# message(str_interp('filtered out ${length(nonmeasured_pids)} protein groups\
#  that were unique to the total MS data, e.g.\n ${sample(specungenes,10)}'))

countstable <- data.table::fread(countfile)

#carry out individual processing of the data sources
ids <- fread('ids.txt')%>%set_colnames(c('feature_id','gene_name'))%>%distinct
countstable <- left_join(countstable,ids)
multidgenes <- countstable$gene_name%>%table%>%keep(~ . > 1)
countstable <- filter(countstable,!gene_name %in% multidgenes) 

    
	countstable%>% { set_rownames(as.matrix(.[-1]),.[[1]]) }


countsmat <- countstable %>% select(gene_name,everything())%>%select(-feature_id)%>%{ set_rownames(as.matrix(.[-1]),.[[1]])} %>% vst
#we first need to produce some plots demonstraitng that the cdata are homoskedastic


#we first need to produce some plots demonstraitng that the cdata are homoskedastic


mstable=data.table::fread(msfile)

message('Taking genes for which every timepoint has at least some information')
#for each gene, take the protein with most signal
mstable_comp <- mstable
mstable_comp %<>% group_by(Protein_IDs)%>%filter(!any(frac_time_missing))


n_filtered_protids <- n_distinct(mstable$Protein_IDs) - n_distinct(mstable_comp$Protein_IDs)
n_filtered_gidids <- n_distinct(mstable$gene_name) - n_distinct(mstable_comp$gene_name)

message(str_interp('filtered out ${n_filtered_protids} protein groups\
for ${n_filtered_gidids} genes which leaves us with ${n_distinct(mstable_comp$gene_name)} gene ids'))


mstable_gene <- 
  mstable_comp%>%semi_join(., 
    group_by(.,gene_name,Protein_IDs)%>%
    summarize(msig=median(signal,na.rm=TRUE))%>%
    ungroup%>%
    arrange(desc(msig))%>%
    distinct(gene_name,.keep_all=TRUE)
  )


msmatrix<-mstable_gene%>%
  ungroup%>%
  select(gene_name,dataset,signal)%>%
  filter(!is.na(gene_name))%>%
  spread(dataset,signal)%>%
  { set_rownames(as.matrix(.[,-1]),.[[1]]) }%>%
  log2

countsmatrix <- countstable %>% { set_rownames(as.matrix(.[,-1]),.[[1]]) }
countmatrix<-DESeq2::vst(countsmatrix)

#print a pot showing homoskedasticity to the 
ms_meanvar_plotname <- basename(msfile)
pdf(file.path('plots','mean_variance_plots',ms_meanvar_plotname)%T>%message)
vsn::meanSdPlot(msmatrix)
dev.off()
















#join the data sources together
exprtable <- bind_rows(.id='assay',list(
	countstable,
	mstable
))

exprtable

#turn them into a matrix
exprmat = exprtable%>%
	filter(!is.na(gene_name))%>%
	group_by(gene_name,dataset)%>%
	summarize(signal=sum(na.omit(signal)))%>%
	mutate(signal = ifelse(signal==0,NA,signal))%>%#not sure about this bit - missing data is tricky
	select(gene_name,dataset,signal)%>%
	spread(dataset,signal)%>%
	{ set_rownames(as.matrix(.[-1]),.[[1]]) }

designmat <- exprtable%>%
	select(dataset,time,fraction)%>%
	filter(dataset%in%colnames(monomat))%>%
	distinct(dataset,.keep_all=TRUE)%>%
	select(-dataset)%>%
	select_if(~ n_distinct(.)>1 )

design = model.matrix( ~ time + assay + time:assay , designmat )

limmafit = limma::lmFit(expr_monomat,design=design)
bayesfit = limma::eBayes(fit,trend=TRUE, robust=TRUE)
coefs<-bayesfit$p.value%>%colnames
allcoefs<-coefs%>%setNames(.,.)%>%map(~topTable(bayesfit,number=nrow(lmonomat),coef=.,confint=0.95)%>%{cbind(gene=rownames(.),.)})%>%
  bind_rows(.id='coefficient')%>%
  as_data_frame


#this assertion requires the full model
assert_that(n_distinct(allcoefs$coefficient)==1+(ntps-1)+(nassays-1)+((ntps-1)*(nassays-1)))


#this will need to change if/when we bring in additional model coefficients
predictdf <- allcoefs%>%
  # filter(gene=='Rps3')%>%
  group_by(gene)%>%
  mutate(intcept     = coefficient=='(Intercept)')%>%
  mutate(coefficient = ifelse(intcept,'timeE13',coefficient))%>%
  mutate(signal = ifelse(intcept,logFC,logFC+logFC[which(intcept)]))%>%
  mutate(upper = ifelse(intcept,CI.R,CI.R+logFC[which(intcept)]))%>%
  mutate(lower = ifelse(intcept,CI.L,CI.L+logFC[which(intcept)]))%>%
  transmute(time=coefficient,gene_name_simp=gene,signal,upper,lower)%>%
  mutate(time = as.factor(str_replace(time,'time','')))

datadf<-
  allexprdatadf %>% 
  filter(fraction=='total')%>%
  select(gene_name_simp,time,signal)%>%
  mutate(signal=log2(signal))%>%
  mutate(time=as.factor(time))

samplerbps<-predictdf$gene_name_simp%>%sample(20)
plotdatadf<-datadf%>%ungroup%>%filter(gene_name_simp%in%samplerbps)
predictpredictdf<-predictdf%>%ungroup%>%filter(gene_name_simp%in%samplerbps)
plotdatadf$gene_name_simp


#creating trajectory plots from fold changes
traj.plots<-
  ggplot(data=plotdatadf,aes(x=as.numeric(time),y=signal)) +
  geom_point(size=I(3)) +
  geom_ribbon(data=predictpredictdf%>%plotdatafilt,aes(x=as.numeric(time),y=signal,ymax=upper,ymin=lower),
    alpha=I(0.5)) +
  scale_y_continuous(name='log2 iBAQ') +
  ggtitle('cyto proteomic trajectory - 95% CI')+
  facet_wrap(~gene_name_simp,ncol=4,scale='free')+
  theme_bw()

traj.plots%>%ggsave(file='~/projects/cortexomics/figures/all_ms_trajplots.pdf')



