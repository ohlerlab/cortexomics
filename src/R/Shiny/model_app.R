# install.packages(c('magrittr','stringr','ggpubr','data.table','assertthat','tidyverse','dplyr','here','conflicted','ggplot2'))
#install.packages('tidyverse')

suppressMessages(library(magrittr))
suppressMessages(library(shiny))
suppressMessages(library(stringr))
suppressMessages(library(ggpubr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(here))
library(ggplot2)
library(rsconnect)
conflict_prefer("filter",'dplyr')
conflict_prefer("last",'dplyr')
# conflict_prefer("setdiff",'dplyr')
conflict_prefer("setdiff", "BiocGenerics")
if(file.exists('trajplotobjects.Rdata')){
  datafile <- 'trajplotobjects.Rdata'
}else{
  # datafile <- here('data/trajplotobjects.Rdata')
    datafile <- '/fast/groups/ag_ohler/work/dharnet_m/cortexomics/data/trajplotobjects.Rdata'

}
if(!exists('exprdf')){
   load(file=datafile)
}


#' Need to think of a way to show if proDA is working well.
#' and ADDITIONALLY, if our dropout inclusive model works better...
#' an example of that would be great.
#' Also reports that work would be great....

#list of limits for particular genes
ulimits = list('Bcl11b'=c(-2,4),'Flna'=c(-4,1),'Nes'=c(-4,1),'Satb2'=c(-1,5))


#we'll use the voom for all the counts, so as to 
exprdatacount<-mscountvoom$E%>%
  as.data.frame%>%
  rownames_to_column('uprotein_id')%>%
  filter(uprotein_id%in%best_uprotein_ids)%>%
  gather(dataset,signal,-uprotein_id)%>%
  separate(dataset,c('time','assay','rep'))%>%
  filter(assay!='MS')
pid2upid <- best_uprotein_ids%>%setNames(.,str_extract(.,'[^_]+'))
if(!'uprotein_id'%in%colnames(exprdatacount)) exprdatacount$uprotein_id <- pid2upid[exprdatacount$protein_id]
#Now we need to add in the mass spec data, with confidence intervals
ms_conf_df <- 
    left_join(
      ((postmeanmat[best_uprotein_ids,,drop=F])-(1.96*postprecsdmat[best_uprotein_ids,,drop=F]))%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,CI.L,-uprotein_id)%>%separate(dataset,c('time','assay')),
      ((postmeanmat[best_uprotein_ids,,drop=F])+(1.96*postprecsdmat[best_uprotein_ids,,drop=F]))%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,CI.R,-uprotein_id)%>%separate(dataset,c('time','assay'))
    )%>%
    left_join(
      ((postmeanmat[best_uprotein_ids,,drop=F]))%>%as.data.frame%>%rownames_to_column('uprotein_id')%>%gather(dataset,signal,-uprotein_id)%>%separate(dataset,c('time','assay')),
      )
ms_conf_df$estimate<-NA
#now stick these together
exprdata<-exprdatacount%>%bind_rows(ms_conf_df)
# exprdatacount$protein_id%>%n_distinct
# ms_conf_df$uprotein_id%>%n_distinct
#this manually derives the TE which is... probably bad. But actually needed as data
exprdata%<>%arrange(uprotein_id)%>%bind_rows(exprdata%>%filter(assay=='ribo')%>%mutate(assay='TE',signal = signal - (exprdata%>%filter(assay=='total')%>%.$signal)))
ntimes <- seq_along(tps)%>%setNames(tps)
#also add gene name
exprdata%<>%mutate(ntime=ntimes[time])%>%safe_left_join(allow_missing=T,metainfo%>%distinct(uprotein_id,gene_name))
#add in the linear predictions from limma
models2plot = prediction_df%>%mutate(ntime=ntimes[time])%>%rename('estimate':=logFC)%>%safe_left_join(allow_missing=T,metainfo%>%distinct(uprotein_id,gene_name))%>%split(.,.$assay)%>%.[c('total','ribo','MS')]
#now add stan models

models2plot <- c(
  models2plot,
  stanmodels
)
stanmodels[[1]]%>%head
names(models2plot)
models2plot[[names(stanmodels)[[1]]]]%>%head
models2plot%>%map(head)
#
geneids<-'Flna'
# exprdata%>%filter(gene_name=='Flna')
# models2plot[['MS']]%>%filter(gene_name=='Flna')
# txtplot(
#   models2plot[['MS']]%>%select(time,assay,uprotein_id,estimate)%>%.[match(best_uprotein_ids,.$uprotein_id),]%>%.$estimate,
#   ms_conf_df%>%select(time,assay,uprotein_id,signal)%>%.[match(best_uprotein_ids,.$uprotein_id),]%>%.$signal
# )
# correctionmodels<-c('MS',names(stanmodels))
correctionmodels<-c('MS')
for(corectionmodel in correctionmodels){
  coruids <- models2plot[[corectionmodel]]$uprotein_id%>%unique
  modalests <- models2plot[[corectionmodel]]%>%select(ntime,assay,uprotein_id,estimate)%>%.[match(coruids,.$uprotein_id),]%>%.$estimate
  datasigs  <- (ms_conf_df%>%select(time,assay,uprotein_id,signal)%>%.[match(coruids,.$uprotein_id),]%>%.$signal)
  datarange <- range(modalests-datasigs)
  stopifnot(diff(datarange)/datarange[1] < 0.0001)
  models2plot[[corectionmodel]]$estimate <- models2plot[[corectionmodel]]$estimate - median(datarange)
  models2plot[[corectionmodel]]$CI.L <- models2plot[[corectionmodel]]$CI.L - median(datarange)
  models2plot[[corectionmodel]]$CI.R <- models2plot[[corectionmodel]]$CI.R - median(datarange)  
}

#
stanmodels[[1]]%>%filter(gene_name=='Flna')
#
prettytpnames<-c('E12.5','E14','E15.5','E17','P0')
assaytitles <- c('RNA-seq','Ribo-seq','Mass Spec')%>%setNames(c(''))
#function to make a particular plot
#Should optinally do linerange, and ribbons
#These can then be tiled.
#our exprmodel needs to be scaled to the same level
expr_model_plot <- function(exprpoint=NULL,exprmodel=NULL,modelname,assaynames=prettyassaynames,tpnames=prettytpnames){
  #checking inputs
  stopifnot(c('signal','gene_name','assay','ntime')%in%colnames(exprpoint))
  stopifnot(c('estimate','CI.R','CI.L','gene_name','assay','ntime')%in%colnames(exprmodel))
  #plot our data
  ggplot(data = exprpoint,aes(x=ntime,y=signal))+
    geom_point()+
    geom_linerange(data=exprpoint,aes(ntime,y=signal,ymin=CI.L,ymax=CI.R))+
    geom_ribbon(data=exprmodel,aes(x=ntime,y=estimate,ymin=CI.L,ymax=CI.R),fill='darkgreen',alpha=I(0.5))+
    geom_line(data=exprmodel,aes(x=ntime,y=estimate,ymin=CI.L,ymax=CI.R),linetype=2,fill='darkgreen',alpha=I(0.5))+
    scale_x_continuous(name='Stage',labels=tpnames)+
    theme_bw()+
    scale_y_continuous(name=str_interp('Fold Change ( ${unique(exprpoint$gene_name)} )'))+
    ggtitle(modelname)
}
id_search <- function(gnames2plot,metainfo=metainfo,idvect=idvect){
  # gnames2plot='Flna'
  qs(gnames2plot,'S>=1[2,]')
  #
  #convert from upprotein ids etc if needed
  is_id <- gnames2plot%in%names(idvect)
  gnames2plot[is_id] <- idvect[gnames2plot[is_id]] 
  #we want to conver them if they are actually uids
  plotable_genes<-metainfo$gene_name
  nonfoundgenes <- gnames2plot%>%setdiff(plotable_genes)
  if(length(nonfoundgenes)==0){
    attr(gnames2plot,'notfound')<-FALSE
    return(gnames2plot)
  }else{
    nonfoundgenes <- paste(sep=',',nonfoundgenes)
    nonfoundgenecol <- paste0(nonfoundgenes,collapse=',')
    othergenes <- sample(unique(plotable_genes))
    alternativelist <- othergenes[head(order((adist(toupper(nonfoundgenes[1]),toupper(othergenes)))),n=5)]
    alternativelist <- paste0(paste0(alternativelist,' ?\n'),collapse='')
    output <- str_interp('Genes ${nonfoundgenecol} Not Found. By ${nonfoundgenes[1]} Did you mean...\n${alternativelist}\n (Case sensitive matching)')
    attr(output,'notfound')<-TRUE
    return(output)
  }
}
#iterate over models2plot for a single gene.
gene_model_plots <- function(exprdata,models2plot,myylim){
  modelassays <- models2plot%>%map_chr(.%>%.$assay%>%unique)
  #
  stopifnot(n_distinct(exprdata$gene_name)==1L)
  #
  #iterate over models
  plotlist<-lapply(seq_along(models2plot),function(i){
      mname=names(models2plot)[i]
      message(mname)
      expr_model_plot(
        exprpoint=exprdata%>%filter(assay==modelassays[i]),
        exprmodel=models2plot[[i]],
        modelname = names(models2plot)[i]
      )
  })
  plotlist[-1]%<>%lapply(function(plot) plot + scale_y_continuous(name=''))
  plotlist%<>%lapply(function(plot) plot + coord_cartesian(ylim=myylim))
  plotlist
}
geneset_model_plots <- function(geneset,exprdata,models2plot,idvect,myylims=NULL){
  gene_ids <- id_search(geneset,metainfo,idvect)
  stopifnot(  gene_ids %in% metainfo$gene_name)
  stopifnot(all(  gene_ids%in%exprdata$gene_name))
  for(i in seq_along(models2plot))stopifnot(all(  gene_ids%in%models2plot[[i]]$gene_name))
  #deal with missing names
  if(attr(gene_ids,'notfound')){
    print(qplot(1,1,label=notfoundtext,geom='text',size=I(10))+theme_bw())
  }
  #otherwise plot
  geneset_plotlist <- lapply(gene_ids,function(gid){
    gexprdata<-exprdata%>%filter(gene_name==gid)
    gmodels<-models2plot%>%map(.%>%filter(gene_name==gid))
    myylim = myylims[gid] 
    #
    gene_model_plots(gexprdata,gmodels,myylim) 
  })
  geneset_plotlist[-1]%<>%lapply(function(plotlist){
    plotlist%<>%lapply(function(plot) plot + ggtitle(''))
  })
  #
  assert_that(length(geneset_plotlist) < 9,msg='Too many genes!')
  trajectoryplot<-ggarrange(plotlist=geneset_plotlist%>%flatten,ncol=length(models2plot),nrow=length(geneset_plotlist))
  trajectoryplot<-annotate_figure(trajectoryplot)
  trajectoryplot
}
#vector that'll convert ids to uprotein_ids
convcol='gene_name'
idvect<-lapply(c('gene_name','gene_id','transcript_id','protein_id'),function(convcol){
  message(convcol)
  distinct(metainfo,gene_name,!!sym(convcol),isbest,hascount)%>%
    arrange(-isbest,-hascount)%>%
    {setNames(.[[1]],.[[2]])}
})%>%purrr::reduce(.f=c)%>%c(.,setNames(metainfo$gene_name,metainfo$gene_name))
#now test the func
plotfile<-'tmp.pdf'
pdf(plotfile,w=3*length(models2plot),h=2*length(models2plot)*length(geneids))
geneset_model_plots(
  # stanmodels[[1]]$gene_name%>%sample(1),
  c('Satb2','Flna',lowlikgenes),
  exprdata,
  models2plot,
  idvect
)
dev.off()
message(normalizePath(plotfile))



