
library(rtracklayer)
library(stringr)
library(magrittr)
library(assertthat)
library(tidyverse)
library(here)
library(data.table)
select =dplyr::select


make_expr_dotplot<-function(expr,scatcols=scattercols){
  assert_that(all(expr %has_name% c('signal','gene_name','assay','time')))
  scale_y_log2_linlabels<-scale_y_continuous(
    name =   'T/Median TPM or LFQ',
    breaks = function(lims) {
      message(lims)
      nb = 7
      step = ceiling((lims[2] - lims[1])/7)
      lims[1] = nb*(floor(lims[1]/nb))
      lims[2] = nb*(ceiling(lims[2]/nb))
      seq(from=lims[1],to=lims[2],by=step);
    },
    labels = function(x) format((2^x),digits = 2)
  )
  
  #don't display a particular assay if it's all zero or na
  expr <- expr%>%group_by(assay)%>%filter(sum(is.na(signal)|(signal==0))<2)
  if(nrow(expr)==0) return(NULL)
  
  exprscaled <- 
    expr%>%
    mutate(signal = ifelse(signal==0,NA,signal))%>%
    group_by(assay)%>%
    mutate(signal = log2(signal/median(na.omit(signal))))%>%
    mutate(n=all(is.na(signal)),signal=ifelse(n,0,signal))%>%
    select(-n)
  
  #remove instances of crazy outliers
  # oldn=exprscaled%>%nrow
  exprscaled%<>%group_by(assay)%>%do({
    tbl<-.
    # browser()
    
    while(
      (nrow(tbl)>3) & 
      (tbl$signal%>%na.omit%>%outliers::grubbs.test(.)%>%.$p%>%`<`(0.05))
      ){
      minsig=min(tbl$signal,na.rm=TRUE)
      if ( ( median(tbl$signal,na.rm=TRUE)-minsig) < 100 ) break 
      tbl <- tbl%>%filter(signal!=minsig )
    }
    tbl
    })
  # hasoutlier = nrow(exprscaled)<oldn 
  #get colors 
  exprscaled%<>%left_join_ov(scatcols,time)
  #et limits
  ylims <- c(floor(min(exprscaled$signal)),ceiling(max(exprscaled$signal)))
  
  exprscaled%>%{ 
    ggplot(data=.,aes(y=signal,x = time,color=I(scol)))+
      geom_point(position='identity',size=I(3))+
      facet_grid( scale = 'free',.~assay )+
      expand_limits(y=ylims)+
      scale_y_log2_linlabels+
      ggplot2::labs( ylab =  'T/T0 (TPM or LFQshares)')+
      # ggtitle(subtitle = ifelse(hasoutlier,'outlier removed',''))+
      theme_bw()+
      theme(text = element_text(size = 16))
  }
  
  
}
# 
# #combine the expression data
expression_table<-bind_rows(
  totalmslfq%>%
      mutate(assay='MS')%>%ungroup%>%
      select(signal,isoform=Protein_IDs,gene_name,assay,time,replicate),
  kallistocounts%>%
    select(signal,isoform=transcript_id,gene_name,assay,time,replicate=rep)
)
totalmslfq%>%filter(gene_name=='Pa2g4')
# 
# # expression_table%>%filter(gene_name=='Pa2g4')
# #now add together isoforms
# expression_table%<>%group_by(assay,time,replicate,gene_name)%>%summarise(signal=sum(signal))
# 
# #changeing names in tables to match
# expression_table%<>%ungroup%>%mutate(time=time%>%str_replace('p',''))

#diagnostics
# totalmslfq%>%ungroup%>%filter(gene_name=='Pa2g4')%>%.$signal%>%format(scientific=FALSE)
# totalmslfq%>%ungroup%>%filter(gene_name=='Pa2g4')%>%.$signal%>%na.omit%>%
# expression_table%>%filter(gene_name=='Sox2',assay=='MS')%>%qplot(data=.,x=time,y=signal)
# expression_table%>%filter(gene_name=='Sox2',assay=='ribo')%>%qplot(data=.,x=time,y=signal)
# expression_table%>%filter(gene_name=='Sox2',assay=='total')%>%qplot(data=.,x=time,y=signal)


#get colors
scattercols<-read_csv('~/projects/cortexomics/stages_colors.tsv',col_names=F)%>%set_colnames(c('time','scol'))%>%mutate(time=time%>%str_replace('\\.',''))
#rename assays
expression_table$assay%<>%recode('total'='RNASeq','MS'='MassSpec','ribo'='RiboSeq')

# make_expr_dotplot(expression_table%>%filter(gene_name%in%'Pa2g4'))

uname = unames[12]
dir.create(showWarnings =FALSE,'~/projects/cortexomics/figures/3d_expr_dotplots')
unames <- expression_table$gene_name%>%unique%>%sample
expression_table$gene_name%>%str_subset('144_NA')

unames = unames%>%setdiff(list.files('~/projects/cortexomics/figures/3d_expr_dotplots/')%>%str_replace('.pdf',''))

mclapply(mc.cores=2,unames,function(uname){
  pdffile = str_interp('~/projects/cortexomics/figures/3d_expr_dotplots/${uname}.pdf')
  p=make_expr_dotplot(expr = expression_table%>%filter(gene_name%in%uname))
  if(!is.null(p)){
    cairo_pdf(pdffile,width = 1+(5*p$data$assay%>%n_distinct))
    print(p+ggtitle(str_interp('Expression Trajectory for ${uname}')));
    dev.off()
  }
})
uname='Flna'
unames%>%str_subset('Flna')

