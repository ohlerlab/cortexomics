
library(nlme)
library(splines)
# save.image('justincaseb4talk.Rdata')

gene2name<-file.path(root,'data/my_gencode.vM12.annotation.gtf')%>%import%>%mcols%>%as_data_frame%>%distinct(gene_id,gene_name)%>%select(feature_id=gene_id,gene_name)


totaltable<-
  Sys.glob(file.path(root,'exploration/intermediate_results/*_results.rds'))%>%
  str_subset('total')%>%
  setNames(.,str_extract(basename(.),regex('.*(?=_results.rds)')))%>%
  map(.%>%readRDS%>%select(feature_id,adj_p_value,log2fc,log2fc_se))%>%bind_rows(.id='set')%>%
  separate(set,c('assay','time'))
totaltable$log2fc[totaltable$time=='Intercept']=0
totaltable$time[totaltable$time=='Intercept']='E13'
totaltable%<>%left_join(gene2name,by='feature_id')%>%select(-feature_id)
totaltable%<>%mutate(lmin = log2fc-1.96*log2fc_se,lmax = log2fc+1.96*log2fc_se)%>%select(-log2fc_se)


riboseqfctable<-
  Sys.glob(file.path(root,'exploration/intermediate_results/*_results.rds'))%>%
  str_subset('ribo')%>%
  setNames(.,str_extract(basename(.),regex('.*(?=_results.rds)')))%>%
  map(.%>%readRDS%>%select(feature_id,adj_p_value,log2fc,log2fc_se))%>%bind_rows(.id='set')%>%
  separate(set,c('assay','time'))
riboseqfctable$log2fc[riboseqfctable$time=='Intercept']=0
riboseqfctable$time[riboseqfctable$time=='Intercept']='E13'
riboseqfctable%<>%left_join(gene2name,by='feature_id')%>%select(-feature_id)
riboseqfctable%<>%mutate(lmin = log2fc-1.96*log2fc_se,lmax = log2fc+1.96*log2fc_se)%>%select(-log2fc_se)

te_fc_table =
 Sys.glob(file.path(root,'exploration/tables/riboseqres_*txt'))%>%
   setNames(.,str_extract(basename(.),regex('.*(?=.txt)')))%>%
  map(.%>%read_tsv%>%
  .[,c(1,4,7)]%>%
  set_colnames(c('feature_id','adj_p_value','log2fc')))%>%
  bind_rows(.id='set')%>%
  separate(set,c('assay','time'))
te_fc_table$assay<-'TE'
te_fc_table%<>%left_join(gene2name,by='feature_id')%>%select(-feature_id)

mstbl<-
  read_tsv(file.path(root,'exploration/tables/ms_limma_coefs.tsv'))%>%
  # filter(!coefficient%>%str_detect('ntercept'))%>%
  select(gene_name=gene,log2fc=logFC,lmin = CI.L, lmax = CI.R,time=coefficient)%>%
  mutate(assay='MassSpec')
mstbl[str_detect(mstbl$time,'ntercept'),]%<>%mutate(lmin = lmin - log2fc,lmax=lmax-log2fc,log2fc=0,time='E13')
mstbl$time%<>%str_replace('time','')%>%str_replace('p','')

totaltable%<>%distinct(time,gene_name,assay,.keep_all=TRUE)
te_fc_table%<>%distinct(time,gene_name,assay,.keep_all=TRUE)
mstbl%<>%distinct(time,gene_name,assay,.keep_all=TRUE)

totaltable$gene_name%>%table%>%table
te_fc_table$gene_name%>%table%>%table
mstbl$gene_name%>%table%>%table


confinttbl<-  bind_rows(mstbl,riboseqfctable,totaltable)
confinttbl%<>%filter(!is.na(gene_name))
confinttbl%<>%mutate(ntime = time%>%
  match(c('E13','E145','E16','E175','P0'))
)
confinttbl$assay%<>%factor(.,levels= rev(unique(.)))

confinttbl%>%write_tsv(file.path(root,'exploration/tables/foldchangetbl.tsv'))

#########Split tables used in othe rplotting code later

#combine the fold change tables
lfc_tbl <-  
  bind_rows(mstbl,riboseqfctable,te_fc_table,totaltable)%>%
  filter(!is.na(gene_name))%>%
  mutate(time=time%>%
    factor(levels=c('E13','E145','E16','E175','P0'))
  )%>%
  {.$assay%<>%factor(.,levels= rev(unique(.)));.}

#format the levels for plots
levels(lfc_tbl$assay)%<>%recode('total'='TotalRNAseq')
levels(lfc_tbl$assay)%<>%recode('ribo'='RiboSeq')

stopifnot(lfc_tbl$assay%>%levels==c('TotalRNAseq','TE','RiboSeq','MassSpec'))

confinttbl <- lfc_tbl%>%
    filter(assay!='TE')%>%
    {.$assay%<>%factor(.,levels= rev(unique(.)));.}%>%
    mutate(ntime = as.integer(time))


l2fctable <- lfc_tbl%>%
    filter(assay!='ribo')%>%
    {.$assay%<>%factor(.,levels= rev(unique(.)));.}%>%
    group_by(gene_name)%>%filter(
      any(
        (adj_p_value[assay=='TE']<0.05)&(abs(log2fc)>log2(1.5))
      ),
    )%>%
    select(-adj_p_value)%>%
    mutate(ntime = as.integer(time))
