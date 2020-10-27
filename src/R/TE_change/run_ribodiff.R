suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,stringr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tibble))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,magrittr))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,assertthat))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,tidyverse))
suppressMessages(library(warn.conflicts = FALSE,quietly=TRUE,data.table))

args <- c(
  countfile='feature_counts/all_feature_counts',
  dispdiff=0,
  outdir= 'ribodiff'
)

args <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
for(i in names(args)) assign(i,args[i])

featurecountsagg <- data.table::fread(countfile)

timepoints = fread('sample_parameter.csv')$time%>%unique
samples = fread('sample_parameter.csv')%>%
  filter(is.na(fraction))%>%
  filter(!str_detect(sample_id,'test'))%>%
  .$sample_id
featurecountsagg = featurecountsagg%>%select(feature_id,one_of(samples))


ribocountsfiles = file.path(paste0(outdir,'/',timepoints[-1],'_ribodiff_counts.txt'))%>%setNames(timepoints[-1])
ribotblfiles = file.path(paste0(outdir,'/riboseqtbl_',timepoints[-1],'.txt'))%>%setNames(timepoints[-1])
riboseqresfiles = file.path(paste0(outdir,'/riboseqres_',timepoints[-1],'.txt'))%>%setNames(timepoints[-1])

timepoints = timepoints
my_contrasts = timepoints[-1] %>% setNames(.,.)


for(tp2 in timepoints[-1]){
  #samples to use
  ribodiffcounts =   
    featurecountsagg%>%
    dplyr::select(matches(paste0('feature_id|E13|',tp2)))
  
  ribodiffcounts = ribodiffcounts[!apply(ribodiffcounts,1,.%>%is.na%>%any),]

  #print counts
  ribodiffcounts%>%
    write_tsv(ribocountsfiles[[tp2]])
  
  #and sample info
  data_frame(Samples = ribodiffcounts%>%colnames%>%.[-1])%>%
    mutate(Data_Type=ifelse(Samples%>%str_detect('ribo'),'Ribo-Seq','RNA-Seq'))%>%
    mutate(Conditions=ifelse(Samples%>%str_detect('E13'),'E13',tp2))%>%
    write_csv(ribotblfiles[[tp2]])
}


#needs conda tho
RiboDiffscript<- file.path("../Applications/RiboDiff/scripts/TE.py")
cmds = map(timepoints[-1], function (tp2){
    paste0("source activate  'ribodiff';\ 
      python ",RiboDiffscript,
      " -d", dispdiff,
      " -e  ",ribotblfiles[[tp2]],
      " -c ",ribocountsfiles[[tp2]],
      " -o ",riboseqresfiles[[tp2]]
      )
  })
cmds[2]%>%paste0(collapse='\n')%>%cat(file='tmp.sh')
system('bash tmp.sh ')

#run the riboseq 

# mclapply(cmds,system)
# assert_that(all(file.exists(riboseqresfiles),msg = 'ribodiff files not created'))

#read in 
ribodiffcolsnms=c('feature_id','disper','p_value','adj_p_value','TE1','TE2','log2fc')
ribodiffcols=cols(
  col_character(),
  col_double(),
  col_double(),
  col_double(),
  col_double(),
  col_double(),
  col_double()
)
ribodiffcols$cols%<>%setNames(ribodiffcolsnms)

ribodiffcontrastobs <- riboseqresfiles%>%map(read_tsv,skip=1,col_names=ribodiffcolsnms,col_types=ribodiffcols)

exprtablessizefactdf<-data_frame(sample_id=sample_list,sizefact=sizeFactors(dds)[sample_list])
ribobasemeans<-ribocountsfiles%>%map(.%>%{quietly(read_tsv)(.)$result}%>%gather(sample_id,count,-feature_id)%>%
    left_join(sizefactdf,by='sample_id')%>%
    group_by(feature_id
      )%>%summarize(base_mean=sum(count/sizefact))
)
ribodiffcontrastobs<- map2(ribodiffcontrastobs,ribobasemeans,left_join)
for(i in seq_along(ribodiffcontrastobs)){ribodiffcontrastobs[[i]]$log2fc_se<-NA}
for(i in seq_along(ribodiffcontrastobs)){ribodiffcontrastobs[[i]]$stat<-NA}

names(ribodiffcontrastobs) = paste0('TE_',names(ribodiffcontrastobs))
my_contrast_objects %<>% {append(.[setdiff(names(.),names(ribodiffcontrastobs))],ribodiffcontrastobs)}




ALPHA=0.05
riboseqresfiles%>%map(.%>%read_tsv%>%filter(padj<ALPHA)%>%.$geneIDs%>%unlist%>%unique)

inputfolder<-'exprdata'
outputfolder <- 'exprdata_filt'
ribodifffolder = 'ribodiff'
exprtablefiles <- Sys.glob(paste0(inputfolder,'/*'))
exprtablefiles %<>% str_subset('tsv$|txt$')
exprtables <- map(exprtablefiles,read_tsv)

ribodifffiltgenes <- Sys.glob(file.path(ribodifffolder,'riboseqres_*'))%>%
  map(.%>%read_tsv%>%.[[1]])

for (exprtablefile in exprtablefiles){
  exprtable<-read_tsv(exprtablefile)
  if('gene' in colnames(exprtable))     gcol = 'gene'
  if('gene_name' in colnames(exprtable))     gcol = 'gene_name'
  passesfilt <- exprtable[[gcol]]%in%ribodifffiltgenes
  assert_that(!all(passesfilt==FALSE))
  message(paste0('Got ',sum(passesfilt),' genes passing filter for ',exprtablefile))
  exprtablefilt <- exprtable[passesfilt,]
  write_tsv(file.path(outputfolder,basename(exprtablefile)))
}

 exprtables%>%map( ~select(.,gene_name))

ribodifffilt

# save.image(file.path('prepimage.R'))


#let's make some toy images. Genes



# colnames(ribodiffres)%<>%str_replace_all('[\\(\\s\\)]','_')
# 

# annotation_gr <- rtracklayer::import("~/bih_cluster/projects/cubit/current/static_data/annotation/GENCODE/M12/GRCm38/gencode.vM12.annotation.gff3")
# 
# ribodiffres<-left_join(
#   ribodiffres,
#   annotation_gr%>%mcols%>%as.data.frame%>%select(gene_id,gene_name)%>%distinct(gene_id,gene_name),
#   by=c('geneIDs'='gene_id')
# )%>%  select(-disper,-pval,-X8)%>%
#   select(gene_name,everything())
# 
# ribodiffressig = ribodiffres %>% filter(padj<0.05,!is.nan(padj))
# 
# ribodiffressig%>%arrange(log2FC_TE_P0_vs_E13_)%>%distinct(gene_name,.keep_all=TRUE)
# ribodiffressig%>%arrange(-log2FC_TE_P0_vs_E13_)%>%distinct(gene_name,.keep_all=TRUE)
# 
# ribodiffres %<>% mutate(l2fc_te = `log2FC_TE(P0 vs E13)`)
# ribodiffres %<>% mutate(log10_mean_te = log10(TEE13+TEP0) / 2)
# 
# ribodiffres$aveLogCPM <-
#   edgeR::aveLogCPM(ribodiffcounts%>%select(-Entry)%>%as.matrix%>%as.integer)%>%
#   setNames(ribodiffcounts$Entry)%>%
#   .[ribodiffres$geneIDs]
#   
# ribodiffres%>%{qplot(data=.,x =aveLogCPM,y=l2fc_te,color=padj<0.05)}
