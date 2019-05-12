#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(GenomicRanges))
suppressMessages(library(DESeq2))
suppressMessages(library(assertthat))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(DESeq2))

message('...done' )

filter<-dplyr::filter
select<-dplyr::select
slice<-dplyr::slice

LOWCOUNTLIM <- 10
setwd('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/pipeline/')

args <- c(
	countfile=here('pipeline','feature_counts/all_feature_counts'),
	msfile=here('pipeline',file.path('ms_tables/ms_LFQ_total_ms_tall.tsv')),
	transformdexprfile=here('pipeline',file.path('exprdata/transformed_data.txt')),
  transformd_scale_cent_exprfile=here('pipeline',file.path('exprdata/cent_scaled_exprdata.txt')),
  designmatrixfile=here('pipeline',file.path('exprdata/designmatrix.txt')),
  normcountstable=here('pipeline','exprdata/allcounts_snorm.tsv')
)

for(i in names(args)) assign(i,args[i])

if(length(base::commandArgs(trailingOnly=TRUE))) args[] <- commandArgs(trailingOnly=TRUE)[1:length(args)]%>%setNames(names(args))
args <- args[!is.na(args)]
message(capture.output(dput(args)))
for(i in names(args)) assign(i,args[i])

# message(str_interp('filtered out ${length(nonmeasured_pids)} protein groups\
#  that were unique to the total MS data, e.g.\n ${sample(specungenes,10)}'))

countstable <- data.table::fread(countfile)%>%select(-dplyr::matches('test'))

countstable%>%gather(library,count,-feature_id)%>%
  {ggplot(.,aes(x=log10(count+1),fill=str_detect(library,'ribo')))+geom_histogram(binwidth = 0.1)+facet_wrap(scale='free_y',~library,ncol=4)+
  theme_bw()
  }%T>%
  ggsave(file='../plots/readcounthistogram.pdf',h=24,w=24)


#carry out individual processing of the data sources
ids <- fread('ids.txt')%>%set_colnames(c('feature_id','gene_name'))%>%distinct
countstable <- left_join(countstable,ids)
multidgenes <- countstable$gene_name%>%table%>%keep(~ . > 1)
countstable <- dplyr::filter(countstable,!gene_name %in% multidgenes) 
countstable %<>% select(gene_name,everything())

    
#convert to matrix
countsmatrix <- countstable%>% { set_rownames(as.matrix(.[-(1:2)]),.[[1]]) }
countsmatrixall <- countsmatrix
countsmatrix %<>% .[,colnames(.)%>%str_detect('ribo|total')]
#filter out stuff with very low counts
lowmediancounts <- countsmatrix %>% apply(1,median) %>%`<`(LOWCOUNTLIM)
countsmatrix <- countsmatrix[!lowmediancounts,]

totcols<-countsmatrix%>%colnames%>%str_subset('total')
ribcols<-countsmatrix%>%colnames%>%str_subset('ribo')
# pdf(file.path('../plots','tmp.pdf')%T>%{normalizePath(.)%>%message})
# countsmatrix%>%apply(1,median)%>%add(0.5)%>%log10%>%hist(breaks=20)
# dev.off()




#and then transform the counts
countsmatrix<-cbind(
  DESeq2::vst(countsmatrix[,ribcols]),
  DESeq2::vst(countsmatrix[,totcols])
)

countsmatrix['Satb2',]

sizefactors<-DESeq2::estimateSizeFactorsForMatrix(countsmatrixall)
write_tsv(enframe(sizefactors,'sample_id','sizefactor'),here('pipeline','sizefactors.tsv'))


countsmatrix_snorm <- countsmatrix %>% {sweep(.,2,STATS = sizefactors[colnames(countsmatrix)],FUN='/')}



countsmatrix_snorm %>%{cbind(gene_name=rownames(.),as_data_frame(.))} %>% write_tsv(normcountstable)

getwd()
mstable=data.table::fread(msfile)
#some formatting differences
mstable$dataset%<>%str_replace('p5','5')
mstable$dataset%<>%str_replace('_rep','_')
mstable$dataset%<>%str_replace('^[^_]+_','')#no need to annotate what signal type it is
mstable$dataset%<>%str_replace('total','MS')#no need to annotate what signal type it is



message('Taking genes for which every timepoint has at least some information')
#for each gene, take the protein with most signal
mstable_comp <- mstable
mstable_comp %<>% group_by(Protein_IDs)%>%filter(!any(frac_time_missing))


n_filtered_protids <- n_distinct(mstable$Protein_IDs) - n_distinct(mstable_comp$Protein_IDs)
n_filtered_gidids <- n_distinct(mstable$gene_name) - n_distinct(mstable_comp$gene_name)

message(str_interp('filtered out ${n_filtered_protids} protein groups\
for ${n_filtered_gidids} genes which leaves us with ${n_distinct(mstable_comp$gene_name)} gene ids'))


mstable_comp




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

#See fi satb2 diffs are appropriate here
mstable_gene%>%filter(gene_name=='Satb2')%>%filter(!fraction_missing)%>%
  group_by(time)%>%summarise(sig=mean(signal))%>%
  {.$sig[5] / .$sig[1]}

msmatrix%>%as.data.frame%>%rownames_to_column('gene_name')%>%filter(gene_name=="Satb2")%>%
  gather%>%
  separate(key,c('time','assay','rep'))%>%
  group_by(time,assay)%>%summarise(sig=mean(as.numeric(value)))%>%
  filter(assay%in%c('ribo','MS'),time %in% c('E13','P0'))%>%arrange(assay)%>%.$sig%>%{.[2]-.[1]}

assert_that(! msmatrix%>%colnames%>%str_detect('E17p5')%>%any)
assert_that(! msmatrix%>%colnames%>%str_detect('rep\\d+')%>%any)

#we first need to produce some plots demonstraitng that the cdata are homoskedastic


#print a pot showing homoskedasticity to the 


ms_meanvar_plotname <- basename(msfile)%>%paste0(.,'.pdf')
# pdf(file.path('../plots','mean_variance_plots',ms_meanvar_plotname)%T>%message)
# vsn::meanSdPlot(msmatrix)
# dev.off()

svg(h=5,w=8,'../plots/Variance_stabilized_signal_ms.svg')
msmatrix%>%{
  standard_deviation = apply(.,1,sd,na.rm=T)
  mean_signal = apply(.,1,mean,na.rm=T)
  qplot(mean_signal,standard_deviation,geom='blank')+geom_hex()+geom_smooth()+
  theme_minimal()+
  ggtitle('Variance Stabilized Data Mass-Spec')
}
dev.off()


svg(h=5,w=8,'../plots/Variance_stabilized_signal_rnaseq.svg'%>%normalizePath%T>%message)
countsmatrix_snorm[,totcols]%>%{
  standard_deviation = apply(.,1,sd,na.rm=T)
  mean_signal = apply(.,1,mean,na.rm=T)
  qplot(mean_signal,standard_deviation,geom='blank')+geom_hex()+geom_smooth()+
  theme_minimal()+
  ggtitle('Variance Stabilized Counts Rnaseq')
}
dev.off()

svg(h=5,w=8,'../plots/Variance_stabilized_signal_ribo.svg'%>%normalizePath%T>%message)
countsmatrix_snorm[,ribcols]%>%{
  standard_deviation = apply(.,1,sd,na.rm=T)
  mean_signal = apply(.,1,mean,na.rm=T)
  qplot(mean_signal,standard_deviation,geom='blank')+geom_hex()+geom_smooth()+
  theme_minimal()+
  ggtitle('Variance Stabilized Counts Riboseq')
}
dev.off()





# count_meanvar_plotname <- basename(countfile)%>%paste0(.,'.pdf')
# pdf(file.path('../plots','mean_variance_plots',count_meanvar_plotname)%>%normalizePath%T>%message)
# vsn::meanSdPlot(countsmatrix)
# dev.off()

#join the data sources together
exprmatrix <- inner_join(
  countsmatrix %>% as_data_frame %>% cbind(gene_name = row.names(countsmatrix), .),
  msmatrix %>% as_data_frame %>% cbind(gene_name = row.names(msmatrix), .)
)%>%  { set_rownames(as.matrix(.[,-1]),.[[1]]) }

#normalize the various 'libraries' with DESeq norm factors
exprmatrix %<>% sweep(.,2,STATS = estimateSizeFactorsForMatrix(.),FUN='/')

#centered and scaled data
cent_scaled_exprmatrix<-
  exprmatrix%>% 
  sweep(.,1,STATS = rowMeans(.),FUN='-')%>%
  sweep(.,1,STATS = rowSds(.),FUN='/')


svg(h=5,w=8,'../plots/Variance_stabilized_signal_all.svg'%>%normalizePath%T>%message)
exprmatrix[,ribcols]%>%{
  standard_deviation = apply(.,1,sd,na.rm=T)
  mean_signal = apply(.,1,mean,na.rm=T)
  qplot(mean_signal,standard_deviation,geom='blank')+geom_hex()+geom_smooth()+
    theme_minimal()+
    ggtitle('Variance Stabilized Expr Sig')
}
dev.off()


#and export
dir.create('exprdata')

exprmatrix %>% as_data_frame %>% cbind(gene_name = row.names(exprmatrix))%>%write_tsv(transformdexprfile)

cent_scaled_exprmatrix %>% as_data_frame %>% cbind(gene_name = row.names(exprmatrix))%>%write_tsv(transformd_scale_cent_exprfile)

designmatrix<-exprmatrix%>%colnames%>%as_data_frame%>%set_colnames('dataset')%>%separate(dataset,c('time','assay','rep'),remove=FALSE)

designmatrix%>%write_tsv(designmatrixfile)




####################Fixing transform
#and then transform the counts
msmatrix<-mstable_gene%>%
  ungroup%>%
  select(gene_name,dataset,signal)%>%
  filter(!is.na(gene_name))%>%
  spread(dataset,signal)%>%
  { set_rownames(as.matrix(.[,-1]),.[[1]]) }

intnames <- intersect(countsmatrix%>%rownames,msmatrix%>%rownames)
message('this is quick can be fixed names intersexction')

#combine variance stabilized counts  with log mas sspec values
exprmatrix2<-cbind(
    (countsmatrix[,ribcols]),
    (countsmatrix[,totcols])
  )%>%
  .[intnames,]%>%cbind(
    log2(msmatrix[intnames,])
  )

exprmatrix2%>%as.data.frame%>%rownames_to_column('gene_name')%>%write_tsv('./exprdata/transformed_data.txt')

#exprmatrix2%<>%{sweep(.,2,STATS = meanlogdeviance(.),FUN='-')}

#exprmatrix2 <- DESeq2::vst(exprmatrix2%>%replace_na(0))

exprmatrix['Satb2',]
exprmatrix2['Satb2',]


library(rtracklayer)
allanno<-Sys.glob(here('pipeline/my_gencode*annotation*cds.gtf'))%>%import
allanno


library(biomaRt)

mousemart<- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

listFilters(mousemart)$description%>%str_subset(regex(ignore_case=T,'Uni'))
listFilters(mousemart)%>%filter(description%>%str_detect(regex(ignore_case=T,'Prot')))

'ensembl_peptide_id'
'with_uniprotswissprot'
listAttributes(mousemart)%>%filter(description%>%str_detect(regex(ignore_case=T,'uni')))
# 2       uniprot_gn  UniProtKB Gene Name ID feature_page
# 3 uniprotswissprot UniProtKB/Swiss-Prot ID feature_page
# 4  uniprotsptrembl     UniProtKB/TrEMBL ID feature_page

#pids-pid table
allprotids <- mstable%>%select(Protein_IDs)%>%distinct%>%mutate(uniprot_id=strsplit(Protein_IDs,';'))%>%unnest
#now collapse to unique 


allanno$protein_id %in% allprotidsbm$ensembl_peptide_id



# BiocManager::install('EnsDb.Mmusculus.v75')



table(allanno$transcript_id %in% names(trs_w_uniprot)) %>%{./sum(.)}


table(allprotids$uniprot_id %in% trs_w_uniprot$uniprot_id%>%str_replace('_MOUSE$',''))

attributePages(humanmart)%>%head

listAttributes(mousemart,page='feature_page')

listDatasets(mousemart,1)
# > attributes =

c("ensembl_gene_id","mmusculus_homolog_ensembl_gene","mmusculus_homolog_perc_id_r1")

# > attributes=c(attributes,"mmusculus_homolog_orthology_type",
# "mmusculus_homolog_subtype", "mmusculus_homolog_perc_id")

# >  orth.mouse = getBM( attributes,filters="with_homolog_mmus",values
# =TRUE,
# mart = human, bmHeader=FALSE)


#Construct a table of swissprotID -> transcript, with the source matched.
#Match our ids as best we can

################################################################################
########Matching protein ID groups to transcript groups
################################################################################


#we'll look at the protein ids for which there are gene names
allprotids <- mstable%>%
  # filter(is.na(gene_name))%>%
  select(Protein_IDs,gene_name)%>%
  distinct%>%
  mutate(uniprot_id=strsplit(Protein_IDs,';'))%>%
  unnest


#lod a bioconductor object
library(EnsDb.Mmusculus.v79)
edb <- EnsDb.Mmusculus.v79
trs_w_uniprot<-transcripts(edb,columns=c('uniprot_id','gene_name'))
#now get gene name, trid an uniprot id
bioc_protiddf<- trs_w_uniprot%>%mcols%>%as_data_frame%>%select(uniprot_id,transcript_id=tx_id,gene_name)
bioc_protiddf%<>%mutate(source='bioc_package')


#alsoget uniprotID-ensembl_peptide links from biomart
bm <- getBM(filters='uniprotswissprot',attributes=c('uniprotswissprot','ensembl_peptide_id'),values=allprotids$uniprot_id,mart=mousemart)
#now join our pids to the gencode piddf, and thence to gencode gids
anno_pid2tr<-mcols(allanno)%>%as.data.frame%>%filter(transcript_type=='protein_coding')%>%select(ensembl_peptide_id=protein_id,transcript_id,gene_name)
bm_protiddf <- bm%>%inner_join(anno_pid2tr)%>%select(uniprot_id=uniprotswissprot,everything())
bm_protiddf%<>%mutate(source = 'bioMart')



#load the swissprot data from this gencode version
gc_protiddf<-fread(here('annotation//gencode.vM12.metadata.SwissProt'),header=F)%>%set_colnames(c('transcript_id','uniprotkb_id','swissprot_id'))
gc_protiddf$transcript_id%<>%str_replace('\\.\\d+$','')
gc_protiddf%<>%gather(protidtype,uniprot_id,-transcript_id)%>%select(-protidtype)
#add gene name id
annotrgnamedf<-mcols(allanno)%>%as.data.frame%>%filter(transcript_type=='protein_coding')%>%select(tr_gene_name=gene_name,transcript_id)
gc_protiddf%<>%left_join(annotrgnamedf)
gc_protiddf%<>%filter(!is.na(tr_gene_name))
gc_protiddf%<>%mutate(source = 'gencode')




#Now join up all our uniID->tr links
gc_protiddf%>%head

gc_protiddf%>%filter(is.na(tr_gene_name))

allpid_tr_df<-bind_rows(
  gc_protiddf%>%select(uniprot_id,transcript_id,gene_name=tr_gene_name,source),
  bioc_protiddf%>%select(uniprot_id,transcript_id,gene_name,source)
  bm_protiddf%>%select(uniprot_id,transcript_id,gene_name,source),
)


allprotids_trs<-allprotids%>%left_join(allpid_tr_df)
allprotids_trs%<>%as_tibble

#some genes have two gene names associated
allprotids_trs%>%group_by(Protein_IDs)%>%filter(!is.na(gene_name))%>%summarise(n_distinct_gnames = n_distinct(gene_name))%>%group_by(n_distinct_gnames)%>%tally

protids_trs<- allprotids_trs%>%
  group_by(Protein_IDs)%>%
  filter(!is.na(gene_name))%>%
  filter(!((n_distinct(gene_name)>1) &(source!='gencode')))%>%
  filter(!n_distinct(gene_name)>1)%>%
  distinct(Protein_IDs,transcript_id)

protids_trs%<>%filter(transcript_id%in%allanno$transcript_id)




allcds <- allanno%>%subset(type=='CDS')%>%split(.,.$transcript_id)%>%.[protids_trs$transcript_id]

reduced_prot_cds <- allcds%>%unlist%>%split(.,protids_trs$Protein_IDs[match(names(.),protids_trs$transcript_id)])%>%reduce%>%unlist

allprotids_trs


bams<-Sys.glob(here('pipeline/star/data/*/*.bam'))%>%
  str_subset(neg=T,regex('transcript'))%>%
  str_subset(neg=T,regex('test'))

testbam<-bams%>%str_subset('E16_ribo_1')







#Get the sums for all the bams
allbamsiglist <- list()
for(bam in bams){
  message(bam)
 
 allbamsiglist %<>% append(list(bamsignals::bamCount(bam, shift=12,reduced_prot_cds , verbose=FALSE,mapqual=200)%>%split(names(reduced_prot_cds))%>%map_dbl(sum)))
}

length(allbamsiglist[[1]])==n_distinct(names(reduced_prot_cds))






library(GenomicAlignments)
summariseOverlaps(reduced_prot_cds%>%split(names(.)),bams[1])

allcountsold<-fread('feature_counts/all_feature_counts')


ms2gid <- mstable%>%select(Protein_IDs,gene_name)%>%distinct%>%tail(n=-1)%>%left_join(allanno%>%mcols%>%as.data.frame%>%distinct(gene_name,gene_id))

ms2gid%>%filter(n_distinct(gene_name)>1)
ms2gid%>%group_by(Protein_IDs)%>%filter(n_distinct(gene_id)==1)



testbamsig[mstable$Protein_IDs]

testfcountmat <- allcountsold%>%.[match(ms2gid$gene_id,allcountsold$feature_id),]

ms2gid%>%head

cor(testfcountmat[[str_subset(colnames(allcountsold),str_replace(basename(bams[1]),'.bam',''))]] ,testbamsig[ms2gid$Protein_IDs],use='complete')

compmat<-data_frame(fcount =testfcountmat[[str_subset(colnames(allcountsold),str_replace(basename(testbam),'.bam',''))]] ,bamsig=testbamsig[ms2gid$Protein_IDs])
compmat%<>%cbind(ms2gid)
compmat%>%  filter(fcount>100)%>%{cor(.$fcount,.$bamsig,use='complete')}

ms2gid%>%head%>%cliplongtext

cliplongtext<-function(df){
  df%>%mutate_if(is.character,funs(str_extract(.,'.{0,10}')))
}
ms2gid%>%cliplongtext%>%head

compmat%>%mutate(set=fcount > bamsig)%>%group_by(set)%>%tally

pdf('tmp.pdf')
qplot(data=compmat,x=log10(fcount),y=log10(bamsig))
dev.off()

compmat%>%filter(fcount>100,fcount/bamsig > 2)
compmat%>%filter(fcount>100,bamsig/fcount > 2)
compmat%>%filter(gene_name=='Orc3')



#plenty of these persist even given a single source
allprotids_trs%>%filter(source=='gencode')%>%group_by(Protein_IDs)%>%summarise(n_distinct_gnames = n_distinct(gene_name))%>%group_by(n_distinct_gnames)%>%tally
#
allprotids_trs%>%filter(source=='gencode')%>%group_by(Protein_IDs)%>%filter( n_distinct(gene_name)>1)%>%distinct(Protein_IDs,gene_name,source)%>%ungroup%>%mutate(Protein_IDs=str_extract(Protein_IDs,'.{10}'))%>%as.data.frame
#and with multiple sources....
allprotids_trs%>%group_by(Protein_IDs)%>%filter(!is.na(gene_name))%>%filter( n_distinct(gene_name)>1)%>%distinct(Protein_IDs,gene_name)%>%
allprotids_trs%>%group_by(Protein_IDs)%>%filter(!is.na(gene_name))%>%filter( n_distinct(gene_name)>1)%>%distinct(Protein_IDs,gene_name)%>%ungroup%>%mutate(Protein_IDs=str_extract(Protein_IDs,'.{0,10}'))%>%as.data.frame
#
allprotids_trs%>%group_by%>%summarise(is.na(allprotids_trs))



#
allprotids_trs



#Now having done this linking
allprotids%<>%left_join(gc_protiddf)

#now, do we ever get conflicting gene names?
#cases where we get a match via gencode pids, but it's different to the existing one
gnameconflictdf<-allprotids_trs%>%filter(source=='gencodem12_metadata')%>%filter(!is.na(tr_gene_name))%>%filter(gene_name!=tr_gene_name)%>%distinct(gene_name,tr_gene_name)%>%mutate(ms_gname_n_in_gencode=!(gene_name%in%allanno$gene_name))
gnameconflictdf%>%head
gnameconflictdf%>%.$ms_gname_n_in_gencode%>%table
gnameconflictdf%>%write_tsv(here('pipeline/gnames_conflict_pids.tsv'))

#How often do we get a gene_name where there wasn't one?
#Pretty often actually.
allprotids%>%filter(is.na(gene_name))%>%group_by(Protein_IDs)%>%summarise(isnewmatch=any(!is.na(tr_gene_name)))%>%.$isnewmatch%>%table

#How often do we get more than one gene name matching in the new gene name matches?
multmatchdf<-allprotids%>%filter(!is.na(tr_gene_name))%>%group_by(Protein_IDs)%>%summarise(hasmultmatch=n_distinct(tr_gene_name)>1)
#somewhat quite often
.$hasmultmatch%>%table

multmatchdf%>%filter(hasmultmatch)%>%write_tsv(here('pipeline/multmatch_pids.tsv'))






#So we'll just use the tr_gene_name (i.e. from gencode pid tables) where possible

#now cases where we DON"T match via gencode pids, but do have a gene name in the ms data
#We can jsut get all the transcripts for teh gene....

gnamematchonly<-allprotids%>%
  filter(is.na(transcript_id))%>%
  filter(!is.na(gene_name))%>%
  select(Protein_IDs,gene_name,uniprot_id)%>%
  inner_join(annotrgnamedf%>%select(gene_name=tr_gene_name,transcript_id))%>%
  left_join(protiddf,by='transcript_id')


gnamematchonly%>%head
table(unique(gnamematchonly$uniprot_id.x) %in% trs_w_uniprot$uniprot_id)
trs_w_uniprot

#okay so we now have multiple pid sets per gene name, sometimes.
#Theory - most genes only have a single pidset which isn't shit

# fdsafd - indorporate the bioconductor ensembl ojbect into our pid table as well.
# fdafdsa - join up the above tablles now, select best ms row for each gene

#get a df that measures this shitness
timemissingsdf<-mstable%>%group_by(Protein_IDs,time)%>%summarise(missing=all(is.na(signal)))%>%summarise(n_missing_times=sum(missing))

#so yeah it looks like most of them have a clear winner
allprotids%>%filter(!is.na(tr_gene_name))%>%distinct(Protein_IDs,tr_gene_name)%>%left_join(timemissingsdf)%>%group_by(tr_gene_name)%>%summarise(nmissingset=paste0(collapse=',',sort(n_missing_times)))%>%
  .$nmissingset%>%table

#but why are some always missing???
timemissingsdf%>%filter(n_missing_times==5)%>%left_join(mstable)%>%as.data.frame%>%head(40)

#Probably those guys are a) super low scores edited in as nas and b) things only present in other fractions
highsignal_pidsets<-timemissingsdf%>%filter(n_missing_times!=5)


#So the strategy is to use our matches where we can, and if not, just use the protein coding transcripts
#based on the gene name

#is there a 1:1 match between gene_names and majority protein IDs?

allprotids%>%group_by(Protein_IDs,gene_name)
allprotids%>%distinct(Protein_IDs,gene_name)%>%nrow
allprotids%>%distinct(Protein_IDs)%>%nrow
allprotids%>%distinct(gene_name)%>%nrow



allprotids%>%filter(is.na(transcript_id))%>%.$gene_name%>%n_distinct



allprotidsmatch<-allprotids%>%mutate(matches = uniprot_id %in% c(protiddf$uniprotkb_id,protiddf$swissprot_id))%>%group_by(Protein_IDs)%>%summarise(matches=any(matches))
allprotidsmatch$matches%>%table %>%{.[2]/sum(.)}#so 85% of genes with 
#also add gene name data
allprotidsmatch%<>%left_join(mstable%>%select(gene_name,Protein_IDs)%>%distinct)

#second table - see if the newer encode version is any better
protiddf2<-fread(here('annotation//gencode.vM21.metadata.SwissProt'),header=F)%>%set_colnames(c('transcript_id','uniprotkb_id','swissprot_id'))
allprotidsmatch2<-allprotids%>%mutate(matches = uniprot_id %in% c(protiddf$uniprotkb_id,protiddf$swissprot_id,protiddf2$uniprotkb_id,protiddf2$swissprot_id))%>%group_by(Protein_IDs)%>%summarise(matches=any(matches))
allprotidsmatch2$matches%>%table %>%{.[2]/sum(.)}#very few extras if using newest gencode


nomatchpids<-allprotidsmatch%>%filter(!matches)%>%.$Protein_IDs

mstable%>%filter(Protein_IDs %in% nomatchpids)%>%head

trids<-allanno%>%subset(gene_name%in%'Ktn1')%>%.$transcript_id%>%unique

protiddf%>%filter(str_replace(transcript_id,'\\.\\d+','') %in% trids)

nomatchgenenames <- allprotidsmatch%>%filter(!matches)%>%.$gene_name

nomatchtrids<-allanno%>%subset(gene_name%in%nomatchgenenames)%>%.$transcript_id%>%unique

#The counts I'm getting with bamsignals is quite different form those I get with feature_counts....

#So let's 



#HOw can the ids have a match in teh gene names and NOT in the swissprot metadata???
#Let's look at an example, a gene which has  a gee