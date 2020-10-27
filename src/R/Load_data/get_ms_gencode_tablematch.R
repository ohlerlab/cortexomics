################################################################################
########Matching protein ID groups to transcript groups
################################################################################


#we'll look at the protein ids for which there are gene names
allprotids <- allms%>%
  # filter(is.na(gene_name))%>%
  select(ms_id,gene_name)%>%
  distinct%>%
  mutate(uniprot_id=strsplit(ms_id,';'))%>%
  unnest
allprotids$uniprot_id%<>%str_replace('\\-\\d+$','')

#lod a bioconductor object
library(EnsDb.Mmusculus.v79)
edb <- EnsDb.Mmusculus.v79
trs_w_uniprot<-mymemoise(transcripts)(edb,columns=c('uniprot_id','gene_name'))
#now get gene name, trid an uniprot id
bioc_protiddf<- trs_w_uniprot%>%mcols%>%as_data_frame%>%select(uniprot_id,transcript_id=tx_id,gene_name)
bioc_protiddf%<>%mutate(source='bioc_package')


#alsoget uniprotID-ensembl_peptide links from biomart
mousemart<- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
getBM<-mymemoise(getBM)
bm <- getBM(filters='uniprotswissprot',attributes=c('uniprotswissprot','ensembl_peptide_id'),values=allprotids$uniprot_id,mart=mousemart)

#now join our pids to the gencode piddf, and thence to gencode gids
anno_pid2tr<-mcols(gtf_gr)%>%as.data.frame%>%filter(transcript_type=='protein_coding')%>%select(ensembl_peptide_id=protein_id,transcript_id,gene_name)
bm_protiddf <- bm%>%inner_join(anno_pid2tr)%>%select(uniprot_id=uniprotswissprot,everything())
bm_protiddf%<>%mutate(source = 'bioMart')



#load the swissprot data from this gencode version
gc_protiddf<-fread(here('annotation//gencode.vM12.metadata.SwissProt'),header=F)%>%set_colnames(c('transcript_id','uniprotkb_id','swissprot_id'))
gc_protiddf$transcript_id%<>%str_replace('\\.\\d+$','')
gc_protiddf%<>%gather(protidtype,uniprot_id,-transcript_id)%>%select(-protidtype)
#add gene name id
annotrgnamedf<-mcols(gtf_gr)%>%as.data.frame%>%filter(transcript_type=='protein_coding')%>%select(tr_gene_name=gene_name,transcript_id)
gc_protiddf%<>%left_join(annotrgnamedf)
gc_protiddf%<>%filter(!is.na(tr_gene_name))
gc_protiddf%<>%mutate(source = 'gencode')




#Now join up all our uniID->tr links


gc_protiddf%>%filter(is.na(tr_gene_name))

allpid_tr_df<-bind_rows(
  gc_protiddf%>%select(uniprot_id,transcript_id,gene_name=tr_gene_name,source),
  bioc_protiddf%>%select(uniprot_id,transcript_id,gene_name,source),
  bm_protiddf%>%select(uniprot_id,transcript_id,gene_name,source)
)


allprotids_trs<-allprotids%>%left_join(allpid_tr_df)
allprotids_trs%<>%as_tibble
allprotids_trs%>%head

#some genes have two gene names associated
allprotids_trs%>%group_by(ms_id)%>%filter(!is.na(gene_name))%>%summarise(n_distinct_gnames = n_distinct(gene_name))%>%group_by(n_distinct_gnames)%>%tally

protids_trs<- allprotids_trs%>%
  group_by(ms_id)%>%
  filter(!is.na(gene_name))%>%
  filter(!((n_distinct(gene_name)>1) &(source!='gencode')))%>%
  filter(!n_distinct(gene_name)>1)%>%
  distinct(ms_id,transcript_id)

protids_trs%<>%filter(transcript_id%in%cds$transcript_id)

allprotids_trs%>%ungroup%>%sample_n(10)
gc_protiddf%>%ungroup%>%sample_n(10)

ms_id2protein_id <- protids_trs%>%safe_left_join(as.data.frame(gtf_gr)%>%distinct(transcript_id,protein_id))

# allms%>%filter(gene_name=='Satb2')
satb2ids <- cds%>%mcols%>%data.frame%>%filter(gene_name=='Satb2')%>%pluck('protein_id')%>%unique
# gc_protiddf%>%head
# gc_protiddf%>%filter(tr_gene_name=='Satb2')%>%distinct


# allcds <- gtf_gr%>%subset(type=='CDS')%>%split(.,.$transcript_id)%>%.[protids_trs$transcript_id]

# reduced_prot_cds <- allcds%>%unlist%>%split(.,protids_trs$ms_id[match(names(.),protids_trs$transcript_id)])%>%reduce%>%unlist


# bams<-Sys.glob(here('pipeline/star/data/*/*.bam'))%>%
#   str_subset(neg=T,regex('transcript'))%>%
#   str_subset(neg=T,regex('test'))

# testbam<-bams%>%str_subset('E16_ribo_1')







# #Get the sums for all the bams
# allbamsiglist <- list()
# for(bam in bams){
#   message(bam) 
#  allbamsiglist %<>% append(list(bamsignals::bamCount(bam, shift=12,reduced_prot_cds , verbose=FALSE,mapqual=200)%>%split(names(reduced_prot_cds))%>%map_dbl(sum)))
# }

# length(allbamsiglist[[1]])==n_distinct(names(reduced_prot_cds))






# library(GenomicAlignments)
# summarizeOverlaps(reduced_prot_cds%>%split(names(.)),bams[1])

# allcountsold<-fread('feature_counts/all_feature_counts')


# ms2gid <- allms%>%select(ms_id,gene_name)%>%distinct%>%tail(n=-1)%>%left_join(gtf_gr%>%mcols%>%as.data.frame%>%distinct(gene_name,gene_id))

# ms2gid%>%filter(n_distinct(gene_name)>1)
# ms2gid%>%group_by(ms_id)%>%filter(n_distinct(gene_id)==1)



# testbamsig[allms$ms_id]

# testfcountmat <- allcountsold%>%.[match(ms2gid$gene_id,allcountsold$feature_id),]

# ms2gid%>%head

# cor(testfcountmat[[str_subset(colnames(allcountsold),str_replace(basename(bams[1]),'.bam',''))]] ,testbamsig[ms2gid$ms_id],use='complete')

# compmat<-data_frame(fcount =testfcountmat[[str_subset(colnames(allcountsold),str_replace(basename(testbam),'.bam',''))]] ,bamsig=testbamsig[ms2gid$ms_id])
# compmat%<>%cbind(ms2gid)
# compmat%>%  filter(fcount>100)%>%{cor(.$fcount,.$bamsig,use='complete')}

# ms2gid%>%head%>%cliplongtext

# cliplongtext<-function(df){
#   df%>%mutate_if(is.character,funs(str_extract(.,'.{0,10}')))
# }
# ms2gid%>%cliplongtext%>%head

# compmat%>%mutate(set=fcount > bamsig)%>%group_by(set)%>%tally

# pdf('tmp.pdf')
# qplot(data=compmat,x=log10(fcount),y=log10(bamsig))
# dev.off()

# compmat%>%filter(fcount>100,fcount/bamsig > 2)
# compmat%>%filter(fcount>100,bamsig/fcount > 2)
# compmat%>%filter(gene_name=='Orc3')





# #plenty of these persist even given a single source
# allprotids_trs%>%filter(source=='gencode')%>%group_by(ms_id)%>%summarise(n_distinct_gnames = n_distinct(gene_name))%>%group_by(n_distinct_gnames)%>%tally
# #
# allprotids_trs%>%filter(source=='gencode')%>%group_by(ms_id)%>%filter( n_distinct(gene_name)>1)%>%distinct(ms_id,gene_name,source)%>%ungroup%>%mutate(ms_id=str_extract(ms_id,'.{10}'))%>%as.data.frame
# #and with multiple sources....
# allprotids_trs%>%group_by(ms_id)%>%filter(!is.na(gene_name))%>%filter( n_distinct(gene_name)>1)%>%distinct(ms_id,gene_name)%>%
# allprotids_trs%>%group_by(ms_id)%>%filter(!is.na(gene_name))%>%filter( n_distinct(gene_name)>1)%>%distinct(ms_id,gene_name)%>%ungroup%>%mutate(ms_id=str_extract(ms_id,'.{0,10}'))%>%as.data.frame
# #
# allprotids_trs%>%group_by%>%summarise(is.na(allprotids_trs))



#




# #Now having done this linking
# allprotids%<>%left_join(gc_protiddf)

# #now, do we ever get conflicting gene names?
# #cases where we get a match via gencode pids, but it's different to the existing one
# gnameconflictdf<-allprotids_trs%>%filter(source=='gencodem12_metadata')%>%filter(!is.na(tr_gene_name))%>%filter(gene_name!=tr_gene_name)%>%distinct(gene_name,tr_gene_name)%>%mutate(ms_gname_n_in_gencode=!(gene_name%in%gtf_gr$gene_name))
# gnameconflictdf%>%head
# gnameconflictdf%>%.$ms_gname_n_in_gencode%>%table
# gnameconflictdf%>%write_tsv(here('pipeline/gnames_conflict_pids.tsv'))

# #How often do we get a gene_name where there wasn't one?
# #Pretty often actually.
# allprotids%>%filter(is.na(gene_name))%>%group_by(ms_id)%>%summarise(isnewmatch=any(!is.na(tr_gene_name)))%>%.$isnewmatch%>%table

# #How often do we get more than one gene name matching in the new gene name matches?
# multmatchdf<-allprotids%>%filter(!is.na(tr_gene_name))%>%group_by(ms_id)%>%summarise(hasmultmatch=n_distinct(tr_gene_name)>1)
# #somewhat quite often
# multmatchdf%>%filter(hasmultmatch)%>%write_tsv(here('pipeline/multmatch_pids.tsv'))






# #So we'll just use the tr_gene_name (i.e. from gencode pid tables) where possible

# #now cases where we DON"T match via gencode pids, but do have a gene name in the ms data
# #We can jsut get all the transcripts for teh gene....
# allprotids%>%colnames

# gnamematchonly<-allprotids%>%
#   filter(is.na(transcript_id))%>%
#   filter(!is.na(gene_name))%>%
#   select(ms_id,gene_name,uniprot_id)%>%
#   inner_join(annotrgnamedf%>%select(gene_name=tr_gene_name,transcript_id))%>%
#   left_join(protiddf,by='transcript_id')


# gnamematchonly%>%head
# table(unique(gnamematchonly$uniprot_id.x) %in% trs_w_uniprot$uniprot_id)
# trs_w_uniprot

# #okay so we now have multiple pid sets per gene name, sometimes.
# #Theory - most genes only have a single pidset which isn't shit

# # fdsafd - indorporate the bioconductor ensembl ojbect into our pid table as well.
# # fdafdsa - join up the above tablles now, select best ms row for each gene

# #get a df that measures this shitness
# timemissingsdf<-allms%>%group_by(ms_id,time)%>%summarise(missing=all(is.na(signal)))%>%summarise(n_missing_times=sum(missing))

# #so yeah it looks like most of them have a clear winner
# allprotids%>%filter(!is.na(tr_gene_name))%>%distinct(ms_id,tr_gene_name)%>%left_join(timemissingsdf)%>%group_by(tr_gene_name)%>%summarise(nmissingset=paste0(collapse=',',sort(n_missing_times)))%>%
#   .$nmissingset%>%table

# #but why are some always missing???
# timemissingsdf%>%filter(n_missing_times==5)%>%left_join(allms)%>%as.data.frame%>%head(40)

# #Probably those guys are a) super low scores edited in as nas and b) things only present in other fractions
# highsignal_pidsets<-timemissingsdf%>%filter(n_missing_times!=5)


# #So the strategy is to use our matches where we can, and if not, just use the protein coding transcripts
# #based on the gene name

# #is there a 1:1 match between gene_names and majority protein IDs?

# allprotids%>%group_by(ms_id,gene_name)
# allprotids%>%distinct(ms_id,gene_name)%>%nrow
# allprotids%>%distinct(ms_id)%>%nrow
# allprotids%>%distinct(gene_name)%>%nrow



# allprotids%>%filter(is.na(transcript_id))%>%.$gene_name%>%n_distinct



# allprotidsmatch<-allprotids%>%mutate(matches = uniprot_id %in% c(protiddf$uniprotkb_id,protiddf$swissprot_id))%>%group_by(ms_id)%>%summarise(matches=any(matches))
# allprotidsmatch$matches%>%table %>%{.[2]/sum(.)}#so 85% of genes with 
# #also add gene name data
# allprotidsmatch%<>%left_join(allms%>%select(gene_name,ms_id)%>%distinct)

# #second table - see if the newer encode version is any better
# protiddf2<-fread(here('annotation//gencode.vM21.metadata.SwissProt'),header=F)%>%set_colnames(c('transcript_id','uniprotkb_id','swissprot_id'))
# allprotidsmatch2<-allprotids%>%mutate(matches = uniprot_id %in% c(protiddf$uniprotkb_id,protiddf$swissprot_id,protiddf2$uniprotkb_id,protiddf2$swissprot_id))%>%group_by(ms_id)%>%summarise(matches=any(matches))
# allprotidsmatch2$matches%>%table %>%{.[2]/sum(.)}#very few extras if using newest gencode


# nomatchpids<-allprotidsmatch%>%filter(!matches)%>%.$ms_id

# mstable%>%filter(ms_id %in% nomatchpids)%>%head

# trids<-gtf_gr%>%subset(gene_name%in%'Ktn1')%>%.$transcript_id%>%unique

# protiddf%>%filter(str_replace(transcript_id,'\\.\\d+','') %in% trids)

# nomatchgenenames <- allprotidsmatch%>%filter(!matches)%>%.$gene_name

# nomatchtrids<-gtf_gr%>%subset(gene_name%in%nomatchgenenames)%>%.$transcript_id%>%unique




# ###Choosing transcripts...


# #The counts I'm getting with bamsignals is quite different form those I get with feature_counts....
# trs_redundant_cds <- cds%>%as.data.frame%>%subset(type=='CDS')%>%group_by(gene_id,transcript_id)%>%select(start,end)%>%nest%>%group_by(gene_id)%>% filter(n_distinct(data)<n_distinct(transcript_id))
# assert_that(nrow(trs_redundant_cds)==0)

# trs_redundant_cds_width<- cds%>%
#   as.data.frame%>%
#   subset(type=='CDS')%>%
#   group_by(gene_id,transcript_id)%>%
#   select(width)%>%
#   summarise(width=sum(width))%>%
#   group_by(gene_id)%>%
#   filter(n_distinct(width)<n_distinct(transcript_id))%>%
#   group_by(gene_id,width)%>%filter(n()==2)

#HOw can the ids have a match in teh gene names and NOT in the swissprot metadata???
#Let's look at an example, a gene which has  a gee