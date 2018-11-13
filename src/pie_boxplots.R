library(svglite)
library(readr)
library(Biostrings)
library(Rsamtools)
library(rtracklayer)
library(GenomicFeatures)
library(stringr)
library(data.table)
library(assertthat)

filter<-dplyr::filter
slice<-dplyr::slice
filter<-dplyr::filter
filter<-dplyr::filter


#pie charts of orf types
satannfiles <- Sys.glob('./../pipeline/SaTAnn/*/*Final_ORFs*')
satannorfs <- 
	# Sys.glob('SaTAnn/*/*Final_ORFs*')%>%
	satannfiles%>%
	setNames(.,basename(dirname(.)))%>%
	mclapply(load_objs)

#now filter out the weird GRanges columns, for now, and aggregate into one data table 
all_orfs <- satannorfs%>%map(.%>%.$ORFs_tx)

for (i in seq_along(all_orfs))all_orfs[[i]]$Protein%<>%as.character
orfs_dt <- all_orfs%>%map(.%>%{
	issimplelist <- mcols(.)%>%vapply(is,TRUE,'atomic')
	.[,issimplelist]%>%GR2DT
})%>%bind_rows(.id='sample')


ccdsorftypenames <- c('CCDS_ORFS','uORFs','non CCDS_ORFs','NMD','dORFs')
nccdsorftypenames <- c('Processed _transcript','other','lincRNA','Antisense','Retained Intron','Pseudogene')

sample_i = 'OD5P'

orfs_dt%<>%group_by(sample)%>%
	mutate(class = case_when(
		compatible_biotype %in%c('nonsense_mediated_decay') ~ 'NMD',
		ORF_category_Tx_compatible %in%c('dORF','overl_dORF') ~ 'dORFs',
		ORF_category_Tx_compatible %in%c('uORF') ~ 'uORF',
		ORF_category_Tx_compatible %in%c('novel') ~ 'non CCDS_ORFs',
		TRUE  ~ 'CCDS_ORFs',
))

orfs_dt%<>%group_by(sample)%>%
	mutate(nonCCDS_type = case_when(
		class != 'non CCDS_ORFs' ~ 'NA',
		transcript_biotype %in% c('lincRNA') ~ 'lincRNA',
		transcript_biotype %in%c('antisense') ~ 'Antisense',
		transcript_biotype %in%c('retained_intron') ~ 'Retained Intron',
		transcript_biotype %in%c('processed_transcript') ~ 'Processed Transcript',
		transcript_biotype %>%str_detect('pseudogene') ~ 'Pseudogene',
		TRUE  ~ 'other'
))
orfs_dt$nonCCDS_type%<>%replace(.,.=='NA',NA)
orfs_dt$nonCCDS_type%>%table
sample_i <- 'OD5P'
sampdt <- orfs_dt%>%filter(sample==sample_i)

pieboxfile<-str_interp('../plots/merged_summary2/${sample_i}_pie_box.pdf')
pieboxfile%>%dirname%>%dir.create

#make charts
pdf(pieboxfile,h=3.5,w=10)
par(oma=c(0,0,0,5))
piecols <- rainbow(sampdt$class%>%n_distinct) pieorder<-c(3, 4, 2L, 1L, 5L)
piedata<-sampdt$class%>%table%>%.[pieorder]
pie(
	piedata,
	piedata%>%{paste0(names(.),':',format(.,big.mark=','))},
	col=piecols[pieorder]
	# col=c('red','yellow','green','blue')
	)
legend('topright',legend=sampdt$class%>%unique,fill=piecols,cex=0.5)
#bxplots of length for classes
sampdt_nccds <- sampdt%>%filter(!is.na(nonCCDS_type))
#now the non ccds types
sampdt_nccds%>%ggplot(aes(x=nonCCDS_type,y=width,fill=nonCCDS_type))+geom_boxplot()+scale_y_log10(name='length')+theme_minimal()+theme(axis.text.x=element_text(vjust=0.5,angle=45))
#and the piechart
piecols <- rainbow(sampdt_nccds$nonCCDS_type%>%n_distinct)
pie(
	sampdt_nccds$nonCCDS_type%>%table,
	sampdt_nccds$nonCCDS_type%>%table%>%{paste0(names(.),':',format(.,big.mark=','))},
	col=piecols
	# col=c('red','yellow','green','blue')
	)
legend('topright',legend=sampdt_nccds$nonCCDS_type%>%unique,fill=piecols,cex=0.5)
#bxplots of length for nonCCDS_typees
sampdt_nccds%>%ggplot(aes(x=nonCCDS_type,y=width,fill=nonCCDS_type))+geom_boxplot()+scale_y_log10(name='length')+theme_minimal()+theme(axis.text.x=element_text(vjust=0.5,angle=45))
dev.off()
normalizePath(pieboxfile)


#
orfs_dt[1,]%>%t
orfs_dt$ORF_category_Gen%>%table
orfs_dt$compatible_biotype%>%table
orfs_dt$ORF_category_Tx%>%table
orfs_dt$ORF_category_Tx_compatible%>%table
orfs_dt$transcript_biotype%>%table
orfs_dt$gene_biotype%>%table
orfs_dt$gene_biotype%>%table


