library(tidyverse)
library(data.table)
library(magrittr)
library(readxl)
add_descrip_lines <- function(musicres,descrips){
	descriplines = map_df(descrips,function(descrip){
		musicres%>%
		head(1)%>%
		as.data.frame%>%{.[,-1]=rep('',ncol(.)-1);.}%>%
		{.[1,1]=descrip;.}
	})
	musicres %<>% rbind(descriplines)
	musicres
}

#Table 1
tpconv = function(x) x %>% str_replace_all('E13','E12.5')%>%
	  str_replace_all('E145','E14')%>%
	  str_replace_all('E16','E15.5')%>%
	  str_replace_all('E175','E17')
dftpconv = function(x) x  %>% {colnames(.)%<>%tpconv;.}%>%
	{
		if('time'%in%colnames(.)).$time%<>%tpconv
		if('contrast'%in%colnames(.)).$contrast%<>%tpconv
		if('sample'%in%colnames(.)).$contrast%<>%tpconv
		.
	}

message('S1 - writing count and MS raw data...')
s1counts = readRDS('data/tx_countdata.rds')%>%
	.$counts%>%
	dftpconv%>%as.data.frame%>%
	rownames_to_column('gene_id')%>%
	select(gene_id,matches('ribo'))
s1tpms = readRDS('data/tx_countdata.rds')%>%
	.$abundance%>%
	dftpconv%>%as.data.frame%>%
	rownames_to_column('gene_id')%>%
	select(gene_id,matches('ribo'))
#
mspeptidefile = here('ext_data/MS_Data_New/Ages_Brain_PEP_summ.txt')
msgenefile = here('ext_data/MS_Data_New/Ages_Brain_PG_summ.txt')
#
stopifnot(s1counts%>%nrow%>%identical(22373L))
stopifnot(highcountgenes%>%length%>%identical(12228L))
list(
	count_data = s1counts,
	tpm_data = s1counts,
	peptide_ms_data = fread(mspeptidefile),
	protein_ms_data = fread(msgenefile),
	matched_ms_data = readRDS('data/sel_ms_mat.rds')%>%as.data.frame%>%rownames_to_column('gene_id'),
	highcountgenes = tibble(gene_id=highcountgenes,gene_name=highcountgnms)
)%>%
map(dftpconv)%>%
openxlsx::write.xlsx( file = "tables/S1.xlsx")


#S2
message('S2 - writing fold changes...')
allxtailtimes <- Sys.glob('pipeline/xtail/*')%>%str_extract('(?<=xtail_).*(?=.txt)')
allxtail <- Sys.glob('pipeline/xtail/*')%>%
	setNames(allxtailtimes)%>%
	lapply(fread)%>%
	map_df(.id='time',.%>%
		{colnames(.)%<>%
			str_replace_all('E13','T1')%>%
			str_replace_all('(E|P)\\d+','T2');
		.})
allxtail$gene_id = gnm2gid[[allxtail$gene_name]]

list(
	xtail = allxtail,
	startstop_eff_lfc = fread(here('tables/ribo_position_effect.tsv')),
	limma_count_lfc = readRDS(here('data/countcontr_df.rds')),
	proDA_ms_lfc_proDA = readRDS(here('data/contrdf.rds')),
	limma_count_stepwise_lfc = readRDS(here('data/stepcountcontrdf.rds')),
	proDA_ms_stepwise_lfc = readRDS(here('data/stepcontrdf.rds'))
)%>%
map(dftpconv)%>%
map(~filter(.,gene_id%in%highcountgenes))%>%
openxlsx::write.xlsx( file = "tables/S2.xlsx")


#S3
codonstats<-read_tsv('tables/tRNA_stat_df')
codonstats%<>%filter(fraction=='total')%>%
	select(-fraction,-common,-abundance_enrich,-availability_enrich)
list(
	isodecoder_data = fread('tables/isodecoder_data.tsv'),
	codon_level_data =  codonstats
)%>%
map(dftpconv)%>%
openxlsx::write.xlsx( file = "tables/S3.xlsx")


#S4
cluster_genes = read_tsv('tables/gene_clusters.tsv')%>%arrange(cluster)
cluster_genes$gene_id = gnm2gid[[cluster_genes$gene_name]]
cluster_genes%<>%filter(gene_id %in% highcountgenes)
cluster_goterms <- read_tsv('tables/cluster_go.tsv')
#output clustering reuslts
list(
	cluster_genes = cluster_genes ,
	cluster_goterms = cluster_goterms
)%>%
openxlsx::write.xlsx( file = "tables/S4.xlsx")


#S5
musicres <- fread('tables/musicdata.tsv')
musicres %<>% add_descrip_lines(
	c('# time - collection stage',
		'#assay - riboseq or rnaseq',
	'#ftset - the flashtag set from telley et al used as a signature',
	'#music_prop - estimated proportion according to MuSiC')
)
list(
	 music_proportion_data = musicres
)%>%
openxlsx::write.xlsx( file = "tables/S5.xlsx")
normalizePath("tables/S5.xlsx")%>%message