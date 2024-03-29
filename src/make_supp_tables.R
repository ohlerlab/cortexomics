library(tidyverse)
library(data.table)
library(magrittr)
library(readxl)
source('src/Preprocess/load_annotation.R')
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
	  str_replace_all('E175','E17')%>%
	  replace_na('all')
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
if(!exists('highcountgenes')) load('data/1_integrate_countdata.R')
#
stopifnot(s1counts%>%nrow%>%identical(22373L))
stopifnot(highcountgenes%>%length%>%identical(12228L))
list(
	count_data = s1counts,
	tpm_data = s1tpms,
	peptide_ms_data = fread(mspeptidefile),
	protein_ms_data = fread(msgenefile),
	matched_ms_data = readRDS('data/sel_ms_mat.rds')%>%as.data.frame%>%rownames_to_column('gene_id'),
	highcountgenes = tibble(gene_id=highcountgenes,gene_name=highcountgnms)
)%>%
map(dftpconv)%>%
openxlsx::write.xlsx( file = "tables/S1.xlsx")

highcountgenes <- readxl::read_xlsx('tables/S1.xlsx','highcountgenes')%>%.$gene_id

#S2
message('S2 - writing fold changes...')
allxtailtimes <- Sys.glob('pipeline/xtail/*')%>%
	str_extract('(?<=xtail_).*(?=.txt)')%>%
	str_subset(neg=T,'_v_')
allxtail <- Sys.glob('pipeline/xtail/*')%>%
	str_subset(neg=T,'_v_')%>%
	setNames(allxtailtimes)%>%
	lapply(fread)%>%
	map_df(.id='time',.%>%
		{colnames(.)%<>%
			str_replace_all('E13','T1')%>%
			str_replace_all('(E|P)\\d+','T2');
		.})
stepallxtail <- Sys.glob('pipeline/xtail/*')%>%
	str_subset('_v_')%>%
	setNames(allxtailtimes)%>%
	lapply(fread)%>%
	map_df(.id='time',.%>%
		{
			tp1 = colnames(.)%>%str_extract('.*_log2TE')%>%
				na.omit%>%head(1)%>%str_replace('_log2TE','')
			x=.
				colnames(x)%<>%
				str_replace_all(tp1,'T1')%>%
				str_replace_all('(E|P)\\d+','T2');
		x})
allxtail$gene_id = gnm2gidv[allxtail$gene_name]
stepallxtail$gene_id = gnm2gidv[stepallxtail$gene_name]
S2=list(
	xtail = allxtail,
	xtail_stepwise = stepallxtail,
	startstop_eff_lfc = fread(here('tables/ribo_position_effect.tsv')),
	limma_count_lfc = readRDS(here('data/countcontr_df.rds'))%>%
		mutate(time=replace_na(time,'initial')),
	proDA_ms_lfc_proDA = readRDS(here('data/contrdf.rds')),
	limma_count_stepwise_lfc = readRDS(here('data/stepcountcontrdf.rds'))%>%
		mutate(time=replace_na(time,'initial')),
	proDA_ms_stepwise_lfc = readRDS(here('data/stepcontrdf.rds'))
)%>%
# map(dftpconv)%>%
map(~filter(.,gene_id%in%highcountgenes))%>%
writexl::write_xlsx( path = "tables/S2.xlsx")

#numbers of start up/down
readxl::read_xlsx(  "tables/S2.xlsx",2,col_types=c(time='text'))%>%
	mutate(st_sig = str_pvalue<0.05,st_up=str_lfc>0,st_down=str_lfc<0)%>%
	filter(st_sig)%>%group_by(st_up,st_down)%>%tally
readxl::read_xlsx(  "tables/S2.xlsx",2,col_types=c(time='text'))%>%
	mutate(st_sig = str_pvalue>0.05,st_up=str_lfc>0,st_down=str_lfc<0)
readxl::read_xlsx(  "tables/S2.xlsx",2,col_types=c(time='text'))%>%
	mutate(end_sig = end_pvalue<0.05,end_up=end_lfc>0,end_down=end_lfc<0)%>%
	filter(end_sig)%>%group_by(end_up,end_down)%>%tally

allcodsigmean_isomerge%>%
	filter(time=='E13',fraction=='Total')%>%


#S3
trnacodonstats<-read_tsv('tables/tRNA_decoder_data.tsv')
codonoccs <- read_tsv('tables/codonoccs_final.tsv')
trnacodonstats$rep%<>%str_replace('rep','')%>%as.numeric
codonstats <- codonoccs%>%
	left_join(trnacodonstats%>%filter(fraction=='Total'))
codonstats%<>%select(-common,-iscodon,-Ct,-sample,-w_cAI,-fraction)
s3cols = c("time", "rep", "codon", "p_site_occ", "a_site_occ",
"anticodon", "abundance", "weightedusage", "availability", "AA", 
"freq", "tAI")
stopifnot(setequal(codonstats%>%colnames,s3cols))
GENETIC_CODEnostop = GENETIC_CODE[GENETIC_CODE!='*']
#all non stop codons have occupancies
stopifnot(codonoccs%>%filter(time=='E13')%>%.$codon%>%
	DNAStringSet%>%translate%>%as.character%>%setequal(GENETIC_CODEnostop))
#53 non-stop codons have tRNA measurements
stopifnot(trnacodonstats%>%filter(time=='E13',is.finite(abundance))%>%
	filter(codon%in%names(GENETIC_CODEnostop))%>%.$codon%>%n_distinct%>%`==`(50))
stopifnot(codonstats%>%filter(time=='E13',is.finite(abundance))%>%
	filter(codon%in%names(GENETIC_CODEnostop))%>%.$codon%>%n_distinct%>%`==`(50))
#
fread('tables/isodecoder_data.tsv')%>%addcodon%>%filter(time=='E13',codon=='AAT')
list(
	isodecoder_data = fread('tables/isodecoder_data.tsv'),
	codon_level_data =  codonstats,
	psite_offsets = fread('ext_data/offsets_rustvar.tsv')
)%>%
map(dftpconv)%>%
openxlsx::write.xlsx( file = "tables/S3.xlsx")
read_xlsx("tables/S2.xlsx",'limma_count_lfc')%>%head%>%as.data.frame

#S4
cluster_goterms <- read_tsv('tables/cluster_go.tsv')
cluster_genes = read_tsv('tables/gene_clusters.tsv')%>%arrange(cluster)
cluster_genes$gene_id = gnm2gidv[cluster_genes$gene_name]
cluster_genes%<>%filter(gene_id %in% highcountgenes)
#output clustering reuslts
list(
	cluster_genes = cluster_genes ,
	cluster_goterms = cluster_goterms,
	traj_model_table = read_tsv('tables/traj_model_df.tsv'),
	pihalf_est_table = read_tsv('tables/est_pihalf_v_mcshane.tsv')
)%>%
openxlsx::write.xlsx( file = "tables/S4.xlsx")