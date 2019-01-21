#This script is going to take in a bunch of files from various other papers and output a GMT file
message('loading libraries')


CHUNKSIZE=1e5

message('...done')

#load arguments
argv <- c(
	geneids = 'ids.txt',
	outputfile = 'gene_set_enrichment/allgenesets.gmt'
)
argv <- commandArgs(trailingOnly=TRUE)[1:length(argv)]%>%setNames(names(argv))
for(i in names(argv)) assign(i,argv[i])
geneids

outputfolder<-dirname(outputfile)
outputfolder <- paste0(outputfolder,'/')
outputfolder%>%dir.create(showWarn=F,rec=TRUE)

ids<-fread(geneids)

file.remove(outputfile)

add_to_gmt<-function(catname,catdescription,catgenes,gmt=outputfile){
	cat(
		paste0(catname,'\t',catdescription,'\t',paste0(collapse='\t',catgenes),'\n'),
		file = gmt,
		append=TRUE
	)
}



####Molyneux cell specificity data
molyneuxfile<-'../ext_data/molyneux_etal_2015_tS3.xlsx'
molyneuxdata <- readxl::read_excel(molyneuxfile,sheet=1,skip=1)
molyneuxdata[1,]%>%t

#lots of data in this table, fold changes an pvalues and some clusteirng I can't remember
#the point of. 
neurites <- '../ext_data/neurites_zappulo_etal_2017.csv'
neurites%<>%fread(skip=2)


#match by name or gene id
neuriteidmatch <- match(neurites$gene_id, ids$gene_id) %>% 
	ifelse(!is.na(.),.,match(neurites$gene_name,ids$gene_name))
stopifnot(mean(is.na(neuriteidmatch))<0.06)
neurites$gene_id <- ids$gene_id[neuriteidmatch]
neurites%<>%filter(!is.na(gene_id))

neurite_ribo_gids<-neurites%>%filter(RiboSeq_padj_Neurite_Soma<0.05)%>%.$gene_id
neurite_ribo_gids<-neurites%>%filter(RiboSeq_padj_Neurite_Soma<0.05)%>%.$gene_id
neurite_te_gids<-neurites%>%filter(`Prot log2FC >1, RNA log2FC <0`==1|`Prot log2FC >1, 1> RNA log2FC >0`==1)%>%.$gene_id
neurite_te_gids%>%colnames

add_to_gmt('Neurites_Riboseq_Zappulo',
	'Genes with a Padj<0.05 for riboseq in neurite/soma Zapullo etal 2017',
	neurite_ribo_gids,
)
add_to_gmt('Neurites_Protl2fc_gt_RNA_Zappulo',
	'Genes with a Prot l2fc of >1, and RNA lf2c of < 1 in Zappulo',
	neurite_te_gids
)
file.info(outputfile)

##5' TOPs
ftoptrs <- fread('cat fpTOP_scan/fimo.gff  | cut -f1   | grep -v "#"',header=F) %>%set_colnames('transcript_id')
ftoptrs$transcript_id %<>% str_replace('\\(.*$','')
gtrmap <-'rsemref/gene_transcript_map.tsv' %>% fread(header=F) %>% set_colnames(c('gene_id','transcript_id'))
fptopgids <- ftoptrs %>% left_join(gtrmap)%>%.$gene_id%>%unique

add_to_gmt('Gids with 5\' TOP',
	'Genes with at least one transcript whose TSS(+/- 100bp) matches the TOP PWM at p<0.05',
	fptopgids
)
file.info(outputfile)


#p-bodies
p_bodies_prot_tab<-'../ext_data/p_bodies_hubstenberger_etal_2017/prot_enrich.csv'%>%
	fread%>%
	set_colnames(c("gene_name", "UniProtKB entry", "log2fc", 
"pval", 
"pre-sorted protein relative abundance (normalized to total spectra)", 
"sorted protein relative abundance (normalized to total spectra)", 
"DDX6 immunoprecipitation with RNAse inhibitor, Mascot score (Ayache et al., 2015)", 
"DDX6 immunoprecipitation with RNAse inhibitor, peptide number (Ayache et al., 2015)", 
"DDX6 immunoprecipitation with RNAse, Mascot score (Ayache et al., 2015)", 
"DDX6 immunoprecipitation with RNAse, peptide number (Ayache et al., 2015)"
)
)
p_bodies_prot_tab%<>%filter(pval<0.05)
p_bodies_prot_tab$gene_id<-ids$gene_id[match(p_bodies_prot_tab$gene_name,toupper(ids$gene_name))]
stopifnot(p_bodies_prot_tab$gene_id%>%is.na%>%mean%>%`<`(0.05))
p_bodies_prot_tab%<>%filter(!is.na(gene_id))

add_to_gmt('P_body_proteins_Hubstenberger_etal_2017',
	'Proteins that are enriched in human p-bodies',
	p_bodies_prot_tab$gene_id
)

p_bodies_rna_tab<-'../ext_data/p_bodies_hubstenberger_etal_2017/RNA_enrich.csv'%>%
	fread%>%
	set_colnames(
		c("Ensembl Gene ID", "gene_name", "log2fc", "pval", "padj", 
		"sorted P-body replicate 1 (count per million reads (CPM))", 
		"sorted P-body replicate 2 (count per million reads (CPM))", 
		"sorted P-body replicate 3 (count per million reads (CPM))", 
		"sorted P-body replicate average (count per million reads (CPM))", 
		"pre-sorted fraction replicate 1 (count per million reads (CPM))", 
		"pre-sorted fraction replicate 2 (count per million reads (CPM))", 
		"pre-sorted fraction replicate 3 (count per million reads (CPM))", 
		"pre-sorted fraction replicate average (count per million reads (CPM))", 
		"sorted P-body replicate 1 (mapped read number)", "sorted P-body replicate 2 (mapped read number)", 
		"sorted P-body replicate 3 (mapped read number)", "pre-sorted fraction replicate 1 (mapped read number)", 
		"pre-sorted fraction replicate 2 (mapped read number)", "pre-sorted fraction replicate 3 (mapped read number)"
))


setdiff(p_bodies_rna_tab$gene_name,toupper(ids$gene_name))
intersect(p_bodies_rna_tab$gene_name,toupper(ids$gene_name))

p_bodies_rna_tab%<>%filter(padj<0.05)
p_bodies_rna_tab$gene_id<-ids$gene_id[match(p_bodies_rna_tab$gene_name,toupper(ids$gene_name))]
# stopifnot(p_bodies_rna_tab$gene_id%>%is.na%>%mean%>%`<`(0.05))
p_bodies_rna_tab%<>%filter(!is.na(gene_id))

add_to_gmt('P_body_RNAs_Hubstenberger_etal_2017',
	'RNAs that are enriched in human p-bodies',
	p_bodies_rna_tab$gene_id
)


#Let's look at poly a genes
altpagenes <- fread('../ext_data/ctag_paperclip_hwang_et_al_2017.csv')
altpagenes%<>%left_join(ids)

add_to_gmt('The few genes identified as having alternative polyadenylation in hwang et al 2017',
	'',
	altpagenes$gene_id
)

molyneuxdata<-readxl::read_xlsx('../ext_data/molyneux_etal_2015_tS3.xlsx',skip=1)%>%
	mutate(tissue=
		case_when(
				cluster%in%c(5,0,15,10) ~ 'CPN', 
				cluster %in% c(16,7,18,11,3) ~ 'SCPN',
 				cluster %in% c(6,9,2,12,17) ~ 'CThPN', 
 				TRUE ~ 'celltype_indep')
	)
molyneuxdata%<>%rename(gene_name=gene_id)
molyneuxdata%<>%left_join(ids)

add_to_gmt('From Molyneux et al - genes in clusters specific to CPN ( 5,0,15,10)',
	'These genes are in clusters specific to corticothalamic projection neurons, i followd their figure 2 but excluded things that looked specific to 2 cell types',
	molyneuxdata%>%filter(tissue=='CPN')%>%.$gene_id
)

add_to_gmt('From Molyneux et al - genes in clusters specific to SCPN ( 16,7,18,11,3 )',
	'These genes are in clusters specific to subcerebral projection neurons, i followd their figure 2 but excluded things that looked specific to 2 cell types',
	molyneuxdata%>%filter(tissue=='SCPN')%>%.$gene_id
)

add_to_gmt('From Molyneux et al - genes in clusters specific to CthPN ( 6,9,2,12,17 )',
	'These genes are in clusters specific to corticospinal motor neurons, i followd their figure 2 but excluded things that looked specific to 2 cell types',
	molyneuxdata%>%filter(tissue=='CThPN')%>%.$gene_id
)



#let's quickly get some clusters from the yuzwa data
# library(biomaRt)

# humanmart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")



# listAttributes(humanmart,page='feature_page')
# attributePages(humanmart)%>%head

# listDatasets(mousemart,1)
# # > attributes =

# c("ensembl_gene_id","mmusculus_homolog_ensembl_gene","mmusculus_homolo
# g_perc_id_r1")

# # > attributes=c(attributes,"mmusculus_homolog_orthology_type",
# "mmusculus_homolog_subtype", "mmusculus_homolog_perc_id")

# # >  orth.mouse = getBM( attributes,filters="with_homolog_mmus",values
# # =TRUE,
# mart = human, bmHeader=FALSE)

# # > dim(orth.mouse)


# humanmart = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# listDatasets(mousemart)

# library(biomaRt)
# genes <- df$ensembl_gene_id


# allttrs<-listAttributes(mousemart,page='feature_page')%>%filter(description%>%str_detect(regex(ignore.case=T,'ID')))%>%.[[1]]
# # listFilters(mousemart)%>%filter(description%>%str_detect(regex(ignore.case=T,'ID')))
# # %>%h

# mmart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

# symbol <- getBM(filters = "ensembl_gene_id",
#                 attributes = 'ensembl_gene_id',
#                 values = neurite_te_gids%>%head(1), 
#                 mart = mmart)

# mousemart = useMart("ENSEMBL_MART_MOUSE",dataset='m129s1svimj_gene_ensembl')


# 'UniProtKB entry'
# humanmart%>%listAttributes
# allttrs<-listAttributes(humanmart,page='feature_page')%>%filter(description%>%str_detect(regex(ignore.case=T,'uniprot')))%>%.[[1]]

# listMarts%>%args
# uniProt <- useMart("unimart")

#  listMarts(host="www.biomart.org")

#  ensembl_entrez_tab<-org.Hs.eg.db::org.Hs.egENSEMBL%>%as.list%>%stack%>%set_colnames(c('ensembl','entrez'))
#  uniprot2entreztab<-org.Hs.eg.db::org.Hs.egUNIPROT%>%as.list%>%stack%>%set_colnames(c('uniprot','entrez')) 

#  p_bodies_prot_tab$UniProtK%>%str_replace('_HUMAN','')%>%is_in(uniprot2entreztab$uniprot)

# up <-  UniProt.ws(taxId=9606)
# ensemblgenes <- select(up,keys=p_bodies_prot_tab[[1]],columns='ENSEMBL',keytype='UNIGENE')


rmarkdown::render(
			"fast/groups/ag_ohler/work/dharnet_m/cortexomics/src/gage_rep.Rmd",
		knit_root_dir='fast/groups/ag_ohler/work/dharnet_m/cortexomics/pipeline',
		params=list(signalfile='xtail/xtail_P0.txt',gene_set_file='gene_set_enrichment/allgenesets.gmt',sigcol=NULL)
)