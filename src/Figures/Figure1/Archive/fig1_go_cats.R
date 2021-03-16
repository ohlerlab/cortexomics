# foldchangetblall_spread%>%head
# 
library(topGO)
source("src/R/Functions/go_term_funcs.R")


# signnum=1

# foldchangecatdf$sigstatus%>%table
# itime = 'P0'







countexprdata <- readRDS(here('data/fig1countexprdata_w_high.rds'))
# highcountgenes <- fData(countexprdata)%>%as.data.frame%>%filter(highcount)%>%.$gene_id

allTEchangelists<-lapply(c(-1,1),function(changedir){
	translsigvect <- 
		foldchangecatdf %>% 
		group_by(gene_id)%>%
		ungroup%>%
		filter(gene_id %in% highcountgenes)%>%
		mutate(sig = 
			(translational_xtail_adj.P.Val%>%`<`(0.05) ) & 
			(abs(translational_xtail_logFC)>0.32) & 
			(sign(translational_xtail_logFC)%in%changedir)
		)%>%
		group_by(gene_id)%>%
		summarise(any(sig))
	translsigvect = translsigvect[[2]]%>%setNames(translsigvect[[1]])
	# translsigvect%>%as.numeric%>%as.factor%>%setNames(names(translsigvect))%>%table
	translsigvect
})

allTEchangedf <- allTEchangelists%>%setNames(c('down','up'))%>%lapply(enframe,'gene_id','sig')%>%bind_rows(.id='direction')
allTEchangedf %<>% spread(direction,sig)
#add in few falses so this matches our 'highcount' set
allTEchangedf <- data.frame(gene_id=highcountgenes)%>%left_join(allTEchangedf)%>%mutate_at(vars(down,up),list(~ replace_na(.,0)))
stopifnot(setequal(allTEchangedf$gene_id,highcountgenes))
allTEchangedf%<>%safe_left_join(ids_nrgname%>%distinct(gene_id,gene_name))

stopifnot(1==(allTEchangedf%>%filter(gene_name=='Satb2')%>%.$up))
stopifnot(1==(allTEchangedf%>%filter(gene_name=='Flna')%>%.$down))
stopifnot(0==(allTEchangedf%>%filter(gene_name=='Nes')%>%.$down))

allTEchangedf%>%group_by(gene_id)%>%filter(any(up==1) & any(down==1))


myontology='MF'
my_stat='Elim'
my_group='All_translational'
dircol='up'
my_background='all'
minnodesize=20

# for(myontology in c('BP','MF','CC')) {
for(myontology in c('MF')) {
  for(dircol in c('down')){
    go_data <- new("topGOdata",
                   ontology = myontology,
                   allGenes = allTEchangedf[[dircol]]%>%as.numeric%>%as.factor%>%setNames(allTEchangedf$gene_id),
                   nodeSize = minnodesize,
                   # annotationFun = annFUN.db
                   annotationFun = annFUN.org,
                   mapping = "org.Mm.eg",
                   ID = "ensembl"
                  )

    results <- runTest(go_data, algorithm = "elim", statistic = "fisher")

    results.tab <- GenTable(object = go_data, elimFisher = results,topNodes = 100)
    results.tab%<>%mutate(Enrichment = Significant / Expected )
    results.tab%<>%mutate(elimFisher = as.numeric(elimFisher) )
    results.tab%<>%dplyr::mutate(gene_ratio = Significant/Annotated)

    goplotfile=str_interp('plots/figures/figure1/ggoplot_minsize${minnodesize}_${myontology}_${my_stat}_${dircol}_vs_${my_background}.pdf')
    pdf(goplotfile,w=6)
    print(plot_go_enrich(results.tab,sort_var = 'elimFisher',str_interp('minsize_${minnodesize}_${myontology}_${my_stat}_${dircol}_vs_${my_background}')))
    dev.off()
    message(normalizePath(goplotfile))

  }
}


# results.tab%>%filter(Term%>%str_detect('protein stab'))

# results.tab%>%filter(Term%>%str_detect('transl'))

# results.tab%>%filter(Term%>%str_detect('hemo'))


# go_data%>%str

# results.tab%>%filter(Term%>%str_detect('protein ubiq'))

# results.tab%>%filter(Term%>%str_detect('protein autoph'))

# results.tab%>%filter(Term%>%str_detect('of GTPase'))

# results.tab%>%filter(Term%>%str_detect('of GTPase'))
# results.tab%>%filter(Term%>%str_detect(''))

# results.tab%>%filter(Term%>%str_detect('exocytos'))




# ###TODO, fix on matt's results





# proteinstabgoid <- 
# hemogoid <- 'GO:0060216'

# assayData(countexprdata)$exprs[ms_id2protein_id%>%filter(gene_id%in%genesInTerm(go_data)[[hemogoid]])%>%.$protein_id,]


# fData(countexprdata)[ms_id2protein_id%>%filter(gene_id%in%genesInTerm(go_data)[[hemogoid]])%>%.$protein_id,]%>%filter(is_gid_highest)

# fData(countexprdata)[ms_id2protein_id%>%filter(gene_id%in%genesInTerm(go_data)[[hemogoid]])%>%.$protein_id,]%>%filter(is_gid_highest)%>%.$gene_name

# fData(countexprdata)[ms_id2protein_id%>%filter(gene_id%in%genesInTerm(go_data)[[hemogoid]])%>%.$protein_id,]%>%filter(is_gid_highest)%>%mutate(sig = gene_id %in% translsigvect[gene_id])

# hemogene_ids = fData(countexprdata)[ms_id2protein_id%>%filter(gene_id%in%genesInTerm(go_data)[[hemogoid]])%>%.$protein_id,]%>%filter(is_gid_highest)%>%filter(translsigvect[gene_id])
# hemogene_ids$gene_name

# translsigvect[hemogene_ids]


# GTOGO%>%head

# GTOGO%>%filter(go_id==hemogoid )


# hemogene_ids %>%mutate( otherdb_agrees = hemogene_ids$gene_id %in% (GTOGO%>%filter(go_id==hemogoid )%>%.$ensembl_gene_id))




# ################################################################################
# ########Here's what i based the old analysis on
# ################################################################################
# set.seed(1234)
# require(org.Mm.eg.db)
# require(DBI)
# require(topGO)

# # select a random list of gene symbol
# x <- unique(unlist(as.list(org.Mm.egENSEMBL)))
# names(x)=x
# # format  this list for topGO
# genesOfInterest <- names(translsigvect)[translsigvect]
# geneList = x
# geneList[!geneList %in% genesOfInterest] <- 0
# geneList[geneList %in% genesOfInterest] <- 1
# geneList = factor(geneList)
# table(geneList)

# # Create topGO object
# GOdata_BP = NULL
# GOdata_BP <-
#   new(
#     "topGOdata",
#     ontology = "BP",
#     allGenes = geneList,
#     description = "Test",
#     nodeSize = 50,
#     annot = annFUN.org,
#     mapping = "org.Mm.eg.db",
#     ID = "ENSEMBL"
#   )

# #?runTest
# results <- runTest(GOdata_BP, algorithm = "elim", statistic = "fisher")
# results.tab <- GenTable(object = GOdata_BP, elimFisher = results,topNodes = 100)
# results.tab%<>%mutate(Enrichment = Significant / Expected )
# results.tab%<>%mutate(elimFisher = as.numeric(elimFisher) )
# results.tab%<>%dplyr::mutate(gene_ratio = Significant/Annotated)
# results.tab%>%filter(Term%>%str_detect('protein stab'))


# pdf('tmp.pdf',w=12)
# plot_go_enrich(results.tab,sort_var = 'elimFisher',
# 	paste0(''))
# dev.off()
# message(normalizePath('tmp.pdf'))

