foldchangetblall_spread%>%head
# 



signnum=1

foldchangecatdf$sigstatus%>%table
itime = 'P0'




changedir = 1
foldchangecatdf%>%head
translsigvect <- 
	foldchangecatdf %>% 
	group_by(gene_id)%>%
	# semi_join(ribo_TE_tbl%>%filter(highcount)%>%left_join(ms_id2protein_id))%>%
	mutate(sig = 
		(translational_xtail_adj.P.Val%>%`<`(0.05) ) & 
		(abs(translational_xtail_logFC)>0.32) & 
		(sign(translational_xtail_logFC)%in%changedir)
	)%>%
	summarise(any(sig))

translsigvect = translsigvect[[2]]%>%setNames(translsigvect[[1]])
translsigvect%>%as.numeric%>%as.factor%>%setNames(names(translsigvect))%>%table

myontology='BP'
my_stat='Elim'
my_group='All_translational'
my_dir='up'
my_background='all'

# pdf(str_interp('plots/figures/figure1/'))

go_data <- new("topGOdata",
               ontology = "BP",
               allGenes = translsigvect%>%as.numeric%>%as.factor%>%setNames(names(translsigvect)),
               nodeSize = 50,
               annotationFun = annFUN.org,
               mapping = "org.Mm.eg",
               ID = "ensembl")

results <- runTest(go_data, algorithm = "elim", statistic = "fisher")
results.tab <- GenTable(object = go_data, elimFisher = results,topNodes = 200)
results.tab%<>%mutate(Enrichment = Significant / Expected )
results.tab%<>%mutate(elimFisher = as.numeric(elimFisher) )
results.tab%<>%dplyr::mutate(gene_ratio = Significant/Annotated)


results.tab%>%arrange(elimFisher)%>%head

results.tab%>%filter(Term%>%str_detect('protein stab'))

results.tab%>%filter(Term%>%str_detect('transl'))

results.tab%>%filter(Term%>%str_detect('protein ubiq'))

results.tab%>%filter(Term%>%str_detect('protein autoph'))

results.tab%>%filter(Term%>%str_detect('of GTPase'))

results.tab%>%filter(Term%>%str_detect('of GTPase'))
results.tab%>%filter(Term%>%str_detect(''))

results.tab%>%filter(Term%>%str_detect('exocytos'))

pdf('tmp.pdf',w=12)
plot_go_enrich(results.tab,sort_var = 'elimFisher',
	paste0(''))
dev.off()
message(normalizePath('tmp.pdf'))
