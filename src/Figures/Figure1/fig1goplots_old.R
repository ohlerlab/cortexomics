	# res_file <-
	#   file.path(modeldir,"intermediate_results",
	#             paste0(my_contrast_name, "_results.rds"))

	# isdeseq<-is.null(parameters$alt_res_file)

	# if(!isdeseq){

	#   alt_res_file = parameters$alt_res_file
	#   contr = unlist(my_contrasts[[my_contrast_name]])
	#   contrastisresfile = (length(contr)==1) && (is.character(contr)) && (file.exists(here(contr)))
	#   if(contrastisresfile) alt_res_file = contr
	#   res <- read_tsv(here::here(alt_res_file))
	#   saveRDS(res,res_file)
	#   modelcoeffs=NA
	# }else{
	#   res <- readRDS(file = res_file)
	#   modelcoeffs <- res@elementMetadata@listData$description[2]%>%str_extract('[^:]*$')%>%str_split(',')%>%.[[1]]
	# }



run_go_ensembl <- function(geneList,my_ontology){
set.seed(1234)
require(org.Mm.eg.db)
require(DBI)
require(topGO)

table(geneList)

# Create topGO object

GOdata <-
  new(
    "topGOdata",
    ontology = my_ontology,
    allGenes = geneList,
    description = "Test",
    annot = annFUN.org,
    mapping = "org.Mm.eg.db",
    ID = "ENSEMBL"
  )
GOdata
}


foldchangecatdf%>%head(1)%>%t

allsigstats=foldchangecatdf$sigstatus%>%unique

f1cats <- c('Translational Only','Concurrent Change','Compensating Change')
paramcombs <- expand.grid(
	# my_cat=c('Translational Only','Concurrent Change','Compensating Change'),
	my_cat='All Translational',
	mydir = c('up','down','updown'),
	# itime = 'P0',
	my_ontology = c('BP','CC','MF'),
	background = c('vs_all','vs_changing')[1],
	my_stat = 'Fisher.elim',stringsAsFactors=FALSE)

paramcombs%<>%filter((!mydir=='updown') | (my_cat%in%c('Translational Only','All Translational')))

background_filterfuncs = list(
			.%>%filter(gene_id %in% highcountgenes),
			.%>%filter(any(!sigstatus %in% 'No Sig Change'))
)%>%setNames(c('vs_all','vs_changing'))

mydir='down'
my_cat='All Translational'


plotfiles <- mclapply(mc.cores=10,1:nrow(paramcombs),function(i){
# plotfiles <- mclapply(mc.cores=20,1,function(i){

	# attach(paramcombs[i,])
	with(paramcombs[i,],{
		updownsign = case_when(mydir=='up' ~ 1,mydir=='down' ~ -1, TRUE ~ c(-1,1))
		my_catnm = my_cat
		my_cat = if (my_cat=='All Translational')  f1cats else  my_cat

		stopifnot('No Sig Change' %in% allsigstats)
		stopifnot(my_cat %in% allsigstats)

		backfilt = background_filterfuncs[[background]]

		godatadf <- foldchangecatdf%>%
			group_by(gene_id) %>%
			backfilt%>%
			mutate(padj = as.numeric(translational_xtail_adj.P.Val))%>%
			mutate(log2fc = as.numeric(translational_xtail_logFC))%>%
			summarise(is_sig = any( (sigstatus%in%my_cat) & (sign(log2fc)%in%updownsign)))

		godata<-godatadf%>%	{setNames(as.factor(as.numeric(.$is_sig)),.$gene_id)}
		
		stopifnot((c(0,1)%in%names(table(godata))))
		
		goob<-quietly(run_go_ensembl)(godata,my_ontology)$result

		goob2plot <- quietly(run_go_enrich_int)(goob,as.character(my_ontology),foldchangecatdf$gene_id%>%setNames(.,.))$result

		titlestring <- 		paste0(as.character(my_ontology),'_',my_stat,'_',str_replace_all(my_catnm,' ','_'),'_',mydir,'_',background)
		p = plot_go_enrich(go_table=goob2plot$table[['elim']],
				sort_var='Fisher.elim',
				title_str=titlestring)


		stopifnot(length(my_catnm)==1)
		plotfile <- paste0('plots/figures/figure1/goplot_',as.character(my_ontology),'_',my_stat,'_',str_replace_all(my_catnm,' ','_'),'_',mydir,'_',background,'.pdf')
		plotfile <- paste0('plots/figures/figure1/goplot_',as.character(my_ontology),'_',my_stat,'_',str_replace_all(my_catnm,' ','_'),'_',mydir,'_',background,'.pdf')
		pdf(plotfile,w=6,h=7)
		print(p)
		dev.off()
		message(normalizePath(plotfile))

	})
})

