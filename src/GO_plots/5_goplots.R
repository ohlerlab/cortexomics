
source("src/R/Functions/go_term_funcs.R")

foldchangecatdf <- readRDS(here('data/foldchangecatdf.rds'))

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

allsigstats=foldchangecatdf$sigstatus%>%unique

background_filterfuncs = list(
			.%>%filter(gene_id %in% highcountgenes),
			.%>%filter(any(!sigstatus %in% 'No Sig Change'))
)%>%setNames(c('vs_all','vs_changing'))

mydir='up'
my_cat = "All Translational"
my_ontology='MF'
background='vs_all'
my_stat='Fisher.elim'
f1cats <- c('Translational Only','Concurrent Change','Compensating Change')
my_catnm <- 'All Translational'
my_cat = if (my_cat=='All Translational')  f1cats else  my_cat

for(mydir in c('up','down')){

		updownsign = case_when(mydir=='up' ~ 1,mydir=='down' ~ -1, TRUE ~ c(-1,1))
		# my_catnm = my_cat
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
		  results     <- runTest(goob, algorithm = 'elim', statistic = "fisher")
		  results.tab <- GenTable(object = goob, elimFisher = results,topNodes = 100)
		  results.tab%<>%mutate(Enrichment = Significant / Expected )
		  results.tab%<>%mutate(elimFisher = elimFisher%>%str_replace("<","")%>%as.numeric )
		  results.tab%<>%mutate(elimFisher = as.numeric(elimFisher) )
		  results.tab

		 p = results.tab%>%arrange(-elimFisher)%>%tail%>%
		 	mutate(Term=as_factor(Term))%>%
		 	ggplot(aes(x=Term,y=-log10(elimFisher),fill=Enrichment))+
		 	stat_identity(geom='bar')+
		 	geom_text(aes(label=GO.ID,y=4,color=I('white')))+
			ggplot2::scale_fill_continuous(
			      # limits = coltranslims,
			      # breaks = round(seq(from=min(coltranslims),to=max(coltranslims),len=n_colbreaks)),
			      # labels = col_label_format,
			      high = '#231f20',low='#d1d1d1'
			    )+
		 	coord_flip()+
		 	theme(axis.text.x=element_text(vjust=-1))+
		 	theme_bw()
		stopifnot(length(my_catnm)==1)
		plotfile <- paste0('plots/GO_plots/ggoplot_',as.character(my_ontology),'_',my_stat,'_',str_replace_all(my_catnm,' ','_'),'_',mydir,'_',background,'.pdf')
		pdf(plotfile,w=12,h=4)
		print(p)
		dev.off()
		message(normalizePath(plotfile))
		
}

if(F){

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

stop()


paramcombs= paramcombs[8,]
plotfiles <- mclapply(mc.cores=10,1:nrow(paramcombs),function(i){
	with(paramcombs[i,],{
		browser()
		paramcombs
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
		  results     <- runTest(goob, algorithm = 'elim', statistic = "fisher")
		  results.tab <- GenTable(object = goob, elimFisher = results,topNodes = 100)
		  results.tab%<>%mutate(Enrichment = Significant / Expected )
		  results.tab%<>%mutate(elimFisher = elimFisher%>%str_replace("<","")%>%as.numeric )
		  results.tab%<>%mutate(elimFisher = as.numeric(elimFisher) )
		  results.tab

		# titlestring <- 		paste0(as.character(my_ontology),'_',my_stat,'_',str_replace_all(my_catnm,' ','_'),'_',mydir,'_',background)
		# p = plot_go_enrich(go_table=goob2plot$table[['elim']],
		# 		sort_var='Fisher.elim',
		# 		title_str=titlestring)
		
		  # #now plot
		  # go_plot <-
		  #   ggplot(got_table,aes(size = log2(Significant) , y = Term, color = -log10(p_value), x = gene_ratio)) + 
		  #   geom_point() + 
		  #   ggplot2::theme_bw() + 
		  #   ggplot2::labs(size = "# significant  genes", y = "go term", color = "p-value", x = "gene-ratio", title=title_str) + 
		  #   #apply transformed color scale
		  #   ggplot2::scale_color_continuous(
		  #     limits = coltranslims,
		  #     breaks = round(seq(from=min(coltranslims),to=max(coltranslims),len=n_colbreaks)),
		  #     labels = col_label_format,high = 'darkgreen',low='lightgreen'
		  #   )+
		  #   #apply transformed size scale
		  #   ggplot2::scale_size_continuous(
		  #     limits = sizetranslims,
		  #     breaks = round(seq(from=min(sizetranslims),to=max(sizetranslims),len=n_sizebreaks)),
		  #     labels = size_label_format
		  #   )+
		  #   theme(axis.text=element_text(size=7), axis.title=element_text(size=14,face="bold"))
		  # #and return


		 p = results.tab%>%arrange(-elimFisher)%>%tail%>%
		 	mutate(Term=as_factor(Term))%>%
		 	ggplot(aes(x=Term,y=-log10(elimFisher),fill=Enrichment))+
		 	stat_identity(geom='bar')+
		 	geom_text(aes(label=GO.ID,y=4,color=I('white')))+
		 	scale_x_discrete('')+
		 	coord_flip()+
		 	theme(axis.text.x=element_text(vjust=-1))+
		 	theme_bw()
		stopifnot(length(my_catnm)==1)
		plotfile <- paste0('plots/GO_plots/ggoplot_',as.character(my_ontology),'_',my_stat,'_',str_replace_all(my_catnm,' ','_'),'_',mydir,'_',background,'.pdf')
		pdf(plotfile,w=12,h=4)
		print(p)
		dev.off()
		message(normalizePath(plotfile))

	})
})


}
