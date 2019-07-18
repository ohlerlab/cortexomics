#!/usr/bin/env Rscript
message('loading libraries')
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(assertthat))
message('...done')

defaultargs <- c(
  ribodifffolder = './ribodiff',
  limmafolder = './exprdata/',
  outputfolder = 'ribodiff_limma_comp'
)
args <- coalesce(
  commandArgs(trailingOnly=TRUE)[1:length(defaultargs)]%>%setNames(names(defaultargs)),
  defaultargs
)
for(i in names(args)) assign(i,args[i])

#read in ribodiff
ribodiffresfiles <- Sys.glob(file.path(ribodifffolder,'/riboseqres_*.txt'))

#read in 
ribodiffcolsnms=c('feature_id','disper','p_value','adj_p_value','TE1','TE2','log2fc')
ribodiffcols=cols(
  col_character(),
  col_double(),
  col_double(),
  col_double(),
  col_double(),
  col_double(),
  col_double()
)
ribodiffcols$cols%<>%setNames(ribodiffcolsnms)

ribodiffcontrastobs <- ribodiffresfiles%>%map(read_tsv,skip=1,col_names=ribodiffcolsnms,col_types=ribodiffcols)

#read in the limma fcs
limmafcs<-fread(file.path(limmafolder,'limma_fold_changes.txtfull.txt'))

#compare them, timetpoint by timepoint....
timepoints = ribodiffresfiles%>%str_extract('(?<=_)(E|P)\\d+')
names(ribodiffcontrastobs) <- timepoints

for(timepoint in timepoints){
	ribodiffdf<-ribodiffcontrastobs%>%.[[timepoint]]%>%select(gene_id=feature_id,padj=adj_p_value,log2fc)
	limmacol <- limmafcs$coefficient%>%unique%>%str_subset(timepoint)%>%str_subset(':assayribo')
	limmadf <- limmafcs%>%filter(coefficient==limmacol)%>%select(gene_name=gene,log2fc=logFC,padj=adj.P.Val)


	ids<-fread('ids.txt')%>%distinct


	ribodiffdf$gene_id%>%n_distinct
	limmadf$gene_name%>%n_distinct

	limmadf$gene_name%>%n_distinct
	compdf<-limmadf%>%left_join(ids)%>%left_join(ribodiffdf)%>%left_join(ribodiffdf,by='gene_id')%>%
		set_colnames(c('gene_name','limma_fc','limma_pval','gene_id','ribodiff_pval','ribodiff_fc'))
	alpha=0.05

	compdf%<>%mutate(limma_ribodiff_sig = paste0('Rdiff:',limma_pval<alpha,'_','Limma:',ribodiff_pval<alpha))

	compdf$limma_ribodiff_sig%>%table

	ribodiffdf$ribodiff_pval
	compdf%>%head

	correlation = round(cor(compdf$ribodiff_fc,compdf$limma_fc),2)

	comp_plot<-qplot(data=compdf,x=ribodiff_fc,y=limma_fc,color=limma_ribodiff_sig)+
		geom_point()+
		scale_x_continuous('Ribodiff TE Change')+
		scale_y_continuous('Limma_Fold_Change')+
		ggtitle(paste0('Fold Change Correlation at ',timepoint,'\n between limma and Ribodiff:',correlation))

	normalizePath(plotfile)%>%message
	outdir <- 'ribodiff_limma_comp'
	outdir %>% dir.create
	
	plotfile<-paste0(outdir,'/',timepoint,'compscatter.pdf')
	ggsave(plot=comp_plot,file=plotfile)


}