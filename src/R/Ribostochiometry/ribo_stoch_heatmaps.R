suppressMessages(library(limma))
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(here))
library(gplots)
library(RColorBrewer)


# source(here('src/R/Load_data/make_expression_table.R'))

root<-here()

vlookup <- function(query,dicttable,key,vals){
	dict = dicttable%>%ungroup%>%distinct_(key,vals)
	stopifnot(!anyDuplicated(dict))
	data_frame(tmp=query)%>%
		left_join(dict,by=c('tmp'=key))%>%
		.[[vals]]
}

str_split_fast = function(x,sep=';') x %>% {str_split(.,sep,n = str_count(.,';')%>%max%>%add(1))}
#' setup, eval = TRUE
sep_element_in<-function(colonlist,ridssplit,sep=';'){
	assert_that(is.character(colonlist))
	assert_that(is.character(ridssplit))
	values<-colonlist%>%str_split_fast

	inds <- rep(seq_along(colonlist),lengths(values))

	values<-flatten_chr(values)

	data_frame(inds,values)%>%
		mutate(match = values %in% ridssplit)%>%
		group_by(inds)%>%
		summarize(match=any(match))%>%
		.$match

}


#load the iBAQ ms data
ms_tall <- Sys.glob(here('pipeline/ms_tables/ms_iBAQ*'))%>% map_df(fread)

message('looking up protein ID annotations')
#'get info on the ribosomal subu,nts from mats table
rids <- read_tsv(here('ext_data/riboprotids.tsv'))
lridssplit<-rids%>%filter(!`Mito-RP`,`RPL_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist
sridssplit<-rids%>%filter(!`Mito-RP`,`RPS_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist


ridssplit<-c(lridssplit,sridssplit)
rids%>%filter(!`Mito-RP`)%>%nrow
#get info on proteins so we can catagorize them
library(biomaRt)
# mart <- useMart(biomart = "ENSEMBL_MART_MOUSE")
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
listAttributes(mart)%>%filter(description%>%str_detect('Uni'))%>%filter(page!='homologs')
# ribogoterm <- "GO:0005840"
# lribogoterm <- 'GO:0015935'
# sribogoterm <- 'GO:0015934'
# riboprotids <- biomaRt::getBM(attributes = c("uniprot_gn"), 
#                  filters=('go'),values=ribogoterm,
#                  mart = mart)
# sriboprotids <- biomaRt::getBM(attributes = c("uniprot_gn"), 
#                  filters=('go'),values=sribogoterm,
#                  mart = mart)%>%.[[1]]
# lriboprotids <- biomaRt::getBM(attributes = c("uniprot_gn"), 
#                  filters=('go'),values=lribogoterm,
#                  mart = mart)%>%.[[1]]
transreggoterm <- "GO:0006417"
transregprotids <- biomaRt::getBM(attributes = c("uniprot_gn_id"), 
                 filters=('go'),values=transreggoterm,
                 mart = mart)
translationgoterm <- "GO:0006417"
translationprotids <- biomaRt::getBM(attributes = c("uniprot_gn_id"), 
                 filters=('go'),values=translationgoterm,
                 mart = mart)[[1]]
ebp1pid = ms_tall$Protein_IDs[match('Pa2g4',ms_tall$gene_name)]



# ms_tall_trans%<>%filter(!gene_name%>%str_detect('Hbs1l'))


#' Thi sis a title with inline R code `r foo`

#' First we load the list of protein IDs, handpicked by Matt, using only the small
#' or large subunits - no mitochondrial riboproteins  



#define ambigous protein groups as those which have elements that appear in more than one protein group
allpgroups <- ms_tall$Protein_IDs%>%unique
multids<-allpgroups%>%unique%>%str_split_fast(';')%>%unlist%>%table%>%keep(~ . > 1)%>%names
all_ambig_pgroups<-allpgroups%>%sep_element_in(multids)
library(data.table)




# ms_tall_trans<-ms_tall%>%mutate(ambig = sep_element_in(Protein_IDs,multids))


#' We pick out the iBAQ data, since LFQ isn't recommended for these kinds of enriched datasets
#' We also include Pa2g4.  
#+ select ribo data  , eval = TRUE, cache=TRUE
n_groups = n_distinct(ms_tall$Protein_IDs)
n_protid = ms_tall$Protein_IDs%>%str_split_fast(';')%>%unlist%>%n_distinct


message('Labeling protein IDs with annotations')
pids = ms_tall%>%ungroup%>%distinct(Protein_IDs)

pids%<>%mutate(pcat = case_when(
 (Protein_IDs==ebp1pid) ~ "Ebp1",
 (Protein_IDs=='P68040') ~ "P68040",
  sep_element_in(Protein_IDs,sridssplit) ~ "Rps",
  sep_element_in(Protein_IDs,lridssplit) ~ "Rpl",
  sep_element_in(Protein_IDs,translationprotids) ~ "translation-associated",
  TRUE ~ "other"
))
stopifnot("Ebp1" %in% pids$pcat)
pids$pcat%>%table

ms_tall_trans <- ms_tall%>%
	ungroup%>%
	filter(sigtype=='iBAQ')%>%
	inner_join(pids%>%filter(!pcat=='other'),by='Protein_IDs')


#besides those, do we have unique quantities per gene_name?
#i.e. every colon seperated gene name appears once?
stopifnot(
	ms_tall_trans$gene_name%>%
	na.omit%>%unique%>%
	str_split_fast(';')%>%unlist%>%
	table%>%`==`(1)
)

#+ collect nums, eval =TRUE, cache = FALSE
n_transids = n_distinct(ridssplit)
n_ribogroups = n_distinct(ms_tall_trans$Protein_IDs)
#then we can just keep the rp name and shorten them to the first element otherwise
ms_tall_trans%<>%mutate(gene_name_simp=str_replace(gene_name,regex('(.*?)([Rr]p[sl][^;]+)(.*)'),'\\2'))
ms_tall_trans$gene_name_simp%<>%str_replace(';.*$','')
#unique correspondance between simplified gene names and protein IDS
ms_tall_trans%<>%mutate(ambig = Protein_IDs %in% all_ambig_pgroups)
stopifnot(ms_tall_trans%>%filter(!ambig)%>%distinct(Protein_IDs,gene_name_simp)%>%map_lgl(.%>%anyDuplicated%>%`==`(1))%>%not%>%all)


message('normalizing iBAQ')

rowcol = "Protein_IDs"
sigcol = "signal"
colcols = c('time','fraction','replicate')
#get matrix of all data for size factors
datamat<-ms_tall_trans%>%
	filter(!ambig)%>%
	select(!!rowcol,!!sigcol,!!!colcols)%>%
	unite_('column',colcols)%>%
	spread_('column',sigcol)%>%
	{set_rownames(as.data.frame(.[,-1]),.[[1]])}
#calculate them
sizefactors = DESeq2::estimateSizeFactorsForMatrix(datamat)
sizefactors=sizefactors%>%stack%>%set_colnames(c('sizefactor','tmp'))%>%
	separate(tmp,into=c('time','fraction','replicate'))

sizefactors$replicate %<>% as.character
ms_tall_trans$replicate %<>% as.character
ms_tall$replicate %<>% as.character

ibnormms_tall <- ms_tall%>%filter(sigtype=='iBAQ')%>%left_join(sizefactors,by=c('time','fraction','replicate'))%>%mutate(sigtype='norm_iBAQ',signal=signal*sizefactor)

#create normalized signal column in ms data
ms_tall_trans%<>%left_join(sizefactors)%>%
	mutate(normsignal = signal / sizefactor)%>%
	select(-sizefactor)

message('getting stochiometry matrices')
#now matrix for our 
sigmat = ms_tall_trans%>%
	filter(!ambig)%>%
	select(Protein_IDs,normsignal,time,fraction,replicate)%>%
	group_by(Protein_IDs,time,fraction)%>%
	summarise(l2normsignal = log2(median(normsignal,na.rm=TRUE)))%>%
	unite_('column',colcols[1:2])%>%
	spread_('column','l2normsignal')

#names of our datasets and protein ids
matpids <- sigmat[[1]]
dsetnames = colnames(sigmat)[2:ncol(sigmat)]

#get stoch matrix for each dataset
stochmats <- lapply(2:ncol(sigmat),function(j) {
		outer(sigmat[[j]],sigmat[[j]], FUN = '-' )%>%
			set_rownames(matpids)%>%
			set_colnames(matpids)
	})%>%setNames(dsetnames)

#plot the stochiometry heatmap
catcolors = data_frame(color=c('#000000',"#FFF200","#00AEEF","#BF1E2E","#E6E6FA"),pcat = c("translation-associated","Rps","Rpl","Ebp1","P68040"))
#colors for fold changes
colors = c(seq(-15,-log2(1.25),length=100),seq(-log2(1.25),log2(1.25),length=100),seq(log2(1.25),15,length=100))
my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)
#plot our heatmaps
hmouts <- list()
for(dset in dsetnames){
	# file.path(root,paste0('plots/ribostochheatmaps/stochiometry_heatmaps.',dset,'.pdf'))%>%dirname%>%dir.create
	hmapfile <- file.path(root,paste0('plots/ribostochheatmaps/stochiometry_heatmaps_newcol.',dset,'.pdf'))
	hmapfile%>%dirname%>%dir.create(showWarnings=FALSE)
	pdf(normalizePath(hmapfile) %T>% message,w = 12, h = 12)
	stochmat = stochmats[[dset]]
	stochmat = stochmat[apply(stochmat,1,.%>%is.na%>%all%>%not),apply(stochmat,2,.%>%is.na%>%all%>%not)]

	#look up the color for our protein catagories
	protcolors <-
		stochmat%>%
		rownames%>%
		vlookup(ms_tall_trans,'Protein_IDs','pcat')%>%
		vlookup(catcolors,'pcat','color')

	par(lend = 1)
	hmouts<-append(hmouts,list(heatmap.2(stochmat, 
		col=colorpanel(75,'purple','black','orange'),
		trace="none",
		keysize=1,
		margins=c(8,6),
		# scale="row",
		# dendrogram="none",
		# Colv = FALSE,
		# Rowv = FALSE,
		# cexRow=0.5 + 1/log10(dim(mymat)[1]),
		# cexCol=1.25,
		main=paste0(dset, " Stochiometry Matrix \n Rps and Translation associated"),
		RowSideColors=protcolors,
		ColSideColors=protcolors,
		key.xlab="Log2(norm_iBAQ ratio)"
	)))
	legend(x=0,y=0.85, legend=catcolors$pcat,fill=catcolors$color,cex=0.7)
	dev.off()
}
stop()

lapply(hmouts,function(hmout){
cutree(as.hclust(hmouts[[1]]$colDendrogram),3)%>%enframe('Protein_IDs','cluster')%>%left_join(ms_tall_trans%>%distinct(Protein_IDs,pcat))%>%
	group_by(cluster)%>%filter(any(pcat=='Ebp1'))%>%
	filter(!pcat %in% c('Ebp1','Rpl','Rps'))
})%>%bind_rows

#figure out which of the translation associated proteins is clustering with the 
#RBps



rids$Protein_IDs%>%str_subset('Q3V1Z5|E9Q070|Q9D8M4|Q9D823')
ms_tall%>%filter(Protein_IDs%>%str_detect('Q3V1Z5|E9Q070|Q9D8M4|Q9D823'))%>%distinct(Protein_IDs,.keep_all=T)
ms_tall$gene_name

pdf(file.path(root,'plots/ribostochdist_new.pdf')%T>%message)
	stochmats[[1]]%>%hist(breaks=42,xlim=c(-5,5))
dev.off()

#now also produce scatterplots for this data.
# onegeneprotgroups = protiddt%>%
#let's replicate Matt's scatter plots
ibnormms_tall <- ms_tall%>%filter(sigtype=='iBAQ')%>%left_join(sizefactors)%>%mutate(sigtype='norm_iBAQ',signal=signal*sizefactor)


scatterdf<-ibnormms_tall%>%
	ungroup%>%
	filter(sigtype=='norm_iBAQ')%>%
	left_join(pids,by='Protein_IDs')%>%
	filter(Protein_IDs%>%str_detect(';'))
lfqe13<- scatterdf%>%   
    filter(time=='E13')%>%
    group_by(pcat,Protein_IDs,fraction)%>%
    summarize(Signal_E13=median(na.omit(signal)))
lfqlater<-scatterdf%>%
    filter(!time=='E13')%>%
    group_by(Protein_IDs,fraction,time)%>%
    summarize(Signal=median(na.omit(signal)))
scatterdf<-left_join(by=c('fraction','Protein_IDs'),lfqe13,lfqlater)
 
ggdf<-scatterdf%>%
  # filter(time=='P0')%>%
  filter(fraction %in% c('poly','cyto','80S'))%>%
  mutate(label=ifelse(pcat=='Ebp1','Ebp1',''))%>%
  filter(Protein_IDs%>%str_detect(';'))

scatterdf$pcat%>%table

inggdf<-ggdf

protcatcols = setNames(c(scattercoldf$col,gethex('grey'),gethex('black')),c('RPL +','RPS +','Ebp1','other','Translation Associated'))

ggdf$Protein%>%table
inggdf$pcat%>%table
    # facet_grid(fraction~time)+
scatterplots<-ggdf%>%
  group_by(fraction,time)%>%
  # group_slice(1)%>%
  # 	sample_frac(0.1)%>%
  do({
    # browser();
    list(
    ggplot(.,aes(label=label,color = pcat,y=Signal_E13,x=Signal))+

    scale_color_manual(values=c('grey',catcolors$color)%>%setNames(c('other',catcolors$pcat)))+
    geom_point(alpha=I(0.2),data=filter(.,pcat=='Other'),size=I(3))+
    geom_point(alpha=I(1),data=filter(.,pcat!='Other'),size=I(3))+
    geom_point(alpha=I(1),data=filter(.,pcat=='Ebp1'),size=I(3))+
    geom_text(show.legend=FALSE,nudge_x = 0.5,nudge_y=-0.5)+
    scale_x_log10(name = paste0('LFQ ',.$time[1]),lim=c(1e5,1e11))+
    scale_y_log10(name = paste0('LFQ ','E13'    ),lim=c(1e5,1e11))+
    geom_abline(slope=1,linetype='dashed')+
    ggtitle(paste0('Fraction: ',.$fraction[1]))+
    guides(text=FALSE)+
    theme_bw()+theme(aspect.ratio=1)+ 
    theme(plot.title = element_text(hjust = 0.5))
  )%>%data_frame(plot=.)})

scatterplotfile<-file.path(root,str_interp('plots/ms_fraction_scatterplots_norm_iBAQ_newcols.pdf'))
pdf(scatterplotfile,useDingbats=FALSE,w=24,h=18)
do.call(gridExtra::grid.arrange,c(scatterplots$plot[1:2],ncol=4))
dev.off()
scatterplotfile%>%message