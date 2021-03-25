suppressMessages(library(limma))
library(biomaRt)
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(assertthat))
suppressMessages(library(here))
library(gplots)
library(RColorBrewer)
library(hashmap)
library(ggrepel)
library(ggpubr)
library(dplyr)
# source(here('src/R/Load_data/make_expression_table.R'))

root<-here()

safe_hashmap<-setRefClass("Safe_Rcpp_Hashmap",
      contains="Rcpp_Hashmap",
      inheritPackage=TRUE
)
setMethod('[[','Safe_Rcpp_Hashmap',function (x, i, j, default,...){
	hashmapname = ''
    .local <- function (x, i, j, ..., exact = TRUE)
    {
        x$`[[`(i)
    }
    out <- .local(x, i, j, ...)
    if(missing(default)){
    	if(any(is.na(out))){
	    	keymissingtxt = as.character(i[is.na(out)])%>%head%>%{ifelse(nchar(.)>15,paste0(substr(.,0,13),'...'),.)}%>%paste0(collapse='...')
    		stop(paste0('Keys missing from safe hashmap: ',hashmapname,':',keymissingtxt))    		
    	}
    } else if(length(default)==1){
    	out[is.na(out)] <- default
    }else{
	    out[is.na(out)] <- default[is.na(out)]
    }
    return(out)
})


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

manual_tl_df<-fread(here('ext_data/transl_prots_2check.txt'))
manual_tl_df%<>%filter(!gene_name=='Pa2g4')
# manual_tl_df[[2]]%in%mstall$Protein_IDs
translsetlist <- list(
	'taset_manual'=manual_tl_df[[2]]%>%str_split_fast%>%unlist,
	'taset_go'=translationprotids
)
translsetname='taset_manual'

{

for(mssigcol in c('LFQ','iBAQ')){
for(translsetname in names(translsetlist)){
	translset <- translsetlist[[translsetname]]

#load the iBAQ ms data
Sys.glob(here(str_interp('pipeline/ms_tables/*')))
mstallfiles <- Sys.glob(here(str_interp('pipeline/ms_tables/ms_${mssigcol}*_ms_*')))%>%
	str_subset(neg=T,'__')%>%
	str_subset(neg=T,'cyto')
stopifnot(mstallfiles%>%str_subset('total')%>%fread%>%nrow%>%`>`(0))
mstall <- mstallfiles%>%map(fread)%>%map(function(.){.$signal <- log2(.$signal) ;. })%>%bind_rows
mstall$signal%<>%{exp(.*log(2))}
mstall$fraction%>%unique

ebp1pid = mstall$Protein_IDs[match('Pa2g4',mstall$gene_name)]

#'get info on the ribosomal subu,nts from mats table
rids <- read_tsv(here('ext_data/riboprotids.tsv'))
lridssplit<-rids%>%filter(!`Mito-RP`,`RPL_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist
sridssplit<-rids%>%filter(!`Mito-RP`,`RPS_+`)%>%.$`Protein_IDs`%>%str_split_fast%>%unlist
ridssplit<-c(lridssplit,sridssplit)
rids%>%filter(!`Mito-RP`)%>%nrow
#get info on proteins so we can catagorize them



#' Thi sis a title with inline R code `r foo`

#' First we load the list of protein IDs, handpicked by Matt, using only the small
#' or large subunits - no mitochondrial riboproteins  
#define ambigous protein groups as those which have elements that appear in more than one protein group
allpgroups <- mstall$Protein_IDs%>%unique
multids<-allpgroups%>%unique%>%str_split_fast(';')%>%unlist%>%table%>%keep(~ . > 1)%>%names
all_ambig_pgroups<-allpgroups%>%sep_element_in(multids)
library(data.table)

#' We pick out the iBAQ data, since LFQ isn't recommended for these kinds of enriched datasets
#' We also include Pa2g4.  
#+ select ribo data  , eval = TRUE, cache=TRUE
n_groups = n_distinct(mstall$Protein_IDs)
n_protid = mstall$Protein_IDs%>%str_split_fast(';')%>%unlist%>%n_distinct


message('Labeling protein IDs with annotations')
pids = mstall%>%ungroup%>%distinct(Protein_IDs)

pids%<>%mutate(pcat = case_when(
 (Protein_IDs==ebp1pid) ~ "Ebp1",
 # (Protein_IDs=='P68040') ~ "P68040",
 # (Protein_IDs=='Q99LG4') ~ 'Ttc5' ,
  sep_element_in(Protein_IDs,sridssplit) ~ "Rps",
  sep_element_in(Protein_IDs,lridssplit) ~ "Rpl",
  sep_element_in(Protein_IDs,translset) ~ "translation-associated",
  TRUE ~ "other"
))
stopifnot("Ebp1" %in% pids$pcat)
pids$pcat%>%table

mstall$fraction%>%unique

stopifnot('total'%>%is_in(mstall$fraction%>%unique))

mstall_trans <- mstall%>%
	ungroup%>%
	filter(sigtype==mssigcol)%>%
	inner_join(pids%>%filter(!pcat=='other'),by='Protein_IDs')


#besides those, do we have unique quantities per gene_name?
#i.e. every colon seperated gene name appears once?
stopifnot(
	mstall_trans$gene_name%>%
	na.omit%>%unique%>%
	str_split_fast(';')%>%unlist%>%
	table%>%`==`(1)
)

#+ collect nums, eval =TRUE, cache = FALSE
n_transids = n_distinct(ridssplit)
n_ribogroups = n_distinct(mstall_trans$Protein_IDs)
#then we can just keep the rp name and shorten them to the first element otherwise
mstall_trans%<>%mutate(gene_name_simp=str_replace(gene_name,regex('(.*?)([Rr]p[sl][^;]+)(.*)'),'\\2'))
mstall_trans$gene_name_simp%<>%str_replace(';.*$','')
#unique correspondance between simplified gene names and protein IDS
mstall_trans%<>%mutate(ambig = Protein_IDs %in% all_ambig_pgroups)
stopifnot(mstall_trans%>%filter(!ambig)%>%distinct(Protein_IDs,gene_name_simp)%>%map_lgl(.%>%anyDuplicated%>%`==`(1))%>%not%>%all)


message('normalizing iBAQ')

rowcol = "Protein_IDs"
sigcol = "signal"
colcols = c('time','fraction','replicate')
#get matrix of all data for size factors
datamat<-mstall_trans%>%
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
mstall_trans$replicate %<>% as.character
mstall$replicate %<>% as.character

ibnormmstall <- mstall%>%filter(sigtype==sigcol)%>%left_join(sizefactors,by=c('time','fraction','replicate'))%>%mutate(sigtype='norm_iBAQ',signal=signal*sizefactor)

#create normalized signal column in ms data
mstall_trans%<>%left_join(sizefactors)%>%
	mutate(normsignal = signal / sizefactor)%>%
	select(-sizefactor)

message('getting stochiometry matrices')

stopifnot('total'%>%is_in(mstall_trans$fraction%>%unique))



stopifnot('Q99LG4'%in%mstall$Protein_IDs)

#now matrix for our 
sigmat = mstall_trans%>%
	filter(!ambig)%>%
	select(Protein_IDs,normsignal,time,fraction,replicate)%>%
	group_by(Protein_IDs,time,fraction)%>%
	summarise(l2normsignal = log2(median(normsignal,na.rm=TRUE)))%>%
	unite_('column',colcols[1:2])%>%
	spread_('column','l2normsignal')

#names of our datasets and protein ids
matpids <- sigmat[[1]]
dsetnames = colnames(sigmat)[2:ncol(sigmat)]
stopifnot(length(dsetnames)==15)

#get stoch matrix for each dataset
stochmats <- lapply(2:ncol(sigmat),function(j) {
		outer(sigmat[[j]],sigmat[[j]], FUN = '-' )%>%
			set_rownames(matpids)%>%
			set_colnames(matpids)
	})%>%setNames(dsetnames)




#plot the stochiometry heatmap
# catcolors = data_frame(color=c('#000000',"#FFF200","#00AEEF","#BF1E2E","green"),pcat = c("translation-associated","Rps","Rpl","Ebp1","P68040"))
catcolors = data_frame(color=c('#000000',"#FFF200","#00AEEF","#BF1E2E","green"),pcat = c("translation-associated","Rps","Rpl","Ebp1","Ttc5"))
#colors for fold changes
colors = c(seq(-15,-log2(1.25),length=100),seq(-log2(1.25),log2(1.25),length=100),seq(log2(1.25),15,length=100))
my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)

# #plot our heatmaps
# hmouts <- list()
# for(dset in dsetnames){
# 	# file.path(root,paste0('plots/ribostochheatmaps/stochiometry_heatmaps.',dset,'.pdf'))%>%dirname%>%dir.create
# 	hmapfile <- file.path(root,paste0('plots/ribostochheatmaps/stochiometry_heatmaps_newcol.',dset,'.pdf'))
# 	hmapfile%>%dirname%>%dir.create(showWarnings=FALSE)
# 	pdf(normalizePath(hmapfile) %T>% message,w = 12, h = 12)
# 	stochmat = stochmats[[dset]]
# 	stochmat = stochmat[apply(stochmat,1,.%>%is.na%>%all%>%not),apply(stochmat,2,.%>%is.na%>%all%>%not)]

# 	#look up the color for our protein catagories
# 	protcolors <-
# 		stochmat%>%
# 		rownames%>%
# 		vlookup(mstall_trans,'Protein_IDs','pcat')%>%
# 		vlookup(catcolors,'pcat','color')

# 	par(lend = 1)
# 	hmouts<-append(hmouts,list(heatmap.2(stochmat, 
# 		col=colorpanel(75,'purple','black','orange'),
# 		trace="none",
# 		keysize=1,
# 		margins=c(8,6),
# 		# scale="row",
# 		# dendrogram="none",
# 		# Colv = FALSE,
# 		# Rowv = FALSE,
# 		# cexRow=0.5 + 1/log10(dim(mymat)[1]),
# 		# cexCol=1.25,
# 		main=paste0(dset, " Stochiometry Matrix \n Rps and Translation associated"),
# 		RowSideColors=protcolors,
# 		ColSideColors=protcolors,
# 		key.xlab="Log2(norm_iBAQ ratio)"
# 	)))
# 	legend(x=0,y=0.85, legend=catcolors$pcat,fill=catcolors$color,cex=0.7)
# 	dev.off()
# }


# BiocManager::install('hashmap')
# library(hashmap)

# lapply(hmouts,function(hmout){
# cutree(as.hclust(hmouts[[1]]$colDendrogram),3)%>%enframe('Protein_IDs','cluster')%>%left_join(mstall_trans%>%distinct(Protein_IDs,pcat))%>%
# 	group_by(cluster)%>%filter(any(pcat=='Ebp1'))%>%
# 	filter(!pcat %in% c('Ebp1','Rpl','Rps'))
# })%>%bind_rows

#figure out which of the translation associated proteins is clustering with the 
#RBps


# library(hashmap)
# library(concatenate)
msid_cat <- distinct(mstall_trans,Protein_IDs,pcat)%>%{safe_hashmap(keys=.[[1]],values=.[[2]])}
gname <- mstall_trans%>%distinct(Protein_IDs,gene_name)%>%filter(!is.na(gene_name))%>%{safe_hashmap(.[[1]],.[[2]])}
#

dsetnames[[1]]->dset
ribostochlists<-lapply(dsetnames,function(dset){
	#now get the median of all rps
	ribomsids<-rownames(stochmats[[dset]])%>%{.[msid_cat[[.]]%>%is_in(c('Rpl','Rps'))]}
	dsetribomed<-sigmat%>%filter(Protein_IDs%in%ribomsids)%>%ungroup%>%summarise_at(vars(-1),list(~ median(.,na.rm=T)))%>%t%>%
		safe_hashmap(keys=rownames(.),values=.)
	stopifnot(dset%in%colnames(sigmat))
	relvals <- sigmat[[dset]] - dsetribomed[[dset]]
	#now plot the stochiometery relative to the median for all ribosomal proteins
	# plotfile<- file.path(root,paste0('plots/ribostochbars',dset,'.pdf'))
	# pdf(plotfile)
	pos <- position_jitter(width = 0.3, seed = 2)
	# plotfile<-'tmp.pdf'
	# pdf('tmp.pdf')	
	ggdf <- data.frame(relativeLFQ =relvals,pcat = msid_cat[[sigmat$Protein_IDs]])%>%
		mutate(gname = sigmat$Protein_IDs%>%gname[[.,default=.]] ) %>%
		mutate(label = ifelse(
			(abs(relativeLFQ)>2.5) & (pcat%in%c('Rps','Rpl'))|
				((translsetname=='taset_manual') & (pcat%in%c('translation-associated'))),
			gname,
			'')
		) %>%
		mutate(pcat = factor(pcat,c('translation-associated','Ebp1','Rps','Rpl')))%>%
		mutate(pcat = recode(pcat,'translation-associated'='translation\nassociated'))%>%
		filter(!is.na(pcat))
	# stopifnot('Q99LG4'%in%sigmat$Protein_IDs)
	if(translsetname=='taset_manual' ) stopifnot(msid_cat[['Q99LG4']]=='translation-associated')
	plot<-ggdf%>%
		ggplot(.,aes(y=relativeLFQ,color=pcat,fill=pcat,x=pcat,label=label))+
		geom_jitter(position=pos,alpha=I(0.5))+
		geom_text_repel(position=pos)+
		geom_boxplot(outlier.shape=NA,alpha=I(0.5))+
		# stat_identity(geom='bar')+
		scale_color_discrete(name='category',catcolors%>%{setNames(.$color,.$pcat)})+
		scale_fill_discrete(name='category',catcolors%>%{setNames(.$color,.$pcat)})+
		scale_x_discrete('Category')+
		scale_y_continuous(str_interp('log2(${mssigcol}) - log2(Median Rp(l/s) ${mssigcol})'),limits=c(-15,5))+
		ggtitle('Protein Stochiometry - ',dset)+
		theme(axis.text.x=element_text(angle=45))+
		theme_bw()
		list(plot,ggdf)
})
ribostochplots <- map(ribostochlists,1)
ribostochdfs <- map(ribostochlists,2)

ribostochdfs%<>%setNames(dsetnames)%>%bind_rows(.id='dataset')

tablefile<- here(paste0('tables/','ribostoch_',mssigcol,'_',translsetname,'.tsv'))

left_join(
	ribostochdfs%>%group_by(dataset)%>%summarise(ebp1_vs_Rpl = 2^(relativeLFQ[pcat=='Ebp1'] - median(na.rm=T,relativeLFQ[pcat=='Rpl']))),
	ribostochdfs%>%group_by(dataset)%>%summarise(ebp1_vs_Rps = 2^(relativeLFQ[pcat=='Ebp1'] - median(na.rm=T,relativeLFQ[pcat=='Rps'])))
)%>%write_tsv(tablefile)

message(normalizePath(tablefile))
#now plot
#plot without labels

ribostochplotsnorepel<-ribostochplots%>%map(function(x){x$layers[[2]]=NULL;x})
plotfile<- here(paste0('plots/','ribostochjitters_',mssigcol,'_',translsetname,'.pdf'))
#
pdf(plotfile,w=12,h=20)
print(ggarrange(plotlist=ribostochplotsnorepel,nrow=5,ncol=3))
dev.off()
message(normalizePath(plotfile))
#
library(ggrepel)
plotfile<- here(paste0('plots/','ribostochjitters_label_',mssigcol,'_',translsetname,'_','.pdf'))
pdf(plotfile,w=24,h=40)
print(ggarrange(plotlist=ribostochplots,nrow=5,ncol=3))
dev.off()
normalizePath(plotfile)
message(normalizePath(plotfile))

}
}


}

stop()

save.image(here('data/ribo_stoch_heatmaps.R'))
stop()


# rids$Protein_IDs%>%str_subset('Q3V1Z5|E9Q070|Q9D8M4|Q9D823')
# mstall%>%filter(Protein_IDs%>%str_detect('Q3V1Z5|E9Q070|Q9D8M4|Q9D823'))%>%distinct(Protein_IDs,.keep_all=T)
# mstall$gene_name

# pdf(file.path(root,'plots/ribostochdist_new.pdf')%T>%message)
# 	stochmats[[1]]%>%hist(breaks=42,xlim=c(-5,5))
# dev.off()

# #now also produce scatterplots for this data.
# # onegeneprotgroups = protiddt%>%
# #let's replicate Matt's scatter plots
# ibnormmstall <- mstall%>%filter(sigtype==sigcol)%>%left_join(sizefactors)%>%mutate(sigtype='norm_iBAQ',signal=signal*sizefactor)


# scatterdf<-ibnormmstall%>%
# 	ungroup%>%
# 	filter(sigtype=='norm_iBAQ')%>%
# 	left_join(pids,by='Protein_IDs')%>%
# 	filter(Protein_IDs%>%str_detect(';'))
# lfqe13<- scatterdf%>%   
#     filter(time=='E13')%>%
#     group_by(pcat,Protein_IDs,fraction)%>%
#     summarize(Signal_E13=median(na.omit(signal)))
# lfqlater<-scatterdf%>%
#     filter(!time=='E13')%>%
#     group_by(Protein_IDs,fraction,time)%>%
#     summarize(Signal=median(na.omit(signal)))
# scatterdf<-left_join(by=c('fraction','Protein_IDs'),lfqe13,lfqlater)
 
# ggdf<-scatterdf%>%
#   # filter(time=='P0')%>%
#   filter(fraction %in% c('poly','total','80S'))%>%
#   mutate(label=ifelse(pcat=='Ebp1','Ebp1',''))%>%
#   filter(Protein_IDs%>%str_detect(';'))

# scatterdf$pcat%>%table

# inggdf<-ggdf

# protcatcols = setNames(c(scattercoldf$col,gethex('grey'),gethex('black')),c('RPL +','RPS +','Ebp1','other','Translation Associated'))

# ggdf$Protein%>%table
# inggdf$pcat%>%table
#     # facet_grid(fraction~time)+
# scatterplots<-ggdf%>%
#   group_by(fraction,time)%>%
#   # group_slice(1)%>%
#   # 	sample_frac(0.1)%>%
#   do({
#     # browser();
#     list(
#     ggplot(.,aes(label=label,color = pcat,y=Signal_E13,x=Signal))+

#     scale_color_manual(values=c('grey',catcolors$color)%>%setNames(c('other',catcolors$pcat)))+
#     geom_point(alpha=I(0.2),data=filter(.,pcat=='Other'),size=I(3))+
#     geom_point(alpha=I(1),data=filter(.,pcat!='Other'),size=I(3))+
#     geom_point(alpha=I(1),data=filter(.,pcat=='Ebp1'),size=I(3))+
#     geom_text(show.legend=FALSE,nudge_x = 0.5,nudge_y=-0.5)+
#     scale_x_log10(name = paste0(mssigcol,' ',.$time[1]),lim=c(1e5,1e11))+
#     scale_y_log10(name = paste0(mssigcol,' ','E13'    ),lim=c(1e5,1e11))+
#     geom_abline(slope=1,linetype='dashed')+
#     ggtitle(paste0('Fraction: ',.$fraction[1]))+
#     guides(text=FALSE)+
#     theme_bw()+theme(aspect.ratio=1)+ 
#     theme(plot.title = element_text(hjust = 0.5))
#   )%>%data_frame(plot=.)})

# scatterplotfile<-file.path(root,str_interp('plots/ms_fraction_scatterplots_norm_iBAQ_newcols.pdf'))
# pdf(scatterplotfile,useDingbats=FALSE,w=24,h=18)
# do.call(gridExtra::grid.arrange,c(scatterplots$plot[1:2],ncol=4))
# dev.off()
# scatterplotfile%>%message

# gname[['']]