# #+ set opts, echo=FALSE, message = FALSE, warning = FALSE, eval = FALSE 
# knitr::opts_chunk$set(eval=FALSE, echo=FALSE,warning=FALSE,eval=FALSE)

# #+ spin,eval=FALSE, echo=FALSE 
# scriptname='/home/dharnett/cortexomics/exploration/ribo_stochiometry.R'
# spindir = dirname(scriptname)
# oldwd = setwd(spindir)
# setwd(spindir)
# on.exit(setwd(oldwd))
# message(paste0("Spinning document in:", spindir))
# d = knitr::spin(scriptname,report=FALSE)
# rmarkdown::render(d, 
# 	output_dir = spindir,
# 	knit_root_dir = spindir
# )

# stop('Spinning Successful!')

#' # Assessing ribosomal stiochiometry
#' We' like to assess the stochiometry of the variouus ribosomal fractions in our data.  
#' To do this we need to  
#' #. Get our data  
#' #. Select the iBAQ  
#' #. Normalize the iBAQ (we do this in the style of DESeq)  
#' #. Do some basic plots to make sure our normalization 
#' #. Get the median ratio of everything to everything else  
#' 	+ foo
#' 	+ bar
# dicttable=ms_tall_trans
# vals='pcat'
# key='Protein_IDs'
# query=matpids
#get stoch matrix for each dataset
library(gplots)
library(RColorBrewer)

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

#' 


library(biomaRt)
mart <- useMart(biomart = "ENSEMBL_MART_MOUSE")
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
listAttributes(mart)%>%filter(description%>%str_detect('Uni'))%>%filter(page!='homologs')
ribogoterm <- "GO:0005840"
riboprotids <- biomaRt::getBM(attributes = c("uniprot_gn"), 
                 filters=('go'),values=ribogoterm,
                 mart = mart)
transreggoterm <- "GO:0006417"
transregprotids <- biomaRt::getBM(attributes = c("uniprot_gn"), 
                 filters=('go'),values=transreggoterm,
                 mart = mart)
translationgoterm <- "GO:0006412"
translationprotids <- biomaRt::getBM(attributes = c("uniprot_gn"), 
                 filters=('go'),values=translationgoterm,
                 mart = mart)[[1]]

ebp1pid = ms_tall$Protein_IDs[match('Pa2g4',ms_tall$gene_name)]



#' Thi sis a title with inline R code `r foo`

#' First we load the list of protein IDs, handpicked by Matt, using only the small
#' or large subunits - no mitochondrial riboproteins  


#define ambigous protein groups as those which have elements that appear in more than one protein group
allpgroups <- ms_tall$Protein_IDs%>%unique
multids<-allpgroups%>%unique%>%str_split_fast(';')%>%unlist%>%table%>%keep(~ . > 1)%>%names
all_ambig_pgroups<-allpgroups%>%sep_element_in(multids)

# ms_tall_trans<-ms_tall%>%mutate(ambig = sep_element_in(Protein_IDs,multids))


#' We pick out the iBAQ data, since LFQ isn't recommended for these kinds of enriched datasets
#' We also include Pa2g4.  
#+ select ribo data  , eval = TRUE, cache=TRUE
n_groups = n_distinct(ms_tall$Protein_IDs)
n_protid = ms_tall$Protein_IDs%>%str_split_fast(';')%>%unlist%>%n_distinct

pids = ms_tall%>%ungroup%>%distinct(Protein_IDs)
pids%<>%mutate(pcat = case_when(
  identical(Protein_IDs,EXTRAPROTEINS) ~ ebp1pid,
  sep_element_in(Protein_IDs,ridssplit) ~ "Ribosomal",
  sep_element_in(Protein_IDs,translationprotids) ~ "Translation Associated",
  TRUE ~ "other"
))


ms_tall_trans <- ms_tall%>%
	ungroup%>%
	filter(sigtype=='iBAQ')%>%
	inner_join(pids%>%filter(!pcat=='other'),by='Protein_IDs')

#+ collect nums, eval =TRUE, cache = FALSE
n_transids = n_distinct(ridssplit)
n_ribogroups = n_distinct(ms_tall_trans$Protein_IDs)


#' ## Riboprotein - Numbers
#' In our total data, there are `r n_groups` protein groups with `r n_protid` protein IDs
#' of which `r n_distinct(multids)` are shared between more than one group,

#' Of the `r n_riboids` riboproteins specified in matt's table, we find
#' `r n_ribogroups` Included in protein groups in our data. Of these `r n_noovribogroups`
#' contain no protein ids shared with other groups and are therefore ambigous
#' riboprotein groups for which we can unambigously assign Mass spec signal. Of these, 


#` We select only row names from the ms_data which have gene names annotated
# nogenename_ribo_pids <- ms_tall_trans%>%filter(gene_name%>%is.na)%>%.$Protein_IDs%>%unique
# data_frame(nogenename_ribo_pids)%>%
# 	write_tsv('~/projects/cortexomics/nogenename_ribo_pids.tsv')
# ms_tall_trans %<>% filter(!is.na(gene_name))

#besides those, do we have unique quantities per gene_name?
#i.e. every colon seperated gene name appears once?
stopifnot(
	ms_tall_trans$gene_name%>%
	na.omit%>%unique%>%
	str_split_fast(';')%>%unlist%>%
	table%>%`==`(1)
)

#then we can just keep the rp name and shorten them to the first element otherwise
ms_tall_trans%<>%mutate(gene_name_simp=str_replace(gene_name,regex('(.*?)([Rr]p[sl][^;]+)(.*)'),'\\2'))
ms_tall_trans$gene_name_simp%<>%str_replace(';.*$','')


#unique correspondance between simplified gene names and protein IDS
ms_tall_trans%<>%mutate(ambig = Protein_IDs %in% all_ambig_pgroups)
stopifnot(ms_tall_trans%>%filter(!ambig)%>%distinct(Protein_IDs,gene_name_simp)%>%map_lgl(.%>%anyDuplicated%>%`==`(1))%>%not%>%all)




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
#create normalized signal column in ms data
ms_tall_trans%<>%left_join(sizefactors)%>%
	mutate(normsignal = signal / sizefactor)%>%
	select(-sizefactor)

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

colnames(mymat) <- c(rep("treatment_1", 3), rep("treatment_2", 3), rep("treatment_3", 3), rep("treatment_4", 3))
rownames(mymat) <- paste("gene", 1:dim(mymat)[1], sep="_")
mymat
mydf <- data.frame(gene=paste("gene", 1:dim(mymat)[1], sep="_"), category=c(rep("CATEGORY_1", 10), rep("CATEGORY_2", 10), rep("CATEGORY_3", 10), rep("CATEGORY_4", 10), rep("CATEGORY_5", 10)))
mydf

#get stoch matrix for each dataset
stochmats <- lapply(2:ncol(sigmat),function(j) {
		outer(sigmat[[j]],sigmat[[j]], FUN = '-' )%>%
			set_rownames(matpids)%>%
			set_colnames(matpids)
	})%>%setNames(dsetnames)

#plot the stochiometry heatmap
catcolors = data_frame(color=c("Red","Black"),pcat = c("Ribosomal","Translation Associated"))

hmapfile <- file.path(root,'plots/stochiometry_heatmaps.pdf')
pdf(hmapfile,w = 12, h = 12)
for(dset in dsetnames){
	stochmat = stochmats[[dset]]
	stochmat = stochmat[apply(stochmat,1,.%>%is.na%>%all%>%not),apply(stochmat,2,.%>%is.na%>%all%>%not)]

	protcolors <-
		stochmat%>%
		rownames%>%
		vlookup(ms_tall_trans,'Protein_IDs','pcat')%>%
		vlookup(catcolors,'pcat','color')

	par(lend = 1)
	heatmap.2(stochmat, 
		# col=greenred(75),
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
	)
	legend(x=0,y=0.85, legend=catcolors$pcat,fill=catcolors$color,cex=0.7)
}
dev.off()

normalizePath(hmapfile) %>% message
match(rownames(stochmat))

# dev.off()

# #normalize those ibaqs babay
# ibaqnormfacts<-ms_tall_trans%>%group_by(time,Protein_IDs)%>%summarize(m=median(na.omit(signal)))%>%
# 	group_by(time)%>%summarize(medexpr=median(na.omit(m)))%>%
# 	mutate(normfact = medexpr - median(medexpr))

# ms_tall_trans%<>%left_join(ibaqnormfacts)%>%mutate(signal = signal+ normfact)



ms_tall_trans$fraction%>%unique

#fit the most basic model where there are no timepoint specific factors
exprmat = ms_tall_trans%>%
	split(.,.$fraction)%>%
	map(.%>%
	# .[[2]]%>%
	select(gene_name_simp,dataset,signal)%>%
	spread(dataset,signal)%>%
	{ set_rownames(as.matrix(.[-1]),.[[1]]) }
	)

# #let's look at the median counts and the number of NAs
# exprmat%>%map(.%>%apply(1,.%>%is.na%>%sum))
# exprmat%>%map(.%>%apply(1,.%>%median))

#select 1? Or just 
monomat <- exprmat[['cyto']]
monomat <- monomat[monomat%>%apply(1,.%>%is.na%>%sum%>%`<`(1)),]
monomat <- monomat[monomat%>%apply(1,median)%>%`>`(10),]

refratios <- 
refexpr = apply(datamat,1,.%>%na.omit%>%log2%>%mean%>%`^`())

refratios = sweep(datamat,MARGIN=1,FUN='/',refexpr)

#take the log
lmonomat = log2(monomat)


means <- 	monomat%>%apply(1,mean,na.rm=TRUE)
sds <- 	monomat%>%apply(1,sd,na.rm=TRUE)


refexpr = apply(datamat,1,function(x){browser()})
refratios = sweep(datamat,MARGIN=1,FUN='/',refexpr)



qplot(means,sds,label=rownames(monomat))+geom_text()


# lmeans <- 	lmonomat%>%apply(1,mean,na.rm=TRUE)
# lsds <- 	lmonomat%>%apply(1,sd,na.rm=TRUE)

# qplot(lmeans,lsds,label=rownames(monomat))+geom_text()

#okay ten we've reliminated all that messy mean-variance.
#Now let's make our model matrix
designmat <- ms_tall_trans%>%
	select(dataset,time,fraction)%>%
	filter(dataset%in%colnames(monomat))%>%
	distinct(dataset,.keep_all=TRUE)%>%
	select(-dataset)%>%
	select_if(~ n_distinct(.)>1 )

design = model.matrix( ~ time , designmat )


fit = limma::lmFit(lmonomat,design=design)
bayesfit = limma::eBayes(fit,trend=TRUE, robust=TRUE)

coefs<-bayesfit$p.value%>%colnames
allcoefs<-coefs%>%setNames(.,.)%>%map(~topTable(bayesfit,number=nrow(lmonomat),coef=.,confint=0.95)%>%{cbind(gene=rownames(.),.)})%>%
	bind_rows(.id='coefficient')%>%
	as_data_frame

#this will need to change if/when we bring in additional model coefficients
predictdf <- allcoefs%>%
	# filter(gene=='Rps3')%>%
	group_by(gene)%>%
	mutate(intcept     = coefficient=='(Intercept)')%>%
	mutate(coefficient = ifelse(intcept,'timeE13',coefficient))%>%
	mutate(signal = ifelse(intcept,logFC,logFC+logFC[which(intcept)]))%>%
	mutate(upper = ifelse(intcept,CI.R,CI.R+logFC[which(intcept)]))%>%
	mutate(lower = ifelse(intcept,CI.L,CI.L+logFC[which(intcept)]))%>%
	transmute(time=coefficient,gene_name_simp=gene,signal,upper,lower)%>%
	mutate(time = as.factor(str_replace(time,'time','')))

datadf<-ms_tall_trans %>% 
	filter(fraction=='cyto')%>%
	select(gene_name_simp,time,signal)%>%
	mutate(signal=log2(signal))%>%
	mutate(time=as.factor(time))

samplerbps<-predictdf$gene_name_simp%>%sample(20)
plotdatadf<-datadf%>%ungroup%>%filter(gene_name_simp%in%samplerbps)
predictpredictdf<-predictdf%>%ungroup%>%filter(gene_name_simp%in%samplerbps)
plotdatadf$gene_name_simp

traj.plots<-
	ggplot(data=plotdatadf,aes(x=as.numeric(time),y=signal)) +
	geom_point(size=I(3)) +
	geom_ribbon(data=predictpredictdf%>%plotdatafilt,aes(x=as.numeric(time),y=signal,ymax=upper,ymin=lower),
		alpha=I(0.5)) +
	scale_y_continuous(name='log2 iBAQ') +
	ggtitle('cyto proteomic trajectory - 95% CI')+
	facet_wrap(~gene_name_simp,ncol=4,scale='free')+
	theme_bw()

traj.plots%>%ggsave(file='~/projects/cortexomics/figures/all_ribo_80s_trajplots.pdf')













# ##looking at ambiguities



# P14869;S4R1N1;D3YVM5

# ms_data%>%filter(Protein_IDs=='')



# ms_tall_trans%>%select(gene_name_simp,dataset,signal)%>%spread(dataset,signal)

# ms_tall_trans%>%filter(Protein_IDs=='A2A547;P84099')%>%.$gene_name
# ms_tall_trans%>%filter(Protein_IDs=='Q9CZM2;E9QAZ2;B8JKK2')%>%.$gene_name
# ms_tall_trans%>%filter(gene_name=='Rpl19')%>%.$Protein_IDs
# ms_tall_trans%>%select(Protein_IDs,gene_name)
# ms_tall_trans%>%select(Protein_IDs,gene_name)%>%distinct
# ms_tall_trans%>%select(gene_name)%>%distinct
# ms_tall_trans%>%select(Protein_IDs)%>%distinct

# #do we have multiple protein IDs corresponding to this rp gene?
# allpgroups      <- ms_tall$Protein_IDs%>%unique
# allpids         <- allpgroups%>%str_split_fast(';')%>%unlist%>%unique
# allpgroups_ribo <- allpgroups[sep_element_in(allpgroups,ridssplit)]
# allpids_wribo    <- allpgroups_ribo%>%str_split_fast(';')%>%unlist%>%unique
# allpids_ribo    <- ridssplit
# allpgene_names <- ms_tall$gene_name%>%str_extract('[Rr]p[sl]\\d+')%>%.[!is.na(.)]%>%unique
# multids<-allpgroups%>%unique%>%str_split_fast(';')%>%unlist%>%table%>%keep(~ . > 1)%>%names


# #ribosomal ids which appear in more than one group
# multids_ribo<-allpgroups_ribo%>%unique%>%str_split_fast(';')%>%unlist%>%table%>%keep(~ . > 1)%>%names

# all_ids_tall<-ms_data_all%>%
# 	select(gene_names=gene_name.x,pids=Protein_IDs)%>%
# 	mutate(row=seq_len(n()))%>%
# 	# filter(!is.na(gene_names))%>%
# 	mutate(genename=str_split_fast(gene_names,';'))%>%
# 	mutate(genenamenum=lengths(genename))%>%
# 	unnest%>%
# 	mutate(pid = str_split_fast(pids,';'))%>%
# 	mutate(pidnamenum=lengths(pid))%>%
# 	unnest%>%
# 	mutate(isribogenename=genename%in%allpgene_names)%>%
# 	mutate(isribopid=pid%in%ridssplit)

# #rows and pid groupss correspond perfectly
# all_ids_tall%>%select(pids,row)%>%distinct
# all_ids_tall%>%select(row)%>%distinct
# all_ids_tall%>%select(pids)%>%distinct

# #so it looks like often when ribo pid is in multiple rows, that's because it
# #has a second with with NA
# all_ids_tall%>%filter(pid%in%multids_ribo)%>%group_by(pid)%>%group_slice(10)%>%as.data.frame
# #always?
# all_ids_tall%>%filter(pid%in%multids_ribo)%>%group_by(pid)%>%filter(sum(!is.na(unique(gene_names)))>1)%>%as.data.frame%>%
# arrange(pid)
# all_ids_tall%>%filter(pid%in%multids_ribo)%>%group_by(pid)%>%filter(sum(!is.na(unique(genename)))>1)%>%as.data.frame%>%
# arrange(pid)

#  P14869;S4R1N1;D3YVM5 2176           1    Rplp0          3 D3YVM5

# all_ids_tall%>%filter(pid%in%multids_ribo)%>%group_by(pid)%>%group_slice(10)%>%as.data.frame

# #here is an example of 
# sharedwithriboids<-all_ids_tall%>%group_by(pids)%>%mutate(pgroup_has_Rplp0=any(pid=='D3YVM5'))%>%
# 	group_by(pid)%>%
# 	filter(any(pgroup_has_Rplp0))%>%as.data.frame

# all_ids_tall%>%group_by(pids)%>%mutate(pgroup_has_Rplp0=any(pid=='E9Q1N8'))%>%
# 	group_by(pid)%>%
# 	filter(any(pgroup_has_Rplp0))%>%as.data.frame

# all_ids_tall%>%group_by(pids)%>%mutate(pgroup_has_Rplp0=any(pid==multids_ribo))%>%
# 	group_by(pid)%>%
# 	filter(any(pgroup_has_Rplp0))%>%as.data.frame

# sharedwithriboids%>%arrange(pid)%>%as.data.frame
# sharedwithriboids%>%group_by(pid)%>%filter(n_distinct(pids)>1)%>%arrange(pid)%>%as.data.frame


# #yes okay, so now what about the other proteins in teh groups with our 
# #ribosomal proteins?
# #let's see the pids that appear in groups 
# sharedwithriboids<-all_ids_tall%>%group_by(pids)%>%mutate(pgroup_has_ribo=any(isribopid))%>%
# 	group_by(pid)%>%
# 	filter(any(!isribopid))%>%group_by(pid)%>%filter(any(pgroup_has_ribo))
# sharedwithriboids%>%group_by(pid)%>%filter(n_distinct(pids)>1)%>%arrange(pid)%>%as.data.frame




# all_ids_tall%>%filter(row==454)%>%as.data.frame

# %>%group_by(pid)%>%filter(sum(!is.na(unique(genename)))>1)%>%as.data.frame%>%
# arrange(pid)


# #get the protein groups that have ribopids that are in more than one group
# all_ids_tall%>%group_by(pids)%>%filter(any(pid%in%multids_ribo))%>%as.data.frame


# all_ids_tall%>%group_by(pids)%>%as.data.frame%>%
# 	group_by(pids)%>%summarize(n=n_distinct(gene_names))%>%.$n%>%table
# all_ids_tall%>%group_by(pids)%>%as.data.frame%>%
# 	group_by(pids)%>%summarize(n=n_distinct(row))%>%.$n%>%table


# #in these cases, could we just sum up by gene name (or catenated genenames)
# all_ids_tall%>%group_by(pids)%>%filter(any(pid%in%multids_ribo))%>%as.data.frame%>%
# 	group_by(pids)%>%summarize(n=n_distinct(gene_names))%>%.$n
# all_ids_tall%>%group_by(pids)%>%filter(any(pid%in%multids_ribo))%>%as.data.frame%>%
# 	group_by(pids)%>%summarize(n=n_distinct(row))%>%.$n


# #do our gene names sometimes appear in more than one gene name group?
# #some gene names appear in more than one gene_name group - yes
# all_ids_tall%>%distinct(genename,gene_names)%>%group_by(genename)%>%filter(n()>1)%>%
# 	arrange(genename)
# #but non of these are ribosomal proteins
# all_ids_tall%>%filter(isribogenename)%>%distinct(genename,gene_names)%>%group_by(genename)%>%filter(n()>1)%>%
# 	arrange(genename)
# #protein IDs are unique per row of our data
# ms_data_all$Protein_IDs%>%duplicated%>%any

# #so we can get a unique gene_name row for each of our ribosomal proteins
# #using the gene_name column....

# #what if we want to use the pid column instead?
# all_ids_tall%>%distinct(pid,Protein_IDs)%>%group_by(pid)%>%filter(n()>1)%>%
# 	arrange(pid)
# all_ids_tall%>%filter(isribopid)%>%distinct(pid,Protein_IDs)%>%group_by(pid)%>%filter(n()>1)
# 	arrange(pid)
# #yes - each of our pids always appears on a single row

# #of course each gene name group often corresponds to more than one protein id
# all_ids_tall%>%distinct(gene_names,pids)%>%group_by(gene_names)%>%filter(n()>1)%>%
# 	arrange(gene_names)%>%as.data.frame

# #ribosomal genes in more than one protein group?
# all_ids_tall%>%distinct(gene_names,pids)%>%group_by(gene_names)%>%filter(n()>1)%>%
# 	arrange(gene_names)%>%as.data.frame


# #some ribosomal pids appear in more than one group? - yes
# allpgroups%>%str_split_fast(';')%>%unlist%>%table%>%{.[names(.)%in%allpids_ribo]}%>%{.[.>1]}


# #if our ribosomal protein ids appear in more than one group, does this group
# #always have a single gene ID attached to it?
# #

# okay so when a gene_name is present, it seema ll of those pids ar
# #always present for only that gene_name, or without an assigned gene_name
# stopifnot(all_ids_tall%>%filter(!is.na(gene_names))%>%distinct(pid,gene_names)%>%group_by(pid)%>%tally%>%.$n%>%`==`(1))
# #do 
# stopifnot(all_ids_tall%>%distinct(pid,genename)%>%group_by(pid)%>%tally%>%.$n%>%`==`(1))
# all_ids_tall%>%distinct(pid,genename)%>%group_by(pid)%>%tally%>%.$n%>%`==`(1)


# #some groups with ribosomal protein ids in them have other ids that are also in other gorups


# # ribogenename<-'Rpl19'
# # ribopids = ms_tall_trans%>%filter(gene_name%>%str_detect(ribogenename))%>%.$Protein_IDs%>%unique%>%str_split_fast(';')%>%unlist
# # allpgroups%>%str_subset(ribopids[1])%>%n_distinct	
# # allpgroups%>%str_subset(ribopids[2])%>%n_distinct

# #what about in general? Do we have overlapping protein groups?



# #okay so we have 39 protein ids which ovvur in more than one group
# multids_ribo%>%n_distinct

# #what groups have overlapping members?
# str_subset(allpgroups,multids[2])
# all_ids_tall%>%filter(pid%in%multids)%>%group_by(pid)%>%summarize(nag=n() - sum(is.na(gene_names)))%>%.$nag%>%table
# multgroups <- allpgroups[sep_element_in(allpgroups,multids)]
# multgroups_ribo <- allpgroups[sep_element_in(allpgroups,multids)]
# multgroups_ribo%>%n_distinct


# #the groups
# #all ids
# #how many ids appear in more than one group globally
# #how many appear in 


# #now fit the more complex model - that accounts for changes in stochiometry over time

# #use the results of this second model to get maximum likelihood estimates of the level of each protein over time

# #choose a reference, the protein which is closest to the median level say

# #now make our stochiometry matrix.

# #now figure out which proteins have significant timepoint effects

# #and annotate these as significantly changing.

# ms_tall_trans%>%group_slice(1:2)
    

# colonlist<-ms_data_all$Protein_IDs

# colonlistmatch(c('a;b;c','b;c','a;d'),c('a'))