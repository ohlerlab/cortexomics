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

#' setup, eval = TRUE
sep_element_in<-function(colonlist,ridssplit,sep=';'){
	assert_that(is.character(colonlist))
	assert_that(is.character(ridssplit))
	values<-colonlist%>%str_split(sep)

	inds <- rep(seq_along(colonlist),lengths(values))

	values<-flatten_chr(values)

	data_frame(inds,values)%>%
		mutate(match = values %in% ridssplit)%>%
		group_by(inds)%>%
		summarize(match=any(match))%>%
		.$match

}



ms_tall_ib <- ms_tall%>%
	ungroup%>%
	filter(sigtype=='LFQ')%>%
	filter(fraction=='total')%>%
	filter(Protein_IDs %in% unique(Protein_IDs))

# stopifnot(
# 	ms_tall_ib$gene_name%>%
# 	na.omit%>%unique%>%
# 	str_split(';')%>%unlist%>%
# 	table%>%`==`(1)
# )

#for now, we'll just keep those rows which have riboproteins on them,
# n_groups = n_distinct(ms_tall$Protein_IDs)
# n_protid = ms_tall%>%str_split(';')%>%unlist%>%n_distinct
# allpgroups <- ms_tall$Protein_IDs%>%unique
# multids<-allpgroups%>%unique%>%str_split(';')%>%unlist%>%table%>%keep(~ . > 1)%>%names
# all_ambig_pgroups<-allpgroups%>%sep_element_in(multids)
# ms_tall_ib%<>%mutate(ambig = sep_element_in(Protein_IDs,multids))

#then we can just keep the rp name and shorten them to the first element otherwise
# ms_tall_ib%<>%mutate(gene_name_simp=str_replace(gene_name,regex('(.*?)([Rr]p[sl][^;]+)(.*)'),'\\2'))
ms_tall_ib%<>%mutate(gene_name_simp=gene_name)
# ms_tall_ib$gene_name_simp%<>%str_replace(';.*$','')
#unique correspondance between simplified gene names and protein IDS
# stopifnot(ms_tall_ib%>%filter(!ambig)%>%distinct(Protein_IDs,gene_name_simp)%>%map_lgl(.%>%anyDuplicated%>%`==`(1))%>%not%>%all)


#looking at numbers of missing genes.
missing<-ms_tall_ib%>%
	filter(!is.na(gene_name_simp))%>%
	group_by(gene_name_simp,time,dataset)%>%
	summarize(signal=sum(na.omit(signal)))%>%
	mutate(signal = ifelse(signal==0,NA,signal))%>%#sum signal for genes with multiple rows
	select(gene_name_simp,dataset,time,signal)

missing%>%	group_by(gene_name_simp,time)%>%summarise(nanum=sum(is.na(signal)))%>%
	group_by(gene_name_simp)%>%summarise(miss=any(nanum>1))%>%
	filter(miss)


#fit the most basic model where there are no timepoint specific factors
monomat = ms_tall_ib%>%
	filter(!is.na(gene_name_simp))%>%
	group_by(gene_name_simp,dataset)%>%
	summarize(signal=sum(na.omit(signal)))%>%mutate(signal = ifelse(signal==0,NA,signal))%>%
	select(gene_name_simp,dataset,signal)%>%
	spread(dataset,signal)%>%
	{ set_rownames(as.matrix(.[-1]),.[[1]]) }



# #let's look at the median counts and the number of NAs
# exprmat%>%map(.%>%apply(1,.%>%is.na%>%sum))
# exprmat%>%map(.%>%apply(1,.%>%median))

monomat%>%apply(1,.%>%is.na%>%sum%>%`<=`(1))%>%sum

#select 1? Or just 
monomat <- monomat[monomat%>%apply(1,.%>%is.na%>%sum%>%`<=`(1)),]
monomat <- monomat[monomat%>%apply(1,median)%>%`>`(10),]

#take the log
lmonomat = log2(monomat)

means <- 	lmonomat%>%apply(1,mean,na.rm=TRUE)
sds <- 	lmonomat%>%apply(1,sd,na.rm=TRUE)


#okay ten we've reliminated all that messy mean-variance.

#Now let's make our model matrix
designmat <- ms_tall_ib%>%
	select(dataset,time,fraction)%>%
	filter(dataset%in%colnames(monomat))%>%
	distinct(dataset,.keep_all=TRUE)%>%
	select(-dataset)%>%
	select_if(~ n_distinct(.)>1 )

design = model.matrix( ~ time , designmat )

library(limma)
lmonomat%>%nrow
fit = limma::lmFit(lmonomat,design=design)
bayesfit = limma::eBayes(fit,trend=TRUE, robust=TRUE)
coefs<-bayesfit$p.value%>%colnames
allcoefs<-coefs%>%setNames(.,.)%>%map(~topTable(bayesfit,number=nrow(lmonomat),coef=.,confint=0.95)%>%{cbind(gene=rownames(.),.)})%>%
	bind_rows(.id='coefficient')%>%
	as_data_frame

allcoefs%>%write_tsv('ms_limma_coefs.tsv')

allcoefs%>%sample_n(10)

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

datadf<-ms_tall_ib %>% 
	filter(fraction=='total')%>%
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

traj.plots%>%ggsave(file='~/projects/cortexomics/figures/all_ms_trajplots.pdf')













# ##looking at ambiguities



# P14869;S4R1N1;D3YVM5

# ms_data%>%filter(Protein_IDs=='')



# ms_tall_ib%>%select(gene_name_simp,dataset,signal)%>%spread(dataset,signal)

# ms_tall_ib%>%filter(Protein_IDs=='A2A547;P84099')%>%.$gene_name
# ms_tall_ib%>%filter(Protein_IDs=='Q9CZM2;E9QAZ2;B8JKK2')%>%.$gene_name
# ms_tall_ib%>%filter(gene_name=='Rpl19')%>%.$Protein_IDs
# ms_tall_ib%>%select(Protein_IDs,gene_name)
# ms_tall_ib%>%select(Protein_IDs,gene_name)%>%distinct
# ms_tall_ib%>%select(gene_name)%>%distinct
# ms_tall_ib%>%select(Protein_IDs)%>%distinct

# #do we have multiple protein IDs corresponding to this rp gene?
# allpgroups      <- ms_tall$Protein_IDs%>%unique
# allpids         <- allpgroups%>%str_split(';')%>%unlist%>%unique
# allpgroups_ribo <- allpgroups[sep_element_in(allpgroups,ridssplit)]
# allpids_wribo    <- allpgroups_ribo%>%str_split(';')%>%unlist%>%unique
# allpids_ribo    <- ridssplit
# allpgene_names <- ms_tall$gene_name%>%str_extract('[Rr]p[sl]\\d+')%>%.[!is.na(.)]%>%unique
# multids<-allpgroups%>%unique%>%str_split(';')%>%unlist%>%table%>%keep(~ . > 1)%>%names


# #ribosomal ids which appear in more than one group
# multids_ribo<-allpgroups_ribo%>%unique%>%str_split(';')%>%unlist%>%table%>%keep(~ . > 1)%>%names

# all_ids_tall<-ms_data_all%>%
# 	select(gene_names=gene_name.x,pids=Protein_IDs)%>%
# 	mutate(row=seq_len(n()))%>%
# 	# filter(!is.na(gene_names))%>%
# 	mutate(genename=str_split(gene_names,';'))%>%
# 	mutate(genenamenum=lengths(genename))%>%
# 	unnest%>%
# 	mutate(pid = str_split(pids,';'))%>%
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
# allpgroups%>%str_split(';')%>%unlist%>%table%>%{.[names(.)%in%allpids_ribo]}%>%{.[.>1]}


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
# # ribopids = ms_tall_ib%>%filter(gene_name%>%str_detect(ribogenename))%>%.$Protein_IDs%>%unique%>%str_split(';')%>%unlist
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

# ms_tall_ib%>%group_slice(1:2)
    

# colonlist<-ms_data_all$Protein_IDs

# colonlistmatch(c('a;b;c','b;c','a;d'),c('a'))