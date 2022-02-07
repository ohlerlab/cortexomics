titanic_wide <- data.frame(Titanic)
head(titanic_wide)

if(!'ggalluvial'%in%installed.packages()) remotes::install_github("corybrunson/ggalluvial", build_vignettes = FALSE)
library(ggalluvial)
#>   Class    Sex   Age Survived Freq
#> 1   1st   Male Child       No    0
#> 2   2nd   Male Child       No    0
#> 3   3rd   Male Child       No   35
#> 4  Crew   Male Child       No    0
#> 5   1st Female Child       No    0
#> 6   2nd Female Child       No    0


# 'plots/figures/figureS1'%>%dir.create
plotfile%>%dirname%>%dir.create(rec=TRUE)

# LOWCOUNTLIM <- 32
# allsegcounts%<>%separate(sample,into=c('time','assay','rep'))%>%group_by(protein_id,assay,time)%>%mutate(hascounts = sum(total)>LOWCOUNTLIM)
# allsegcounts%<>%unite(sample,time,assay,rep)

#in annotation
#in salmon (filterd because too short, etc)

#detect in MS
#detect in rnaseq
#detect in riboseq

#in ms table
#matched to gene
#in selected ms_ids




stepstepprotcontrdf$gene_name%>%n_distinct


{
count_detection_df <-
	tx_countdata$counts%>%
	as.data.frame%>%
	rownames_to_column('gene_id')%>%
	pivot_longer(cols=-gene_id,names_to='dataset',values_to='count')%>%
	# group_by(protein_id)%>%group_slice(1:1000)%>%
	filter(dataset %>%str_detect('ribo|total'))%>%
	separate(dataset,into=c('time','assay','rep'))%>%
	mutate(hascounts = count!= 0)%>%
	distinct(gene_id,time,assay,hascounts)%>%
	group_by(gene_id,time,assay)%>%
	summarise(hascount = any(hascounts))


count_detection_df$gene_id%>%inclusiontable(tx_countdata$counts%>%rownames)

count_detection_df$gene_name=gid2gnm[[count_detection_df$gene_id]]

# sampprotids <- sample(best_uprotein_ids,10)
# sampprotids%<>%str_replace('_\\d+$','')

count_detection_df$hascount%>%table

}


{
proteinmsdata<-readRDS('data/proteinmsdata.rds')

msdetectdf = proteinmsdata%>%
	filter(!is.na(gene_name))%>%
 	select(gene_name,matches('iBAQ.*input'))%>%
 	gather(dataset,signal,,-gene_name)%>%
 	mutate(gene_id=gnm2gid[[gene_name]])%>%
 	separate(dataset,c('sigtype','time','fraction','rep'))%>%
 	group_by(gene_name,gene_id,time)%>%
 	summarise(hascount=any(!is.na(signal)))%>%
 	mutate(time=str_replace(time,'p',''))

msdetectdf$gene_id = gnm2gid[[msdetectdf$gene_name]]
timepoints <- count_detection_df$time%>%unique
allids = count_detection_df%>%distinct(gene_id,gene_name)


msdetectdf%>%group_by(.,gene_id,time,hascount)%>%filter(n()>1)%>%group_slice(1)


msdetection_df<-bind_rows(
	msdetectdf,
	data.frame(gene_id=allids$gene_id%>%unique%>%setdiff(msdetectdf$gene_id)%>%rep(each=5),time=timepoints,hascount=FALSE,assay='MS')
)
msdetection_df$assay='MS'
msdetection_df%>%group_by(.,gene_id,time,hascount)%>%filter(n()>1)%>%group_slice(1)


countriverdfall <- count_detection_df%>%
	bind_rows(msdetection_df)

countriverdfall%>%filter(is.na(gene_id))
countriverdfall%>%filter(is.na(gene_name))

countriverdfall$gene_name=gid2gnm[[countriverdfall$gene_id]]

stopifnot(countriverdfall$assay%>%is.na%>%`!`%>%all)

countriverdfall%>%group_by(.,gene_id,time,assay,hascount)%>%filter(n()>1)%>%group_slice(1)
countriverdfall%>%distinct(.,gene_id,time,assay,hascount)%>%nrow
countriverdfall%>%nrow

stopifnot(countriverdfall%>%{nrow(distinct(.,gene_id,time,assay,hascount))==nrow(.)  })


countriverdf<-countriverdfall%>%	group_by(assay,gene_id)%>%filter(!all(hascount),any(hascount))%>%
	mutate(hascount = ifelse(hascount,'Detected','Not-Detected'))%>%
	# ungroup%>%
	spread(time,hascount)%>%
	mutate(Annotated='Detected')

gcols <- c("Annotated",unique(count_detection_df$time))
stopifnot(gcols%in%colnames(countriverdf))

countriverdf%<>%group_by(Annotated,E13,E145,E16,E175,P0,assay)%>%tally

# countriverdf%>%head
countriverdflong <- to_lodes_form(data.frame(countriverdf),
                              key = "Time",
                              axes = 1:6)
countriverdflong%>%head
countriverdf%>%head
countriverdflong%>%group_by(alluvium)%>%group_slice(3)
countriverdflong$Freq = countriverdflong$n
is_alluvia_form(as.data.frame(countriverdf), axes = 1:6, silent = TRUE)

detecttbl_toprint<-countriverdfall%>%group_by(time,assay,hascount)%>%tally%>%filter(hascount)%>%
	select(n_gene_ids_detected=n)%>%
	bind_rows(data.frame(time='Annotation',assay='Annotation',n_gene_ids_detected=allids$gene_id%>%n_distinct),.)


countriverdfall%>%group_by(time,assay,hascount)%>%tally%>%filter(hascount)%>%
	select(n_gene_ids_detected=n)


(here('tables/manuscript/n_detections.tsv'))%>%dirname%>%dir.create
detecttbl_toprint%>%write_tsv(here('tables/manuscript/n_detections.tsv'))


countriverdfall%>%group_by(time,gene_id)%>%summarise(n=n_distinct(assay))%>%filter(n!=3)
countriverdfall%>%group_by(time,gene_id)%>%filter(hascount[assay=='ribo'],!hascount[assay=='total'])%>%distinct(gene_id,time)%>%group_by(time)%>%tally
countriverdfall%>%group_by(time,gene_id)%>%filter(!hascount[assay=='ribo'],hascount[assay=='total'])%>%distinct(gene_id,time)%>%group_by(time)%>%tally


gainlossdfall <- countriverdfall%>%
	group_by(gene_id,time)%>%
	mutate(hascount = hascount & hascount[assay=='ribo'])%>%
	group_by(assay,gene_id)%>%
	mutate(lost = 	(!hascount) & (hascount[time=='E13']))%>%
	mutate(gained = (hascount) & (!(hascount[time=='E13'])))
gainlossdf <- gainlossdfall%>%	group_by(assay,time)%>%
	summarise(still_present = sum(hascount) - sum(gained), lost = sum(lost), gained = sum(gained))
gainlossdf$time[gainlossdf$time!='E13']%<>%paste0(' vs. E13')
E13pres = gainlossdf%>%filter(time=='E13')%>%.$still_present

gainlosspcrows <- data.frame(time=rep('P0 vs E13 (%)',3),
	still_present = gainlossdf%>%filter(time=='P0 vs. E13')%>%{.$still_present/E13pres}%>%round(2),
	lost = gainlossdf%>%filter(time=='P0 vs. E13')%>%{.$lost/E13pres}%>%round(2),
	gained = gainlossdf%>%filter(time=='P0 vs. E13')%>%{.$gained/E13pres}%>%round(2),
	assay=unique(gainlossdf$assay)
)

gainlossdf%>%bind_rows(gainlosspcrows)%>%unite(tmp,still_present,lost,gained)%>%spread(time,tmp)%>%
	mutate(E13 = str_replace(E13,'_0_0',''))%>%
	identity%>%as.data.frame%>%
	.[,colnames(.)%>%str_detect('%')%>%order]%>%
	mutate_at(vars(matches('vs')),list(~str_replace_all(.,'_',',')))%>%
	# mutate_at(vars(matches('vs')),list(~str_replace_all(.,'_',' ')))%>%
	arrange(assay=='MS',assay=='ribo')%>%
	write_tsv('tables/manuscript/assay_time_gainloss.tsv')

'tables/manuscript/assay_time_gainloss.tsv'%>%normalizePath

}
# # countriverdf%>%tail
# # #For each protein ID get a factor 


# titanic_long <- to_lodes_form(data.frame(Titanic),
#                               key = "Demographic",
#                               axes = 1:3)
# titanic_long%>%head
# # (countriverdf%>%.[,-c(1:2)]%>%colSums%>%n_distinct%>%`==`(1)%>%stopifnot)
# countriverdflong


# # countriverdflong%<>%set_colnames(c('assay','stratum','alluvium','Time','Freq'))

# countriverdflong <- to_lodes_form(data.frame(table(assay=countriverdf$assay,E13=countriverdf$E13,E16=countriverdf$E16,P0=countriverdf$P0)),key='Time',axes=2:4)

# countriverdflong%>%head
# countriverdflong$Time
# countriverdflong$Time

{
stopifnot(gainlossdfall%>%group_by(gene_id)%>%tally%>%.$n%>%table%>%n_distinct%>%`==`(1))

filterfuncs <- list(
	'_' = identity,
	'_detectchange_'=	.%>%filter(any(hascount))%>%filter(!all(hascount)),
	'_isdetected_'=.%>%filter(any(hascount))
)
imap(filterfuncs,function(filterfunc,filtername){
	countriverdflong <- gainlossdfall%>%
		group_by(assay,gene_id)%>%
		select(-lost,-gained)%>%
		# group_by(assay,gene_id)%>%§
		filterfunc%>%
		mutate(hascount = ifelse(hascount,'Detected','Non-Detected'))%>%
		spread(time,hascount)%>%
		group_by(assay,E13,E145,E16,E175,P0)%>%
		tally(name = 'Freq')%>%
		to_lodes_form(key = 'Time',axes = 2:6)

	countriverdflong

	plotfile <- paste0('plots/figures/figureS1/detection_riverplot',filtername,'.pdf')
	pdf(plotfile,w=24*.6,h=24*.6)
	p=ggplot(data = countriverdflong,
	       aes(x = Time, stratum = stratum, alluvium = alluvium,
	           y = Freq, label = stratum,fill=stratum,color=stratum)) +
	  # geom_alluvium() +
	  geom_flow(aes(fill = stratum)) +
	  geom_stratum() + 
	  geom_text(stat='stratum',color='Black') +
	  theme_minimal() +
	  facet_grid(assay ~ . )+
	  ggtitle("")
	 print(p)
	dev.off()
	plotfile%>%normalizePath%>%message
})
}





