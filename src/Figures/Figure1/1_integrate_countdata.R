################################################################################
########This (hopefully final) version of the script uses deepshape prime
########to quantify Riboseq data
################################################################################
base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source(here("src/Figures/Figure0/0_load_annotation.R"))
}


################################################################################
########Load dpprime data, write as salmon files
################################################################################

dpfiles = Sys.glob(here('pipeline/deepshapeprime/*/run*'))%>%
	data.frame(file=.)%>%
	mutate(sample=basename(dirname(file)))%>%
	mutate(name=basename(file))%>%
	mutate(number=str_extract(name,'\\d+'))%>%
	arrange(desc(as.numeric(number)))%>%
	group_by(sample)%>%slice(1)%>%{setNames(.$file,.$sample)}

dpexprdata = dpfiles%>%map_df(.id='sample',fread)
dpexprdata%<>%set_colnames(c('sample','transcript_id','TPM','diff','count'))
dpexprdata%<>%mutate(oldtrid=transcript_id,transcript_id=str_replace(transcript_id,'\\.\\d+',''))
stopifnot(dpexprdata$transcript_id%>%unique%>%setdiff(alltrs)%>%length%>%`==`(273L))
dpexprdata%<>%filter(transcript_id%in%alltrs)
dpexprdata%<>%mutate(gene_id = trid2gid[[transcript_id]])

notrimtr2g = dpexprdata%>%distinct(oldtrid,gene_id)
dpexprdata_s = dpexprdata%>%split(.,.$sample)%>%head(1)
dpexprdata%<>%tibble

dpexprdata%<>%safe_left_join(cdsgrl%>%width%>%sum%>%enframe('transcript_id','length'))

library(tximport)

dpoutfiles = dpexprdata%>%split(.,.$sample)%>%imap_chr(function(dpexprdata_s,sample){
	fakesalmonfile=str_interp('pipeline/deepshapeprime/fakesalmonfiles/${sample}/quant.sf')
	message(fakesalmonfile)
	fakesalmonfile%>%dirname%>%dir.create(showWarn=F,rec=TRUE)
	data.frame(
		Name=dpexprdata_s$transcript_id,
		length=dpexprdata_s$length,
		EffectiveLength=dpexprdata_s$length,
		TPM=dpexprdata_s$TPM,
		NumReads=dpexprdata_s$count
	)%>%write_tsv(fakesalmonfile)
	invisible()
	fakesalmonfile
})

################################################################################
########Import it all to get tr length scaled counts
################################################################################
{	
library(tximport)
library(tidyverse)

rnasalmonfiles = Sys.glob(here('pipeline/salmon/data/*total*/quant.sf'))

allquantfiles = c(rnasalmonfiles,dpoutfiles)
names(allquantfiles) <- allquantfiles%>%dirname%>%basename

dptrs = dpoutfiles[[1]]%>%fread%>%.$Name
salmontrs = allquantfiles[[1]]%>%fread%>%.$Name%>%str_extract('ENSMUST\\w+')
inclusiontable(dptrs,salmontrs)
trs = intersect(dptrs,salmontrs)
tx2genetbl = cds%>%mcols%>%as.data.frame%>%select(transcript_id,gene_id)
tx_countdata = tximport(files=allquantfiles,
	ignoreTxVersion=TRUE,
	tx2gene=tx2genetbl,
	type='salmon',
	countsFromAbundance='scaledTPM',
	importer=function(file){
		read_tsv(file,col_types=cols())%>%
			mutate(Name=str_extract(Name,'ENSMUST\\w+'))%>%
			filter(Name%in%trs)%>%
			arrange(match(Name,trs))
})

randomround = function(x)floor(x)+rbinom(length(x),1,x%%1)
for (i in 1:ncol(tx_countdata$counts)){ tx_countdata$counts[,i]%<>%randomround}

tx_countdata$counts%>%as.data.frame%>%
	rownames_to_column('gene_id')%>%
	mutate(gene_name = gid2gnm[[gene_id]])%>%
	select(-gene_id)%>%select(gene_name,everything())%>%
	mutate_at(vars(-gene_name),list(randomround))%>%
	write_tsv('data/tx_scaled_countData.tsv')

}
# alldpgeneids = tx2genetbl
# allgids = tx2genetbl
# tx_countdata$abundance%>%rownames%>%inclusiontable(allgids)

# stopifnot(tx_countdata$abundance%>%rownames%>%setequal(allgids))

################################################################################
########Save as exprset object
################################################################################
ribosamples = str_subset(colnames(tx_countdata$counts),'ribo')
rnasamples = str_subset(colnames(tx_countdata$counts),'total')
mainsamples = c(ribosamples,rnasamples)

ishighcount = tx_countdata$counts[,ribosamples]%>%
	as.data.frame%>%
	rownames_to_column('gene_id')%>%
	pivot_longer(-gene_id)%>%
	separate(name,c('time','assay','rep'))%>%
	filter(assay=='ribo')%>%
	group_by(gene_id,time)%>%
	summarise(ishighcount=sum(value)>=32)%>%
	{setNames(.$ishighcount,.$gene_id)}

ishighcount = ishighcount[rownames(tx_countdata$counts)]

highcountgenes = rownames(tx_countdata$counts)[ishighcount]
highcountgnms = gid2gnm[[highcountgenes]]

allcountmat <- tx_countdata$counts[,mainsamples]
allcountdesign = colnames(allcountmat)%>%data.frame(sample=.)%>%
	separate(sample,into=c('time','assay','rep'),remove=F)
allcountdesign = allcountdesign%>%arrange(assay=='ribo')%>%
	mutate(assay=as_factor(assay))%>%
	as.data.frame%>%set_rownames(.$sample)
allcountmat <- allcountmat[,rownames(allcountdesign)]

{

featuredata = data.frame(
	gene_id = allcountmat%>%rownames,
	gene_name = allcountmat%>%rownames%>%{gid2gnm[[.]]},
	ishighcount = ishighcount
)%>%set_rownames(rownames(allcountmat))

dim(allcountmat)
dim(allcountdesign)
dim(featuredata)
colnames(allcountmat)

colnames(allcountmat)==rownames(allcountdesign)

countexprdata <- ExpressionSet(
	allcountmat,
	AnnotatedDataFrame(allcountdesign),
	AnnotatedDataFrame(featuredata)
)

countexprdata%>%saveRDS(here('data/countexprset.rds'))

}

################################################################################
########Also get isoform level 
################################################################################
tx2genemap=data.frame(trid2gid$keys(),trid2gid$values())
iso_tx_countdata = tximport(files=allquantfiles,
	txOut=TRUE,
	ignoreTxVersion=TRUE,
	tx2gene=tx2genemap,
	type='salmon',
	countsFromAbundance='scaledTPM',,
	importer=function(file){
		read_tsv(file,col_types=cols())%>%
			mutate(Name=str_extract(Name,'ENSMUST\\w+'))%>%
			filter(Name%in%alltrs)%>%
			arrange(match(Name,alltrs))
})

stopifnot(iso_tx_countdata$abundance%>%rownames%>%setequal(alltrs))

################################################################################
########Now run limma on the main samples
################################################################################
library(limma)
library(edgeR)
allcountdge <- DGEList(tx_countdata$counts)
allcountkeep = filterByExpr(allcountdge)
allcountkeep[] = TRUE #NB - sicne we filter above with highcountgenes, we can just do this
allcountdge = allcountdge[allcountkeep,]
allcountdge <- calcNormFactors(allcountdge)

monosamples = str_subset(colnames(allcountdge),'80S')
rnasamples = str_subset(colnames(allcountdge),'total')

ribodesign = model.matrix(~rep+time,allcountdesign[ribosamples,])
ribovoom = voom(allcountdge[,ribosamples],design=ribodesign)
rnadesign = model.matrix(~rep+time,allcountdesign[rnasamples,])
rnavoom = voom(allcountdge[,rnasamples],design=rnadesign)
fulldesign = model.matrix(~rep:time+assay*time,
		allcountdesign[c(rnasamples,ribosamples),]
	)
allvoom = voom(allcountdge[,c(rnasamples,ribosamples)],design=fulldesign)

stopifnot(colnames(allvoom$E)==c(rnasamples,ribosamples))
allvoom$weights <- cbind(ribovoom$weights,rnavoom$weights)
allcountebayes <- eBayes(lmFit(allvoom))

timecoefs = allcountebayes$coef%>%colnames%>%str_subset('^time')
tetimecoefs = allcountebayes$coef%>%colnames%>%str_subset('assayribo:time')
allcontrasts = c(timecoefs,tetimecoefs)%>%setNames(.,.)

alllimmares = 
	map_df(.id='contrast',allcontrasts,function(contrasti){
		topTable(contrasts.fit(allcountebayes,coef=c(contrasti)),n=Inf,confint=.95)
	})
#unmangle the gene ids
alllimmares%<>%rownames_to_column('gene_id')%>%
	mutate(gene_id=str_replace(gene_id,'\\..*',''))
alllimmares%<>%mutate(gene_name=gid2gnm[[gene_id]])


limmatechangesumtable <- alllimmares%>%
	dplyr::rename('adj_p_value' := adj.P.Val,'log2fc':=logFC)%>%
	group_by(gene_id,gene_name)%>%
  mutate(sig = (adj_p_value < 0.05)& (abs(log2fc)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (log2fc > 0))),
    down = as.numeric(any(sig & (log2fc < 0)))
  )

limmatechangesumtable%>%write_tsv(here('tables/limmaTEchange.tsv'))

besttrs = iso_tx_countdata$abundance%>%
	as.data.frame%>%
	rownames_to_column('transcript_id')%>%
	left_join(tx2genemap%>%set_colnames(c('transcript_id','gene_id')))%>%
	pivot_longer(-one_of('gene_id','transcript_id'))%>%
	separate(name,c('time','assay','rep'))%>%
	group_by(gene_id,transcript_id)%>%summarise(value=median(value))%>%
	group_by(gene_id)%>%slice(which.max(value))%>%.$transcript_id

#get the designs describing the average at a given timepoint/assay
tpavdesign = allcountebayes$design[rownames(allcountebayes$design)%>%str_subset('_2$'),]
rownames(tpavdesign) %<>% str_replace('_2','')
tpavdesign[,colnames(tpavdesign)%>%str_detect('rep2:')]%<>%multiply_by(0.5)
tpavdesign%<>%t
datagroup = colnames(tpavdesign)[1]
ribcols <- colnames(tpavdesign)%>%str_subset('ribo')
totcols <- colnames(tpavdesign)%>%str_subset('total')
tedesign = tpavdesign[,ribcols] - tpavdesign[,totcols]
colnames(tedesign)%<>%str_replace('ribo','TE')
tpavdesign%<>%cbind(tpavdesign, tedesign)

countpred_df<-	lapply(colnames(tpavdesign)%>%setNames(.,.),function(datagroup){
		message(datagroup)
		contrlm_object = contrasts.fit(allcountebayes,tpavdesign[,datagroup,drop=F])
		topTable(eBayes(contrlm_object),coef=1,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
#
countpred_df$gene_name = gid2gnm[[countpred_df$gene_id]]
countpred_df%>%saveRDS('data/countpred_df.rds')
tx_countdata%>%map_if(is.matrix,~.[,mainsamples])%>%saveRDS('data/tx_countdata.rds')
allvoom%>%saveRDS('data/allvoom.rds')


# contrnames = as.list(contrnames)
# contrnames = c(contrnames,lapply(tps[-1],function(tp){
	# contrnames%>%str_subset(tp)
# })%>%setNames(paste0('ribo_',tps[-1])))
################################################################################
########Contrasts df, (including riboseq, rather than TE)
################################################################################

	
{
tps = allcountdesign$time%>%unique
contrnames = allcountebayes$design%>%colnames%>%str_subset(neg=T,'I|rep')
contrnamesassay = contrnames%>%str_detect('assayribo')%>%ifelse('TE','all')
contrnamestime = contrnames%>%str_extract('(?<=time)\\w+')

names(contrnames) = paste(sep='_',contrnamesassay,contrnamestime)
names(contrnames)%<>%str_replace('_NA','')

contr = contrnames[1]

stopifnot(contr %in% colnames(allcountebayes$design))

countcontr_df<-	lapply(contrnames,function(contr){
		topTable(allcountebayes,coef=contr,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
countcontr_df%<>%separate(contrast,c('assay','time'))

countcontr_df%>%filter(time=='E145',gene_id=='ENSMUSG00000000794')


ribocountebayes <- eBayes(lmFit(ribovoom))

ribocontr_names <- colnames(ribocountebayes$design)%>%
	str_subset(neg=T,'Inter|rep|ribo')%>%
	setNames(paste0('ribo_',tps[-1]))
ribocontr_df<-	lapply(ribocontr_names,function(contr){
		topTable(ribocountebayes,coef=contr,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
ribocontr_df%<>%separate(contrast,c('assay','time'))

countcontr_df <- bind_rows(countcontr_df,ribocontr_df)

countcontr_df %>% saveRDS(here('data/countcontr_df.rds'))

}

################################################################################
########Stepwise version
################################################################################
{
desmat = allvoom$design
tpcontrnames = desmat%>%colnames%>%str_subset('^time')
for(i in 2:length(tpcontrnames)){
	haslatr = desmat[,tpcontrnames[i]]==1
	for(j in (i-1):1){
		desmat[haslatr,tpcontrnames[j]]<- 1
	}
}
tpcontrnames = desmat%>%colnames%>%str_subset('^assayribo:time|(time.*:assayribo)')
for(i in 2:length(tpcontrnames)){
	haslatr = desmat[,tpcontrnames[i]]==1
	for(j in (i-1):1){
		desmat[haslatr,tpcontrnames[j]]<- 1
	}
}
allvoomstep = allvoom
allvoomstep$design = desmat

stepallebayes = eBayes(lmFit(allvoomstep))

{
contrnames = stepallebayes$design%>%colnames%>%str_subset(neg=T,'I|rep')
names(contrnames) = stepallebayes$design%>%colnames%>%
	str_subset(neg=T,'I|rep')%>%str_replace('assayribo','TE')%>%
	str_replace(':','_')%>%
	str_replace('time','')%>%str_replace('(?=^[^T])','all_')

stepcountcontrdf<-	lapply(contrnames,function(contr){
		topTable(stepallebayes,coef=contr,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
stepcountcontrdf%<>%separate(contrast,c('assay','time'))


ribocountebayes <- eBayes(lmFit(ribovoom))
ribocontr_names <- colnames(ribocountebayes$design)%>%
	str_subset(neg=T,'Inter|rep|ribo')%>%
	setNames(paste0('ribo_',tps[-1]))
ribocontr_df<-	lapply(ribocontr_names,function(contr){
		topTable(ribocountebayes,coef=contr,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
ribocontr_df%<>%separate(contrast,c('assay','time'))

stepcountcontrdf <- bind_rows(stepcountcontrdf,ribocontr_df)
eBayes(lmFit(allvoom))

stepcountcontrdf %>% saveRDS(here('data/stepcountcontrdf.rds'))

}
}

allvoom %>% saveRDS(here('data/allvoom.rds'))
iso_tx_countdata %>% saveRDS(here('data/iso_tx_countdata.rds'))

highcountgnms %>% saveRDS(here('data/highcountgnms.rds'))




save.image('data/1_integrate_countdata.R')
# load('data/1_integrate_countdata.R')

#next
# src/R/TE_change/run_xtail.R
# 


################################################################################
########Run Limma (different variances)
################################################################################
satb2trid = 'ENSMUST00000177424'
iso_tx_countdata$counts[satb2trid,]

{
tx_countdata_highcounts = tx_countdata%>%map_if(is.matrix,~ .[highcountgenes,])

monosamples = tx_countdata$abundance%>%colnames%>% str_subset('80S')
polysamples = tx_countdata$abundance%>%colnames%>% str_subset('Poly')
fracsamples = c(monosamples,polysamples)

ribofrac_design = data.frame(sample=fracsamples)%>%	
	mutate(fraction = str_extract(sample,'80S|Poly'))%>%
	mutate(time = sample%>%str_extract('(?<=80S|Poly).*?(?=_)'))%>%
	mutate(rep = sample%>%str_extract('[0-9]+$'))%>%
	arrange(fraction=='80S',time=='P0',time=='E16')%>%
	set_rownames(.$sample)
ribofrac_design$fraction%<>%as_factor

fracvoom = voom(allcountdge[,ribofrac_design$sample],design=model.matrix(~rep:time+fraction*time, ribofrac_design))

stopifnot(colnames(fracvoom$E)==ribofrac_design$sample)
fr_allcountebayes <- eBayes(lmFit(fracvoom))

fr_timecoefs = fr_allcountebayes$coef%>%colnames%>%str_subset('^time')
fr_tetimecoefs = fr_allcountebayes$coef%>%colnames%>%str_subset('assayribo:time')
fr_allcontrasts = c('fraction80S',fr_timecoefs,fr_tetimecoefs)%>%setNames(.,.)

fr_alllimmares = 
	map_df(.id='contrast',fr_allcontrasts,function(contrasti){
		topTable(contrasts.fit(fr_allcountebayes,coef=c(contrasti)),n=Inf,confint=.95)
	})
#unmangle the gene ids
fr_alllimmares%<>%rownames_to_column('gene_id')%>%mutate(gene_id=str_replace(gene_id,'\\..*',''))
fr_alllimmares%<>%mutate(gene_name=gid2gnm[[gene_id]])

fr_alllimmares %>% saveRDS(here('data/fr_alllimmares.rds'))


fr_limmatechangesumtable <- fr_alllimmares%>%
	dplyr::rename('adj_p_value' := adj.P.Val,'log2fc':=logFC)%>%
	group_by(gene_id,gene_name)%>%
  mutate(sig = (adj_p_value < 0.05)& (abs(log2fc)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (log2fc > 0))),
    down = as.numeric(any(sig & (log2fc < 0)))
  )

fr_limmatechangesumtable%>%write_tsv(here('tables/fr_limmaTEchange.tsv'))
stop()

#get the designs describing the average at a given timepoint/assay
f_tpavdesign = fr_allcountebayes$design[rownames(fr_allcountebayes$design)%>%str_subset('_2$'),]
rownames(f_tpavdesign) %<>% str_replace('_2','')
f_tpavdesign[,colnames(f_tpavdesign)%>%str_detect('rep2:')]%<>%multiply_by(0.5)
f_tpavdesign%<>%t
datagroup = colnames(f_tpavdesign)[1]
monocols <- colnames(f_tpavdesign)%>%str_subset('80S')
polycols <- colnames(f_tpavdesign)%>%str_subset('Poly')
mfracdesign = f_tpavdesign[,monocols] - f_tpavdesign[,polycols]
colnames(mfracdesign)%<>%str_replace('RP.*80S','MonoFrac_')
f_tpavdesign%<>%cbind(f_tpavdesign, mfracdesign)
#
fr_countpred_df<-	lapply(colnames(f_tpavdesign)%>%setNames(.,.),function(datagroup){
		message(datagroup)
		topTable(eBayes(contrasts.fit(fr_allcountebayes,f_tpavdesign[,datagroup,drop=F])),coef=1,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
#
fr_countpred_df$gene_name = gid2gnm[[fr_countpred_df$gene_id]]
fr_countpred_df%>%saveRDS('data/fr_countpred_df.rds')
tx_countdata%>%map_if(is.matrix,~.[,fracsamples])%>%saveRDS('data/fr_tx_countdata.rds')
allvoom%>%saveRDS('data/allvoom.rds')

}
	
{

ftps = ribofrac_design$time%>%unique
f_contrnames = fr_allcountebayes$design%>%colnames%>%str_subset(neg=T,'I|rep')
names(f_contrnames) = fr_allcountebayes$design%>%colnames%>%str_subset(neg=T,'I|rep')%>%
	str_replace('(.*)fraction80S','MonoFrac_\\1')%>%
	str_replace(':','_')%>%
	str_replace('time','')%>%
	str_replace('_$','')%>%
	str_replace('^(E16|P0)$','time_\\1')%T>%print
	
contr = f_contrnames[1]

stopifnot(contr %in% colnames(fr_allcountebayes$design))

fr_countcontr_df<-	lapply(f_contrnames,function(contr){
		topTable(fr_allcountebayes,coef=contr,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
fr_countcontr_df%<>%separate(contrast,c('assay','time'))
fr_countcontr_df%>%distinct(time,assay)

monovoom = voom(allcountdge[,monosamples],design=model.matrix(~time, ribofrac_design[monosamples,]))
monocountebayes <- eBayes(lmFit(monovoom))

monocontrasts = colnames(monocountebayes$design)%>%setNames(.,str_replace(.,'time','mono_')%>%str_replace('\\(Intercept\\)','mono_E13'))

monocontr_df<-	lapply(monocontrasts,function(contr){
		topTable(monocountebayes,coef=contr,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
monocontr_df%<>%separate(contrast,c('assay','time'))

fr_countcontr_df <- bind_rows(fr_countcontr_df,monocontr_df)
fr_countcontr_df$contrast%>%unique
fr_countcontr_df %>% saveRDS(here('data/fr_countcontr_df.rds'))
fr_countcontr_df%>%distinct(time,assay)
fr_countcontr_df$gene_name = gid2gnm[[fr_countcontr_df$gene_id]]

}
################################################################################
########Stepwise
################################################################################
	
{
desmat = fracvoom$design
tpcontrnames = desmat%>%colnames%>%str_subset('^time')
for(i in 2:length(tpcontrnames)){
	haslatr = desmat[,tpcontrnames[i]]==1
	for(j in (i-1):1){
		desmat[haslatr,tpcontrnames[j]]<- 1
	}
}
tpcontrnames = desmat%>%colnames%>%str_subset('^time.*:fraction80S')
for(i in 2:length(tpcontrnames)){
	haslatr = desmat[,tpcontrnames[i]]==1
	for(j in (i-1):1){
		desmat[haslatr,tpcontrnames[j]]<- 1
	}
}
fracvoomstep = fracvoom
fracvoomstep$design = desmat

fr_stepallebayes = eBayes(lmFit(fracvoomstep))
}
{
f_contrnames = fr_stepallebayes$design%>%colnames%>%str_subset(neg=T,'I|rep')
names(f_contrnames) = fr_allcountebayes$design%>%colnames%>%str_subset(neg=T,'I|rep')%>%
	str_replace('(.*)fraction80S','MonoFrac_\\1')%>%
	str_replace(':','_')%>%
	str_replace('time','')%>%
	str_replace('_$','')%>%
	str_replace('^(E16|P0)$','time_\\1')%T>%print

fr_stepcountcontrdf<-	lapply(f_contrnames,function(contr){
		topTable(fr_stepallebayes,coef=contr,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
fr_stepcountcontrdf%<>%separate(contrast,c('assay','time'))
fr_stepcountcontrdf%>%distinct(assay,time)

fr_stepcountcontrdf %>% saveRDS(here('data/fr_stepcountcontrdf.rds'))
}

