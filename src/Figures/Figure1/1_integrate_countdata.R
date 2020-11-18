################################################################################
########This (hopefully final) version of the script uses deepshape prime
########to quantify Riboseq data
################################################################################
base::source(here::here('src/R/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("/fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/Figures/Figure0/0_load_annotation.R")
}


################################################################################
########Load dpprime data, write as salmon files
################################################################################

dpfiles = Sys.glob(here('pipeline/deepshapeprime/*ribo*/run*'))%>%
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
	#now import	
	# tximport(files=fakesalmonfile,type='salmon',tx2gene=dpexprdata_s%>%select(transcript_id,gene_id))
})

# test_that({
# 	all(dpexprdata$transcript_id %in% salmontrs)
# 	all(names(cdsgrl) %in% salmontrs)
# })
################################################################################
########Import it all to get tr length scaled counts
################################################################################
{	
library(tximport)
library(tidyverse)

rnasalmonfiles = Sys.glob(here('pipeline/salmon/data/*/quant.sf'))

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
# test_that({
	# txIdCol
	# inclusiontable(raw[[txIdCol]],tx2genetbl[[1]])
	# inclusiontable(txId,tx2genetbl[[1]])
# 
# })
randomround = function(x)floor(x)+rbinom(length(x),1,x%%1)
for (i in 1:ncol(tx_countdata$counts)){ tx_countdata$counts[,i]%<>%randomround}
# tx_countdata$counts%<>%as.data.frame%>%mutate_all(randomround)%>%as.matrix%>%set_rownames(rownames(tx_countdata$counts))

tx_countdata$counts%>%as.data.frame%>%
	rownames_to_column('gene_id')%>%
	mutate(gene_name = gid2gnm[[gene_id]])%>%
	select(-gene_id)%>%select(gene_name,everything())%>%
	mutate_at(vars(-gene_name),list(randomround))%>%
	write_tsv('data/tx_scaled_countData.tsv')

}
alldpgeneids = tx2genetbl
tx_countdata$abundance%>%rownames%>%inclusiontable(allgids)

stopifnot(tx_countdata$abundance%>%rownames%>%setequal(allgids))


################################################################################
########Save as exprset object
################################################################################
ribosamples = str_subset(colnames(tx_countdata$counts),'ribo')
rnasamples = str_subset(colnames(tx_countdata$counts),'total')


ishighcount = tx_countdata$counts%>%
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

allcountmat <- tx_countdata$counts
allcountdesign = colnames(allcountmat)%>%data.frame(sample=.)%>%separate(sample,into=c('time','assay','rep'),remove=F)
allcountdesign = allcountdesign%>%arrange(assay=='ribo')%>%mutate(assay=as_factor(assay))%>%as.data.frame%>%set_rownames(.$sample)

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

countexprdata <- ExpressionSet(
	allcountmat,
	AnnotatedDataFrame(allcountdesign),
	AnnotatedDataFrame(featuredata)
)

countexprdata%>%saveRDS(here('pipeline/exprdata/countexprset.rds'))

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
########Run Limma (different variances)
################################################################################
tx_countdata_highcounts = tx_countdata%>%map_if(is.matrix,~ .[highcountgenes,])


library(limma)
library(edgeR)
allcountdge <- DGEList(tx_countdata_highcounts$counts)
allcountkeep = filterByExpr(allcountdge)
allcountkeep[] = TRUE #NB - sicne we filter above with highcountgenes, we can just do this
allcountdge = allcountdge[allcountkeep,]
allcountdge <- calcNormFactors(allcountdge)

ribosamples = str_subset(colnames(allcountdge),'ribo')
rnasamples = str_subset(colnames(allcountdge),'total')

ribovoom = voom(allcountdge[,ribosamples],design=model.matrix(~rep+time,allcountdesign[ribosamples,]))
rnavoom = voom(allcountdge[,rnasamples],design=model.matrix(~rep+time,allcountdesign[rnasamples,]))
allvoom = voom(allcountdge[,c(rnasamples,ribosamples)],design=model.matrix(~rep:time+assay*time,allcountdesign[c(rnasamples,ribosamples),]))

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
alllimmares%<>%rownames_to_column('gene_id')%>%mutate(gene_id=str_replace(gene_id,'\\..*',''))
alllimmares%<>%mutate(gene_name=gid2gnm[[gene_id]])


limmatechangesumtable <- alllimmares%>%
	rename('adj_p_value' := adj.P.Val,'log2fc':=logFC)%>%
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
		topTable(eBayes(contrasts.fit(allcountebayes,tpavdesign[,datagroup,drop=F])),coef=1,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
#
countpred_df%>%saveRDS('data/countpred_df.rds')
tx_countdata%>%saveRDS('data/tx_countdata.rds')
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
names(contrnames) = allcountebayes$design%>%colnames%>%str_subset(neg=T,'I|rep')%>%str_replace('assayribo','TE')%>%str_replace(':','_')%>%str_replace('time','')%>%str_replace('(?=^[^T])','all_')
contr = contrnames[1]

stopifnot(contr %in% colnames(allcountebayes$design))

countcontr_df<-	lapply(contrnames,function(contr){
		topTable(allcountebayes,coef=contr,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
countcontr_df%<>%separate(contrast,c('assay','time'))


ribocountebayes <- eBayes(lmFit(allvoom[,11:20]))

ribocontr_df<-	lapply(colnames(ribocountebayes$design)%>%str_subset(neg=T,'Inter|rep|ribo')%>%setNames(paste0('ribo_',tps[-1])),function(contr){
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
tpcontrnames = desmat%>%colnames%>%str_subset('^assayribo:time')
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
names(contrnames) = stepallebayes$design%>%colnames%>%str_subset(neg=T,'I|rep')%>%str_replace('assayribo','TE')%>%str_replace(':','_')%>%str_replace('time','')%>%str_replace('(?=^[^T])','all_')

stepcountcontrdf<-	lapply(contrnames,function(contr){
		topTable(stepallebayes,coef=contr,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
stepcountcontrdf%<>%separate(contrast,c('assay','time'))


ribocountebayes <- eBayes(lmFit(allvoom[,11:20]))

ribocontr_df<-	lapply(colnames(ribocountebayes$design)%>%str_subset(neg=T,'Inter|rep|ribo')%>%setNames(paste0('ribo_',tps[-1])),function(contr){
		topTable(ribocountebayes,coef=contr,number=Inf,confint=.95)%>%
		as.data.frame%>%rownames_to_column('gene_id')
	})%>%bind_rows(.id='contrast')	
ribocontr_df%<>%separate(contrast,c('assay','time'))

stepcountcontrdf <- bind_rows(stepcountcontrdf,ribocontr_df)
eBayes(lmFit(allvoom))

stepcountcontrdf %>% saveRDS(here('data/stepcountcontrdf.rds'))
}
}
countcontr_df%>%group_by(gene_id)%>%group_slice(1)
stepcountcontrdf%>%group_by(gene_id)%>%group_slice(1)

allvoom %>% saveRDS(here('data/allvoom.rds'))


#data.frame(trid=besttrs)%>%write_tsv('pipeline/scikitribotrs.txt',col_names=F)

#system('grep -f pipeline/scikitribotrs.txt annotation/gencode.vM12.annotation.gtf > pipeline/scikitribo.gtf')
#'pipeline/scikitribo.gtf'%>%{rtracklayer::import(.)}
#oanno <- 'annotation/gencode.vM12.annotation.gtf'
#oanno %<>% rtracklayer::import(.)
#oanno[oanno$transcript_id%>%trimids%>%is_in(besttrs)]%>%{rtracklayer::export(.,'pipeline/scikitribo.gtf')}

if(F)test_that({
	  alllimmate = alllimmares%>%filter(contrast%>%str_detect('assayribo:time'))%>%mutate(time=str_replace(contrast,'assayribo:time',''))
		alllimmate%>%head
	  alllimmate%<>%rename('adj_p_value' := adj.P.Val,'log2fc':=logFC)

  alllimmate%>%mutate(sig = (adj_p_value < 0.05)& (abs(log2fc)>log2(1.5)))%>%
    mutate(change = ifelse(!sig,0,ifelse(log2fc>0,1,-1)))%>%
    select(gene_name,time,change)%>%
    group_by(gene_name)%>%
    summarise(techangepat = paste(change,collapse=','))%>%
    .$techangepat%>%table%>%sort


	alllimmares%>%filter(gene_name=='Satb2')%>%as.data.frame
})



# save.image('data/1_integrate_countdata.R')
# load('data/1_integrate_countdata.R')

#next
# /fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/R/TE_change/run_xtail.R
# 
