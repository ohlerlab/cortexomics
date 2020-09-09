################################################################################
########This (hopefully final) version of the script uses deepshape prime
########to quantify Riboseq data
################################################################################
source(here::here('src/R/Rprofile.R'))

#get exons
gtf = here('pipeline/my_gencode.vM12.annotation.gtf')
if(!exists('gtf_gr')) gtf_gr<-rtracklayer::import(con=gtf,format='gtf')
exons <- gtf_gr%>%subset(type=='exon')
cds <- gtf_gr%>%subset(type=='CDS')

 #make gene names non redundant
ids_nrgname <- mcols(cds)%>%as.data.frame%>%distinct(gene_name,gene_id)%>%group_by(gene_name)%>%
	mutate(new_gene_name=paste0(gene_name[1],c('',paste0('_',2:n())))[1:n()])%>%
	ungroup%>%
	select(gene_name=new_gene_name,gene_id)%>%
	left_join(cds%>%mcols%>%as.data.frame%>%distinct(transcript_id,gene_id,protein_id))

stopifnot(ids_nrgname%>%group_by(gene_id)%>%filter(n_distinct(gene_name)>1)%>%nrow%>%identical(0L))
cdsgrl = cds%>%split(.,.$transcript_id)


################################################################################
########Load dpprime data, write as salmon files
################################################################################


trid2gid = cds%>%mcols%>%as.data.frame%>%select(transcript_id,gene_id)%>%{safe_hashmap(.[[1]],.[[2]])}
gid2gnm = ids_nrgname%>%distinct(gene_id,gene_name)%>%{safe_hashmap(.[[2]],.[[1]])}
gnm2gid = ids_nrgname%>%distinct(gene_id,gene_name)%>%{safe_hashmap(.[[1]],.[[2]])}

dpfiles = Sys.glob(here('pipeline/deepshapeprime/*/run*'))%>%
	data.frame(file=.)%>%
	mutate(sample=basename(dirname(file)))%>%
	mutate(name=basename(file))%>%
	mutate(number=str_extract(name,'\\d+'))%>%
	arrange(desc(as.numeric(number)))%>%
	group_by(sample)%>%slice(1)%>%{setNames(.$file,.$sample)}

dpexprdata = dpfiles%>%map_df(.id='sample',fread)
dpexprdata%<>%set_colnames(c('sample','transcript_id','TPM','diff','count'))

dpexprdata%<>%mutate(oldtrid=transcript_id,transcript_id=str_replace(transcript_id,'\\.\\d+',''),gene_id = trid2gid[[transcript_id]])
notrimtr2g = dpexprdata%>%distinct(oldtrid,gene_id)
dpexprdata_s = dpexprdata%>%split(.,.$sample)%>%head(1)
dpexprdata%>%group_by(sample)%>%summarise(sum(count))

dpexprdata%<>%safe_left_join(cdsgrl%>%width%>%sum%>%enframe('transcript_id','length'))

library(tximport)

dpexprdata%>%split(.,.$sample)%>%imap(function(dpexprdata_s,sample){
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
	#now import	
	# tximport(files=fakesalmonfile,type='salmon',tx2gene=dpexprdata_s%>%select(transcript_id,gene_id))
})

################################################################################
########Import it all to get tr length scaled counts
################################################################################
	
library(tximport)
library(tidyverse)

allquantfiles = c(Sys.glob(here('pipeline/deepshapeprime/fakesalmonfiles/*ribo*/quant.sf')),Sys.glob(here('pipeline/salmon/data/*/quant.sf')))
names(allquantfiles) <- allquantfiles%>%dirname%>%basename
allquantfiles%>%head(2)%>%map(readLines,3)

dptrs = allquantfiles[[1]]%>%fread%>%.$Name
salmontrs = allquantfiles[[11]]%>%fread%>%.$Name%>%str_extract('ENSMUST\\w+')
# inclusiontable(dptrs,salmontrs)
trs = intersect(dptrs,salmontrs)

tx_countdata = tximport(files=allquantfiles,
	tx2gene=dpexprdata%>%distinct(transcript_id,gene_id),
	type='salmon',
	countsFromAbundance='scaledTPM',
	importer=function(file){
		read_tsv(file)%>%
			mutate(Name=str_extract(Name,'ENSMUST\\w+'))%>%
			filter(Name%in%trs)%>%arrange(match(Name,trs))
})

randomround = function(x)floor(x)+rbinom(length(x),1,x%%1)


tx_countdata$counts%>%as.data.frame%>%rownames_to_column('gene_id')%>%
	mutate(gene_name = gid2gnm[[gene_id]])%>%
	select(-gene_id)%>%select(gene_name,everything())%>%
	mutate_at(vars(-gene_name),list(randomround))%>%
	write_tsv('data/tx_scaled_countData.tsv')

ishighcount = tx_countdata$counts%>%apply(1,max)%>%`>`(32)
highcountgenes = rownames(tx_countdata$counts)[ishighcount]
highcountgnms = gid2gnm[[highcountgenes]]


################################################################################
########Save as exprset object
################################################################################

allcountmat <- tx_countdata$counts
allcountdesign = colnames(allcountmat)%>%data.frame(sample=.)%>%separate(sample,into=c('time','assay','rep'),remove=F)
allcountdesign = allcountdesign%>%arrange(assay=='ribo')%>%mutate(assay=as_factor(assay))%>%as.data.frame%>%set_rownames(colnames(allcountmat))

{
conflict_prefer("rowMedians", "matrixStats")
conflict_prefer("which", "BiocGenerics")

featuredata = data.frame(
	gene_id = allcountmat%>%rownames,
	gene_name = allcountmat%>%rownames%>%{gid2gnm[[.]]}
)%>%set_rownames(rownames(allcountmat))

countexprdata <- ExpressionSet(
	allcountmat,
	AnnotatedDataFrame(allcountdesign),
	AnnotatedDataFrame(featuredata)
)

countexprdata%>%saveRDS(here('pipeline/exprdata/countexprset.rds'))

}




################################################################################
########Run Limma (different variances)
################################################################################
tx_countdata_highcounts = tx_countdata%>%map_if(is.matrix,~ .[highcountgenes,])

mscountonlyvoom <- limma::voom(tx_countdata)

library(limma)
library(edgeR)
allcountdge <- DGEList(tx_countdata_highcounts$counts)
allcountkeep = filterByExpr(allcountdge)
allcountdge = allcountdge[allcountkeep,]
allcountdge <- calcNormFactors(allcountdge)

ribosamples = str_subset(colnames(allcountdge),'ribo')
rnasamples = str_subset(colnames(allcountdge),'total')

ribovoom = voom(allcountdge[,ribosamples],design=model.matrix(~time,allcountdesign[ribosamples,]))
rnavoom = voom(allcountdge[,rnasamples],design=model.matrix(~time,allcountdesign[rnasamples,]))
allvoom = voom(allcountdge[,c(rnasamples,ribosamples)],design=model.matrix(~assay*time,allcountdesign[c(rnasamples,ribosamples),]))

stopifnot(colnames(allvoom$E)==c(rnasamples,ribosamples))
allvoom$weights <- cbind(ribovoom$weights,rnavoom$weights)
allcountebayes <- eBayes(lmFit(allvoom))

timecoefs = allcountebayes$coef%>%colnames%>%str_subset('^time')
tetimecoefs = allcountebayes$coef%>%colnames%>%str_subset('assayribo:time')
allcontrasts = c(timecoefs,tetimecoefs)%>%setNames(.,.)

alllimmares = 
	map_df(.id='contrast',allcontrasts,function(contrasti){
		topTable(contrasts.fit(allcountebayes,coef=c(contrasti)),n=Inf)
	})

#unmangle the gene ids
alllimmares%<>%rownames_to_column('gene_id')%>%mutate(gene_id=str_replace(gene_id,'\\..*',''))




#next
# /fast/work/groups/ag_ohler/dharnet_m/cortexomics/src/R/TE_change/run_xtail.R
# 