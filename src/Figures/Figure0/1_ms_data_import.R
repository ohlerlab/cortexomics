{
library(Biostrings)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicFeatures)
mspeptidefile = here('ext_data/MS_Data_New/Ages_Brain_PEP_summ.txt')
msgenefile = here('ext_data/MS_Data_New/Ages_Brain_PG_summ.txt')



gid2gnm = ids_nrgname%>%distinct(gene_id,gene_name)%>%
	{safe_hashmap(.[[2]],.[[1]])}
peptidemsdata = fread(mspeptidefile)%>%
	dplyr::rename('gene_name':=Gene.names)
proteinmsdata = fread(msgenefile)%>%
	dplyr::rename('gene_name':=Gene.names)

cdsgrl = cdsgrl[str_order_grl(cdsgrl)]
proteinsequences = GenomicFeatures::extractTranscriptSeqs(cdsgrl,x=fafileob)
aaproteinsequences = translate(proteinsequences)
}

trlist = ids_nrgname%>%distinct(gene_name,transcript_id)%>%
	{split(.[[2]],.[[1]])}

igene='Satb2'

if(!file.exists(here('data/pepmatch.rds'))){
	pepmatch = mclapply(AAStringSet(peptidemsdata$Sequence),function(pattern){
		match = vmatchPattern(pattern,aaproteinsequences)
		ismatch = match %>%elementNROWS%>%`>`(0)%>%which
		names(match)[ismatch]
	})
	saveRDS(pepmatch,here('data/pepmatch.rds'))
}else{
	pepmatch<-readRDS(here('data/pepmatch.rds'))
}

gnamecompdf <- data.frame(name=1:nrow(peptidemsdata),
	msgname = peptidemsdata$gene_name)%>%
	left_join(pepmatch%>%enframe%>%unnest(value),allow_dups=TRUE)%>%
	mutate(gene_name = trid2gnm[[value]])


gnamecompdf$tr_id <- gnamecompdf$value

#so take these unambiguous matches and let's use them to update cases wher ethe gene name doesn't match
#our data
gnamematches = gnamecompdf%>%group_by(msgname)%>%distinct(msgname,gene_name)%>%
	filter(n()==1)
gnamematchdict = gnamematches%>%{setNames(.$gene_name,.$msgname)}

proteinmsdata%<>%mutate(gene_name = ifelse(gene_name%in%cds$gene_name,
	gene_name,
	gnamematchdict[gene_name]))
peptidemsdata%<>%mutate(gene_name = ifelse(gene_name%in%cds$gene_name,
	gene_name,
	gnamematchdict[gene_name]))

protpepdf <- proteinmsdata%>%{.$Peptide.IDs}%>%str_split(';')%>%
	setNames(proteinmsdata$gene_name)%>%IntegerList%>%stack
protpepdf$value <- (protpepdf$value %>%match(peptidemsdata$id))
protgnmtrids <- setNames(pepmatch[protpepdf$value],protpepdf$name)%>%
	enframe('gene_name','tr_id')%>%unnest(tr_id)
table(protgnmtrids$gene_name==trid2gnm[[protgnmtrids$tr_id]])
protgnmtrids <- protgnmtrids%>%semi_join(tibble(tr_id=fmcols(cdsgrl,transcript_id),gene_name=fmcols(cdsgrl,gene_name)))


trid2gnm <- ids_nrgname%>%select(transcript_id,gene_name)%>%{hashmap(.[[1]],.[[2]])}
proteinmsdata$gene_name%>%inclusiontable(trid2gnm[[alltrs]])

gnm2gid <- ids_nrgname%>%distinct(gene_name,gene_id)%>%{hashmap(.[[1]],.[[2]])}
proteinmsdata%<>%tibble
proteinmsdata$g_id <- gnm2gid[[proteinmsdata$gene_name]]


proteinmsdata%>%saveRDS('data/proteinmsdata.rds')
protgnmtrids%>%saveRDS('data/protgnmtrids.rds')

#as(names(cdsgrl),'CharacterList')[pepmatch[1:2]]


