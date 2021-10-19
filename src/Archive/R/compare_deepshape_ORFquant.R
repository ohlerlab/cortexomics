################################################################################
########This script tries to figure out why ORFquant doesn't identify a lot of the
########Genes that have MS in my data
################################################################################
source('src/Rprofile.R')
#get our expr data and metadata
load('data/integrate_exprdata2.Rdata')
#id conversions
trid2id = metainfo%>%
	distinct(transcript_id,gene_id)%>%
	{safe_hashmap(.[['transcript_id']],.[['gene_id']])}

trimids = function(str) str_replace(str,'\\.\\d+','')

#load various runs of ORFquant

#load deepshape
#src/R/load_deepshape_cov.R

#load the deepshape results
deepshapedf = fread(here('pipeline/deepshape/DeepShapePrimeOutputs/runlog_199.txt'))%>%head(-1)%>%
	select(transcript_id=1,TPM=2)%>%
	mutate(transcript_id=str_replace(transcript_id,'\\.\\d+',''))
deepshapegenes = deepshapedf%>%filter(TPM>1)%>%.$transcript_id%>%trid2id[[.]]

#load ribotaper
ribotaperdf = fread('pipeline/ribotaper/orfs_found')
ribotapergnms = ribotaperdf$gene_id

genes_w_mass_spec = metainfo%>%filter(has_ms)%>%.$gene_id

inclusiontable(genes_w_mass_spec,deepshapegenes)

inclusiontable(genes_w_mass_spec,ribotapergnms)


