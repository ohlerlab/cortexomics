################################################################################
################################################################################
base::source(here::here('src/Rprofile.R'))
if(!exists("cdsgrl")) {
	base::source("src/Figures/load_annotation.R")
}
if(!exists('here')) base::source(here::here('src/Rprofile.R'))
intersect <- BiocGenerics::intersect
if(!exists('codonprofiles')) load('data/codon_coverage.Rdata')
if((!exists('allcodsigmean_isomerge'))|(!'availability'%in%colnames(allcodsigmean_isomerge))) base::source(here('src/Figures/Figure2/3_tRNA_array_analysis.R'))
library(rlang)
if('occupancy' in colnames(codonoccs)) codonoccs%<>%rename('dwell_time':=occupancy)


N_ST_CODONS = 10

	#Now merge with the tRNA data
tRNA_occ_df<-trna_ab_df%>%
	left_join(codonoccs)
dttrna_df = tRNA_occ_df%>%filter(fraction=='total')
stopifnot(codondf%>%group_by(codon,time)%>%tally%>%.$n%>%`==`(1))


#get the te changing gnes
allxtail = Sys.glob('pipeline/xtail/*')%>%map_df(.id='time',fread)%>%group_by(gene_name)
allxtail$gene_id = gnm2gid[[ allxtail$gene_name]]
techangedf <- allxtail%>%group_by(gene_id,gene_name)%>%
  mutate(sig = (adj_p_value < 0.05)& (abs(log2fc)>log2(1.25)))%>%
  summarise(
    up = as.numeric(any(sig & (log2fc > 0))),
    down = as.numeric(any(sig & (log2fc < 0)))
  )
techangegenes = techangedf%>%filter(up==1|down==1)%>%.$gene_name
teupgenes = techangedf%>%filter(up==1)%>%.$gene_name
tedowngenes = techangedf%>%filter(down==1)%>%.$gene_name





codondf = codondf%>%
	group_by(time)%>%
	mutate(abundance=abundance-median(na.rm=T,abundance))%>%
	group_by(codon,time)%>%select(abundance,dwell_time)




iso_tx_countdata<-readRDS('data/iso_tx_countdata.rds')

totalcols = colnames(iso_tx_countdata[[1]])%>%str_subset('total')
ribocols = colnames(iso_tx_countdata[[1]])%>%str_subset('ribo')

isos_te = iso_tx_countdata$abundance[,ribocols] / iso_tx_countdata$abundance[,totalcols]

highribo = iso_tx_countdata$abundance[,ribocols]%>%rowMins%>%`>`(1)
hightot = iso_tx_countdata$abundance[,totalcols]%>%rowMins%>%`>`(1)
isos_te = isos_te[highribo & hightot,]

isos_techange = rowMeans(log2(isos_te[,7:8])) - rowMeans(log2(isos_te[,1:2]))

isos_techange %<>%enframe('tr_id','techange')

cdsgrl %<>% sort_grl_st
cdsstartseq = extractTranscriptSeqs(x=fafile,cdsgrl%>%resize_grl(N_ST_CODONS*3,'start'))

codondftall = map_df(.id='pos',c(1:N_ST_CODONS)%>%setNames(.,.),function(x){
	cdsstartseq%>%substr((x-1)*3+1,x*3)%>%enframe('tr_id','codon')
	# c((x-1)*3+1,x*3)
})%>%
mutate(pos = paste0('codon_',pos))
codondf = codondftall%>%spread(pos,codon)




#First jsut try to model the te with these codons
leavoutcod=1

cods2use = 2:N_ST_CODONS
cods2use%<>%setdiff(leavoutcod)

codonform = paste0('2^techange ~ ',paste0(collapse='+','codon_',cods2use))

codonglm = codondf%>%left_join(isos_techange)%>%glm(data=.,formula=as.formula(codonform))
codonglm%>%aov

#do codons have similiar effects in similiar places?
codonglm$coef%>%enframe('codeff','value')%>%mutate(pos = codeff%>%str_extract('\\d+'),codon=codeff%>%str_extract('...$'))%>%
	mutate(cododn = rank(codon))%>%
	{txtplot(.$cododn,.$value)}

#change between E13 and E175


#now calcuate a 'slewness score'
tr_ramp_df = codondftall%>%
	mutate(pos = pos%>%str_extract('\\d+'))%>%
	filter(pos%in%2:5)%>%
	left_join(dttrna_df)%>%
	group_by(time,tr_id)%>%
	summarise(slowness=log2(mean(2^abundance)),startocc=mean(dwell_time))

tr_ramp_df_alltp = tr_ramp_df%>%
	group_by(tr_id)%>%
	summarise(slowness=log2(mean(2^slowness)),startocc=mean(startocc))


tr_ramp_df_alltp%>%left_join(isos_techange)%>%{quicktest(.$slowness,.$techange)}


techangetrs = ids_nrgname%>%filter(gene_name%in%tedowngenes)%>%.$transcript_id

tr_ramp_df_alltp%>%left_join(
	isos_techange%>%filter(tr_id%in%techangetrs)
)%>%{quicktest(.$slowness,.$techange)}

	#.[,c('tr_id',sort(as.numeric(colnames(.))[-1]))]


#test with proper te calc
#test with TM proteins sep.