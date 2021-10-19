#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")
source('src/Rprofile.R')
library(BSgenome.Scerevisiae.UCSC.sacCer3)
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")
ygenome = BSgenome.Scerevisiae.UCSC.sacCer3
writeXStringSet(getSeq(ygenome),filepath='yeast_test/saccer3.fa')

library(GenomicFeatures)
yanno = import('yeast_test/Yeast.sacCer3.sgdGene.gtf')
ydb = GenomicFeatures::makeTxDbFromGRanges(yanno)
trs = exonsBy((ydb),use.names=T)
extractTranscriptSeqs(x=ygenome,trs)%>%writeXStringSet('yeast_test/yeast_transcripts.fa')
bamfile = 'yeast_test/weinberg_yeast_riboAligned.out.bam'

# RiboseQC::prepare_annotation_files('riboseqcanno',gtf_file='yeast_test/Yeast.sacCer3.sgdGene.gtf',genome_seq='yeast_test/saccer3.fa')
# ranalysis=RiboseQC_analysis('/fast/work/groups/ag_ohler/dharnet_m/cortexomics/riboseqcanno/Yeast.sacCer3.sgdGene.gtf_Rannot',
# 	bamfile)

offsets=fread('yeast_test/weinberg_yeast_riboAligned.out.bam_P_sites_calcs')%>%rename('offset':=cutoff,'length'=read_length)
#load the yeast proteomics data #https://www.mcponline.org/content/15/4/1309/tab-figures-data
proteomicsdata = readxl::read_xlsx(here('ext_data/Lawless_etal_2016_supp1.xlsx'),sheet=1)%>%select(Protein,Name,`Final Quant`)

#load the scikit ribo data #https://www.sciencedirect.com/science/article/pii/S2405471217305495#bbib42

#I'll need to 
#I'll need to check if my glm approach gets similiar results to scikit ribo
#I"ll need to see if the spec if better than scikit ribo
#I'll need to run scikit ribo on this.
#For which I need salmon

bamfile = 'yeast_test/weinberg_yeast_riboAligned.out.sort.bam'
library(RiboseQC)

cdsgrl = yanno%>%subset(type=='CDS')%>%split(.,.$transcript_id)%>%sort_grl_st
exonsgrl = yanno%>%subset(type=='exon')%>%split(.,.$transcript_id)%>%sort_grl_st


# get_exp_cds <- function(cdsgrl,exonsgrl,expn){
# 	#create the expanded cds objects, first using exanded exon objects
# 	stopifnot(all(names(exonsgrl)==names(cdsgrl)))
# 	exonsexp <- exonsgrl%>%
# 		resize_grl(sum(width(.))+expn,'end',check=FALSE)%>%
# 		resize_grl(sum(width(.))+expn,'start',check=FALSE)%>%
# 		.[!any(is_out_of_bounds(.))]
# 	#expcds
# 	expcds <- pmapToTranscripts(
# 		cdsgrl,
# 		exonsexp
# 	)%>%
# 	unlist%>%
# 	resize(width(.)+expn,'end')%>%
# 	resize(width(.)+expn,'start')%>%
# 	{spl_mapFromTranscripts(.,exons_grl = exonsexp)}
# }


cdsexonsgrl <- exonsgrl%>%.[(fmcols(cdsgrl,transcript_id))]
# names(cdsexonsgrl)=names(cdsgrl)
# expcdsgrl <- get_exp_cds(cdsgrl,cdsexonsgrl,0)
# expcdsgrl <- expcdsgrl%>%split(.,names(.))
expcdsgrl <- cdsgrl
expcdsgrl <- expcdsgrl%>%dropSeqlevels('chrM',pruning='coarse')

get_cds_vects<-function(bamfile,cdsgrl,offsets,expn=0){
    #read psites from a bam file using our offsets
    psites = get_genomic_psites(bamfile,cdsgrl%>%unlist,offsets,mapqthresh=20,comps=c('chrM'='chrM'))
    #now convert to a cocverage vector
    psitecovrles = mapToTranscripts(psites,cdsgrl)%>%coverage
}


if(!file.exists('data/yeast_psitecovrles.rds')) {
	psitecovrles = 
		# Sys.glob('yeast_test/*ribo*sort*bam')%>%
		bamfile%>%
		setNames(.,basename(dirname(.)))%>%
		lapply(FUN=get_cds_vects,cdsgrl=expcdsgrl[],offsets=offsets,expn=0)
	psitecovrles%>%saveRDS('data/yeast_psitecovrles.rds')

}else{
	psitecovrles<-readRDS('data/yeast_psitecovrles.rds')	
}
densities = psitecovrles[[1]]%>%mean
tpms = 1e6*densities/sum(densities)

if(!file.exists(here('data/spec_vals.rds'))){
	spec_vals = psitecovrles[[1]]%>%lapply(.%>%as.vector%>%ftestvect)%>%map(1)
	saveRDS(spec_vals,here('data/spec_vals.rds'))
}else{
	spec_vals<-readRDS(here('data/spec_vals.rds'))
}


stop('')


if(!file.exists('data/yeast_codmat.rds')) {
	bestcdsseq = expcdsgrl%>%extractTranscriptSeqs(x=ygenome,.)
	codposdf = lapply(bestcdsseq,function(cdsseq){
		codonmat = codons(cdsseq)%>%{cbind(pos = .@ranges@start,as.data.frame(.))}%>%
			identity
	})
	codposdf%<>%bind_rows(.id='protein_id')
	codposdf%>%saveRDS('data/yeast_codposdf.rds')
}else{
	codposdf<-readRDS('data/yeast_codposdf.rds')	
}

psitemats = psitecovrles[[1]]%>%lapply(.%>%matrix(byrow=T,ncol=3))

codon_cov_df = cbind(codposdf,psitemats%>%map_df(as.data.frame))

codon_cov_df <- codon_cov_df%>%set_colnames(c('protein_id','pos','codon','readsf1','readsf2','readsf3'))

codon_cov_df%<>%mutate(reads = readsf1 + readsf2 + readsf3)

codon_cov_df%<>%group_by(protein_id)%>%mutate(a_codon=lead(codon))

sitedf = codon_cov_df%>%select(-reads)%>%pivot_longer(matches('readsf'),names_to='phase',values_to='count',names_prefix='readsf')

#Rnafold data
rnafoldfile = 'yeast_test/yeast_sacCer3_rnafold.txt'
rnafolddf = readLines(rnafoldfile)%>%
	str_split_fixed('\t',n=2)%>%
	{setNames(.[,2],.[,1])}%>%
	strsplit(' ')%>%map(as.numeric)%>%
	enframe('protein_id','rnafold')%>%
	unnest(rnafold)%>%
	group_by(protein_id)%>%
	mutate(pos=1:n())

sitedf%<>%left_join(rnafolddf)

stopcodons <-  (GENETIC_CODE=='*')%>%.[.]%>%names

sections = (1:length(psitecovrles[[1]])) %>% {split(sample(.),ceiling(./50))}
section=sections[[1]]

glmfits = mclapply(sections%>%sample(10),function(section){
	# 
	subsitedf = sitedf%>%filter(protein_id %in% names(expcdsgrl)[section])
	subsitedf%<>%rename('gene':=protein_id)
	#
	library(MASS)
	message('fit')
	codglmfit = glm.nb(data=subsitedf,formula= count ~ 0 + gene)
	message('done')
	codglmfit
}
)


require("MatrixModels")
require("Matrix")
require('MASS') 

mtheta = glmfits%>%map_dbl(~.$theta)%>%mean
message('final fit')
glmfit = glm4(count ~ 0 + gene + codon+a_codon+phase + rnafold,
	data=sitedf%>%rename('gene':=protein_id),
	family=negative.binomial(mtheta),MXITER=400,doFit=T, sparse=T, verbose=T)
glmfits = list(glmfit)%>%setNames(c('wb_ribo'))
message('done')

glmfit_p_codon <-  map_df(.id='sample',glmfits, ~ coef(.)%>%enframe)%>%
	filter(name%>%str_detect('^codon'))%>%
	mutate(p_codon = name%>%str_replace('codon',''))%>%
	filter(translate(DNAStringSet(p_codon))!='*')%>%
	select(sample,p_codon,codon_dt_glm=value)
glmfit_a_codon <-  coef(glmfit)%>%enframe%>%
	filter(name%>%str_detect('a_codon'))%>%
	mutate(a_codon = name%>%str_replace('a_codon',''))%>%
	filter(translate(DNAStringSet(a_codon))!='*')%>%
	select(a_codon,codon_dt_glm=value)

left_join(glmfit_a_codon,
	fread('ext_data/weinberg_etal_2016_S2.tsv')%>%rename('a_codon':=Codon)
)%>%{quicktest(.$codon_dt_glm,log2(.$`RiboDensity at A-site`))}

##Codon level correspondences
glmfitgeffs <-  coef(glmfit)%>%enframe%>%
	filter(name%>%str_detect('gene'))%>%
	mutate(gene = name%>%str_replace('gene',''))%>%
	select(transcript_id=gene,gene_flux_glm=value)

quantdf

if(!file.exists(here('data/spec_vals.rds'))){
	spec_vals = psitecovrles[[1]]%>%lapply(.%>%as.vector%>%ftestvect)%>%map(1)
	saveRDS(spec_vals,here('data/spec_vals.rds'))
}else{
	spec_vals<-readRDS(here('data/spec_vals.rds'))
}
spec_valssq = spec_vals%>%map_dbl(sqrt)
spec_valssq_pm = 1e6 * (spec_valssq/sum(spec_valssq,na.rm=TRUE))
densities = psitecovrles[[1]]%>%mean
tpms = 1e6*(densities/sum(densities))

quantdf = proteomicsdata%>%mutate(transcript_id=str_split(Protein,'_'))%>%unnest(transcript_id)%>%
	left_join(enframe(tpms,'transcript_id','tpm'))%>%
	left_join(enframe(spec_valssq_pm,'transcript_id','spec'))%>%
	left_join(glmfitgeffs)%>%
	left_join(fread('yeast_test/weinberg_yeast_rna_salm/quant.sf')%>%rename('transcript_id':=Name))%>%
	group_by(Protein)%>%summarise(gene_flux_glm=log(sum(exp(gene_flux_glm))),TPM=sum(TPM),tpm=sum(tpm),spec=sum(spec),`Final Quant`=`Final Quant`[1])

quantdf$prot = quantdf$`Final Quant`

wbsalmondf = fread('yeast_test/weinberg_yeast_rna_salm/quant.sf')%>%rename('gene':=Name)

quicktest(log(quantdf$tpm),quantdf$gene_flux_glm)

quicktest(log(quantdf$TPM),log(quantdf$prot))
quicktest(log(quantdf$tpm),log(quantdf$prot))

################################################################################
########Gene level quant
################################################################################

trimrlelist <- function(rlelist,n,rev=F){
	triminds = list(1:n)%>%map(multiply_by,-1)%>%rep(length(rlelist))%>%IntegerList
	revfun = if(rev)revElements else identity
	rlelist%>%
	revfun%>%
	{.[triminds]}%>%
	revfun
}

wbsalmondf = fread('yeast_test/weinberg_yeast_rna_salm/quant.sf')%>%rename('gene':=Name)
wbribosalmondf = fread('yeast_test/weinberg_yeast_ribo_salm/quant.sf')%>%rename('gene':=Name)
wbsktedf = fread('yeast_test/weinberg_yeast_ribo_skout/genesTE.csv')
wbskbetadf = fread('yeast_test/weinberg_yeast_ribo_skout/genesBetas.csv')

wbsalmondf <- wbsalmondf%>%inner_join(wbsktedf)%>%
	left_join(wbribosalmondf%>%select(gene,rTPM=TPM))%>%
	left_join(wbskbetadf)%>%
	left_join(enframe(spec_valssq_pm,'gene','spec'))%>%
	# mutate(log_sTPM = (log(TPM) - mean(log(TPM)))/sd(log(TPM)))%>%
	# mutate(skPA= log_sTPM + log(2^(log2_TE)))%>%
	mutate(skPA= 2^(log2(TPM) + log2_TE))%>%
	left_join(enframe(tpms,'gene','tpm'))

skquantdf<-proteomicsdata%>%inner_join(wbsalmondf,by=c("Protein"='gene'))%>%
	filter(rTPM>1)%>%
	group_by(Protein)%>%
	filter(n()==1)%>%
	summarise_at(vars(spec,skPA,tpm,TPM,rTPM,skPA,`Final Quant`,log2_TE),list(sum))
	# mutate(log2_TE = skPA - TPM)


#now plot
scattertitle = skquantdf%>%{cor.test(use='complete',(log2(.[['TPM']])),log2(.[['Final Quant']]))}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n${.$estimate}(${.$conf.low} - ${.$conf.high})')}
plotfile<- here(paste0('plots/','skriboTPM_vprot','.pdf'))
pdf(plotfile)
	ggplot(data=skquantdf,aes(x=log2(TPM),y=log2(`Final Quant`)))+
	geom_point()+
	scale_x_continuous(paste0('log2(RNAseq TPM)'))+
	scale_y_continuous(paste0('Lawless et al SRM proteomcs'))+
	ggtitle(paste0('Riboseq TPMs - Scikitribo data'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)

scattertitle = skquantdf%>%{cor.test(use='complete',log2(.[['rTPM']]),log2(.[['Final Quant']]))}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n${.$estimate}(${.$conf.low} - ${.$conf.high})')}
plotfile<- here(paste0('plots/','skribo_riboTPM_vprot','.pdf'))
pdf(plotfile)
	ggplot(data=skquantdf,aes(x=log2(rTPM),y=log2(`Final Quant`)))+
	geom_point()+
	scale_x_continuous(paste0('log2(Riboseq TPM)'))+
	scale_y_continuous(paste0('Lawless et al SRM proteomcs'))+
	ggtitle(paste0('Riboseq TPMs - Scikitribo data'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)

scattertitle = skquantdf%>%{cor.test(use='complete',(log2(.[['skPA']])),log2(.[['Final Quant']]))}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n${.$estimate}(${.$conf.low} - ${.$conf.high})')}
plotfile<- here(paste0('plots/','skribo_PA_vprot','.pdf'))
pdf(plotfile)
	ggplot(data=skquantdf,aes(x=log2(skPA),y=log2(`Final Quant`)))+
	geom_point()+
	scale_x_continuous(paste0('log2(SkRibo TE) + log2(RnaseqTPM)'))+
	scale_y_continuous(paste0('Lawless et al SRM proteomcs'))+
	ggtitle(paste0('Scikit_Ribo PAs - Scikitribo data'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)

scattertitle = skquantdf%>%{cor.test(use='complete',.[['log2_TE']],log2(.[['Final Quant']])) }%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n${.$estimate}(${.$conf.low} - ${.$conf.high})')}
plotfile<- here(paste0('plots/','skribo_TE_vprot','.pdf'))
pdf(plotfile)
	ggplot(data=skquantdf,aes(x=log2_TE,y=log2(`Final Quant`)))+
	geom_point()+
	scale_x_continuous(paste0('log2(SkRibo TE)'))+
	scale_y_continuous(paste0('Lawless et al SRM proteomcs'))+
	ggtitle(paste0('Scikit Ribo TEs - Scikitribo data'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)

#
scattertitle = skquantdf%>%{cor.test(use='complete',(log2(.[['spec']])),log2(.[['Final Quant']]))}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n${.$estimate}(${.$conf.low} - ${.$conf.high})')}
plotfile<- here(paste0('plots/','yeast_spec_vprot','.pdf'))
pdf(plotfile)
	ggplot(data=skquantdf,aes(x=log2(spec),y=log2(`Final Quant`)))+
	geom_point()+
	scale_x_continuous(paste0('log2(sqrt(spectral coef(riboseq)))'))+
	scale_y_continuous(paste0('Lawless et al SRM proteomcs'))+
	ggtitle(paste0('Scikit_Ribo PAs - Scikitribo data'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)




proteomicsdata%>%inner_join(wbsalmondf,by=c("Protein"='gene'))%>%
	filter(rTPM>1)%>%
	{quicktest(log2(.$TPM),log2(.$`Final Quant`))}

proteomicsdata%>%inner_join(wbsalmondf,by=c("Protein"='gene'))%>%
	filter(rTPM>1)%>%
	{quicktest(log2(.$`Final Quant`),log2(.$rTPM))}

proteomicsdata%>%inner_join(wbsalmondf,by=c("Protein"='gene'))%>%
	filter(rTPM>1)%>%
	mutate(skPA = TPM + log2_TE)%>%
	{quicktest(log2(.$`Final Quant`),log2(.$skPA))}

proteomicsdata%>%inner_join(wbsalmondf,by=c("Protein"='gene'))%>%
	{quicktest(log2(.$`Final Quant`),.$beta)}

proteomicsdata%>%inner_join(wbsalmondf,by=c("Protein"='gene'))%>%
	{quicktest(log2(.$`Final Quant`),.$logTPM_scaled_x)}


proteomicsdata%>%inner_join(wbsalmondf,by=c("Protein"='gene'))%>%
	{quicktest(log2(.$`tpm`),.$logTPM_scaled_x)}


proteomicsdata%>%inner_join(wbsalmondf,by=c("Protein"='gene'))%>%
	filter(rTPM>1)%>%
	{quicktest(log2(.$`Final Quant`),.$skPA)}

proteomicsdata%>%inner_join(wbsalmondf,by=c("Protein"='gene'))%>%
	{quicktest(log2(.$`Final Quant`),.$log2_TE)}

proteomicsdata%>%inner_join(wbsalmondf,by=c("Protein"='gene'))%>%{lm(data=.,log(`Final Quant`) ~ log(TPM) + log2_TE)}



proteomicsdata%>%inner_join(wbsalmondf,by=c("Protein"='gene'))%>%
	{quicktest(log2(.$`log2_TE`),log2(.$TPM))}

fread('ext_data/weinberg_etal_2016_S2.tsv')%>%colnames

fread('ext_data/weinberg_etal_2016_S2.tsv')%>%rename('codon':=Codon)%>%left_join(
	fread('yeast_test/weinberg_yeast_ribo_skout/codons.csv')
)%>%
	{quicktest(log2(.$`RiboDensity at A-site`),log2(.$codon_dwell_time))}



################################################################################
########Plot output of dermkit ribo
################################################################################
	

#A site correspondance
left_join(glmfit_a_codon%>%rename('codon'=a_codon),
	fread('yeast_test/weinberg_yeast_ribo_skout/codons.csv')%>%rename('codon':='codon')
)%>%{quicktest(.$codon_dt_glm,.$`codon_dwell_time`)}

#A site correspondance with original weinberg
codondensitydf <- glmfit_a_codon%>%rename('codon'=a_codon)%>%
	left_join(
		fread('ext_data/weinberg_etal_2016_S2.tsv')%>%rename('codon':=Codon)
	)%>%left_join(
		fread('yeast_test/weinberg_yeast_ribo_skout/codons.csv')%>%rename('codon':='codon')
	)

#now plot
scattertitle = codondensitydf%>%{cor.test(use='complete',(.[['codon_dwell_time']]),.[['RiboDensity at A-site']])}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
{str_interp('pearsons rho =\n${.$estimate}(${.$conf.low} - ${.$conf.high})')}
plotfile<- here(paste0('plots/','skr_cod_vs_wbRibodens','.pdf'))
pdf(plotfile)
	ggplot(data=codondensitydf,aes(x=codon_dwell_time,y=`RiboDensity at A-site`))+
	geom_point()+
	scale_x_continuous(paste0('codon dwell time - scikit-ribo'))+
	scale_y_continuous(paste0('codon dwell time - weinberg et al'))+
	ggtitle(paste0('Scikit codon DT vs Weinberg'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)

#now plot
scattertitle = codondensitydf%>%{cor.test(use='complete',(.[['codon_dt_glm']]),.[['RiboDensity at A-site']])}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
{str_interp('pearsons rho =\n${.$estimate}(${.$conf.low} - ${.$conf.high})')}
plotfile<- here(paste0('plots/','glm_cod_vs_wbRibodens','.pdf'))
pdf(plotfile)
	ggplot(data=codondensitydf,aes(x=codon_dt_glm,y=`RiboDensity at A-site`))+
	geom_point()+
	scale_x_continuous(paste0('codon dwell time - GLM'))+
	scale_y_continuous(paste0('codon dwell time - weinberg et al'))+
	ggtitle(paste0('GLM codon DT vs Weinberg'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)

#now plot
scattertitle = codondensitydf%>%{cor.test(use='complete',(.[['codon_dt_glm']]),.[['codon_dwell_time']])}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
{str_interp('pearsons rho =\n${.$estimate}(${.$conf.low} - ${.$conf.high})')}
plotfile<- here(paste0('plots/','glm_cod_vs_skrDT','.pdf'))
pdf(plotfile)
	ggplot(data=codondensitydf,aes(x=codon_dt_glm,y=`codon_dwell_time`))+
	geom_point()+
	scale_x_continuous(paste0('codon dwell time - GLM'))+
	scale_y_continuous(paste0('codon dwell time - ScikitRibo'))+
	ggtitle(paste0('GLM codon DT vs Scikit'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)


#now plot
scattertitle = codondensitydf%>%{cor.test(use='complete',(.[['codon_dt_glm']]),.[['RiboDensity at A-site']])}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n ${.$conf.low} - ${.$conf.high}')}
plotfile<- here(paste0('plots/','dkr_cod_vs_wbRibodens','.pdf'))
pdf(plotfile)
	ggplot(data=codondensitydf,aes(x=codon_dt_glm,y=`RiboDensity at A-site`))+
	geom_point()+
	scale_x_continuous(paste0('codon density - GLM'))+
	scale_y_continuous(paste0('codon density - weinberg et al'))+
	ggtitle(paste0('GLM gene term vs Riboseq Density'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)

#now plot
scattertitle = codondensitydf%>%{cor.test(use='complete',(.[['codon_dt_glm']]),.[['codon_dwell_time']])}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n ${.$conf.low} - ${.$conf.high}')}
plotfile<- here(paste0('plots/','dkr_cod_vs_skr_cod','.pdf'))
pdf(plotfile)
	ggplot(data=codondensitydf,aes(x=codon_dt_glm,y=`codon_dwell_time`))+
	geom_point()+
	scale_x_continuous(paste0('codon density - GLM'))+
	scale_y_continuous(paste0('codon density - weinberg et al'))+
	ggtitle(paste0('GLM gene term vs Riboseq Density'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)

#now plot
scattertitle = quantdf%>%filter(!tpm==0)%>%{cor.test(use='complete',log(.[['tpm']]),.[['gene_flux_glm']])}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
	{str_interp('pearsons rho =\n ${.$conf.low} - ${.$conf.high}')}
plotfile<- here(paste0('plots/','dkr_flux_vs_ribodensity','.pdf'))
pdf(plotfile)
	ggplot(data=quantdf,aes(x=log(tpm),y=gene_flux_glm))+
	geom_point()+
	scale_x_continuous(paste0('RiboDensity'))+
	scale_y_continuous(paste0('DKR_flux'))+
	ggtitle(paste0('GLM gene term vs Riboseq Density'),sub=scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)


quicktest<-function(x,y){
	require(txtplot)
	complete = is.finite(x) & is.finite(y)
	message(str_interp('{sum(!complete)} (round(mean(!complete)*100,2) %) missing values'))
	txtplot(x[complete],y[complete])
	cor.test(x[complete],y[complete])
}

quicktest(quantdf$tpm,quantdf$spec)

quicktest(log(quantdf$tpm),log(quantdf$`Final Quant`))

quicktest(log(skquantdf$spec),log(skquantdf$`Final Quant`))
quicktest(log(skquantdf$rTPM),log(skquantdf$`Final Quant`))



lowabquantdf <- quantdf%>%filter(tpm<100)
quicktest(log(lowabquantdf$tpm),log(lowabquantdf$`Final Quant`))
quicktest(log(lowabquantdf$spec),log(lowabquantdf$`Final Quant`))

qdfcomp = quantdf%>%filter(is.finite(tpm),is.finite(spec),is.finite(`Final Quant`))

qdfcomp%>%gather(predictor,value,-Protein,-`Final Quant`)%>%
	lm(data=.,`Final Quant` ~ value * predictor )%>%anova

psitecovrles[[1]][list(1:30)%>%multiply_by(-1)%>%rep(10)%>%IntegerList]%>%as.matrix%>%colSums%>%txtplot
rlelist<-psitecovrles[[1]]


