require("MatrixModels")
require("Matrix")
require('MASS') 
source('src/R/Rprofile.R')
################################################################################
########Get codon sequences and quickly see if dwell time calc works for me
################################################################################
psitecov <- psitecovrles[[1]]

is3nt = psitecov%>%runLength%>%sum%>%`%%`(3)%>%`==`(0)

codposdfspl=codposdf%>%split(.,.$protein_id)


stopcodons <-  (GENETIC_CODE=='*')%>%.[.]%>%names

cdsgrl <- cds%>%split(.,.$protein_id)
cdsexonsgrl <- exons%>%split(.$transcript_id)%>%.[(fmcols(cdsgrl,transcript_id))]
names(cdsexonsgrl)=names(cdsgrl)
# expcdsgrl <- get_exp_cds(cdsgrl,cdsexonsgrl,0)
# expcdsgrl <- expcdsgrl%>%split(.,names(.))
expcdsgrl<-cdsgrl
expcdsgrl <- expcdsgrl%>%dropSeqlevels('chrM',pruning='coarse')

if(!file.exists('data/psitecovrles.rds')) {	
	psitecovrles = Sys.glob('pipeline/star/data/*/*ribo*.bam')%>%
		setNames(.,basename(dirname(.)))%>%
		lapply(FUN=get_cds_vects,cdsgrl=expcdsgrl[names(bestcds)],offsets=offsets)

	psitecovrles%>%saveRDS('data/psitecovrles.rds')

}else{
	psitecovrles<-readRDS('data/psitecovrles.rds')	
}


if(!file.exists('data/codmat.rds')) {

	bestcdsseq = cdsgrl%>%sort_grl_st%>%extractTranscriptSeqs(x=FaFile(REF),.)
	vmatchPattern(bestcdsseq,pattern=DNAString('N'))%>%elementNROWS%>%is_in(0)%>%table

	codposdf = lapply(bestcdsseq,function(cdsseq){
		codonmat = codons(cdsseq)%>%{cbind(pos = .@ranges@start,as.data.frame(.))}%>%
			identity
	})
	codposdf%<>%bind_rows(.id='protein_id')
	codposdf%>%saveRDS('data/codposdf.rds')
}else{
	codposdf<-readRDS('data/codposdf.rds')	
}

stopifnot(all(names(psitecovrles[[1]]) %in%unique(codposdf$protein_id)))
testid=sample(names(psitecovrles[[1]]),1)
codposn = codposdf%>%subset(protein_id==testid)%>%nrow
stopifnot(codposn == (length(psitecovrles[[1]][[testid]])/3))

get_sitedf<-function(testpsitecovs,codposdfspl,stop_codons=stopcodons){
		sitedf <- testpsitecovs%>%lapply(as.vector)%>%stack%>%set_colnames(c('count','gene'))
		
		testpsitecodons = codposdfspl%>%.[names(testpsitecovs)]%>%bind_rows

		sitedf%<>%group_by(gene)%>%mutate(phase=as_factor(((1:length(count))-1)%%3) )
		sitedf$codon=NA
		sitedf$codon[seq(1,nrow(sitedf),by=3)] = testpsitecodons$x
		sitedf$codon[seq(2,nrow(sitedf),by=3)] = testpsitecodons$x
		sitedf$codon[seq(3,nrow(sitedf),by=3)] = testpsitecodons$x
		acod = testpsitecodons%>%group_by(protein_id)%>%mutate(acod=lead(x))%>%.$acod
		sitedf$a_codon = NA
		sitedf$a_codon[seq(1,nrow(sitedf),by=3)] = acod
		sitedf$a_codon[seq(2,nrow(sitedf),by=3)] = acod
		sitedf$a_codon[seq(3,nrow(sitedf),by=3)] = acod
		sitedf <- sitedf%>%filter(!codon %in% stopcodons,!a_codon %in% stop_codons)
		sitedf

}


sections = (1:1000) %>% {split(sample(.),ceiling(./10))}
section=sections[[1]]

codoncounts = lapply(sections,function(section){
	# 
	tnt_entrop_cds <- pmes[is3nt]%>%order(decreasing=TRUE)%>%.[section]

	testpsitecovs = c(psitecov[is3nt][tnt_entrop_cds])

	sitedf = get_sitedf(testpsitecovs,codposdf%>%split(.,.$protein_id))

	sitedf%>%group_by(a_codon)%>%summarise(count=sum(count))

})

glmfits = lapply(sections,function(section){
	# 
	tnt_entrop_cds <- pmes[is3nt]%>%order(decreasing=TRUE)%>%.[section]

	testpsitecovs = c(psitecov[is3nt][tnt_entrop_cds])

	sitedf = get_sitedf(testpsitecovs,codposdf%>%split(.,.$protein_id))

	library(MASS)
	message('fit')
	codglmfit = glm.nb(data=sitedf,formula= count ~ 0 + gene + codon+a_codon+phase)
	message('done')
	codglmfit
})
mtheta = glmfits%>%map_dbl(~.$theta)%>%mean


#look at how consisten values are
glmfits%>%setNames(seq_along(.))%>%map_df(~.$coefficients%>%enframe)%>%
	filter(!name%>%str_detect('geneENS'))%>%
	group_by(name)%>%summarise(val=mean(value),stder=sd(value))

glmfits <- psitecovrles%>%map(.x=.,.f=function(psitecovsfit){
	#
	sitedf = get_sitedf(psitecovsfit,codposdf%>%split(.,.$protein_id)%>%.[names(psitecovsfit)])
	message('fit')
	glmfit = glm4(count ~ 0 + gene + codon+a_codon+phase, data=sitedf,family=negative.binomial(mtheta),MXITER=400,doFit=T, sparse=T, verbose=T)
	
})

codposdf%>%split(.,.$protein_id)%>%.[names(psitecovsfit)]%>%bind_rows%>%group_by(protein_id)%>%tally%>%left_join(enframe(sum(runLength(psitecovrles[[1]])/3),'protein_id','len'))%>%filter(n!=(len))
testpid = 'ENSMUSP00000020231'

psitecovrles[[1]][testpid]%>%runLengths%>%sum
codposdf%>%filter(protein_id==testpid)%>%nrow%>%multiply_by(3)

glmfit = glmfits[1]
quickcodsigs = sitedf%>%group_by(codon)%>%summarise(mcount=mean(count))

glmfitcodon%>%left_join(quickcodsigs)

glmfits%<>%setNames(bams%>%dirname%>%basename)

glmfit_p_codon <-  map_df(.id='sample',glmfits, ~ coef(.)%>%enframe)%>%
	filter(name%>%str_detect('^codon'))%>%
	mutate(p_codon = name%>%str_replace('codon',''))%>%
	filter(translate(DNAStringSet(p_codon))!='*')%>%
	select(sample,p_codon,codon_dt_glm=value)

#now plot
plotfile<- here(paste0('plots/','glm_codondtvar','.pdf'))
pdf(plotfile,w=5,h=5)
glmfit_p_codon%>%group_by(sample)%>%summarise(codsd = sd(codon_dt_glm))%>%.$codsd%>%{data.frame(v=.,x=seq_along(.))}%>%{qplot(data=.,y=v,x=x,geom='point',main='variance in codon dwell times - GLM model')}
dev.off()
normalizePath(plotfile)


################################################################################
########single fit
################################################################################
	



glmfit_p_codon <-  coef(glmfit)%>%enframe%>%
	filter(name%>%str_detect('^codon'))%>%
	mutate(p_codon = name%>%str_replace('codon',''))%>%
	filter(translate(DNAStringSet(p_codon))!='*')%>%
	select(p_codon,codon_dt_glm=value)

glmfit_a_codon <-  glmfits%>%map_df(.id='sample',.%>%coef(glmfit)%>%enframe)%>%
	filter(name%>%str_detect('a_codon'))%>%
	mutate(a_codon = name%>%str_replace('a_codon',''))%>%
	filter(translate(DNAStringSet(a_codon))!='*')%>%
	select(sample,a_codon,codon_dt_glm=value)


glmfit_a_codon%>%group_by(sample)%>%summarise(sd(codon_dt_glm))%>%


glmfit_p_codon <-  glmfits%>%map_df(.id='sample',.%>%coef(glmfit)%>%enframe)%>%
	filter(name%>%str_detect('^codon'))%>%
	mutate(p_codon = name%>%str_replace('codon',''))%>%
	filter(translate(DNAStringSet(p_codon))!='*')%>%
	select(sample,p_codon,codon_dt_glm=value)

glmfit_p_codon%>%group_by(sample)%>%summarise(sd(codon_dt_glm))



glmfitgeffs <-  coef(glmfit)%>%enframe%>%
	filter(name%>%str_detect('gene'))%>%
	mutate(gene = name%>%str_replace('gene',''))%>%
	select(protein_id=gene,gene_flux_glm=value)
	
glmfitgeffs%>%
	left_join(enframe(mean(psitecov),'protein_id','density'))%>%
	mutate(ldens = log(density))%>%
	{cor.test.finite(.$gene_flux_glm,.$ldens)}

################################################################################
########Do our gene fluxes work better than densities?
################################################################################

macountdf = allcountmat[,'E13_ribo_1']%>%enframe('protein_id','count')%>%
	filter(protein_id %in% pid2msid$keys())%>%
	mutate(ms_id = pid2msid[[protein_id]])%>%
	group_by(ms_id)%>%
	slice(which.max(count))%>%
	left_join(matchedms_mat[,'E13_MS_1']%>%enframe('ms_id','E13_MS_1'))%>%
	left_join(ibaqMatl2[,1]%>%enframe('ms_id','ibaq_E13_1'))%>%
	left_join(enframe(sum(width(cds%>%split(.,.$protein_id))),'protein_id','width'))%>%
	left_join(glmfitgeffs)
	



macountdf%>%filter(ms_id%in%ids2comp)%>%{cor.test.finite(.$ibaq_E13_1, log2(.$count) )}
macountdf%>%filter(ms_id%in%ids2comp)%>%{cor.test.finite(.$ibaq_E13_1, log2(.$count/.$width) )}
macountdf%>%filter(ms_id%in%ids2comp)%>%filter(gene_flux_glm> -6.5)%>%{cor.test.finite(.$ibaq_E13_1, .$gene_flux_glm )}

macountdf%>%filter(ms_id%in%ids2comp)%>%{cor.test.finite(log2(.$count/.$width), .$gene_flux_glm )}


################################################################################
########Do these correlate with tRNA etc?
################################################################################
	
codonprofiles <- readRDS('data/codonprofiles.rds')

glmfitcodon%>%left_join(codonoccs%>%filter(time=='E13'))%>%
	{cor.test.finite(.$codon_dt_glm,.$occupancy)}

allcodsigmean_isorepmerge <- allcodsigmean_isomerge%>%group_by(time,fraction,codon)%>%summarise_at(vars(abundance,balance,weightedusage),list(mean))

#Psite dt vs tRNA
glmfit_p_codon%>%rename('codon':=p_codon)%>%left_join(allcodsigmean_isorepmerge%>%filter(time=='E13',fraction=='Total'))%>%
	{cor.test.finite(.$codon_dt_glm,.$balance)}

glmfit_p_codon%>%rename('codon':=p_codon)%>%
	left_join(allcodsigmean_isorepmerge%>%filter(time=='E13',fraction=='Total'))%>%
	{cor.test.finite(.$codon_dt_glm,.$abundance)}

#Asite vs tRNA
glmfit_a_codon%>%rename('codon':=a_codon)%>%left_join(allcodsigmean_isorepmerge%>%filter(time=='E13',fraction=='Total'))%>%
	{cor.test.finite(.$codon_dt_glm,.$balance)}

glmfit_a_codon%>%rename('codon':=a_codon)%>%left_join(allcodsigmean_isorepmerge%>%filter(time=='E13',fraction=='Total'))%>%
	{cor.test.finite(.$codon_dt_glm,.$abundance)}


#Asite vs tRNA
glmfit_a_codon%>%rename('codon':=a_codon)%>%left_join(allcodsigmean_isorepmerge%>%filter(time=='E13',fraction=='80S'))%>%
	{cor.test.finite(.$codon_dt_glm,.$balance)}
glmfit_a_codon%>%rename('codon':=a_codon)%>%left_join(allcodsigmean_isorepmerge%>%filter(time=='E13',fraction=='80S'))%>%
	{cor.test.finite(.$codon_dt_glm,.$abundance)}


glmfitcodon%>%mutate(AA = DNAStringSet(codon)%>%translate%>%as.character)%>%
	group_by(AA)%>%filter(n()>1)%>%ungroup%>%
	mutate(.,cod_dt_glm_aares = lm(.$codon_dt_glm ~ 0+.$AA)$residuals)%>%
	left_join(allcodsigmean_isomerge%>%filter(time=='E13',fraction=='Total'))%>%
	{cor.test.finite(.$cod_dt_glm_aares,.$abundance)}

glmfitcodon%>%mutate(AA = DNAStringSet(codon)%>%translate%>%as.character)%>%
	group_by(AA)%>%filter(n()>1)%>%ungroup%>%
	mutate(.,cod_dt_glm_aares = lm(.$codon_dt_glm ~ 0+.$AA)$residuals)%>%
	left_join(allcodsigmean_isomerge%>%filter(time=='E13',fraction=='Total'))%>%
	{cor.test.finite(.$cod_dt_glm_aares,.$balance)}

codonoccs%>%left_join(allcodsigmean_isomerge%>%filter(time=='E13',fraction=='Total'))%>%
	{cor.test.finite(.$occupancy,.$abundance)}

sample = samples[bamind]

glmfits%>%map(~.$theta)

#now plot
plotfile<- here(paste0('plots/','tRNA_ab_v_GLM_dt','.pdf'))
pdf(plotfile)
scattertitle =glmfitcodon%>%left_join(allcodsigmean_isorepmerge)%>% {cor.test(.$codon_dt_glm,.$abundance)%>%tidy%>%mutate_if(is.numeric,round,5)%>%
	{str_interp('pearsons rho =\n ${.$conf.low} - ${.$conf.high};p =${.$p.value} ')}}
glmfitcodon%>%left_join(allcodsigmean_isorepmerge%>%filter(time=='E13',fraction=='Total'))%>%
	# filter(codon_dt_glm<0.5)%>%
	ggplot(.,aes(label=codon,x=codon_dt_glm,y=abundance))+
	geom_text()+
	# scale_color_discrete(name='colorname',colorvals)+
	scale_x_continuous(paste0(sample,' GLM calculated Codon Dwell Time'))+
	scale_y_continuous(paste0('tRNA Abundance'))+
	ggtitle(scattertitle)+
	theme_bw()
dev.off()
normalizePath(plotfile)


#now plot
plotfile<- here(paste0('plots/','tRNA_bal_v_GLM_dt','.pdf'))
pdf(plotfile)
glmfitcodon%>%left_join(allcodsigmean_isorepmerge%>%filter(time=='E13',fraction=='Total'))%>%
	# filter(codon_dt_glm<0.5)%>%
	ggplot(.,aes(label=codon,x=codon_dt_glm,y=balance))+
	geom_text()+
	# scale_color_discrete(name='colorname',colorvals)+
	scale_x_continuous(paste0(sample,' GLM calculated Codon Dwell Time'))+
	scale_y_continuous(paste0('tRNA Abundance'))+
	# ggtitle(paste0('title'))+
	theme_bw()
dev.off()
normalizePath(plotfile)
