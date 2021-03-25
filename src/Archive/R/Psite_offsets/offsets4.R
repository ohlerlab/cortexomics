base::source('src/R/Rprofile.R')
library(GenomicFiles)
library(GenomicFeatures)
argv <- c(
	bam = here('pipeline/star/data/E175_ribo_1/E175_ribo_1.bam')%T>%{stopifnot(file.exists(.))},
	gtf = here('pipeline/my_gencode.vM12.annotation.gtf'),
	REF = here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa'),
	outfolder = here('pipeline/seqshift_reads/data/E175_ribo_1/')
)
argv[] <- commandArgs(trailing=TRUE)
for (nm in names(argv)) assign(nm,argv[[nm]])
dropemptyseqs = .%>%keepSeqlevels(.,unique(seqnames(.)))
fmcols <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)[start(grl@partitioning)]
}
fmcols_List <- function(grl,...){
  with(grl@unlistData@elementMetadata,...)%>%split(grl@partitioning)
}
windsize = 50
#load annotation choose bam
library(rtracklayer)
library(GenomicFiles)
library(GenomicAlignments)

metainfo = read_tsv('data/metainfo.tsv')
bestcds = metainfo%>%filter(isbest)%>%.$protein_id%>%cdsspl[.]
bam = here('pipeline/star/data/E175_ribo_1/E175_ribo_1.bam')%T>%{stopifnot(file.exists(.))}
#get exons, cds
if(!exists('gtf_gr')) gtf_gr<-mimport(con=gtf,format='gtf')
exons <- gtf_gr%>%subset(type=='exon')%>%split(.,.$transcript_id)
cds <- gtf_gr%>%subset(type=='CDS')
cdsspl <- cds%>%split(.,.$protein_id)
startcods <- gtf_gr%>%subset(type=='start_codon')
exonsexp <- exons[fmcols(bestcds,transcript_id)]%>%
	resize_grl(sum(width(.))+windsize,'end',check=FALSE)%>%
	resize_grl(sum(width(.))+windsize,'start',check=FALSE)%>%
	{.[!any(is_out_of_bounds(.))]}
#windows for metaplots
metaplotwinds <- pmapToTranscripts (bestcds%>%resize_grl(1),exonsexp[fmcols(bestcds,transcript_id)])%>%
  unlist%>%
  resize(windsize+1,fix='end',ignore.strand=TRUE)%>%
  resize(windsize+windsize,'start',ignore.strand=TRUE)%>%
  {spl_mapFromTranscripts(.,exons_grl = exonsexp[seqnames(.)])}  
metaplotwinds%<>%split(.,.$xHits)

#get reads that overlap these
if(!exists('reads_trbak'))reads_trbak <- bam%>%BamFile(yieldSize=NA)%>%
	readGAlignments(param=ScanBamParam(which=c(metaplotwinds)%>%unlist))
reads_tr<-reads_trbak
reads_tr <- GRanges(reads_tr,length=qwidth(reads_tr))
reads_tr%<>%subsetByOverlaps(metaplotwinds)

disttable = reads_tr%>%
	mapToTranscripts(metaplotwinds)%>%
	resize(1,ignore.strand=TRUE)%>%
	.[.$xHits%>%duplicated%>%`!`]%>%
	{tibble(dist = windsize+1-start(.),readn=.$xHits)}
reads_tr$dist[disttable$readn]= disttable$dist
nf2n <- .%>%as.character%>%as.numeric
commonlengths = c(25:31)

reads_tr%<>%subset(length %in% commonlengths)
reads_tr$phase = ((3 - reads_tr$dist)%%3)
reads_tr$coddist = ((reads_tr$dist) + (reads_tr$phase))
reads_tr%<>%subset(dist >=6)
reads_tr%<>%subset(dist <= nf2n(length)-6)

reads_tr$dist%<>%as.factor
reads_tr$length%<>%as.factor
reads_tr$phase%<>%as.factor
reads_tr$coddist%<>%as.factor
reads_tr%<>%subset(!is.na(dist))

reads_tr$coddist%>%table%>%sort

get_read_flankseq <- function(trainreads,seq,nbp=2,trim=TRUE){

	stopifnot('length' %in% colnames(mcols(trainreads)))
	stopifnot(all(seqnames(trainreads)%in%seqinfo(seq)@seqnames))


	fp<-trainreads[TRUE]%>%resize(nbp,'start')%>%resize(nbp*2,'end')
	fp_inbounds <- !is_out_of_bounds(fp,seqinfo(seq))

	tp <- trainreads%>%resize(nbp,'end')%>%resize(nbp*2,'start')
	tp_inbounds <- !is_out_of_bounds(tp,seqinfo(seq))
	
	inbounds <- fp_inbounds & tp_inbounds
	inbounds%>%table
	
	if(trim){
		fp<-fp[inbounds]
		tp<-tp[inbounds]
	}else{
		assert_that(all(fp_inbounds))
		assert_that(all(tp_inbounds))
	}
	startseq <- getSeq(seq,fp)
	endseq <- getSeq(seq,tp)


	seqmat <- cbind( dnaseq2onehot(startseq,'fp.'),dnaseq2onehot(endseq,'tp.'))
	seqmat %>% as.data.frame%>%mutate_all(as.factor)
}

#remove then add seqcols
seqcols = paste0(colnames(mcols(reads_tr))%>%str_subset('^(fp|tp)\\.\\d'))
mcols(reads_tr)%<>%.[,colnames(.)%>%setdiff(seqcols)]
mcols(reads_tr) %<>% cbind(get_read_flankseq(reads_tr,FaFile(REF)))

# #verify the dist and seq are corret
# teststarts = bestcds%>%resize_grl(1)
# readfps = reads_tr%>%resize(1)
# reads_tr[readfps%>%overlapsAny(teststarts)]
# startreads = reads_tr[start(reads_tr)%in%(start(teststarts)%>%unlist)]
# start(startreads)%in%unlist(start(teststarts))%>%head
# startreads%>%sample(10)%>%getSeq(x=FaFile(REF),.)
# startreads%>%get_read_flankseq(FaFile(REF))%>%apply(2,mean)
# startreads$dist%>%table%>%sort


####
seqcols = paste0(colnames(mcols(reads_tr))%>%str_subset('^(fp|tp)\\.\\d'))
rf_formula = as.formula(paste0('coddist ~ length + phase +', paste0(seqcols,collapse=' + ')))
seqshiftmodel_allvars <- ranger::ranger(formula= rf_formula ,
		data=mcols(reads_tr)%>%as.data.frame%>%head(10e3),
		probability=FALSE,num.threads=4
)
predict(seqshiftmodel_allvars,mcols(reads_tr)%>%as.data.frame%>%head(10e3))$prediction%>%table
reads_tr$codshift = predict(seqshiftmodel_allvars,mcols(reads_tr)%>%as.data.frame)$prediction%>%nf2n
reads_tr$seqshift = reads_tr$codshift - nf2n(reads_tr$phase)
reads_tr%>%mcols%>%as.data.frame%>%group_by(length)%>%summarise(seqshift%>%table%>%sort%>%tail(1)%>%names)
reads_tr%>%resize(1)%>%mapToTranscripts(metaplotwinds)%>%coverage%>%as.matrix%>%colSums%>%.[(windsize-12):(windsize+1+12)]%>%txtplot
reads_tr%>%strandshift(9)%>%resize(1)%>%mapToTranscripts(metaplotwinds)%>%coverage%>%as.matrix%>%colSums%>%.[(windsize-12):(windsize+1+12)]%>%txtplot
reads_tr%>%strandshift(.$seqshift)%>%resize(1)%>%mapToTranscripts(metaplotwinds)%>%coverage%>%as.matrix%>%colSums%>%.[(windsize-12):(windsize+1+12)]%>%txtplot



#now test our model
{
is3nt = bestcds%>%width%>%sum%>%mod(3)%>%`==`(0)

cds4test = bestcds[is3nt]%>%head(100)

bcdspsites = psites%>%subsetByOverlaps(cds4test)

lastex = cds4test%>%sort_grl_st%>%.[[1]]%>%tail(1)

if(!exists('bestcdsreadsbak'))bestcdsreadsbak  <- bam%>%BamFile(yieldSize=NA)%>%readGAlignments(param=ScanBamParam(which=cds4test%>%unlist))
bestcdsreads=bestcdsreadsbak

expcds = bestcds%>%sort_grl_st%>%resize_grl(.,sum(width(.))+42,'end')%>%unlist
ov = findOverlaps(resize(as(bestcdsreads,'GRanges'),1),expcds,select='arbitrary')
# ov = ov%>%as.data.table%>%group_by(queryHits)%>%sample_n(1)
ov=ov%>%as.data.table

isneg = (strand(bestcdsreads)=='-')%>%as.vector
phaseoffsets = (3-expcds$phase[ov[[1]]])%%3
mcols(bestcdsreads)$phase=NA
mcols(bestcdsreads)$phase[!isneg] = ((start(bestcdsreads)[!isneg] - start(expcds)[ov[[1]][!isneg]])+phaseoffsets[!isneg])%%3
mcols(bestcdsreads)$phase[isneg] = (end(expcds)[ov[[1]][isneg]] - end(bestcdsreads)[isneg] +phaseoffsets[isneg])%%3
stopifnot((mcols(bestcdsreads)$phase%>%is.na%>%mean)<0.05)
mcols(bestcdsreads)$length=qwidth(bestcdsreads)
bestcdsreads%<>%subset(!is.na(phase))
ov = ov[!is.na(ov[[1]]),]
ov = ov[mcols(bestcdsreads)$length%in%25:31,]
bestcdsreads%<>%subset(length%in%25:31)
invisible(bestcdsreads)
}

{
bestcdsreads%<>%as("GRanges")
mcols(mcols(bestcdsreads) %<>% cbind(get_read_flankseq(bestcdsreads,FaFile(REF))))
bestcdsreads$codshift = predict(seqshiftmodel_allvars,mcols(bestcdsreads)%>%as.data.frame)$prediction%>%nf2n
bestcdsreads$seqshift = bestcdsreads$codshift - nf2n(bestcdsreads$phase)
}

stop()
bestcdsreads%>%strandshift(9)%>%resize(1)%>%mapToTranscripts(metaplotwinds)%>%coverage%>%as.matrix%>%colSums%>%log1p%>%txtplot
bestcdsreads%>%strandshift(.$seqshift)%>%resize(1)%>%mapToTranscripts(bestcds)%>%coverage%>%.[GRanges(names(.),IRanges(1,15))]%>%as.matrix%>%colSums%>%matrix(byrow=TRUE,ncol=3)
bestcdsreads%>%resize(1)%>%mapToTranscripts(bestcds)%>%shift(bestcdsreads$seqshift[.$xHits])%>%coverage%>%.[GRanges(names(.),IRanges(1,15))]%>%as.matrix%>%colSums%>%matrix(byrow=TRUE,ncol=3)
bestcdsreads%>%strandshift(nf2n(-.$phase))%>%resize(1)%>%mapToTranscripts(bestcds)%>%coverage%>%.[GRanges(names(.),IRanges(1,42+30))]%>%as.matrix%>%colSums%>%matrix(byrow=TRUE,ncol=3)
bestcdsreads%>%strandshift(nf2n(-.$phase))%>%mapToTranscripts(expcds%>%split(.,.$proteinid))%>%resize(1,ignore.strand=TRUE)%>%coverage%>%.[GRanges(names(.),IRanges(1,42+30))]%>%as.matrix%>%colSums%>%matrix(byrow=TRUE,ncol=3)

bestcdsreads%>%resize(1)%>%mapToTranscripts(bestcds)%>%shift(bestcdsreads$seqshift[.$xHits])%>%coverage%>%.[GRanges(names(.),IRanges(1,42+30))]%>%as.matrix%>%colSums%>%matrix(byrow=TRUE,ncol=3)

################################################################################
########Now test with glm
################################################################################
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

psitecovseqadj = bestcdsreads%>%
	resize(1)%>%
	mapToTranscripts(bestcds)%>%
	dropemptyseqs%>%
	shift(bestcdsreads$seqshift[.$xHits])%>%
	coverage

psitecovseq = bestcdsreads%>%
	resize(1)%>%
	mapToTranscripts(bestcds)%>%
	dropemptyseqs%>%
	shift(9)%>%
	coverage


codposdf<-readRDS('data/codposdf.rds')	

codposdfspl=codposdf%>%split(.,.$protein_id)
stopcodons <-  (GENETIC_CODE=='*')%>%.[.]%>%names

sitedf = get_sitedf(psitecovseqadj,codposdfspl)
sitedfnoadj = get_sitedf(psitecov[names(psitecovseqadj)],codposdfspl)

txtplot(
	sitedfnoadj%>%group_by(codon)%>%summarise(mcount=sum(count))%>%.$mcount,
	sitedf%>%group_by(codon)%>%summarise(mcount=sum(count))%>%.$mcount
)

sitedf%>%group_by(a_codon)%>%summarise(mcount=sum(count))%>%filter(mcount>7.5e3)

library(MatrixModels)
library(MASS)
mtheta=0.32

glmfit = glm4(count ~ 0 + gene + codon +a_codon, data=sitedf%>%filter(phase==0),family=negative.binomial(mtheta),MXITER=400,doFit=T, sparse=T, verbose=T)
glmfitnoadj = glm4(count ~ 0 + gene + codon + a_codon + phase, data=sitedfnoadj,family=negative.binomial(mtheta),MXITER=400,doFit=T, sparse=T, verbose=T)


glmfit_p_codon <-  list(seqadj = glmfit,noseqadj = glmfitnoadj)%>%
	map_df(.id='sample', ~ coef(.)%>%enframe)%>%
	filter(name%>%str_detect('^codon'))%>%
	mutate(p_codon = name%>%str_replace('codon',''))%>%
	filter(translate(DNAStringSet(p_codon))!='*')%>%
	select(sample,p_codon,codon_dt_glm=value)

glmfit_p_codon%>%spread(sample,codon_dt_glm)%>%{txtplot(.[[2]],.[[3]])}

glmfit_p_codon%>%spread(sample,codon_dt_glm)%>%.[['seqadj']]%>%sd
glmfit_p_codon%>%spread(sample,codon_dt_glm)%>%.[['noseqadj']]%>%sd

glmfit_a_codon <-  list(seqadj = glmfit,noseqadj = glmfitnoadj)%>%
	map_df(.id='sample', ~ coef(.)%>%enframe)%>%
	filter(name%>%str_detect('a_codon'))%>%
	mutate(a_codon= name%>%str_replace('a_codon',''))%>%
	filter(translate(DNAStringSet(a_codon))!='*')%>%
	select(sample,a_codon,codon_dt_glm=value)

glmfit_a_codon%>%spread(sample,codon_dt_glm)%>%{txtplot(.[[2]],.[[3]])}

glmfit_a_codon%>%spread(sample,codon_dt_glm)%>%.[['seqadj']]%>%sd
glmfit_a_codon%>%spread(sample,codon_dt_glm)%>%.[['noseqadj']]%>%sd

##



#Still not QUITE right - this doesn't put everything fully in the right place.
#AAAALMSOT but where is the error coming from?
################################################################################
########Debug
################################################################################
bestcdsreads%>%strandshift(.$seqshift)%>%resize(1)%>%mapToTranscripts(cds4test[[1]])



mappedbcdsreads = bestcdsreads%>%strandshift(nf2n(-.$phase))%>%resize(1)%>%mapToTranscripts(metaplotwinds)%>%{.$mphase = start(.)%%3;.}

bestcdsreads[na.omit(ov[[1]])]%>%strandshift(nf2n(-.$phase))%>%

partimap = bestcdsreads%>%strandshift(nf2n(-.$phase))%>%resize(1)%>%pmapToTranscripts(.,bestcds[names(expcds)[ov[[1]]] ])


ufmap = bestcdsreads%>%strandshift(nf2n(-.$phase))%>%resize(1)%>%mapToTranscripts(.,bestcds)
#use only the matches to the transripts we chose above
ismatch = (names(expcds)[ov[[1]][ufmap$xHits]] == names(bestcds)[ufmap$transcriptsHits])

ismatch%>%table

ufmap = ufmap[ismatch%in%TRUE]
ufmap%>%keepSeqlevels(.,unique(seqnames(.)))%>%coverage%>%.[GRanges(names(.),IRanges(1,42+30))]%>%as.matrix%>%colSums%>%matrix(byrow=TRUE,ncol=3)


dropemptyseqs = keepSeqlevels(.,unique(seqnames(.)))
offreads = ufmap%>%subset(start==6)
bestcdsreads[offreads$xHits]
expcds[ov[[1]][offreads$xHits]]

(end(bestcdsreads[offreads$xHits]) - end(expcds[ov[[1]][offreads$xHits]]))%%3

ulbestcds = bestcds%>%unlist
texphasemap = ulbestcds%>%strandshift(3)%>%resize(1)%>%mapToTranscripts(.,bestcds)
# texphasemap = texphasemap%>%shift(-ulbestcds$phase[texphasemap$xHits])
texphasemap%>%coverage%>%.[GRanges(names(.),IRanges(1,42+30))]%>%as.matrix%>%colSums%>%matrix(byrow=TRUE,ncol=3)
texphasemap%<>%{.$mphase = (-(start(.)-1))%%3;.}

outofphase = (texphasemap$mphase != ulbestcds$phase[texphasemap$xHits])

texphasemap[outofphase]

upstream_dist_till(bestcdsreads,cds4test%>%resize_grl(1,'start')%>%unlist)%>%intersect(30)

bestcdsreads_rin = bestcdsreads[bestcdsreads%in%reads_tr]
bestcdsreads_rin[1]
reads_tr[match(bestcdsreads_rin[1],reads_tr)]
isnaov = ov[[1]]%>%is.na
testind = isnaov%>%which%>%head(1)

testovmerge = bestcdsreads[testind]%>%mergeByOverlaps(cds4test)
end(testovmerge[[1]])
end(testovmerge$cds4test)

bestcds['ENSMUSP00000114649']%>%width%>%sum%>%`%%`(3)

bestcds['ENSMUSP00000114649']%>%mapToTranscripts(unlist(.),.)%>%
	{.$mphase = -(start(.)-1)%%3;.}%>%
	{.$phase = bestcds[['ENSMUSP00000114649']][.$xHits]$phase;.}%>%
	{.$myphase = (start(.)-1)%%3;.}%>%
	{.$backphase = (3-.$phase)%%3;.}%>%
	{.$w = width(bestcds[['ENSMUSP00000114649']][.$xHits]);.}%>%
	{.$mw = width(bestcds[['ENSMUSP00000114649']][.$xHits])%%3;.}

texphasemap['ENSMUSP00000114649']
################################################################################
########
################################################################################
	
#%>%lapply(head,20)%>%as.matrix%>%colSums%>%log1p%>%txtplot


#ensure that our reads have been re-positioned as expected around the start codon
txt_metaplot(alltrainreads%>%apply_psite_offset('manual'),startsites)
txt_metaplot(alltrainreads%>%apply_psite_offset('rfoffset'),startsites)