library(data.table)

base::source(here::here('src/R/Rprofile.R'))
base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))
if(!exists('iso_tx_countdata')) load('data/1_integrate_countdata.R')
# base::source(here::here('src/Figures/Figure0/0_load_annotation.R'))

trlens = exons%>%split(.,.$transcript_id)%>%width%>%sum
cdslens = cds%>%split(.,.$transcript_id)%>%width%>%sum

# dpcovfiles = Sys.glob('pipeline/dpprimecov/*ribo*/dpprime_dist_ribosig.tsv')%>%setNames(.,basename(dirname(.)))
dpcovfiles = Sys.glob('pipeline/deepshapebamdata/*.bam.reformat.absolute')%>%str_subset('ribo|total')
dpcovfiles%<>%setNames(.,basename(.)%>%str_replace('.bam.reformat.absolute',''))
#dpcovfiles = dpcovfiles[1]

# dpcovfiles<-'pipeline/deepshape/tr_pos_ribosig'%>%setNames('E13_ribo_1')
dpcovfile=dpcovfiles[1]

if(!file.exists(here('data/dpcovdata.rds'))){
	dpcovdata <-lapply(dpcovfiles,function(dpcovfile){
		deepshapedata = fread(dpcovfile)
		# deepshapedata = deepshapedata%>%
		# 	mutate(V1=str_replace(V1,'\\.\\d+$',''))%>%
		# 	set_colnames(c('seqnames','start','score'))%>%
		# 	mutate(end=start)%>%
		# 	GRanges
		deepshapedata$V6%<>%trimids
		deepshapedata%<>% select(seqnames=V6,start=V5,end=V5)%>%GRanges
		deepshapedata$score=1		
		# seqnames(deepshapedata)%>%unique %>%is_in(names(cdslens))
		seqlengths(deepshapedata) = cdslens[names(seqlengths(deepshapedata))]
		deepshapecov = deepshapedata%>%coverage(weight='score')
		deepshapecov
	})
	saveRDS(dpcovdata,here('data/dpcovdata.rds'))
}else{
	dpcovdata<-readRDS(here('data/dpcovdata.rds'))
}
# dpcovdata[[1]]=deepshapecov
# names(dpcovdata)[[1]]='E13_ribo_1'
#ENSMUSP00000028342
stopifnot(dpcovdata[['E13_ribo_1']]['ENSMUST00000028342']%>%which.max%>%`==`(359))


dpcovdata[['E16_ribo_1']]


trid2gid = cds%>%mcols%>%as.data.frame%>%select(transcript_id,gene_id)%>%{safe_hashmap(.[[1]],.[[2]])}
tr2genemap=data.frame(transcript_id=trid2gid$keys(),gene_id=trid2gid$values())

ribocovtrs <- iso_tx_countdata$abundance%>%
	as.data.frame%>%
	rownames_to_column('transcript_id')%>%
	left_join(tr2genemap)%>%
	group_by(gene_id)%>%
	pivot_longer(matches('total|ribo'))%>%
	group_by(gene_id,transcript_id)%>%
	summarise(v=median(value))%>%
	filter(transcript_id%in%names(dpcovdata[[1]]))%>%
	slice(which.max(v))%>%
	filter(v>10)%>%
	.$transcript_id

data.frame(trid=ribocovtrs)%>%write_tsv('pipeline/scikitribotrs.txt',col_names=F)

system('grep -f pipeline/scikitribotrs.txt annotation/gencode.vM12.annotation.gtf > pipeline/scikitribo.gtf')

# deepshapecov <- dpcovdata[[1]]

# longdpcovs=dpcovdata[['E16_ribo_1.absolute']]%>%
# 	.[sum(runLength(.))>200]


# longdpcovs%>%
# 	head(10)%>%
# 	.[GRanges(names(.),IRanges(start=1,end=30))]%>%
# 	identity%>%
# 	as.matrix%>%
# 	colSums%>%log2%>%txtplot

# longdpcovs%>%
# 	head(400)%>%
# 	.[GRanges(names(.),IRanges(start=1,end=50))]%>%
# 	identity%>%
# 	as.matrix%>%
# 	colSums%>%matrix(nrow=3)%>%rowSums


# dpcovdata

# readaugsigdf = filt_reads_psite_list[[1]][,]%>%group_by(transcript,psite,cds_start)%>%tally%>%group_by(transcript)%>%mutate(n=n/sum(n))%>%
# 	mutate(pos = psite-cds_start)%>%
# 	group_by(pos)%>%
# 	summarise(signal=mean(n))




if(F){
	library(GenomicAlignments)
	bam = 'pipeline/star_transcript/data/E13_ribo_1/E13_ribo_1.sort.bam'
	mapqthresh=200
	param<-ScanBamParam()
	trreads = GenomicAlignments::readGAlignments(bam,param=param,use.names=TRUE,where=cdsgrl%>%unlist%>%subset(seqnames=='chr1'))
	any(duplicated(names(trreads)))
	seqlevels(trreads) <- seqlevels(trreads)%>%str_extract('ENSMUST\\w+\\.\\d+')
	seqlevels(trreads) <- trimids(seqlevels(trreads))
	sample = basename(dirname(bam))
	tpminds = match(seqnames(trreads),rownames(iso_tx_countdata$abundance))%>%as.numeric
	names(trreads)[duplicated(names(trreads))]
	tibble(grind = 1:length(trreads),read = as.vector(names(trreads)),TPM = iso_tx_countdata$abundance[tpminds,sample])%>%
		group_by(read)%>%mutate(w = TPM/sum(TPM))%>%sample_n(1,weight=TPM)



# dpcovdata%>%saveRDS('data/tpmweightdpcovdata.rds')
msvect = matched_ms_matrix[,'E13_MS_1']
#and our standard TPM vects
tpmvect = tx_countdata$abundance[,'E13_ribo_1']
rnatpmvect = tx_countdata$abundance[,'E13_total_1']

wdpcovdata <- readRDS('data/tpmweightdpcovdata.rds')
#weighted cov density
wdpcovdensvect = wdpcovdata[[1]]%>%mean%>%enframe('trid','dens')%>%mutate(gid=trid2gid[[trid]])%>%group_by(gid)%>%summarise(dens=sum(dens))
wdpcovdensvect = setNames(wdpcovdensvect$dens,wdpcovdensvect$gid)
#non weighted cov density
dpcovdensvect = dpcovdata[[1]]%>%mean%>%enframe('trid','dens')%>%mutate(gid=trid2gid[[trid]])%>%group_by(gid)%>%summarise(dens=sum(dens))
dpcovdensvect = setNames(dpcovdensvect$dens,dpcovdensvect$gid)
#also get cov density from our resampled bam data
shnames = intersect(names(msvect),names(dpcovdensvect))%>%intersect(names(tpmvect))%>%intersect(names(wdpcovdensvect))

quicktest(log2(msvect[shnames]),log2(dpcovdensvect[shnames]))
quicktest(log2(msvect[shnames]),log2(tpmvect[shnames]))
quicktest(log2(msvect[shnames]),log2(wdpcovdensvect[shnames]))
quicktest(log2(msvect[shnames]),log2(rnatpmvect[shnames]))


stop()

singletrtrs = tx2genemap%>%set_colnames(c('trid','gid'))%>%group_by(gid)%>%filter(n()==1)%>%.$trid
# {input.bamreformat}.readtpms
#okay so these don't match up, and should...
iso_tx_countdata$abundance[,'E13_ribo_1',drop=F]%>%
	as.data.frame%>%
	rownames_to_column('trid')%>%
	pivot_longer(-trid)%>%
	left_join(dpcovdata[['E13_ribo_1']]%>%mean%>%
	enframe('trid','covmean'))%>%
	filter(trid%in%singletrtrs)%>%
	{quicktest(.$value,.$covmean)}

	cdslens <- cdsgrl%>%width%>%sum

	

	abundance = iso_tx_countdata$abundance[,'E13_ribo_1',drop=F]%>%
		as.data.frame%>%
		rownames_to_column('trid')%>%set_colnames(c('trid','abundance'))

	fullbamdata = fread('pipeline/deepshapebamdata/E13_ribo_1.bam.reformat.absolute')
	#something about the procedure didn't work...
	fullbamdata$V7%<>%trimids
	fullbamdata$V6%<>%trimids
	fullbamdata%<>%rename('gid'=V7,'trid'=V6)
	# fullbamdata = dtplyr::lazy_dt(fullbamdata)%>%left_join(abundance,by='trid')
	fullbamdata%<>%as_tibble
	
	bamdata = fullbamdata%>%left_join(abundance)%>%group_by(V1)%>%sample_n(1,weight=abundance/sum(abundance))
	bamdata = 
	dpbamcounts = bamdata%>%group_by(gid,trid)%>%tally%>%select(gid,trid,count=n)
	dpbamcounts = bamdata%>%group_by(gid,trid)%>%tally%>%select(gid,trid,count=n)
	dpbamcounts%<>%	
		left_join(cdslens%>%enframe('trid','len'))%>%
		mutate(dens=count/len)

	#okay so the shit in the bam data lines up beautifully... which is kind of odd
	dpbamabundcomp = dpbamcounts%>%	
		left_join(iso_tx_countdata$abundance[,'E13_ribo_1',drop=F]%>%
			as.data.frame%>%
			rownames_to_column('trid')%>%
			pivot_longer(-trid)
		)

	dpbamabundcomp%>%{quicktest(.$value,.$dens)}

	dpbamabundcomp%>%left_join(dpcovdata[[1]]%>%mean%>%enframe('trid','covmean'))%>%
		{quicktest(.$covmean,.$dens)}		


	bamdata.rtpms = fread('pipeline/deepshapebamdata/E13_ribo_1.bam.reformat.absolute.readtpms')
	bamdata.wtpms = fread('pipeline/deepshapebamdata/E13_ribo_1.bam.reformat.absolute.withtpm')
	bamdata.rtpms[1:3,]%>%left_join(bamdata,by=c(V1='V2'))


          print(filestring)
	#okay so the sampling procedure seems to make the correlation worse....



	dpcovdata[['E13_ribo_1']]%>%names%>%head(1)%>%cdslens[.]
	dpcovdata[['E13_ribo_1']]%>%runLength%>%sum%>%head(1)

	jdata = fread('pipeline/dpprimecov/E13_ribo_1/dpprime_dist_ribosig.tsv')

	bamdata%<>%group_by(V1)%>%mutate(m=n())
	
	#it's not exact though, so eh, fuck it.
	#So is our dp coverage object reflictive of that
	#Do these coverage objects look right?


	#Find a big spike somewhere
	# fread('E13_ribo_1.pausepred.csv')%>%arrange(desc(`Z-score`))

	offsets<-read_tsv(here('ext_data/offsets_manual.tsv'))
	bams <- Sys.glob(here('pipeline/star/data/*/*.bam'))%>%str_subset(negate=TRUE,'transcript')

	bam=bams[1]
	testspike = GRanges('chr2:25271475')
	mapqthresh=200
	library(GenomicAlignments)
	riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=mapqthresh,which=testspike)
	reads <- readGAlignments(bam,param=riboparam,use.names=T)
	reads%<>%subset(qwidth%in%offsets$length)

	testspiketr = testspike%>%mapToTranscripts(exonsgrl)
	testspikecds = testspike%>%mapToTranscripts(cdsgrl)
	testtr=seqnames(testspikecds)%>%unlist%>%as.vector


	bamdatatesttrrle = bamdata%>%filter(trid==testtr)%>%group_by(V5)%>%tally%>%left_join(data.frame(V5=1:cdslens[testtr]),.)%>%
		mutate(n=replace_na(n,0))%>%.$n%>%Rle

	stopifnot(which.max(bamdatatesttrrle)==start(testspikecds))
	###oookay so it looksl ike my genomic alignment of psites the bamdata stuff works out juuuust fine.
	#and it also looks like the dp cov object reflects it, but is the wrong size...
	stopifnot(which.max(dpcovdata[[1]][[testtr]])==start(testspikecds))


	testspikecds

	reads%>%apply_psite_offset()
	get_genomic_psites(bam,testspike,	offsets)




inclusiontable(names(tputrs),alltrs)


deepshapecov[startwinds]%>%as.matrix%>%colSums%>%txtplot

cdsgrl['ENSMUST00000000096']%>%unlist%>%mapToTranscripts(exonsgrl['ENSMUST00000000096'])


names(deepshapecov)%>%inclusiontable(names(startwinds))
commonnames = intersect(names(startwinds),names(deepshapecov))

deepshapecov[commonnames][startwinds[commonnames]]

deepshapecovsums = sum(deepshapecov)




	
	get_genomic_psites()
	#Verify that our coverage tracks look right


# #now calculate spec. coefficients.
# deepshapespevals = deepshapecovinorfquant%>%mclapply(.%>%as.vector%>%ftestvect)
# deepshapespecs = deepshapespevals%>%map_dbl(1)%>%sqrt%>%enframe('transcript_id','spec_coef')
# deepshapespec_pM = deepshapespecs%>%mutate(spec_pM = spec_coef / (sum(spec_coef,na.rm=T)) * 1e6)


# deepshapespewave = deepshapecovinorfquant%>%
# 		mclapply(mc.cores=10,function(v)DWPT(as.vector(v)))
# deepshapespewavemn = deepshapespewave%>%map_dbl(mean)

# txtplot(log1p(deepshapespewavemn),log1p(deepshapecovinorfquant%>%mean))

################################################################################
########Just playing with the stats here
################################################################################
	


txtplot(deepshapecovinorfquant%>%runLength%>%sum,(deepshapespevals%>%map_dbl(1)))


txtplot(deepshapecovinorfquant%>%mean,sqrt((deepshapespevals%>%map_dbl(1))))


txtdensity(lm(deepshapecovinorfquant%>%mean~sqrt((deepshapespevals%>%map_dbl(1))))$residuals)

library(outliers)

cooksd <- cooks.distance(lm(deepshapecovinorfquant%>%mean~sqrt((deepshapespevals%>%map_dbl(1)))))

isoutlier = cooksd > 4/length(cooksd)

}
