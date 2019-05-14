# suppressMessages({library(svglite)})
suppressMessages({library(readr)})
suppressMessages({library(Biostrings)})
suppressMessages({library(Rsamtools)})
suppressMessages({library(rtracklayer)})
suppressMessages({library(GenomicFeatures)})
suppressMessages({library(stringr)})
suppressMessages({library(data.table)})
suppressMessages({library(assertthat)})
suppressMessages({library(parallel)})
suppressMessages({library(Biostrings)})
suppressMessages({library(dplyr)})
suppressMessages({library(riboWaltz)})
suppressMessages({library(purrr)})
suppressMessages({library(here)})
suppressMessages({library(GenomicAlignments)})

fp <-function(gr)ifelse(strand(gr)=='-',end(gr),start(gr))
tp <-function(gr)ifelse(strand(gr)=='-',start(gr),end(gr))
strandshift<-function(gr,shift)shift(gr , ifelse(strand(gr)=='-',- shift,shift))

for(fname in lsf.str('package:dplyr')) assign(fname,get(fname,'package:dplyr'))

source(here('/src/R/Rprofile.R'))

argv <- c(
	transcriptbam = 'riboWaltz/RPI8_PolyE16_2/RPI8_PolyE16_2.star_transcript.bam',
	gtf = here('pipeline/my_gencode.vM12.annotation.gtf'),
	outfolder = 'riboWaltz/RPI8_PolyE16_2_RUST/'
)
transcriptbam <- 'riboWaltz/E16_2/E16_2.star_transcript.bam'


argv[] <- commandArgs(trailing=TRUE)

for (nm in names(argv)) assign(nm,argv[[nm]])

RUST <- str_detect(outfolder,'RUST')


bam <- here('pipeline/star/data/E13_ribo_1/E13_ribo_1.bam')%T>%{stopifnot(file.exists(.))}

readlistcolnames <- c("transcript", "end5", "end3", "length", "cds_start", "cds_stop")

if(!is(cds,"GRanges"))cds <-  gtf%>%import%>%subset(type=='CDS')
if(!exists('startcods'))startcods <- gtf%>%import%>%subset(type=='start_codon')

if(!exists('startreadsob')) startreadsob <- bam%>%BamFile(yieldSize=NA)%>%readGAlignments(param=ScanBamParam(which=startcods))
startreadsob%<>% .[njunc(startreadsob)==0]
startreads<-as(startreadsob,"GRanges")
startreads

uniqstartreads<-startreads%>%as("GRanges")%>%unique

FLANKSIZE=4

uniqstartreads$fpseq <- uniqstartreads%>%
	resize(FLANKSIZE,'start')%>%strandshift(-(0.5*FLANKSIZE))%>%
	getSeq(x=FaFile(here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa')),.)

uniqstartreads$tpseq <- uniqstartreads%>%
	resize(FLANKSIZE,'end')%>%strandshift(0.5*FLANKSIZE)%>%
	getSeq(x=FaFile(here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa')),.)

uniqstartreads$tpseq%>%as.matrix%>%apply(2,.%>%table%>%{round(./sum(.),3)})
uniqstartreads$fpseq%>%as.matrix%>%apply(2,table)


scol<-uniqstartreads$fpseq%>%as.matrix%>%.[,1]


breakpointbreaks = function(limits){
	nd2=(limits[2]-limits[1] +1) / 2
	s = seq(limits[1],limits[2])
	n = length(s)
	fhalf<-c(s[1:nd2])
	rhalf<-c(s[-(1:nd2)])
	bpoint<-mean(c(tail(fhalf,1),head(rhalf,1)))
	c(fhalf,bpoint,rhalf)%>%setNames(c(seq(-nd2,-1),'//',seq(1,nd2)))
}



dbrks <- breakpointbreaks(c(1,FLANKSIZE))
labs <- names(dbrks)


pfile<-here(paste0('plots/readprocessing/fpbasecomp.',basename(bam),'.pdf'))%T>%{dir.create(dirname(.))}
pdf(pfile)
ggpubr::ggarrange(ncol=2,
	uniqstartreads$fpseq%>%as.matrix%>%apply(2,.%>%table%>%{round(./sum(.),3)})%>%melt%>%set_colnames(c('Base','Position','Frequency'))%>%
		qplot(data=.,color=Base,x=Position,y=Frequency)+geom_line()+
		ggtitle("Base Composition 5' End of Reads")+
		scale_x_continuous(name="5' Position",breaks=dbrks,labels=labs)+
		theme_minimal(),
	uniqstartreads$tpseq%>%as.matrix%>%apply(2,.%>%table%>%{round(./sum(.),3)})%>%melt%>%set_colnames(c('Base','Position','Frequency'))%>%
		qplot(data=.,color=Base,x= Position,y=Frequency)+geom_line()+
		scale_x_continuous(name="3' Position",breaks=dbrks,labels=labs)+
		ggtitle("Base Composition 3' End of Reads")+
		theme_minimal()
)%>%ggpubr::annotate_figure(top=paste0('Read Composition for ',basename(bam)))
dev.off()
message(pfile)


readlocscore<-match(as(startreads,'GRanges'),uniqstartreads)%>%table%>%enframe('loc','score')

uniqstartreads$score <- readlocscore$score

locidaugs<-uniqstartreads%>%{.$id<-seq_along(.);.}%>%mergeByOverlaps(startcods)%>%{data.frame(locid=.$id,augloc=fp(strandshift(.$startcods,1)))}%>%
	group_by(locid)%>%filter(n()==1)
#?
all(locidaugs$locid==seq_len(nrow(locidaugs)))

uniqstartreads$augloc<-NA
uniqstartreads$augloc[locidaugs$locid]<-locidaugs$augloc

uniqstartreads$augdist <- abs(uniqstartreads$augloc - fp(uniqstartreads) )



uniqstartreads$score%<>%as.numeric
uniqstartreads$width<-NULL
uniqstartreads$length<-uniqstartreads%>%width

#old seperate char cols
#mcols(uniqstartreads) %<>% cbind(uniqstartreads$fpseq%>%as.matrix%>%set_colnames(paste0('fp',seq_len(ncol(.)))))
#mcols(uniqstartreads) %<>% cbind(uniqstartreads$tpseq%>%as.matrix%>%set_colnames(paste0('tp',seq_len(ncol(.)))))

 
string2onehot<-function(scol) vapply(c('A','C','T','G'),function(lev)as.numeric(scol==lev),rep(1,length(scol)))
dnaseq2onehot <- function(mat,pre){
	mat<-as.matrix(mat);
	lapply(1:ncol(mat),function(n) string2onehot(mat[,n])%>%set_colnames(paste0(pre,n,'.',colnames(.))))%>%purrr::reduce(cbind)
}

test<-uniqstartreads%>%head(3)

#create training data from the granges object

traindata<-uniqstartreads%>%
	subset(!is.na(augdist))%>%
	mcols
seqonehot<-cbind(dnaseq2onehot(mat=DNAStringSet(traindata$fpseq),'fp'),dnaseq2onehot(mat=traindata$tpseq,'tp'))

traindata%<>%as.data.frame%>%select(-fpseq,-tpseq)

traindata%<>%cbind(seqonehot)

traindata%<>%subset(!is.na(augdist))%>%{.[rep(1:nrow(.),.$score),]}

readlenfracs <- traindata$length%>%table%>%{./sum(.)}%>%round(3)
readlenfracs <- readlenfracs%>%keep(~ .>0.05)
traindata <- traindata%>%filter(length %in% names(readlenfracs))
traindata%>%colnames%>%{stopifnot('tp2.C'%in%.)}

traindata%>%group_by(length)%>%summarise(dists=list(enframe(tail(sort(table(augdist))),'dist','freq')))%>%unnest%>%as.data.frame



library(randomForest)
traindata$length%>%is.na%>%table



0
datalist<- traindata%>%select(-score,-augloc)%>%split(.,.$augdist)%>%map(~ split(.,sample(c('train','test'),size=nrow(.),p=c(0.9,0.1),rep=TRUE)))
datalist<- list(test=datalist%>%map_df('test'),train=datalist%>%map_df('train'))
#I need to ensure that each offset is in the data 3 times

datalist%>%map(nrow)
sharedoffsets<-datalist%>%map(~.$augdist%>%table%>%enframe)%>%bind_rows(.id='set')%>%filter(n_distinct(set)==2)%>%.$name
datalist%<>%map(~filter(.,augdist %in% sharedoffsets))
sharedoffsets

predvars<-datalist$train%>%select(length,matches('fp|tp'))%>%colnames

rf_fit <- randomForest(formula= augdist ~ . ,
	x = datalist$train%>%select(one_of(predvars)), 
	xtest=datalist$test%>%select(one_of(predvars)),
	y = datalist$train$augdist%>%as_factor, 
	ytest=datalist$test$augdist%>%as_factor,
	keep.forest=TRUE
)

importancedf<-randomForestExplainer::measure_importance(rf_fit)

impvars<-importancedf$variable[importancedf$p_value<0.001]%>%as.character

rf_fit$importance%>%as.data.frame%>%rownames_to_column('variable')%>%mutate(imp=MeanDecreaseGini/max(MeanDecreaseGini))%>%arrange(desc(imp))

importancedf%>%filter(variable%in%impvars)

# imprf_fit <- randomForest(formula= augdist ~ . ,
# 	x = datalist$train%>%select(one_of(impvars)), 
# 	xtest=datalist$test%>%select(one_of(impvars)),
# 	y = datalist$train$augdist%>%as_factor, 
# 	ytest=datalist$test$augdist%>%as_factor,
# 	keep.forest=TRUE
# )

imprf_fit <- ranger::ranger(formula= augdist ~ . ,
	data=datalist$train%>%select(augdist,predvars),
	importance='permutation'
)


###Look at predictions on the uniqe data points
upreddata<-traindata%>%select(predvars)%>%unique

upreddata$predicted_augdist <- predict(rf_fit,newdata=upreddata,type='response')

impupreddata<-traindata%>%select(impvars)%>%unique

impupreddata$predicted_augdist_imp <- predict(imprf_fit,newdata=impupreddata,type='response')

impupreddata$predicted_augdist%>%table

traindata%>%select(one_of(predvars))


#####Change of strategy, just use prior knowledge

imprf_fit_lengthonly <- ranger::ranger(formula= as_factor(augdist) ~ length ,
	data=datalist$train%>%filter(between(augdist,3,max(as.numeric(sharedoffsets))))%>%select(augdist,predvars)
)

lens<-datalist$train%>%distinct(length)%>%arrange(length)
pred<-predict(imprf_fit_lengthonly,lens,type='response')
lengthoffsetdf<-lens%>%mutate(offset=as.numeric(as.character(pred$prediction)))

#Now, we know that we have reasonable length offsets, what if we now try to predict phase with the random forrest?
phasetraindata <- traindata%>%left_join(lengthoffsetdf)
phasetraindata%<>%mutate(phase = ((offset - augdist)%%3))
#phasetraindata$phase[phasetraindata$phase==2]<- -1
#look to make sure this is right
phasetraindata%>%select(length,augdist,offset,phase)%>%sample_n(5)
#this looks right I think?
phasetraindata%>%select(length,augdist,offset,phase)%>%mutate((offset - phase)  - augdist )%>%sample_n(5)
#The 'phase' must be SUBTRACTED From the offset, in order to get the offset that puts you in phase again.

#Now try to predict phase
phasepredvars<-predvars
phasepred <- ranger::ranger(formula = as_factor(phase) ~ . ,data=phasetraindata%>%select(phase,phasepredvars),importance='permutation')

#Look at predictions on the uniqe data points
upreddata<-traindata%>%select(predvars)%>%unique
upreddata$predicted_augdist <- predict(phasepred,data=upreddata,type='response',)$prediction


################################################################################
########Testing how well it worked with a selection of top CDS
################################################################################
	


string2onehot<-function(scol) vapply(c('A','C','T','G'),function(lev)as.numeric(scol==lev),rep(1,length(scol)))

dnaseq2onehot <- function(mat,pre){
	mat<-as.matrix(mat);
	lapply(1:ncol(mat),function(n) string2onehot(mat[,n])%>%set_colnames(paste0(pre,n,'.',colnames(.))))%>%purrr::reduce(cbind)
}


offset_reads_rf <- function(reads,lengthoffsetdf,phasepred,REF='pipeline/my_GRCm38.p5.genome.chr_scaff.fa'){
	#get the sequence info
		# reads <- reads[qwidth(reads)%in%lengthoffsetdf$length]
		qwidth<-qwidth(reads)
		reads <- as(reads,"GRanges")
		reads$qwidth <- qwidth

		ureads<-unique(reads)

		fpseq <- ureads%>%
			as("GRanges")%>%
			resize(FLANKSIZE,'start')%>%strandshift(-(0.5*FLANKSIZE))%>%
			getSeq(x=FaFile(here(REF)),.)

		tpseq <- ureads%>%
			as("GRanges")%>%
			resize(FLANKSIZE,'end')%>%strandshift(0.5*FLANKSIZE)%>%
			getSeq(x=FaFile(here(REF)),.)

		ldf <- data.frame(length=ureads$qwidth)

		ureads$offsets = ldf%>%left_join(lengthoffsetdf)%>%.$offset

		if(phasepred=NULL)

		ureads$phase <- predict(phasepred,data=cbind(ldf,dnaseq2onehot(fpseq,'fp'),dnaseq2onehot(tpseq,'tp')))$prediction%>%as.character%>%as.numeric

		#get the shifts for our reads

		minds <- match(reads,ureads)
		shift <- ureads$offsets[minds] + ureads$phase[minds]

		shift


}


if(!is('exons','GRanges')) exons <- gtf%>%import%>%subset(type=='exon')

cdsstren<-bamsignals::bamCount(bam,cds)
cds$count<-cdsstren

trcounts<-cds$count%>%split(cds$transcript_id)%>%map_dbl(sum)




#gene tr relationshiop df
gtrdf<-exons%>%mcols%>%.[,c('gene_id','transcript_id')]%>%as.data.frame%>%distinct


startcodcount <- startcods%>%.$transcript_id%>%table

#get the tr with the highest signal per gene
toptrs <- gtrdf%>%
	left_join(enframe(trcounts,'transcript_id','count'))%>%
	filter(startcodcount[transcript_id]==1)%>%
	group_by(gene_id)%>%
	slice(which.max(count))%>%
	arrange(desc(count))%>%
	.$transcript_id%>%
	head(1e3)


#now take the top trs with only 1 start codon
topcds <- cds%>%
	subset(transcript_id %in% toptrs)%>%
	identity
	# head(1e3)

topstartcods<-startcods[match(topcds%>%.$transcript_id%>%unique,startcods$transcript_id)]

#get the exons for these
topcdsexons <- exons%>%subset(transcript_id %in% toptrs)

#get reads over them
topcdsreads <- bam%>%BamFile(yieldSize=NA)%>%readGAlignments(param=ScanBamParam(which=topcds))

#reads as a gr
topcdsreadsgr<-topcdsreads%>%as("GRanges")

#get the shifts and store in hte gr
topcdsreadsgr$length<-qwidth(topcdsreads)



#Now calculate the read shifts
topcdsreadsgr$shifts<-offset_reads_rf(topcdsreads,lengthoffsetdf,phasepred)

#now map these to transcript space
cdsread_trmap<-topcdsreadsgr%>%mapToTranscripts(topcdsexons%>%split(.$transcript_id))
cdsread_trmap$shift<-topcdsreadsgr$shifts[cdsread_trmap$xHits]
cdsread_trmap$length<-topcdsreadsgr$length[cdsread_trmap$xHits]

#eliminate those without a shift (wrong read length)
cdsread_trmap%<>%.[!is.na(.$shift)]

#use shift to get the psites/
topcdsreadsgrmap_psites<- cdsread_trmap%>%shift(.$shift)%>%resize(1,'start')

#Now get frame of psites
topstartcodsmapped<-topstartcods%>%pmapToTranscripts(topcdsexons%>%split(.,.$transcript_id)%>%.[topstartcods$transcript_id])
#
topstartcodphases<-topstartcodsmapped%>%start%>%`%%`(3)%>%setNames(seqnames(topstartcodsmapped))

#assess entropy of the phase - should 


get_frame_entropy<-function(gr,topstartcodphases){
	stopifnot(all(seqnames(gr) %in% names(topstartcodphases)))
	adjstart <- start(gr) - topstartcodphases[as.character(gr@seqnames)]
	ptable<-adjstart %>% `%%`(3)%>%table
	list(ptable%>%{./sum(.)}%>%{-sum(.*log(.))},ptable/sum(ptable))
}

stopifnot(topstartcodsmapped%>%get_frame_entropy(topstartcodphases)%>%.[[1]]%>%identical(0))

cdsread_trmap%>%get_frame_entropy(topstartcodphases)
topcdsreadsgrmap_psites%>%get_frame_entropy(topstartcodphases)
topcdsreadsgrmap_psites%>%shift(sample(0:2,size=length(.),rep=T))%>%get_frame_entropy(topstartcodphases)


topcdsmap<-topcds%>%pmapToTranscripts(topcdsexons%>%split(.,.$transcript_id)%>%.[topcds$transcript_id])%>%reduce



library(txtplot)

overalloffsets<-6:28%>%setNames(.,.)
readsizes <- 25:31%>%setNames(.,.)

cdscontent<-lapply(overalloffsets,function(offset){
	lapply(readsizes,function(readsize){
		if((readsize - 6 - offset) < 0) return(NA)
		message(offset)
		message(readsize)
		cdsread_trmap%>%subset(length==readsize)%>%resize(1,'start',ignore.strand=T)%>%shift(offset)%>%countOverlaps(topcdsmap)%>%`>`(0)%>%sum
	})
})


cdscontent%>%map('25')%>%unlist%>%txtplot
cdscontent%>%map('26')%>%unlist%>%txtplot
cdscontent%>%map('27')%>%unlist%>%txtplot
cdscontent%>%map('28')%>%unlist%>%txtplot

cdscontent[[1]]%>%names


topcdsreadsgr






startcds_trmap<-startcods%>%subset(transcript_id %in% trs )%>%mapToTranscripts(topcdsexons%>%split(.,.$transcript_id))%>%resize(1,'center')

topcdsreadsgrmap_psites%>%subsetByOverlaps(startcds_trmap%>%resize(101,'center'))%>%distanceToNearest

#make sure we have one start codon per transcript
stopifnot(startcds_trmap%>%seqnames%>%table%>%table)

w2startcods<-startcds_trmap%>%seqnames%>%table%>%keep(~ .==2)%>%names%>%.[1]

startcods%>%subset(transcript_id==w2startcods)

cdsread_trmap%>%seqinfo
topcdsreadsgrmap_psites%>%seqinfo
startcds_trmap%>%seqinfo


pfile<-'plots/readprocessing/comp_shifts_phaseplot.pdf'%>%here
pdf((pfile))
plot(1)
dev.off()
normalizePath(pfile)


phase_var_importance <- ranger::importance_pvalues(phasepred,data=phasetraindata,method='altmann',formula=as_factor(phase) ~ .)






#####Trying to parallelize

library("randomForest")
library("caret")
names(getModelInfo())
library(parallel)
library(doParallel)
cluster <- makeCluster(10) # convention to leave 1 core for OS
registerDoParallel(cluster)

fitControl <- trainControl(method = "cv",
                           number = 5,
                           allowParallel = TRUE)

rfParam <- expand.grid(
	mtry=floor(c(.3,.3,.8)*ncol(datalist$train)),
	ntree=c(100,500,1000,5000),importance=TRUE,keep.forest=TRUE)

rfParam <- expand.grid(.mtry=10)

m<- train(x=datalist$train%>%select(one_of(predvars)),y=datalist$train$augdist%>%as_factor,method="parRF",
	trControl=fitControl,
	tuneGrid=rfParam,mtry=10)

stopCluster(cluster)
registerDoSEQ()



object.size(rf_fit)/1e6

rf_fit
stop()



fp <-function(gr)ifelse(strand(gr)=='-',end(gr),start(gr))
tp <-function(gr)ifelse(strand(gr)=='-',start(gr),end(gr))
strandshift<-function(gr,shift)shift(gr , ifelse(strand(gr)=='-',- shift,shift))


################################################################################
########Testing how well it worked with a selection of top CDS
################################################################################
	

cdsstren<-bamsignals::bamCount(bam,cds)
cds$count<-cdsstren

trcounts<-cds$count%>%split(cds$transcript_id)%>%map_dbl(sum)




string2onehot<-function(scol) vapply(c('A','C','T','G'),function(lev)as.numeric(scol==lev),rep(1,length(scol)))
dnaseq2onehot <- function(mat,pre){
	mat<-as.matrix(mat);
	lapply(1:ncol(mat),function(n) string2onehot(mat[,n])%>%set_colnames(paste0(pre,n,'.',colnames(.))))%>%purrr::reduce(cbind)
}


offset_reads_rf <- function(reads,lengthoffsetdf,phasepred,REF='pipeline/my_GRCm38.p5.genome.chr_scaff.fa'){
	#get the sequence info
		# reads <- reads[qwidth(reads)%in%lengthoffsetdf$length]
		qwidth<-qwidth(reads)
		reads <- as(reads,"GRanges")
		reads$qwidth <- qwidth

		ureads<-unique(reads)

		fpseq <- ureads%>%
			as("GRanges")%>%
			resize(FLANKSIZE,'start')%>%strandshift(-(0.5*FLANKSIZE))%>%
			getSeq(x=FaFile(here(REF)),.)

		tpseq <- ureads%>%
			as("GRanges")%>%
			resize(FLANKSIZE,'end')%>%strandshift(0.5*FLANKSIZE)%>%
			getSeq(x=FaFile(here(REF)),.)

		ldf <- data.frame(length=ureads$qwidth)

		ureads$offsets = ldf%>%left_join(lengthoffsetdf)%>%.$offset

		ureads$phase <- predict(phasepred,data=cbind(ldf,dnaseq2onehot(fpseq,'fp'),dnaseq2onehot(tpseq,'tp')))$prediction%>%as.character%>%as.numeric

		#get the shifts for our reads

		minds <- match(reads,ureads)
		shift <- ureads$offsets[minds] + ureads$phase[minds]

		shift


}


if(!exists('exons')) exons <- gtf%>%import%>%subset(type=='exon')


#gene tr relationshiop df
gtrdf<-exons%>%mcols%>%.[,c('gene_id','transcript_id')]%>%as.data.frame%>%distinct

#get the tr with the highest signal per gene
toptrs <- gtrdf%>%
	left_join(enframe(trcounts,'transcript_id','count'))%>%
	group_by(gene_id)%>%
	slice(which.max(count))%>%
	arrange(desc(count))%>%
	.$transcript_id

startcodcount <- startcods%>%.$transcript_id%>%table

#now take the top trs with only 1 start codon
topcds <- cds%>%
	subset(transcript_id %in% trs)%>%
	subset(startcodcount[transcript_id]==1)
	head(1e3)

#get the exons for these
topcdsexons <- exons%>%subset(transcript_id %in% trs)

#get reads over them
topcdsreads <- bam%>%BamFile(yieldSize=NA)%>%readGAlignments(param=ScanBamParam(which=topcds))



#reads as a gr
topcdsreadsgr<-topcdsreads%>%as("GRanges")

#get the shifts and store in hte gr
topcdsreadsgr$shifts<-offset_reads_rf(topcdsreads,lengthoffsetdf,phasepred)

#now map these to transcript space
cdsread_trmap<-topcdsreadsgr%>%mapToTranscripts(topcdsexons%>%split(.$transcript_id))
cdsread_trmap$shift<-shifts[cdsread_trmap$xHits]

#eliminate those without a shift (wrong read length)
cdsread_trmap%<>%.[!is.na(.$shift)]

#use shift to get the psites
topcdsreadsgrmap_psites<- cdsread_trmap%>%shift(.$shift)%>%resize(1,'start')

#assess entropy of the phase - should 
topcdsreadsgrmap_psites%>%start%>%`%%`(3)%>%table%>%{./sum(.)}%>%{-sum(.*log(.))}
cdsread_trmap%>%start%>%`%%`(3)%>%table%>%{./sum(.)}%>%{-sum(.*log(.))}


plot_metagene_hm_rmd(res_all, "profiles_fivepr", names(rdata_list)[i], paste0(output_fig_path, "rds/"))




stop()





startcds_trmap<-startcods%>%subset(transcript_id %in% trs )%>%mapToTranscripts(topcdsexons%>%split(.,.$transcript_id))%>%resize(1,'center')

topcdsreadsgrmap_psites%>%subsetByOverlaps(startcds_trmap%>%resize(101,'center'))%>%distanceToNearest

#make sure we have one start codon per transcript
stopifnot(startcds_trmap%>%seqnames%>%table%>%table)

w2startcods<-startcds_trmap%>%seqnames%>%table%>%keep(~ .==2)%>%names%>%.[1]

startcods%>%subset(transcript_id==w2startcods)

cdsread_trmap%>%seqinfo
topcdsreadsgrmap_psites%>%seqinfo
startcds_trmap%>%seqinfo


pfile<-'plots/readprocessing/comp_shifts_phaseplot.pdf'%>%here
pdf((pfile))
plot(1)
dev.off()
normalizePath(pfile)


phase_var_importance <- ranger::importance_pvalues(phasepred,data=phasetraindata,method='altmann',formula=as_factor(phase) ~ .)






#####Trying to parallelize

library("randomForest")
library("caret")
names(getModelInfo())
library(parallel)
library(doParallel)
cluster <- makeCluster(10) # convention to leave 1 core for OS
registerDoParallel(cluster)

fitControl <- trainControl(method = "cv",
                           number = 5,
                           allowParallel = TRUE)

rfParam <- expand.grid(
	mtry=floor(c(.3,.3,.8)*ncol(datalist$train)),
	ntree=c(100,500,1000,5000),importance=TRUE,keep.forest=TRUE)

rfParam <- expand.grid(.mtry=10)

m<- train(x=datalist$train%>%select(one_of(predvars)),y=datalist$train$augdist%>%as_factor,method="parRF",
	trControl=fitControl,
	tuneGrid=rfParam,mtry=10)

stopCluster(cluster)
registerDoSEQ()



object.size(rf_fit)/1e6

rf_fit
stop()


#I should try this with one hot encoding, I think.
library(VennDiagram)
pdf('tmp.pdf')
grid.newpage()
draw.pairwise.venn(area1=267+11969,cross.area= 11969,area2= 3006+11969, category = c("Uniprot Peptides", "OD5P Riboseq-based Canonical Peptides"), lty = rep("blank", 
    2), fill = c("grey", "red"), alpha = rep(0.5, 2), cat.pos = c(0, 
    0), cat.dist = rep(0.025, 2),cex=2,label.cex=2)
dev.off()
normalizePath('tmp.pdf')

