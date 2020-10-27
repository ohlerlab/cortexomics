
#' ---
#' title: "Sequence specific offset determination"
#' author: "Dermot Harnett"
#' ---


#+ setup, include=FALSE, echo=FALSE, eval=T
knitr::opts_chunk$set(root.dir = here::here(),eval=FALSE,cache=FALSE,echo=FALSE,warning = FALSE,message = FALSE,include=FALSE)



# suppressMessages({library(svglite)})
suppressMessages({library(readr)})
suppressMessages({library(Biostrings)})
suppressMessages({library(Rsamtools)})
#suppressMessages({library(psd)})
suppressMessages({library(txtplot)})
suppressMessages({library(rtracklayer)})
suppressMessages({library(stringr)})
suppressMessages({library(data.table)})
suppressMessages({library(assertthat)})
suppressMessages({library(parallel)})
suppressMessages({library(dplyr)})
# suppressMessages({library(riboWaltz)})
suppressMessages({library(purrr)})
suppressMessages({library(here)})
suppressMessages({library(magrittr)})
suppressMessages({library(stringr)})
suppressMessages({library(tidyverse)})
suppressMessages({library(GenomicAlignments)})
suppressMessages({library(GenomicFeatures)})
suppressMessages({library(GenomicFiles)})
suppressMessages({library(bamsignals)})
suppressMessages({library(memoise)})
suppressMessages({library(here)})

source(here('src/R/Rprofile.R'))

STOPWINDOWSTART = -2
STOPWINDOWEND = 2

reduce <- GenomicRanges::reduce

MAPQTHRESH <- 50


 
string2onehot<-function(scol) vapply(c('A','C','T','G'),function(lev)as.numeric(scol==lev),rep(1,length(scol)))
dnaseq2onehot <- function(mat,pre){
	mat<-as.matrix(mat);
	lapply(1:ncol(mat),function(n) string2onehot(mat[,n])%>%set_colnames(paste0(pre,n,'.',colnames(.))))%>%purrr::reduce(cbind)
}



mimport<-mymemoise(function(...)rtracklayer::import(...))
mapToTranscripts<-mymemoise(GenomicFeatures::mapToTranscripts)
pmapToTranscripts<-mymemoise(GenomicFeatures::pmapToTranscripts)
readGAlignments<-mymemoise(GenomicAlignments::readGAlignments)
splitmaptotranscripts <- mymemoise(function(reads,transcripts) reads%>%split(.,ceiling(seq_along(.)/ 50e3))%>%lapply(.%>%mapToTranscripts(transcripts%>%split(.$transcript_id)))%>%Reduce(f=c))

fp <-function(gr)ifelse(strand(gr)=='-',end(gr),start(gr))
tp <-function(gr)ifelse(strand(gr)=='-',start(gr),end(gr))
strandshift<-function(gr,shift)shift(gr , ifelse(strand(gr)=='-',- shift,shift))

for(fname in lsf.str('package:dplyr')) assign(fname,get(fname,'package:dplyr'))



myf <- 2 %>% mymemoise(function(x){message('foo'); x +1})(.)



argv <- c(
	bam = here('pipeline/star/data/E13_ribo_2/E13_ribo_2.bam')%T>%{stopifnot(file.exists(.))},
	gtf = here('pipeline/my_gencode.vM12.annotation.gtf'),
	outfolder = 'riboWaltz/E13_ribo_2.bam/'
)

argv[] <- commandArgs(trailing=TRUE)

for (nm in names(argv)) assign(nm,argv[[nm]])

bam %T>%{stopifnot(file.exists(.))}
REF <- here('pipeline/my_GRCm38.p5.genome.chr_scaff.fa')

#get exons
gtf_gr<-mimport(con=gtf,format='gtf')

if(!is('exons','GRanges')) exons <- gtf_gr%>%subset(type=='exon')
if(!is('cds','GRanges')) cds <- gtf_gr%>%subset(type=='CDS')
if(!exists('startcods')) startcods <- gtf_gr%>%subset(type=='start_codon')

#define compartments
compartments <- rep('nucl',length(seqlevels(cds)))%>%setNames(seqlevels(cds))
circs_in_data <- intersect(DEFAULT_CIRC_SEQS,names(compartments))
compartments[circs_in_data] <- circs_in_data


#' # Improving Periodicity with a sequence specific model. 
#' ## Procedure: 
#' - For a given library, quickly calculate read overlaps for all cds

#get counts over all cds
cdsstren<-mymemoise(bamCount)(bam,cds)
cds$count<-cdsstren
trcounts<-cds$count%>%split(cds$transcript_id)%>%map_dbl(sum)

#gene tr relationshiop df
gtrdf<-exons%>%mcols%>%.[,c('gene_id','transcript_id')]%>%as.data.frame%>%distinct



#' - For a given library, Now take the trancsript with the most overlaps, for the top 1k genes

startcodcount <- startcods%>%.$transcript_id%>%table

simpletrs<-startcodcount[startcodcount==1]%>%enframe%>%select(transcript_id=name)
ccds <- cds%>%mcols%>%as.data.frame%>%select(transcript_id,tag)%>%filter(tag=='CCDS')%>%distinct
stopifnot(nrow(ccds)>1e3)
#get the tr with the highest signal per gene
#get the tr with the highest signal per gene
toptrs <- gtrdf%>%
	left_join(enframe(trcounts,'transcript_id','count'))%>%
	semi_join(simpletrs)%>%
	semi_join(ccds)%>%
	group_by(gene_id)%>%
	slice(which.max(count))%>%
	arrange(desc(count))%>%
	.$transcript_id%>%
	head(1e3)


#' - Ensuring these have only 1 start codon, for simplicity

topcds <- cds%>%
	subset(transcript_id %in% toptrs)%>%
	identity
	# head(1e3)

topstartcods<-startcods[match(topcds%>%.$transcript_id%>%unique,startcods$transcript_id)]

#get the exons for these
topcdsexons <- exons%>%subset(transcript_id %in% toptrs)

#' - Map our reads to these transcripts

#mapped cds
topcdsmap<-topcds%>%pmapToTranscripts(topcdsexons%>%split(.,.$transcript_id)%>%.[topcds$transcript_id])%>%reduce
#also start codons
starcods_trmap<-topstartcods%>%mapToTranscripts(topcdsexons%>%split(.$transcript_id))
starcods_trmap%<>%setNames(as.character(seqnames(.)))


library(GenomicAlignments)
#get reads over them
readtopreads<-function(bam,top_cds=topcds){bam%>%BamFile(yieldSize=NA)%>%readGAlignments(param=ScanBamParam(which=top_cds))}
readtopreads<-mymemoise(readtopreads)
topcdsreadsbak <- bam %>% readtopreads
topcdsreads<-topcdsreadsbak

#reads as a gr
topcdsreadsgr<-topcdsreads%>%as("GRanges")
#get the shifts and store in hte gr
topcdsreadsgr$length<-qwidth(topcdsreads)
#now map these to transcript space


cdsread_trmap<-topcdsreadsgr%>%splitmaptotranscripts(topcdsexons)

# testwidde<-cdsread_trmap%>%subset(width>40)%>%.[1])
cdsread_trmap$length<-topcdsreadsgr$length[cdsread_trmap$xHits]
#select only those which mapped cleanly onto transcripts
cdsread_trmap%<>%trim
cdsread_trmap<-cdsread_trmap%>%subset(width==length)

cdsread_trmap$phase <- ((start(cdsread_trmap) - start(starcods_trmap[seqnames(cdsread_trmap)])) %%3)

non3ttr<-coverage(cdsread_trmap)[topcdsmap]%>%lengths%>%map_dbl(`%%`,3)%>%.[.!=0]%>%names
stopifnot(length(non3ttr)==0)#all the cds should be multiples of 3

#assess entropy of the phase - should 
get_frame_entropy<-function(gr,topstartcodphases){
	stopifnot(all(seqnames(gr) %in% names(topstartcodphases)))
	adjstart <- start(gr) - topstartcodphases[as.character(gr@seqnames)]
	ptable<-adjstart %>% `%%`(3)%>%table
	list(ptable%>%{./sum(.)}%>%{-sum(.*log(.))},ptable/sum(ptable))
}


readsizes <- cdsread_trmap$length%>%table%>%{./sum(.)}%>%keep(. > 0.05)%>%names%>%as.numeric%>%setNames(.,.)
overalloffsets<-seq(6,max(readsizes),by=3)%>%setNames(.,.)

#' - Annotate the reads on these transcripts with compartment as well as their read length, phase (fram relative to start codon), 

#add compartment to the mapped reads
trcomps<-data_frame(chr=as.character(topcdsexons@seqnames),transcript_id=topcdsexons$transcript_id)%>%mutate(compartment=compartments[chr])%>%
	distinct(transcript_id,compartment)%>%{setNames(CharacterList(as.list(.$compartment)),.$transcript_id)}
cdsread_trmap$compartment<-trcomps[cdsread_trmap@seqnames]%>%unlist

#cdsread_trmap%>%split(paste(cdsread_trmap$compartment,cdsread_trmap$length,cdsread_trmap$phase,sep=';'))

#' - Go through each combination of compartment, offset (only multiples of 3, too preserve phase), readlength (using only read lengths accounting for above 5% of the reads) 
#' and get the fraction of the reds on our 1k transripts which map within the CDS


get_offsets<-mymemoise(function(cdsread_trmap,overalloffsets,compartments,readsizes){
	lapply(overalloffsets,function(offset){
		lapply(unique(compartments),function(compartment_i){
			lapply(readsizes,function(length){
				lapply(0:2%>%setNames(.,.),function(phase_i){
					if((length - 6 - offset) < 0) return(data_frame(score=NA))
					cat('.')
					cdsread_trmap%>%
						subset(compartment==compartment_i)%>%
						subset(phase==phase_i)%>%
						subset(length==length)%>%
						resize(1,'start',ignore.strand=T)%>%
						shift(offset)%>%
						countOverlaps(topcdsmap)%>%
						`>`(0)%>%sum%>%data_frame(score=.)
				})%>%bind_rows(.id='phase')
			})%>%bind_rows(.id='length')
		})%>%bind_rows(.id='compartment')
	})%>%bind_rows(.id='offset')
})

offset_cds_scores<-get_offsets(cdsread_trmap,overalloffsets,compartments,readsizes)

offset_cds_scores%<>%mutate_at(vars(everything()),as.numeric)
offset_cds_scores$compartment<-unique(compartments)[offset_cds_scores$compartment]

#' - select the ones that have the best score for each length and phase
#We now have optimal offsets per readlength/phase
bestscores<-offset_cds_scores%>%group_by(length,phase)%>%slice(which.max(score))

#
outfolder%>%dir.create(showWarnings=F,rec=T)
bestscores%>%write_tsv(file.path(outfolder,'cdsmax_offsets.tsv'))


apply_cds_offsets <- function(cdsread_trmap,bestscores){
	mcols(cdsread_trmap)%>%as.data.frame%>%left_join(bestscores,by=c('length','phase'))%>%.$offset
}

cdsread_trmap$cdsshift <- apply_cds_offsets(cdsread_trmap,bestscores)

#confirmed - This gives us maximum cds content per length and phase 
cdsread_trmap%<>%subset(!is.na(.$cdsshift))
cdsread_trmap%<>%intersect(.,trim(.))


topExonseq <- getSeq(x=FaFile(REF),topcdsexons)
topExonseq%<>%split(topcdsexons$transcript_id)

#check - This works to get cds sequences
# topExonseq%>%lapply(.%>%unlist%>%translate%>%str_detect('\\*'))

topExonseq%<>%lapply(.%>%unlist)%>%DNAStringSet

# stopifnot(topExonseq%>%lapply(`[`,1:3)%>%DNAStringSet%>%as.character%>%`==`('ATG'))


###OKay so I want to... make small shifts depending on sequence that will... maximise cds content,
####And also bring my reads into alignment
###I think a reasonable appraoch is to use the stop codon as our alignment point

setpos<- .%>%{strand(.)<-'+';.}

prestops <- topcdsmap%>%setpos%>%resize(3,'end')%>%resize(1,'start')
prestops%<>%setNames(seqnames(.))
#all stop codons!
stopifnot(topExonseq[prestops%>%resize(3)%>%shift(3)]%>%translate%>%`==`('*'))

cdsread_trmap%<>%setpos

####Addstarts and stop annotations to our reads

#get the psites predicted to -+2 of the first nt of our pre-stop codons
cdsread_trmap$stop<-cdsread_trmap%>%resize(1,'start')%>%shift(.$cdsshift)%>%{strand(.)<-'+';.}%>%overlapsAny(resize(prestops,5,'center'))
#get the distance to the 1 site of the prestop
cdsread_trmap$stopdist<-start(cdsread_trmap) + cdsread_trmap$cdsshift- start(prestops[seqnames(cdsread_trmap)])




prestarts <- topcdsmap%>%setpos%>%resize(3,'start')%>%resize(1,'start')
prestarts%<>%setNames(seqnames(.))
#all stop codons!
stopifnot(topExonseq[prestarts%>%resize(3)]%>%translate%>%`==`('M'))
#also starts
cdsread_trmap$startcodpos<-cdsread_trmap%>%resize(1,'start')%>%shift(.$cdsshift)%>%{strand(.)<-'+';.}%>%overlapsAny(resize(prestarts,5,'center'))
#get the distance to the 1 site of the prestop
cdsread_trmap$startdist<-start(cdsread_trmap) - start(prestarts[seqnames(cdsread_trmap)])

prestarts <- topcdsmap%>%setpos%>%resize(3,'start')%>%resize(1,'start')
prestarts%<>%setNames(seqnames(.))
#all stop codons!
stopifnot(topExonseq[prestarts%>%resize(3)]%>%translate%>%`==`('M'))
#also starts
cdsread_trmap$startcodpos<-cdsread_trmap%>%resize(1,'start')%>%shift(.$cdsshift)%>%{strand(.)<-'+';.}%>%overlapsAny(resize(prestarts,5,'center'))
#get the distance to the 1 site of the prestop
cdsread_trmap$startdist<-start(cdsread_trmap) - start(prestarts[seqnames(cdsread_trmap)])





################################################################################
########Train model for sequence shifting
################################################################################

stoppsites<-cdsread_trmap%>%subset(between(stopdist,STOPWINDOWSTART,STOPWINDOWEND))

cdsread_trmap%>%head%>%resize(3)%>%shift(3)%>%shift(-.$stopdist)%>%topExonseq[.]%>%translate

topExonseq[stoppsites%>%shift(- .$cdsshift)%>%resize(2,'start')%>%resize(4,'end')]%>%as.matrix%>%apply(2,table)
topExonseq[stoppsites%>%shift(- .$cdsshift)%>%shift(.$length)%>%resize(2,'start')%>%resize(4,'end')%>%trim%>%subset(width==4)]%>%as.matrix%>%apply(2,table)


midseq <- topExonseq[stoppsites%>%resize(2,'center')%>%resize(4,'start')]

nttable<-.%>%factor(level=c('A','C','T','G'))%>%table
# fpcutpwm<-(apply(as.matrix(startseq),2,nttable) / apply(as.matrix(midseq),2,nttable))%>%log2
# tpcutpwm<-(apply(as.matrix(endseq),2,nttable) / apply(as.matrix(midseq),2,nttable))%>%log2

#Now we should use these pwms to get cut site scores around the edge of our reads.
# startedgeseq <- topExonseq[stoppsites%>%{strand(.)<-'+';.}%>%resize(3,'start')%>%resize(6,'end')%>%trim%>%subset(width=6)]
# endedgeseq <- topExonseq[stoppsites%>%{strand(.)<-'+';.}%>%resize(3,'start')%>%resize(6,'end')%>%trim%>%subset(width=6)]


STOPWINDOWSTART=-2
STOPWINDOWEND=2


#get the psites predicted to -+2 of the first nt of our pre-stop codons
cdsread_trmap$stop<-cdsread_trmap%>%resize(1,'start')%>%shift(.$cdsshift)%>%{strand(.)<-'+';.}%>%overlapsAny(resize(prestops,5,'center'))
#get the distance to the 1 site of the prestop
cdsread_trmap$stopdist<-start(cdsread_trmap) + cdsread_trmap$cdsshift- start(prestops[seqnames(cdsread_trmap)])

stoppsites<-cdsread_trmap%>%subset(between(stopdist,STOPWINDOWSTART,STOPWINDOWEND))



#for now let's not consider base pairs with no strong pwm
library(ranger)
exonseq<-topExonseq

get_seqforrest_traindata <- function(stoppsites,exonseq){

	stopifnot(c('stopdist','length') %in% colnames(mcols(stoppsites)))
	stopifnot(all(seqnames(stoppsites)%in%names(exonseq)))
	stopifnot(is(exonseq,'DNAStringSet'))

	startseq <- exonseq[stoppsites%>%resize(2,'start')%>%resize(4,'end')]
	
	endseq <- exonseq[stoppsites%>%resize(2,'end')%>%resize(4,'start')]

	seqmat <- cbind( dnaseq2onehot(startseq,'fp.'),dnaseq2onehot(endseq,'tp.'))

	seqmat%>%cbind(stopdist=stoppsites$stopdist,length=stoppsites$length)%>%as.data.frame

}

readseqdata <- get_seqforrest_traindata(stoppsites,topExonseq)

#' - Then, taking the reads whose psites now fall within +-2bp of a STOP codon, use a random forrest classifier a. la. the scikit ribo paper
#' to try to learn offsets, such that we maximize the number of Psites landing exactly on stop codons, and with the length, phase, and
#' basepairs +-2bp up/downstream of either cut site as features (n.b. no features selection just yet might improve things)
shiftforrestfit <- ranger::ranger(formula= factor(stopdist) ~ . ,
	data=readseqdata,
	importance='permutation'
)

predictedshifts<-mymemoise(predict)(shiftforrestfit,data=readseqdata)$prediction%>%as.character%>%as.numeric

stoppsites%>%shift( - predictedshifts )%>%shift(3)%>%resize(3)%>%topExonseq[.]%>%translate



#filter set of all reads so we can use them all (ones too close to edges can't get seq)
nofpflank<-topcdsmap%>%start%>%`<`(3)
notpflank<-topcdsmap%>%{seqlengths(.)[as.vector(seqnames(.))] - end(.)}%>%`<`(3)
hasflanks <- (!nofpflank)&(!notpflank)
trswithflanks <- hasflanks%>%keep(~ . ) %>% names

cdsread_trmap%<>%subset(seqnames%in%trswithflanks)
cdsread_trmap%<>%subset(start>2)
cdsread_trmap%<>%{.[( seqlengths(.)[as.vector(seqnames(.))] - end(.) ) > 1] }

startseq <- topExonseq[cdsread_trmap%>%resize(2,'start')%>%resize(4,'end')]
endseq <- topExonseq[cdsread_trmap%>%resize(2,'end')%>%resize(4,'start')]

seqmat <- cbind( dnaseq2onehot(startseq,'fp.'),dnaseq2onehot(endseq,'tp.'))

readseqdata<-seqmat%>%cbind(length=cdsread_trmap$length)%>%as.data.frame

predictedshifts<-predict(shiftforrestfit,data=readseqdata)$prediction%>%as.character%>%as.numeric

#it's the model thats actually sure

cdsread_trmap$seqshift <- predictedshifts



topcdsmap%<>%setNames(seqnames(.))


cdstosamp<-seqnames(cdsread_trmap)%>%unique%>%.[100:200]%>%as.character


cdscov <- function(reads2use,cds2use) coverage(reads2use)[cds2use]%>%unlist%>%as.vector

get_frac_inds<-function(vect,l,r){
	lind <- floor(l*length(vect))
	rind <- floor(r*length(vect))
	vect[lind:rind]
}
# get_frac_inds(1:20,0.3,0.36)
get3bpfrac<-function(x) x %>%
	# add(c(1000000,0,0))%>%
	fft%>%
	abs%>%
	`^`(2)%>%
	{sum(get_frac_inds(.,0.3,0.36))/sum(.) }

#' ## Testing
#' 
#' We can see if our procedure worked by quickly looking at the fft (fast fourier transform) on the top1k genes
#' binned randomly into 5 groups, before and after the application of our sequence specific cutoffs (to ALL reads),
#' rather than just those around the stop codon
get_periodicity_scores<-mymemoise(function(cdstosamp,cdsread_trmap){
	message('.')
	psites2use<-cdsread_trmap%>%
		keepSeqlevels(cdstosamp,'coarse')%>%
		resize(1,'start')%>%
		shift(.$cdsshift)


	cds2use <- topcdsmap[cdstosamp]

	periodicfrac<-cdscov(psites2use,cds2use)%>%get3bpfrac

	data_frame(
		no_seqshift = psites2use%>%cdscov(cds2use)%>%get3bpfrac,
		with_seqshift = psites2use %>% shift(-.$seqshift)%>%cdscov(cds2use)%>%get3bpfrac
	)
})


seqshift_periodicities<- seqnames(cdsread_trmap)%>%
	unique%>%
	as.character%>%
	split(seq_along(.)%%5)%>%
	lapply(F=get_periodicity_scores,cdsread_trmap)

seqshift_periodicities%<>%bind_rows


seqshift_periodicities%>%unlist%>%txtplot


#+ spectral_coefficient_strip_plot, fig.width =4,fig.height=4,out.width=400,out.height=450,dev='pdf',include=TRUE,eval=TRUE
{
	seqshift_periodicities%>%bind_rows%>%gather(set,spectral_coefficient)%>%qplot(data=.,x=set,color=set,y=spectral_coefficient,geom=c('point'))+theme_bw()
}

stop()

saveRDS(.GlobalEnv,'workingenv.rds')

# #Get riboqc_cutoffs
# sample = bam%>%dirname%>%basename 
# riboqccutoffs <- str_interp(here('pipeline/riboqc/data/${sample}/_P_sites_calcs'))%T>%{stopifnot(file.exists(.))}
# riboqcdf <- riboqccutoffs%>%fread%>%select(length=read_length,riboqc_shift=cutoff)


# fcods=20
# lflank=(fcods*3)
# rflank=(fcods*3)+2
# mpoint=((fcods*2)+1)*3
# epoint=((fcods*4)+2)*3
# MIDDLESHIFT=0#WHY do I need this????
# #
# # cdsread_trmap%>%as.data.frame%>%mutate(d=startdist+stopdist)%>%.$d%>%hist
# # cdsread_trmap%>%as.data.frame%>%mutate(d=startdist+stopdist)%>%.$d%>%between(-lflank,rflank)%>%table
# metaplotlabs<-c(paste0('-',lflank/2),'AUG',lflank/2,paste0('mid -',lflank/2),'mid',paste0('mid +',lflank/2),paste0('end -',lflank/2),'stop',paste0('end +',lflank/2))
# metaplotbreaks<-c(-lflank/2,0,lflank/2,mpoint-lflank/2,mpoint,mpoint+lflank/2,epoint-(lflank/2),epoint,epoint+(lflank/2))
# #
# seqshiftfuncs <- list(
# 	Riboqc = .%>% safe_left_join(riboqcdf,by=c('length'))%>% mutate(startdist=startdist+riboqc_shift,stopdist=stopdist+riboqc_shift),
# 	CDSmax_shift = .%>% mutate(startdist=startdist+cdsshift,stopdist=stopdist+cdsshift),
# 	seqshift = .%>% mutate(startdist=startdist+cdsshift,stopdist=stopdist+cdsshift)%>%mutate(startdist=startdist-seqshift,stopdist=stopdist-seqshift))
# #

# #' We can also simply look at a metaplot of the frames, (comparing riboqc and cdsmax shift)
# plist=lapply(seqshiftfuncs,mymemoise(function(seqshiftfunc){
# 	#
# 	p=cdsread_trmap%>%
# 		# sample(100e3)%>%
# 		# shift(.$cdsshift)%>%
# 		mcols()%>%as_tibble	%>%
# 		select(startdist,stopdist,length,cdsshift,seqshift)%>%
# 		# filter(between(abs(startdist-stopdist),-lflank,rflank))%>%
# 		seqshiftfunc%>%
# 		mutate(
# 			mdist = startdist - ((floor(((startdist-stopdist+1)/3)/2))*3),
# 			instart=between(startdist,-lflank,rflank),
# 			inend=between(stopdist,-lflank,rflank),
# 			inmiddle=between(mdist,-lflank,rflank)
# 		)%>%
# 		filter((instart|inend|inmiddle) & (!(instart&inend&inmiddle)))%>%
# 		# filter(instart)%>%
# 		mutate(x = NA,
# 			x=ifelse(instart,startdist,x),
# 			x=ifelse(inmiddle,mdist+mpoint,x),
# 			x=ifelse(inend,stopdist+epoint,x) ) %>%
# 		mutate(length=as.character(length))%>%
# 		bind_rows(.,mutate(.,length='all'))%>%
# 		group_by(x,length)%>%tally%>%
# 		# filter(length=='all')%>%
# 		mutate(phase = factor(x%%3))%>%
# 		# arrange(desc(n))
# 		ggplot(aes(x=x,ymax=n,ymin=0,color=phase))+geom_linerange()+
# 		facet_grid(scale='free',length~.)+
# 		geom_vline(linetype=2,xintercept=c(0,epoint+2),alpha=I(0.3))+
# 		scale_x_continuous(labels=metaplotlabs,breaks=metaplotbreaks)+
# 		# coord_cartesian(xlim=c(-lflank,rflank))+
# 		geom_rect(aes(xmin=0,xmax=epoint,ymin=0,ymax=Inf,alpha=I(0.01)),fill='grey',color='white')+
# 		scale_fill_discrete(guide=F)+
# 		theme_minimal()	
# 	p
# }))

# #+ riboshift_comparison plot, fig.width =12,fig.height=12,out.width=1200,out.height=1250,dev='pdf',include=TRUE,eval=TRUE
# ggpubr::ggarrange(ncol=3,plotlist=plist,labels=names(seqshiftfuncs))


# #' .... All of which seems to indicate our procedure is doing a better job of recovering actual P-site positions. Note that it's p
# #' site positions are consistent between reads - no dramatic difference between 28 and the others - and manage to recover periodicity
# #' in specific read lengths where it was missing. 
# #' 
# #' I remain a little concerned by how the stop codon isn't quite in the CDS for some of these...
# #' We could try choosing our set such that things only have negative shifts - doesn't work as well though (also weirdly ruins periodicity for 26bp)
# #' 




