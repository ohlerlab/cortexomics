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
suppressMessages({library(magrittr)})
suppressMessages({library(Rsamtools)})

source(here("src/Rprofile.R"))

if(!exists('cds'))load(here('data/1_integrate_countdata.R'))
#load arguments
argv <- c(
	transcriptbam = here('pipeline/star_transcript/data/E16_ribo_2/E16_ribo_2.bam'),
	gtf=here('ext_data/gencode.vM12.annotation.gtf'),
	outputtsv='riboWaltz/data/E16_ribo_2/offsets.tsv'
)
argv <- commandArgs(trailingOnly=TRUE)[1:length(argv)]%>%setNames(names(argv))
for(i in names(argv)) assign(i,argv[i])

#I'm lazily having it modify the transcript names for my bams because I don't want
#to regenerate them with the trimmed fasta just yet.
devtools::load_all('Applications/riboWaltz')

RUST=0

roundup <- function(x,n) n*ceiling(x/n)
rounddown <- function(x,n) n*floor(x/n)
number_ticks <- function(limits,n=3){
	out = c(seq(min(rounddown(limits,n)),- n,n),seq(0,max(roundup(limits,n)),by=n))
	out = out[between(out,limits[1],limits[2])]
	out
}

################################################################################
########PEriodiicty with Ribowltz
################################################################################


#transcriptbam<-here('pipeline/star/data/E13_ribo_1/E13_ribo_1.star_transcript.bam')%T>%{stopifnot(file.exists(.))}
sampname<-transcriptbam%>%dirname%>%basename
ribowfolder<-here('pipeline','riboWaltz',sampname)%T>%{dir.create((.))}
file.copy(transcriptbam,ribowfolder)
transcriptbam<-file.path(ribowfolder,basename(transcriptbam))%T>%{stopifnot(file.exists(.))}

message(getwd())
message(argv)
# save.image()
# stop()
# on.exit(save.image())

if(!exists('riboWaltzanno'))riboWaltzanno <- create_annotation(gtf)

i=1
message(i);i=i+1

sampnames <- basename(dirname(transcriptbam))%>%setNames(transcriptbam%>%str_replace('.bam','')%>%basename)

if(!exists('reads_list_raw')) reads_list_raw <- bamtolist(bamfolder = dirname(transcriptbam), name_samples = sampnames,annotation = riboWaltzanno)
reads_list <- reads_list_raw

#do rust if needed
if(RUST) {

	nureads_list<- reads_list
	ureads_list<-reads_list %>% map(~ distinct(.,transcript,end5,end3,.keep_all=TRUE))

}

{
reads_list%<>%setNames(dirname(transcriptbam)%>%basename)
# example_frames <- frame_psite(reads_psite_list, sample = sampnames[[1]], region = "all")#offsets


psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")
reads_psite_list <- psite_info(reads_list, psite_offset)

dir.create(rec=T,dirname(outputtsv))
psite_offset%>%filter(total_percentage>1)%>%transmute(read_length=length,cutoff=corrected_offset_from_5)%>%
	write_tsv(outputtsv)

stopifnot(reads_psite_list[[1]]%>%filter(transcript=='ENSMUST00000000033.9')%>%nrow%>%`==`(0L))
}

fullribocovtrs=reads_psite_list[[1]]$transcript%>%unique%>%trimids%>%intersect(ribocovtrs)

filt_reads_psite_list=	reads_psite_list%>%map(.%>%	filter(trimids(transcript)%in%fullribocovtrs))

filt_reads_psite_list[[1]]%>%
	group_by(transcript)%>%
	mutate()

readaugsigdf = filt_reads_psite_list[[1]][,]%>%group_by(transcript,psite,cds_start)%>%tally%>%group_by(transcript)%>%mutate(n=n/sum(n))%>%
	mutate(pos = psite-cds_start)%>%
	group_by(pos)%>%
	summarise(signal=mean(n))

readstopsigdf = filt_reads_psite_list[[1]][,]%>%group_by(transcript,psite,cds_stop)%>%tally%>%group_by(transcript)%>%mutate(n=n/sum(n))%>%
	mutate(pos = psite-cds_stop)%>%
	group_by(pos)%>%
	summarise(signal=mean(n))

readaugsigdf%>%
	filter(between(pos,-40,100))%>%
	ggplot(data=.,aes(x=pos,y=signal))+
	geom_line()+
	scale_x_continuous(name='position',
					limits= c(-40,100),
					minor_breaks=number_ticks,breaks=partial(number_ticks,n=12))+
	theme_bw()

readstopsigdf%>%
	filter(between(pos,-100,40))%>%
	ggplot(data=.,aes(x=pos,y=signal))+
	geom_line()+
	scale_x_continuous(name='position',
					limits= c(-100,40),
					minor_breaks=number_ticks,breaks=partial(number_ticks,n=12))+
	theme_bw()

#get the psites

#
# example_psite_region <- region_psite(reads_psite_list, riboWaltzanno, sample = sampnames)

# example_frames_stratified <- frame_psite_length(reads_psite_list, sample = sampnames[[1]],
                                                # region = "all", cl = 90)

# example_frames <- frame_psite(reads_psite_list, sample = names(reads_psite_list), region = "all")
# example_frames

# # reads_psite_list[[1]]%>%filter(length%in%c(27:30))%>%split(.$length)

# example_metaprofile <- metaprofile_psite(reads_psite_list, riboWaltzanno, sample = sampnames[[1]],
#                                          utr5l = 20, cdsl = 40, utr3l = 20,
#                                          plot_title = "auto")
# example_metaprofile

lvect=30
lennames = c(lvect)%>%paste0('rl',.)
rl_metaprofiles <- metaprofile_psite(reads_psite_list[[1]]%>%filter(length%in%c(lvect))%>%split(.$length)%>%setNames(lennames), 
		riboWaltzanno, sample = lennames,
                                         utr5l = 20, cdsl = 40, utr3l = 20,
                                         plot_title = "auto")
rl_metaprofiles

# psite_offset%>%filter(length==28)

# filt_metaprofiles <- metaprofile_psite(
# 	reads_psite_list[[1]]%>%filter(length==28)%>%filter(between(psite-cds_start,-3,0))%>%list(filt=.), 
# 		riboWaltzanno, sample = 'filt',
#                                          utr5l = 20, cdsl = 40, utr3l = 20,
#                                          plot_title = "auto")
# filt_metaprofiles

#so, here are the two peaks


{
oldstartreads = reads_psite_list[[1]]%>%filter(length==30)%>%
	#filter(between(psite-cds_start,-3,0))
	identity


oldstartreads%<>%filter(end5>4)%>%left_join(riboWaltzanno)%>%filter(l_tr-end3 > 4)%>%select(-l_tr,-l_utr5,-l_cds,-l_utr3)

fpseq = oldstartreads%>%
	select(seqnames=transcript,start=end5,end=end3)%>%
	{GRanges(.)}%>%
	resize(2)%>%
	resize(4,'end')%>%
	getSeq(x=FaFile('ext_data/gencode.vM12.pc_transcripts_ntrim.fa'))
tpseq = oldstartreads%>%
	select(seqnames=transcript,start=end5,end=end3)%>%
	{GRanges(.)}%>%
	resize(2,'end')%>%
	resize(4,'start')%>%
	getSeq(x=FaFile('ext_data/gencode.vM12.pc_transcripts_ntrim.fa'))

# test_that({
# 	reads_list_raw[[1]]%>%filter(transcript=='ENSMUST00000000033.9')
	
# 	seqinfo(FaFile('ext_data/gencode.vM12.pc_transcripts_ntrim.fa'))%>%.['ENSMUST00000000033.9']

# getSeq('ENSMUST00000000033',x=FaFile('ext_data/gencode.vM12.pc_transcripts_ntrim.fa'))

# getSeq(GRanges('ENSMUST00000000033:1'),x=FaFile('ext_data/gencode.vM12.pc_transcripts_ntrim.fa'))

# })

fpseq=fpseq%>%as.matrix%>%{set_colnames(.,paste0('fpseq',1:ncol(.)))}
tpseq=tpseq%>%as.matrix%>%{set_colnames(.,paste0('tpseq',1:ncol(.)))}

oldstartreads%<>%cbind(fpseq)
oldstartreads%<>%cbind(tpseq)

oldstartreads%<>%mutate(augpos = as.factor(psite-cds_start))
oldstartreads%<>%mutate(phase = as.factor(factnum(augpos) %%3))
# #look at freqs
# startreads
# 	filter((augpos%in%c(-3,0)))%>%
# 	group_by(augpos)%>%summarise(table(fpseq3)%>%enframe('nucl','freq'))%>%
# 	mutate(freq=as.numeric(freq))%>%
# 	mutate(freq=freq/sum(freq))%>%
# 	spread(nucl,freq)

#startreads$oldpsite = startreads$psite
# save.image()
}

trainfrac=0.5
stop()

{
startreads=oldstartreads%>%
	# filter((psite_from_start%>%between(-20,40))|psite_from_stop%>%between(-40,20))
	identity

library(ranger)
message('get sample set')
traintrs = startreads$transcript%>%unique%>%{sample(.,size=floor(length(.)*trainfrac))}
testtrs = startreads$transcript%>%unique%>%setdiff(traintrs)
startreadstrain=startreads%>%filter(transcript%in%traintrs)%>%filter(between(factnum(augpos),-3,3))
startreadstrain$augpos%<>%as.character%<>%as.factor
startreadtest=startreads%>%filter(transcript%in%testtrs)
oldstartreadtest=startreadtest
message('train...')
varstouse = colnames(startreadstrain)%>%str_subset('fpseq\\d|tpseq\\d')
varstouse = colnames(startreadstrain)%>%str_subset('fpseq')
rangermodel = ranger(data=startreadstrain%>%select(one_of(varstouse),augpos), 
	augpos ~ .,
	mtry=function(x)length(varstouse)
)
# seqshiftmodel <- holdoutRF(formula= augpos ~ . ,
				# data=	
			# )
message('predict...')
offsetadj = as.numeric(as.character(predict(rangerpred,data=startreadtest)$prediction))
# offsetadj = factnum(startreadtest$augpos)

startreadtest$psite = startreadtest$psite - factnum(offsetadj)
startreadtest$psite_from_start = startreadtest$psite - startreadtest$cds_start
startreadtest$psite_from_stop = startreadtest$psite - startreadtest$cds_stop
}
startreadstrain$augpos%>%table%>%sort

{
adjsamples=c('noadj','adj')
adjreadlist = list(oldstartreadtest,startreadtest)%>%setNames(adjsamples)

filt_metaprofiles <- metaprofile_psite(
	adjreadlist, 
		riboWaltzanno, sample = adjsamples,
                                         utr5l = 20, cdsl = 40, utr3l = 20,
                                         plot_title = "auto")
ggpubr::ggarrange(plotlist=filt_metaprofiles[-1],ncol=1)
}

frameob = frame_psite(adjreadlist, sample = adjsamples, region = "all")
frameob[[1]]%>%filter(region=='CDS')

codonusageold = codon_usage_psite(list('a'=oldstartreadtest),riboWaltzanno,sample='a',site='asite',fastapath='ext_data/gencode.vM12.pc_transcripts_ntrim.fa',fasta_genome=F)
codonusagenew = codon_usage_psite(list('a'=startreadtest),riboWaltzanno,sample='a',site='asite',fastapath='ext_data/gencode.vM12.pc_transcripts_ntrim.fa',fasta_genome=F)

codonusageold$dt
codonusagenew$dt

codonusageold$dt$raw_value%>%sd
codonusagenew$dt$raw_value%>%sd

bamtime = transcriptbam%>%dirname%>%basename%>%str_extract('[^_]+')


codonusageoldtrnadf = codonusageold$dt%>%left_join(allcodsigmean_isomerge%>%filter(fraction=='Poly',rep=='rep1')%>%filter(time==bamtime))
	# filter(normalized_value<)%>%
scattertitle = codonusageoldtrnadf%>%{cor.test(use='complete',(.[['normalized_value']]),(.[['abundance']]))}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
		{str_interp('pearsons rho =\n${.$estimate}(${.$conf.low} - ${.$conf.high})')}
codonusageoldtrnadf%>%{qplot(data=.,raw_value,abundance,label=codon,geom='blank')+geom_text()+ggtitle(scattertitle)}


codonusageoldtrnadf = codonusagenew$dt%>%left_join(allcodsigmean_isomerge%>%filter(fraction=='Poly',rep=='rep1')%>%filter(time==bamtime))
	# filter(normalized_value<)%>%
scattertitle = codonusageoldtrnadf%>%{cor.test(use='complete',(.[['normalized_value']]),(.[['abundance']]))}%>%tidy%>%mutate_if(is.numeric,round,3)%>%
		{str_interp('pearsons rho =\n${.$estimate}(${.$conf.low} - ${.$conf.high})')}
codonusageoldtrnadf%>%{qplot(data=.,raw_value,abundance,label=codon,geom='blank')+geom_text()+ggtitle(scattertitle)}

oldcodondt <- readRDS('data/psitecovrles.rds')

# }

#okay so directly adjusting with augpos works
#adjusting with 

# message(i);i=i+1
# readlens <- reads_psite_list[[1]]$length%>%table%>%keep(~.>1e3)%>%names%>%sort

# example_metaprofile_i<- rep(NA,length(readlens))%>%setNames(readlens)
# example_metaheatmap_compi <- rep(NA,length(readlens))%>%setNames(readlens)





# example_metaheatmap_compi%>%names

# message(i);i=i+1
# mytheme <- gridExtra::ttheme_default(
#     core = list(fg_params=list(cex = 0.5)),
#     colhead = list(fg_params=list(cex = 0.5)),
#     rowhead = list(fg_params=list(cex = 0.5)))


# ribowaltzpdf <- here('plots','riboWaltz',sampnames[[1]],'_ribowaltplots.pdf')
# ribowaltzpdf%>%dirname%>%dirname%>%dir.create
# ribowaltzpdf%>%dirname%>%dir.create

# pdf(ribowaltzpdf,w=14,h=7)
# rlength_distr(reads_list[1], sample = sampnames, cl = 99)[["plot"]]
# grid::grid.newpage()
# # example_ends_heatmap[["plot"]]
# gridExtra::grid.table(psite_offset, theme = mytheme,rows=NULL)
# # example_psite_region[["plot"]]
# example_frames_stratified[["plot"]]
# example_metaprofile[["plot"]]
# example_frames[["plot"]]

# metaprofile_psite%>%debug
# metaprofile_psite%>%undebug

# for(readlen in readlens%>%keep(~.%in%c(25:30))){
# 	try({
# 	message(paste0('comparison plots for readlength:',readlen))

# 	reads_psite_readlen <- reads_psite_list[[sampnames[[1]]]][length == as.numeric(readlen)]

# 	example_metaprofile_i <- metaprofile_psite(setNames(list(reads_psite_readlen),sampnames[[1]]), riboWaltzanno, sample = sampnames[[1]],
# 	                                            length_range = 25:31, utr5l = 20, cdsl = 60, 
# 	                                            transcripts = reads_psite_readlen$transcript%>%unique,
# 	                                            utr3l = 20, plot_title = "auto")

# 	print(example_metaprofile_i[['plot']])

# 	comparison_dt <- list()
	
	
# 	comparison_dt[[paste0("subsample_",readlen,"nt")]] <- reads_psite_readlen
# 	comparison_dt[["whole_sample"]] <- reads_psite_list[[sampnames[[1]]]]

# 	names_list <- list( paste0("subsample_",readlen,"nt"),"whole_sample" )%>%setNames(c(paste0("Only_",readlen),'All'))

# 	scale_facts <- comparison_dt%>%map_dbl(~ 1e6 / nrow(.))%>%setNames(names(comparison_dt))
# 	# example_metaheatmap_compi <- metaheatmap_psite(comparison_dt, riboWaltzanno, sample = names_list,
# 	                                         # utr5l = 20, cdsl = 40, utr3l = 20, log = F, scale_factors=scale_facts)
# 	# print(example_metaheatmap_compi[['plot']])

# 	example_frames <- frame_psite(reads_psite_list%>%map(~.[length==as.numeric(readlen)]), sample = sampnames[[1]], region = "all")
# 	print(example_frames[["plot"]])

# 	})
# }

# dev.off()
# normalizePath(ribowaltzpdf)


# psite_offset%>%filter(outputtsv


# library(ORFik)

# txdb<-makeTxDbFromGFF(gtf)

# ik_shifts <- detectRibosomeShifts(topcdsreads, txdb, stop = TRUE)

# psite_offset%>%filter(length==28)
# psite_offset%>%filter(length==28)

# psite_offset%>%colnames
# ik_shifts%>%colnames

# ik_psite_offset<-psite_offset
# ik_psite_offset%<>%.[length%in%ik_shifts$fragment_length,]

# ik_psite_offset$offset_from_5 <- - ik_shifts$offsets_start[match(ik_psite_offset$length,ik_shifts$fragment_length)]
# ik_psite_offset$corrected_offset_from_5 <- - ik_shifts$offsets_start[match(ik_psite_offset$length,ik_shifts$fragment_length)]
# ik_psite_offset$offset_from_3 <- - ik_shifts$offsets_stop[match(ik_psite_offset$length,ik_shifts$fragment_length)]
# ik_psite_offset$corrected_offset_from_3 <- - ik_shifts$offsets_stop[match(ik_psite_offset$length,ik_shifts$fragment_length)]


# ik_psite_offset$offset_from_5 <- 17
# ik_psite_offset$corrected_offset_from_5 <- 17



# #Now offset with the ORFik offsets
# reads_psite_list_ik <- psite_info(reads_list, ik_psite_offset)


# ribowaltzpdf <- here('plots','riboWaltz',sampnames[[1]],'orfik_ribowaltplots.pdf')

# pdf(ribowaltzpdf,w=14,h=7)

# example_frames_stratified <- frame_psite_length(reads_psite_list_ik, sample = sampnames[[1]],
#                                                 region = "all", cl = 90)[['plot']]%>%print

# example_frames <- frame_psite(reads_psite_list_ik, sample = sampnames[[1]], region = "all")[['plot']]%>%print



# for(readlen in 25){
# 	try({
# 	message(paste0('comparison plots for readlength:',readlen))

# 	reads_psite_readlen <- reads_psite_list_ik[[sampnames[[1]]]][length == as.numeric(readlen)]

# 	example_metaprofile_i <- metaprofile_psite(setNames(list(reads_psite_readlen),sampnames[[1]]), riboWaltzanno, sample = sampnames[[1]],
# 	                                            length_range = 25:31, utr5l = 20, cdsl = 60, 
# 	                                            transcripts = reads_psite_readlen$transcript%>%unique,
# 	                                            utr3l = 20, plot_title = "auto")

# 	print(example_metaprofile_i[['plot']])

# 	comparison_dt <- list()
	
	
# 	comparison_dt[[paste0("subsample_",readlen,"nt")]] <- reads_psite_readlen
# 	comparison_dt[["whole_sample"]] <- reads_psite_list_ik[[sampnames[[1]]]]

# 	names_list <- list( paste0("subsample_",readlen,"nt"),"whole_sample" )%>%setNames(c(paste0("Only_",readlen),'All'))

# 	scale_facts <- comparison_dt%>%map_dbl(~ 1e6 / nrow(.))%>%setNames(names(comparison_dt))
# 	# example_metaheatmap_compi <- metaheatmap_psite(comparison_dt, riboWaltzanno, sample = names_list,
# 	                                         # utr5l = 20, cdsl = 40, utr3l = 20, log = F, scale_factors=scale_facts)
# 	# print(example_metaheatmap_compi[['plot']])

# 	example_frames <- frame_psite(reads_psite_list_ik%>%map(~.[length==as.numeric(readlen)]), sample = sampnames[[1]], region = "all")
# 	print(example_frames[["plot"]])

# 	})
# }

# dev.off()
# normalizePath(ribowaltzpdf)




#' AImless, tired, stressed
#' okay, so take the read length... 28nts - my hypothesis is I can make the plot better for 28nt by choosing the offsets 
#' differently for different basepairs 
#' also puzzles me why some offsets appear off center...
#' ?Do the two different peaks for say 28 nt seem to hve different 
#' Are the ribowaltz plots both showing the same shift in the relative to start and stop ones?
#' 
#' There is soo much of this.
#' Start and stop offsets tooo, I can imagine choosing one or the other for certain 
#' 