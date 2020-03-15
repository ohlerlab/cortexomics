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

#load arguments
argv <- c(
	transcriptbam = here('pipeline/star/data/E16_ribo_2/E16_ribo_2.star_transcript.bam'),
	gtf=here('annotation/gencode.vM12.annotation.gtf'),
	outputtsv='riboWaltz/data/E16_ribo_2/offsets.tsv'
)
argv <- commandArgs(trailingOnly=TRUE)[1:length(argv)]%>%setNames(names(argv))
for(i in names(argv)) assign(i,argv[i])





RUST=0
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
sampnames <- basename(dirname(transcriptbam))

dirname(transcriptbam)


if(!exists('reads_list_raw')) reads_list_raw <- bamtolist(bamfolder = dirname(transcriptbam), name_samples = sampnames,annotation = riboWaltzanno)

reads_list<-reads_list_raw

#do rust if needed

if(RUST) {

	nureads_list<- reads_list
	ureads_list<-reads_list %>% map(~ distinct(.,transcript,end5,end3,.keep_all=TRUE))

}

reads_list%<>%setNames(dirname(transcriptbam)%>%basename)


reads_psite_list <- psite_info(reads_list, psite_offset)

example_frames <- frame_psite(reads_psite_list, sample = sampnames[[1]], region = "all")#offsets


psite_offset <- psite(reads_list, flanking = 6, extremity = "auto")

psite_offset%>%filter(total_percentage>5)%>%transmute(read_length=length,cutoff=corrected_offset_from_5)%>%
	write_tsv(outputtsv)



# #get the psites

# #
# # example_psite_region <- region_psite(reads_psite_list, riboWaltzanno, sample = sampnames)

# # psite_cds <- psite_per_cds(reads_psite_list, riboWaltzanno)

# example_frames_stratified <- frame_psite_length(reads_psite_list, sample = sampnames[[1]],
#                                                 region = "all", cl = 90)

# example_frames <- frame_psite(reads_psite_list, sample = sampnames[[1]], region = "all")


# example_metaprofile <- metaprofile_psite(reads_psite_list, mm81cdna, sample = sampnames[[1]],
#                                          utr5l = 20, cdsl = 40, utr3l = 20,
#                                          plot_title = "auto")


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