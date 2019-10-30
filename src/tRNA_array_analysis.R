
allread <- Sys.glob('ext_data/Matthew Kraushar J_030719-tRNA-PCR-22-10212019/PCR Array Service Report/*')%>%setNames(.,basename(.))%>%lapply(safely(.%>%readxl::read_xlsx(.,sheet=7)))

allread%<>%map('result')%>%bind_rows(.id='table')

	filter(!`...1`=='Transcript name')%>%

allread%>%mutate(sig = as.numeric(`T-TEST`)<0.05)%>%
		filter(!is.na(sig))%>%
		group_by(table,sig)%>%tally%>%spread(sig,n)%>%as.data.frame



cdstousetrna<-cds%>%split(.,.$protein_id)%>%.[highprotein_ids]%>%unlist

codons <- getSeq(cdstousetrna,x=FaFile(REF))%>%split(cdstousetrna$protein_id)%>%lapply(.%>%unlist%>%xscat)%>%DNAStringSet%>%vmatchPattern('TCT',.)%>%unlist%>%subset((start %% 3) == 1)

codonsgr <- GRanges(names(codons),codons)%>%mapFromTranscripts(cdstousetrna%>%split(.,.$protein_id))%>%subset(strand=='+')

codonsgr_exp <- codonsgr%>%resize(3+9+9,'center')


signalovercodons<-bams%>%str_subset('ribo')%>%mclapply(function(bam){
	riboparam<-ScanBamParam(scanBamFlag(isDuplicate=FALSE,isSecondaryAlignment=FALSE),mapqFilter=MAPQTHRESH,which=codonsgr_exp)
	reads <- readGAlignments(bam,param=riboparam)
	mcols(reads)$cdsshift <- get_cds_offsets(reads,psite_model$offsets,psite_model$compartments)
	#reads <- reads%>%subset(width==27)
	reads <- reads%>%subset(width %in% psite_model$offsets$length)
	mcols(reads)$length <- width(reads)
	reads%<>%subset(!is.na(cdsshift))
	psites <- apply_psite_offset(reads,c('cdsshift'))%>%as("GRanges")
	mcols(psites)$length <- mcols(reads)$length
	startwindpsites <- psites%>%mapToTranscripts(codonsgr_exp)
	start(startwindpsites)
})

codplotdf <- signalovercodons%>%
	lapply(table)%>%
	setNames(basename(bams%>%str_subset('ribo')))%>%
	map_df(.id='bam',enframe,'pos','count')%>%
	group_by(bam)%>%
	mutate(count=count/sum(count))





stagecols <- c(E12.5='#214098',E14='#2AA9DF',E15.5='#F17E22',E17='#D14E28',P0='#ED3124')
stageconv = names(stagecols)%>%setNames(c('E13','E145','E16','E175','P0'))
bamcols <- bams%>%basename%>%str_extract('[^_]+')%>%stageconv[.]%>%stagecols[.]%>%setNames(basename(bams))


plotfile<-'plots/TCT_ribo_cov.pdf'%T>%pdf(h=4,w=8)
qplot(data=codplotdf,color=bam,x=as.numeric(pos)-10,y=count,geom='line')+
	scale_x_continuous(name='Position Relative to TCT')+
	# scale_color_discrete(guide=TRUE)+
	scale_color_manual(values=bamcols)+
	scale_y_continuous(name='Relative P Site Density')+
	theme_bw()
dev.off()
normalizePath(plotfile)%>%message

#TODO - make sure strand spec 