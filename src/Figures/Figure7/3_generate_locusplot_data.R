################################################################################
########This program will generate the data needed for the locuplots for the shiny
########app
################################################################################
	
COMPRESS=T
allgenes = readxl::read_xlsx('tables/S1.xlsx',sheet=6)$gene_name
donegenes = Sys.glob(str_interp('data/Shiny_track_data/*.rds'))%>%basename%>%str_replace('.rds','')
genes2do = allgenes%>%setdiff(donegenes)
# for(igene in c('Satb2','Flna','Nes','Bcl11b','Tle4')[1:2]){
for(igene in genes2do){
	try({
	types = c('UTR','CDS')
	if(!COMPRESS) types = c(types,'gene')


	igeneanno <- fread(str_interp('grep -ie ${igene} ${annogtf}'))%>%
		as.matrix%>%apply(1,paste0,collapse='\t')%>%paste0(collapse='\n')%>%import(text=.,format='gtf')%>%
		subset(gene_name==igene)%>%
		subset(type%in%types)
	stopifnot(length(igeneanno)>0)
	if(COMPRESS){
		igeneanno$toshift <- get_comp_shift(igeneanno,maxwidth=10)
	}else{
		igeneanno$toshift <- 0
	}
	igeneanno$transcript_id%<>%trimids


	exontrack =
	  igeneanno%>%
	  shift(.,-.$toshift)%>%
	  # {.$transcript_id %<>% str_replace('\\.[0-9]+$','');.}%>%
	  subset(type%in%setdiff(types,'gene'))%>%
	  .[,c('type','transcript_id')]%>%
	  {.$feature<-as.character(.$type);.}%>%
	  {.$transcript<-as.character(.$transcript_id);.}%>%
	  Gviz::GeneRegionTrack(.,thinBoxFeature=c("UTR"))
	  
	tps = names(tpcols)

	gviztracks <- list(
		get_cov_track(igeneanno,totbams)%>%nametrack('RNA-Seq (RPM)'),
		get_cov_track(igeneanno,ribobams)%>%nametrack('Ribo-Seq (RPM)'),
		get_peptide_track(igeneanno)%>%nametrack('log2(Normalized Intensity)'),
	 	exontrack%>%nametrack('Transcripts')
	 	# getstart_track(igeneanno)%>%nametrack('AUG')
	 	# getmotiftrack('TGTANATA',igeneanno)%>%nametrack('Pum2 Motifs'),
	 	# peaktrack(zhangetal_pum2clip,igeneanno)%>%nametrack('Zhang et al Clip'),
	 	# peaktrack(pum1clip,igeneanno)%>%nametrack('PUM1 Clip')
	)
	gviztracks%<>%reflecttracks
	dir.create('data/Shiny_track_data/',showWarn=F,rec=T)
	saveRDS(gviztracks,str_interp('data/Shiny_track_data/${igene}.rds'))
	message(normalizePath(str_interp('data/Shiny_track_data/${igene}.rds')))
	})
}

gviztracks<-readRDS('data/Shiny_track_data/Satb2.rds')

	